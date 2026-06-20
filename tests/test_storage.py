from __future__ import annotations

import json

import numpy as np
import pytest
import scipy.sparse as sp

from quantum_lattice_models import (
    EXTERNAL_MATRIX_FAMILY,
    LatticeSpec,
    create_model_spec,
    export_hamiltonian_artifact,
    import_hamiltonian,
    load_hamiltonian,
    metadata_path,
    save_hamiltonian,
)
from quantum_lattice_models.cli import main


def test_dense_npy_round_trip_uses_metadata_sidecar(tmp_path) -> None:
    result = create_model_spec(
        "ssh_model",
        parameters={"n_cells": 3, "t1": 0.25},
    ).build_result()

    path = save_hamiltonian(result, tmp_path / "ssh", format="npy")
    restored = load_hamiltonian(path)

    assert path.suffix == ".npy"
    assert metadata_path(path).exists()
    assert np.allclose(restored.matrix, result.matrix)
    assert restored.model == result.model
    assert restored.basis == result.basis


def test_sparse_npz_round_trip_is_self_contained(tmp_path) -> None:
    result = create_model_spec(
        "tight_binding_chain",
        parameters={"n_sites": 6, "periodic": True},
        representation="sparse",
    ).build_result()

    path = save_hamiltonian(result, tmp_path / "chain.npz")
    restored = load_hamiltonian(path)

    assert sp.issparse(restored.matrix)
    assert np.allclose(restored.matrix.toarray(), result.matrix.toarray())
    assert restored.model == result.model
    with pytest.raises(ValueError, match="NPY export requires a dense matrix"):
        save_hamiltonian(result, tmp_path / "chain.npy")


def test_file_oriented_spectrum_and_export_cli(tmp_path, capsys) -> None:
    model_path = create_model_spec(
        "tight_binding_chain",
        parameters={"n_sites": 4},
        representation="sparse",
    ).save(tmp_path / "chain.json")

    assert main(["spectrum", str(model_path)]) == 0
    assert len(capsys.readouterr().out.strip().splitlines()) == 4

    output = tmp_path / "matrix.npz"
    assert main(["export", str(model_path), "--output", str(output)]) == 0
    assert capsys.readouterr().out.strip() == str(output)
    restored = load_hamiltonian(output)
    assert restored.shape == (4, 4)
    with np.load(output, allow_pickle=False) as archive:
        metadata = json.loads(str(archive["metadata"].item()))
    assert metadata["model"]["family"] == "tight_binding_chain"


def test_import_external_dense_matrix_with_portable_metadata(tmp_path) -> None:
    matrix = np.array([[0.0, 1.0 - 0.25j], [1.0 + 0.25j, 0.5]])
    source = tmp_path / "external.npy"
    np.save(source, matrix, allow_pickle=False)
    metadata = tmp_path / "external.json"
    metadata.write_text(
        json.dumps(
            {
                "basis": "single-particle site basis",
                "basis_dimension": 2,
                "lattice": {
                    "n_sites": 2,
                    "positions": [[0.0, 0.0], [1.0, 0.0]],
                    "bonds": [],
                },
                "local_degrees": [
                    {
                        "index": 0,
                        "site": 0,
                        "kind": "orbital",
                        "local_dimension": 2,
                        "label": "left",
                    },
                    {
                        "index": 1,
                        "site": 1,
                        "kind": "orbital",
                        "local_dimension": 2,
                        "label": "right",
                    },
                ],
                "basis_mappings": [
                    {
                        "local_degree": 0,
                        "basis_index": 0,
                        "role": "single_particle_state",
                    },
                    {
                        "local_degree": 1,
                        "basis_index": 1,
                        "role": "single_particle_state",
                    },
                ],
                "interactions": [
                    {
                        "degrees": [0, 1],
                        "operators": ["create", "annihilate"],
                        "coefficient": {"__complex__": [1.0, -0.25]},
                        "kind": "hopping",
                    }
                ],
                "units": {"energy": "eV"},
                "conventions": {"basis_ordering": "site index"},
                "references": ["doi:10.0000/example"],
                "provenance": {"generator": "external-code"},
                "metadata": {"sample": "A"},
            }
        )
    )

    result = import_hamiltonian(source, metadata_path=metadata)
    portable = save_hamiltonian(result, tmp_path / "imported.npz")
    restored = load_hamiltonian(portable)

    assert np.allclose(restored.matrix, matrix)
    assert restored.model.family == EXTERNAL_MATRIX_FAMILY
    assert restored.model.basis == "single-particle site basis"
    assert restored.model.lattice is not None
    assert restored.model.units == {"energy": "eV"}
    assert len(restored.model.local_degrees) == 2
    assert restored.model.interactions[0].coefficient == 1.0 - 0.25j
    assert restored.model.metadata["matrix_dimension"] == 2
    assert restored.model.provenance["generator"] == "external-code"
    assert restored.metadata["imported"] is True
    with pytest.raises(ValueError, match="cannot reconstruct"):
        restored.model.hamiltonian()


def test_import_external_sparse_npz_and_validation(tmp_path) -> None:
    source = tmp_path / "external-sparse.npz"
    sp.save_npz(source, sp.csr_matrix([[1.0, 0.5], [0.5, -1.0]]))
    metadata = tmp_path / "metadata.json"
    metadata.write_text(json.dumps({"basis": "custom occupation basis"}))

    result = import_hamiltonian(source, metadata_path=metadata)

    assert sp.issparse(result.matrix)
    assert result.representation == "sparse"
    assert np.allclose(result.matrix.toarray(), [[1.0, 0.5], [0.5, -1.0]])

    non_square = tmp_path / "non-square.npy"
    np.save(non_square, np.ones((2, 3)), allow_pickle=False)
    with pytest.raises(ValueError, match="nonempty and square"):
        import_hamiltonian(non_square, metadata_path=metadata)

    non_hermitian = tmp_path / "non-hermitian.npy"
    np.save(non_hermitian, np.array([[0.0, 1.0], [0.0, 0.0]]), allow_pickle=False)
    with pytest.raises(ValueError, match="must be Hermitian"):
        import_hamiltonian(non_hermitian, metadata_path=metadata)
    accepted = import_hamiltonian(
        non_hermitian,
        metadata_path=metadata,
        require_hermitian=False,
    )
    assert accepted.metadata["hermitian"] is False

    mismatched_metadata = tmp_path / "mismatched-metadata.json"
    mismatched_metadata.write_text(
        json.dumps({"basis": "custom occupation basis", "basis_dimension": 3})
    )
    with pytest.raises(ValueError, match="basis_dimension"):
        import_hamiltonian(source, metadata_path=mismatched_metadata)


def test_import_matrix_cli_supports_analysis_and_reexport(tmp_path, capsys) -> None:
    source = tmp_path / "external.npy"
    np.save(source, np.diag([-1.0, 2.0]), allow_pickle=False)
    metadata = tmp_path / "metadata.json"
    metadata.write_text(json.dumps({"basis": "two-state external basis"}))
    imported = tmp_path / "portable.npz"

    assert (
        main(
            [
                "import-matrix",
                str(source),
                "--metadata",
                str(metadata),
                "--output",
                str(imported),
            ]
        )
        == 0
    )
    assert capsys.readouterr().out.strip() == str(imported)

    assert main(["spectrum", str(imported), "--json"]) == 0
    spectrum = json.loads(capsys.readouterr().out)
    assert spectrum["eigenvalues"] == [-1.0, 2.0]

    assert main(["inspect", str(imported)]) == 0
    inspection = json.loads(capsys.readouterr().out)
    assert inspection["family"] == EXTERNAL_MATRIX_FAMILY
    assert inspection["matrix"]["hermitian"] is True

    exported = tmp_path / "reexported.npy"
    assert (
        main(
            [
                "export",
                str(imported),
                "--format",
                "npy",
                "--output",
                str(exported),
            ]
        )
        == 0
    )
    capsys.readouterr()
    assert np.allclose(load_hamiltonian(exported).matrix, np.diag([-1.0, 2.0]))


def test_load_rejects_matrix_metadata_mismatch(tmp_path) -> None:
    result = create_model_spec(
        "tight_binding_chain",
        parameters={"n_sites": 2},
    ).build_result()
    path = save_hamiltonian(result, tmp_path / "matrix.npz")
    with np.load(path, allow_pickle=False) as archive:
        arrays = {name: archive[name] for name in archive.files}
    metadata = json.loads(str(arrays["metadata"].item()))
    metadata["matrix_shape"] = [3, 3]
    arrays["metadata"] = np.array(json.dumps(metadata))
    np.savez_compressed(path, **arrays)

    with pytest.raises(ValueError, match="shape does not match"):
        load_hamiltonian(path)


def test_artifact_exports_are_deterministic_and_round_trippable(tmp_path) -> None:
    spec = create_model_spec(
        "custom_tight_binding",
        lattice=LatticeSpec(
            n_sites=2,
            positions=((0.0, 0.0), (1.0, 0.0)),
            bonds=(),
            units={"position": "angstrom"},
        ),
        parameters={"hopping": 1.0, "onsite": 0.25, "hermitian": True},
        units={"energy": "eV"},
    )
    result = spec.build_result()

    model_path = export_hamiltonian_artifact(spec, tmp_path / "selected-model", artifact="model")[0]
    lattice_path = export_hamiltonian_artifact(
        spec, tmp_path / "selected-lattice", artifact="lattice"
    )[0]
    metadata = export_hamiltonian_artifact(
        result, tmp_path / "selected-metadata", artifact="metadata"
    )[0]
    bundle = export_hamiltonian_artifact(
        result, tmp_path / "bundle", artifact="bundle", format="npz"
    )

    assert model_path.name == "selected-model.json"
    assert lattice_path.name == "selected-lattice.json"
    assert json.loads(model_path.read_text()) == spec.to_dict()
    assert json.loads(lattice_path.read_text()) == spec.lattice.to_dict()
    assert json.loads(metadata.read_text()) == result.to_metadata()
    assert [path.name for path in bundle] == [
        "matrix.npz",
        "model.json",
        "metadata.json",
        "lattice.json",
        "manifest.json",
    ]
    manifest = json.loads((tmp_path / "bundle" / "manifest.json").read_text())
    assert manifest["files"] == [
        "matrix.npz",
        "model.json",
        "metadata.json",
        "lattice.json",
    ]
    restored = load_hamiltonian(tmp_path / "bundle" / "matrix.npz")
    assert np.allclose(restored.matrix, result.matrix)
    assert restored.model == spec


def test_artifact_export_validation(tmp_path) -> None:
    spec = create_model_spec("tight_binding_chain", parameters={"n_sites": 2})

    with pytest.raises(ValueError, match="portable lattice"):
        export_hamiltonian_artifact(spec, tmp_path / "lattice", artifact="lattice")
    with pytest.raises(ValueError, match="constructed Hamiltonian"):
        export_hamiltonian_artifact(spec, tmp_path / "matrix", artifact="matrix")
    with pytest.raises(ValueError, match="Artifact must be"):
        export_hamiltonian_artifact(spec, tmp_path / "unknown", artifact="unknown")


def test_npy_bundle_includes_sidecar_and_removes_stale_artifacts(tmp_path) -> None:
    result = create_model_spec(
        "tight_binding_chain",
        parameters={"n_sites": 2},
    ).build_result()
    bundle_path = tmp_path / "bundle"
    bundle_path.mkdir()
    (bundle_path / "lattice.json").write_text("{}")
    (bundle_path / "unmanaged.txt").write_text("keep")

    outputs = export_hamiltonian_artifact(
        result,
        bundle_path,
        artifact="bundle",
        format="npy",
    )

    assert [path.name for path in outputs] == [
        "matrix.npy",
        "matrix.npy.json",
        "model.json",
        "metadata.json",
        "manifest.json",
    ]
    assert not (bundle_path / "lattice.json").exists()
    assert (bundle_path / "unmanaged.txt").read_text() == "keep"
    assert np.allclose(load_hamiltonian(bundle_path / "matrix.npy").matrix, result.matrix)


def test_cli_selective_artifact_exports_preserve_default_behavior(tmp_path, capsys) -> None:
    model_path = create_model_spec(
        "custom_tight_binding",
        lattice=LatticeSpec(n_sites=2, bonds=()),
        parameters={"hopping": 1.0, "onsite": 0.0, "hermitian": True},
    ).save(tmp_path / "model.json")

    matrix = tmp_path / "matrix.npz"
    assert main(["export", str(model_path), "--output", str(matrix)]) == 0
    assert capsys.readouterr().out.strip() == str(matrix)
    assert load_hamiltonian(matrix).shape == (2, 2)

    lattice = tmp_path / "lattice.json"
    assert (
        main(
            [
                "export",
                str(model_path),
                "--artifact",
                "lattice",
                "--output",
                str(lattice),
            ]
        )
        == 0
    )
    assert capsys.readouterr().out.strip() == str(lattice)
    assert json.loads(lattice.read_text())["n_sites"] == 2

    bundle = tmp_path / "portable.bundle"
    assert (
        main(
            [
                "export",
                str(matrix),
                "--artifact",
                "bundle",
                "--output",
                str(bundle),
            ]
        )
        == 0
    )
    outputs = capsys.readouterr().out.strip().splitlines()
    assert outputs[-1] == str(bundle / "manifest.json")
    assert load_hamiltonian(bundle / "matrix.npz").shape == (2, 2)
