from __future__ import annotations

import json

import numpy as np
import pytest

from quantum_lattice_models import (
    ModelSpec,
    create_model_spec,
    estimate_model_dimension,
    load_hamiltonian,
    save_hamiltonian,
    site_magnetization_z,
    spin_flip_parity_basis,
    symmetry_action_diagnostic,
    transverse_field_ising_parity_sector,
    transverse_field_ising_sparse,
)
from quantum_lattice_models.cli import main
from quantum_lattice_models.reduced import ReducedBasisMapping, reduced_operator


@pytest.mark.parametrize("parity", [-1, 1])
def test_spin_flip_parity_basis_is_an_orthonormal_eigenbasis(parity: int) -> None:
    basis = spin_flip_parity_basis(4, parity)
    transformation = basis.mapping.transformation()
    flip = np.fliplr(np.eye(16))

    assert basis.dimension == 8
    assert np.allclose((transformation.getH() @ transformation).toarray(), np.eye(8))
    assert np.allclose(flip @ transformation.toarray(), parity * transformation.toarray())

    reduced = np.arange(8, dtype=complex)
    assert np.allclose(basis.project(basis.embed(reduced)), reduced)
    restored = ReducedBasisMapping.from_dict(basis.mapping.to_dict())
    assert np.allclose(restored.transformation().toarray(), transformation.toarray())


@pytest.mark.parametrize("periodic", [False, True])
def test_tfim_parity_sectors_match_full_isometric_blocks(periodic: bool) -> None:
    full = transverse_field_ising_sparse(5, j=0.7, h=0.3, periodic=periodic)
    sector_values = []
    for parity in (-1, 1):
        sector = transverse_field_ising_parity_sector(
            5,
            parity,
            j=0.7,
            h=0.3,
            periodic=periodic,
        )
        expected = reduced_operator(full, sector.basis.mapping)
        assert np.allclose(sector.matrix.toarray(), expected.toarray())
        sector_values.extend(np.linalg.eigvalsh(sector.matrix.toarray()))

    assert np.allclose(
        np.sort(sector_values),
        np.linalg.eigvalsh(full.toarray()),
    )


def test_parity_sector_portable_symmetry_observables_and_storage(tmp_path) -> None:
    spec = create_model_spec(
        "transverse_field_ising_parity_sector_sparse",
        parameters={"n_sites": 4, "parity": -1, "j": 0.8, "h": 0.2},
    )
    result = spec.build_result()

    assert spec.symmetry_actions[0].name == "global_spin_flip"
    assert spec.symmetry_actions[0].metadata["selected_eigenvalue"] == -1
    assert result.metadata["sector"]["parity"] == -1
    assert (
        estimate_model_dimension(
            "transverse_field_ising_parity_sector_sparse",
            n_sites=7,
            parity=1,
        )
        == 64
    )

    sector = transverse_field_ising_parity_sector(4, -1, j=0.8, h=0.2)
    _, ground_state = np.linalg.eigh(result.matrix.toarray())
    magnetization = site_magnetization_z(
        ground_state[:, 0],
        4,
        basis=sector.basis,
    )
    assert np.allclose(magnetization, 0.0, atol=1e-10)

    diagnostic = symmetry_action_diagnostic(
        transverse_field_ising_sparse(4, j=0.8, h=0.2),
        spec,
        "global_spin_flip",
    )
    assert diagnostic.metadata["conserved"]

    path = save_hamiltonian(result, tmp_path / "parity.npz")
    restored = load_hamiltonian(path)
    assert restored.model.symmetry_actions == spec.symmetry_actions
    assert restored.metadata["sector"]["mapping"]["components"]

    encoded = spec.to_dict()
    assert ModelSpec.from_dict(encoded) == spec
    assert "symmetry_actions" in encoded


def test_parity_sector_cli_and_validation(tmp_path, capsys) -> None:
    path = tmp_path / "parity.json"
    assert (
        main(
            [
                "create",
                "transverse_field_ising_parity_sector_sparse",
                "--n-sites",
                "5",
                "--parity",
                "1",
                "--output",
                str(path),
            ]
        )
        == 0
    )
    capsys.readouterr()

    assert (
        main(["dry-run", "--model", "transverse_field_ising_parity_sector_sparse", "--json"]) == 0
    )
    report = json.loads(capsys.readouterr().out)
    assert report["dimension"] == 32
    assert report["representation"] == "sparse"

    assert main(["spectrum", str(path), "--json"]) == 0
    assert len(json.loads(capsys.readouterr().out)["eigenvalues"]) == 16

    with pytest.raises(ValueError, match="parity must be"):
        spin_flip_parity_basis(4, 0)
    with pytest.raises(ValueError, match="parity must be"):
        create_model_spec(
            "transverse_field_ising_parity_sector_sparse",
            parameters={"parity": 0},
        )
