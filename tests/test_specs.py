from __future__ import annotations

import json

import numpy as np
import pytest
import scipy.sparse as sp

from quantum_lattice_models import (
    LatticeSpec,
    ModelSpec,
    create_model_spec,
    load_model,
)
from quantum_lattice_models.cli import main
from quantum_lattice_models.lattice import Bond


def test_lattice_spec_round_trip_preserves_complex_bonds() -> None:
    original = LatticeSpec(
        n_sites=3,
        positions=((0.0, 0.0), (1.0, 0.0), (0.5, 0.8)),
        bonds=(Bond(0, 1), Bond(1, 2, 0.25j)),
        sublattice_labels=("A", "B", "A"),
        unit_cells=(0, 0, 1),
        boundary_conditions={"x": "open"},
        metadata={"label": "triangle fragment"},
    )

    restored = LatticeSpec.from_dict(original.to_dict())

    assert restored == original
    assert restored.to_lattice().bonds[1].value == 0.25j


def test_model_spec_json_round_trip_and_dense_sparse_construction(tmp_path) -> None:
    spec = create_model_spec(
        "tight_binding_chain",
        parameters={"n_sites": 5, "hopping": 0.7, "periodic": True},
    )
    path = spec.save(tmp_path / "chain.json")
    restored = load_model(path)

    dense = restored.hamiltonian()
    sparse = restored.hamiltonian(sparse=True)

    assert restored == spec
    assert dense.shape == (5, 5)
    assert sp.issparse(sparse)
    assert np.allclose(sparse.toarray(), dense)
    assert json.loads(path.read_text())["schema_version"] == "1.0"
    assert restored.parameters["onsite"] == 0.0

    onsite_spec = create_model_spec(
        "tight_binding_chain",
        parameters={"n_sites": 3, "onsite": (0.1, 0.2, 0.3)},
    )
    assert np.allclose(np.diag(onsite_spec.hamiltonian()), [0.1, 0.2, 0.3])

    kitaev = create_model_spec(
        "kitaev_chain_bdg",
        parameters={"n_sites": 3, "pairing": 0.25j},
    )
    kitaev_path = kitaev.save(tmp_path / "kitaev.json")
    assert load_model(kitaev_path).parameters["pairing"] == 0.25j


def test_custom_lattice_model_spec_builds_dense_and_sparse() -> None:
    lattice = LatticeSpec(
        n_sites=3,
        bonds=(Bond(0, 1), Bond(1, 2, 0.5j)),
    )
    spec = create_model_spec(
        "custom_tight_binding",
        parameters={"hopping": 0.8, "hermitian": True},
        lattice=lattice,
    )

    dense = spec.hamiltonian()
    sparse = spec.hamiltonian(sparse=True)

    assert dense[0, 1] == -0.8
    assert dense[1, 2] == 0.5j
    assert np.allclose(sparse.toarray(), dense)


def test_model_spec_validation_rejects_invalid_schema_and_representation() -> None:
    with pytest.raises(ValueError, match="schema_version"):
        ModelSpec(family="ssh_model", schema_version="2.0").validate()
    with pytest.raises(ValueError, match="representation"):
        ModelSpec(family="ssh_model", representation="gpu").validate()
    with pytest.raises(ValueError, match="does not provide a sparse builder"):
        ModelSpec(
            family="transverse_field_ising",
            parameters={"n_sites": 2},
            representation="sparse",
        ).validate()
    with pytest.raises(ValueError, match="must have type int"):
        ModelSpec(
            family="ssh_model",
            parameters={"n_cells": "four"},
        ).validate()
    with pytest.raises(ValueError, match="registered basis"):
        ModelSpec(
            family="ssh_model",
            parameters={"n_cells": 2},
            basis="qubit",
        ).validate()


def test_create_inspect_and_validate_cli(tmp_path, capsys) -> None:
    path = tmp_path / "ssh.json"
    assert (
        main(
            [
                "create",
                "ssh_model",
                "--n-cells",
                "4",
                "--t1",
                "0.3",
                "--representation",
                "dense",
                "--output",
                str(path),
            ]
        )
        == 0
    )
    assert path.exists()
    capsys.readouterr()

    assert main(["inspect", str(path)]) == 0
    inspection = json.loads(capsys.readouterr().out)
    assert inspection["family"] == "ssh_model"
    assert inspection["dimension"] == 8

    assert main(["validate", str(path)]) == 0
    assert capsys.readouterr().out.startswith("valid\tssh_model\t1.0")
