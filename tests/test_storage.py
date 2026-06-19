from __future__ import annotations

import json

import numpy as np
import pytest
import scipy.sparse as sp

from quantum_lattice_models import (
    create_model_spec,
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
