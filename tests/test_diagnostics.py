from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from quantum_lattice_models.diagnostics import (
    diagnose_matrix,
    estimate_dense_memory,
    estimate_model_dimension,
    has_particle_hole_symmetric_spectrum,
    is_hermitian,
    matrix_density,
    matrix_storage_bytes,
)
from quantum_lattice_models.models import kitaev_chain_bdg, tight_binding_chain_sparse


def test_model_dimension_and_dense_memory_estimates() -> None:
    assert estimate_model_dimension("transverse_field_ising", n_sites=5) == 32
    assert estimate_model_dimension("heisenberg_ladder", n_rungs=3) == 64
    assert estimate_model_dimension("bose_hubbard_chain", n_sites=3, max_occupancy=2) == 27
    assert estimate_model_dimension("kagome_lattice_tight_binding", n_rows=2, n_cols=4) == 24
    assert estimate_dense_memory(4) == 4 * 4 * np.dtype(np.complex128).itemsize


def test_dense_and_sparse_matrix_diagnostics() -> None:
    dense = np.diag([1.0, 2.0, 3.0])
    sparse = sp.csr_matrix(dense)

    assert matrix_density(dense) == 1 / 3
    assert matrix_density(sparse) == 1 / 3
    assert matrix_storage_bytes(dense) == dense.nbytes
    assert matrix_storage_bytes(sparse) < matrix_storage_bytes(dense)
    assert is_hermitian(dense)
    assert is_hermitian(sparse)

    summary = diagnose_matrix(tight_binding_chain_sparse(n_sites=8))
    assert summary.shape == (8, 8)
    assert summary.sparse
    assert summary.nonzero_entries == 14
    assert summary.hermitian


def test_hermiticity_and_particle_hole_diagnostics() -> None:
    nonhermitian = np.array([[0.0, 1.0], [0.0, 0.0]])
    assert not is_hermitian(nonhermitian)
    assert not has_particle_hole_symmetric_spectrum(nonhermitian)

    kitaev = kitaev_chain_bdg(
        n_sites=5,
        hopping=1.0,
        chemical_potential=0.3,
        pairing=0.4,
    )
    assert has_particle_hole_symmetric_spectrum(kitaev)


def test_diagnostic_validation_failures() -> None:
    with pytest.raises(ValueError, match="positive integer"):
        estimate_dense_memory(0)
    with pytest.raises(ValueError, match="must be an integer"):
        estimate_model_dimension("tight_binding_chain", n_sites=2.5)
    with pytest.raises(ValueError, match="nonempty two-dimensional"):
        matrix_density(np.array([]))
    assert not is_hermitian(np.ones((2, 3)))
