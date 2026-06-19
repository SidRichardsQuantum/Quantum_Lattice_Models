from __future__ import annotations

import numpy as np
import pytest
import scipy.sparse as sp

from quantum_lattice_models.operators import PAULI_X, PAULI_Y, local_operator, two_site_operator
from quantum_lattice_models.spin import (
    SpinField,
    SpinInteraction,
    graph_spin_hamiltonian,
    graph_spin_hamiltonian_sparse,
    heisenberg_chain,
    heisenberg_chain_sparse,
    heisenberg_ladder,
    heisenberg_ladder_sparse,
    j1_j2_heisenberg_chain,
    j1_j2_heisenberg_chain_sparse,
    longitudinal_field_ising,
    longitudinal_field_ising_sparse,
    next_nearest_neighbor_ising,
    next_nearest_neighbor_ising_sparse,
    transverse_field_ising,
    transverse_field_ising_sparse,
    xxz_chain,
    xxz_chain_sparse,
    xy_chain,
    xy_chain_sparse,
)


@pytest.mark.parametrize(
    ("dense_builder", "sparse_builder", "parameters"),
    [
        (transverse_field_ising, transverse_field_ising_sparse, {"n_sites": 3}),
        (longitudinal_field_ising, longitudinal_field_ising_sparse, {"n_sites": 3}),
        (
            next_nearest_neighbor_ising,
            next_nearest_neighbor_ising_sparse,
            {"n_sites": 4, "periodic": True},
        ),
        (heisenberg_chain, heisenberg_chain_sparse, {"n_sites": 3, "periodic": True}),
        (xy_chain, xy_chain_sparse, {"n_sites": 3, "anisotropy": 0.2}),
        (xxz_chain, xxz_chain_sparse, {"n_sites": 3, "anisotropy": 0.7}),
        (j1_j2_heisenberg_chain, j1_j2_heisenberg_chain_sparse, {"n_sites": 4}),
        (heisenberg_ladder, heisenberg_ladder_sparse, {"n_rungs": 2}),
    ],
)
def test_named_sparse_spin_builders_match_dense(dense_builder, sparse_builder, parameters) -> None:
    dense = dense_builder(**parameters)
    sparse = sparse_builder(**parameters)

    assert sp.issparse(sparse)
    assert np.allclose(sparse.toarray(), dense)


def test_graph_spin_builder_supports_mixed_axes_and_site_fields() -> None:
    interactions = (SpinInteraction(0, 1, "X", "Y", 0.25),)
    fields = (SpinField(1, "X", -0.5),)

    sparse = graph_spin_hamiltonian_sparse(2, interactions, fields)
    dense = graph_spin_hamiltonian(2, interactions, fields)
    expected = 0.25 * two_site_operator(PAULI_X, PAULI_Y, 0, 1, 2)
    expected -= 0.5 * local_operator(PAULI_X, 1, 2)

    assert np.allclose(sparse.toarray(), expected)
    assert np.allclose(dense, expected)
    assert len(dense.terms) == 2


def test_graph_spin_builder_rejects_invalid_graph_records() -> None:
    with pytest.raises(ValueError, match="distinct"):
        graph_spin_hamiltonian_sparse(2, [SpinInteraction(0, 0, "X", "X", 1.0)])
    with pytest.raises(ValueError, match="axes"):
        graph_spin_hamiltonian_sparse(2, fields=[SpinField(0, "A", 1.0)])
