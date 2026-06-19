from __future__ import annotations

import numpy as np
import pytest

from quantum_lattice_models.observables import (
    bipartite_entanglement_entropy,
    connected_spin_correlation_matrix,
    correlation_zz,
    expectation,
    inverse_participation_ratio,
    magnetization_z,
    reduced_density_matrix,
    site_magnetization_z,
    spin_correlation_matrix,
    static_spin_structure_factor,
    total_magnetization_z,
)
from quantum_lattice_models.operators import PAULI_Z, local_operator
from quantum_lattice_models.spin import fixed_magnetization_basis


def test_expectation_matches_known_z_value() -> None:
    zero = np.array([1.0, 0.0], dtype=complex)
    assert np.isclose(expectation(zero, PAULI_Z), 1.0)


def test_site_and_total_magnetization_for_product_state() -> None:
    state = np.zeros(8, dtype=complex)
    state[2] = 1.0  # |010>

    assert np.allclose(site_magnetization_z(state, n_sites=3), [1.0, -1.0, 1.0])
    assert np.isclose(total_magnetization_z(state, n_sites=3), 1.0)
    assert np.isclose(magnetization_z(state, n_sites=3), 1.0 / 3.0)
    assert np.isclose(correlation_zz(state, n_sites=3, i=0, j=1), -1.0)
    assert np.isclose(expectation(state, local_operator(PAULI_Z, 0, 3)), 1.0)


def test_bell_state_correlations_structure_factor_and_entropy() -> None:
    bell = np.array([1.0, 0.0, 0.0, 1.0], dtype=complex) / np.sqrt(2.0)

    correlations = spin_correlation_matrix(bell, n_sites=2, axis="Z")
    connected = connected_spin_correlation_matrix(bell, n_sites=2, axis="Z")
    factors = static_spin_structure_factor(
        bell,
        n_sites=2,
        momenta=[0.0, np.pi],
        connected=True,
    )
    density = reduced_density_matrix(bell, n_sites=2, subsystem=[0])

    assert np.allclose(correlations, np.ones((2, 2)))
    assert np.allclose(connected, np.ones((2, 2)))
    assert np.allclose(factors, [2.0, 0.0])
    assert np.allclose(density, np.eye(2) / 2.0)
    assert np.isclose(bipartite_entanglement_entropy(bell, 2, [0]), 1.0)


def test_product_state_has_zero_connected_correlations_and_entropy() -> None:
    state = np.zeros(4, dtype=complex)
    state[0] = 1.0

    assert np.allclose(connected_spin_correlation_matrix(state, 2), 0.0)
    assert np.isclose(bipartite_entanglement_entropy(state, 2, [0]), 0.0)
    assert np.allclose(reduced_density_matrix(state, 2, [0]), [[1.0, 0.0], [0.0, 0.0]])


@pytest.mark.parametrize("axis", ["X", "Y", "Z"])
def test_sector_observables_match_full_space_embedding(axis: str) -> None:
    basis = fixed_magnetization_basis(4, 0)
    reduced = np.arange(1, basis.dimension + 1, dtype=complex)
    reduced /= np.linalg.norm(reduced)
    full = basis.embed(reduced)

    assert np.allclose(
        site_magnetization_z(reduced, 4, basis=basis),
        site_magnetization_z(full, 4),
    )
    assert np.allclose(
        spin_correlation_matrix(reduced, 4, axis=axis, basis=basis),
        spin_correlation_matrix(full, 4, axis=axis),
    )
    assert np.allclose(
        connected_spin_correlation_matrix(reduced, 4, axis=axis, basis=basis),
        connected_spin_correlation_matrix(full, 4, axis=axis),
    )
    assert np.allclose(
        static_spin_structure_factor(
            reduced,
            4,
            [0.0, np.pi / 2.0, np.pi],
            axis=axis,
            basis=basis,
        ),
        static_spin_structure_factor(
            full,
            4,
            [0.0, np.pi / 2.0, np.pi],
            axis=axis,
        ),
    )


def test_sector_reduced_density_matrix_and_entropy_match_embedding() -> None:
    basis = fixed_magnetization_basis(2, 0)
    reduced = np.ones(2, dtype=complex) / np.sqrt(2.0)
    full = basis.embed(reduced)

    sector_density = reduced_density_matrix(reduced, 2, [0], basis=basis)
    full_density = reduced_density_matrix(full, 2, [0])

    assert np.allclose(sector_density, np.eye(2) / 2.0)
    assert np.allclose(sector_density, full_density)
    assert np.isclose(
        bipartite_entanglement_entropy(reduced, 2, [0], basis=basis),
        bipartite_entanglement_entropy(full, 2, [0]),
    )


def test_noncontiguous_subsystem_density_matrix_is_normalized() -> None:
    state = np.zeros(8, dtype=complex)
    state[5] = 3.0  # Unnormalized |101>
    density = reduced_density_matrix(state, 3, [2, 0])

    expected = np.zeros((4, 4), dtype=complex)
    expected[3, 3] = 1.0
    assert np.allclose(density, expected)
    assert np.isclose(np.trace(density), 1.0)


def test_observable_validation_failures() -> None:
    state = np.array([1.0, 0.0])
    basis = fixed_magnetization_basis(2, 0)

    with pytest.raises(ValueError, match="selected spin basis"):
        site_magnetization_z(state, 3)
    with pytest.raises(ValueError, match="basis.n_sites"):
        site_magnetization_z(state, 3, basis=basis)
    with pytest.raises(ValueError, match="axis"):
        spin_correlation_matrix(state, 1, axis="A")
    with pytest.raises(ValueError, match="unique"):
        reduced_density_matrix(state, 1, [0, 0])
    with pytest.raises(ValueError, match="nonzero norm"):
        reduced_density_matrix(np.zeros(2), 1, [0])
    with pytest.raises(ValueError, match="different from one"):
        bipartite_entanglement_entropy(state, 1, [0], base=1.0)


def test_inverse_participation_ratio() -> None:
    localized = np.array([1.0, 0.0, 0.0, 0.0])
    uniform = np.ones(4) / 2.0
    assert np.isclose(inverse_participation_ratio(localized), 1.0)
    assert np.isclose(inverse_participation_ratio(uniform), 0.25)
