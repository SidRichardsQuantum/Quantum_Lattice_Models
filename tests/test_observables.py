from __future__ import annotations

import numpy as np

from quantum_lattice_models.observables import (
    correlation_zz,
    expectation,
    inverse_participation_ratio,
    magnetization_z,
)
from quantum_lattice_models.operators import PAULI_Z, local_operator


def test_expectation_matches_known_z_value() -> None:
    zero = np.array([1.0, 0.0], dtype=complex)
    assert np.isclose(expectation(zero, PAULI_Z), 1.0)


def test_magnetization_and_correlation_for_all_zero_state() -> None:
    state = np.zeros(4, dtype=complex)
    state[0] = 1.0
    assert np.isclose(magnetization_z(state, n_sites=2), 1.0)
    assert np.isclose(correlation_zz(state, n_sites=2, i=0, j=1), 1.0)
    assert np.isclose(expectation(state, local_operator(PAULI_Z, 0, 2)), 1.0)


def test_inverse_participation_ratio() -> None:
    localized = np.array([1.0, 0.0, 0.0, 0.0])
    uniform = np.ones(4) / 2.0
    assert np.isclose(inverse_participation_ratio(localized), 1.0)
    assert np.isclose(inverse_participation_ratio(uniform), 0.25)
