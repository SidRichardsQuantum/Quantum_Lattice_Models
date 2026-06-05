"""Basic observables for small dense systems."""

from __future__ import annotations

import numpy as np

from quantum_lattice_models.operators import PAULI_Z, local_operator, two_site_operator


def expectation(state: np.ndarray, operator: np.ndarray) -> complex:
    """Return ``<state|operator|state>`` for a state vector."""

    vector = np.asarray(state, dtype=complex).reshape(-1)
    matrix = np.asarray(operator, dtype=complex)
    if matrix.shape != (vector.size, vector.size):
        raise ValueError("operator shape must match the state-vector dimension.")
    return complex(np.vdot(vector, matrix @ vector))


def magnetization_z(state: np.ndarray, n_sites: int) -> float:
    """Return average z magnetization per site."""

    values = [expectation(state, local_operator(PAULI_Z, i, n_sites)).real for i in range(n_sites)]
    return float(np.mean(values))


def correlation_zz(state: np.ndarray, n_sites: int, i: int, j: int) -> float:
    """Return the two-point ``<Z_i Z_j>`` correlation."""

    return float(expectation(state, two_site_operator(PAULI_Z, PAULI_Z, i, j, n_sites)).real)


def inverse_participation_ratio(vector: np.ndarray) -> float:
    """Return sum_i |v_i|^4 / (sum_i |v_i|^2)^2."""

    values = np.asarray(vector, dtype=complex).reshape(-1)
    probabilities = np.abs(values) ** 2
    norm = probabilities.sum()
    if norm == 0:
        raise ValueError("vector must have nonzero norm.")
    return float(np.sum(probabilities**2) / norm**2)
