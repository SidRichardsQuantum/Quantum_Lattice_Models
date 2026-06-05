"""Dense operators and Kronecker-product helpers."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

PAULI_X = np.array([[0, 1], [1, 0]], dtype=complex)
PAULI_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
PAULI_Z = np.array([[1, 0], [0, -1]], dtype=complex)
IDENTITY = np.eye(2, dtype=complex)

PAULI_MATRICES: dict[str, np.ndarray] = {
    "I": IDENTITY,
    "X": PAULI_X,
    "Y": PAULI_Y,
    "Z": PAULI_Z,
}


def kron_all(operators: Sequence[np.ndarray]) -> np.ndarray:
    """Return the Kronecker product of a sequence of matrices."""

    if not operators:
        raise ValueError("At least one operator is required.")
    result = np.asarray(operators[0], dtype=complex)
    for operator in operators[1:]:
        result = np.kron(result, np.asarray(operator, dtype=complex))
    return result


def local_operator(operator: np.ndarray, site: int, n_sites: int) -> np.ndarray:
    """Construct a one-site operator acting on a spin chain."""

    _validate_site(site, n_sites)
    factors = [IDENTITY] * n_sites
    factors[site] = np.asarray(operator, dtype=complex)
    return kron_all(factors)


def two_site_operator(
    operator_i: np.ndarray,
    operator_j: np.ndarray,
    i: int,
    j: int,
    n_sites: int,
) -> np.ndarray:
    """Construct a two-site operator acting on sites ``i`` and ``j``."""

    _validate_site(i, n_sites)
    _validate_site(j, n_sites)
    if i == j:
        raise ValueError("Two-site operators require distinct sites.")

    factors = [IDENTITY] * n_sites
    factors[i] = np.asarray(operator_i, dtype=complex)
    factors[j] = np.asarray(operator_j, dtype=complex)
    return kron_all(factors)


def pauli_string_matrix(operators: Sequence[str]) -> np.ndarray:
    """Construct the dense matrix for a Pauli string such as ``('X', 'I', 'Z')``."""

    try:
        factors = [PAULI_MATRICES[label] for label in operators]
    except KeyError as exc:
        raise ValueError(f"Unknown Pauli operator {exc.args[0]!r}.") from exc
    return kron_all(factors)


def _validate_site(site: int, n_sites: int) -> None:
    if n_sites < 1:
        raise ValueError("n_sites must be positive.")
    if not 0 <= site < n_sites:
        raise IndexError(f"site must satisfy 0 <= site < {n_sites}.")
