"""Spectral utilities for dense Hamiltonians."""

from __future__ import annotations

import numpy as np


def eigenvalues(H: np.ndarray) -> np.ndarray:
    """Return sorted eigenvalues, using Hermitian routines when applicable."""

    matrix = np.asarray(H, dtype=complex)
    if _is_hermitian(matrix):
        return np.linalg.eigvalsh(matrix)
    return np.linalg.eigvals(matrix)


def eigensystem(H: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return eigenvalues and eigenvectors as columns."""

    matrix = np.asarray(H, dtype=complex)
    if _is_hermitian(matrix):
        return np.linalg.eigh(matrix)
    return np.linalg.eig(matrix)


def ground_energy(H: np.ndarray) -> float:
    """Return the lowest real part of the spectrum."""

    values = eigenvalues(H)
    return float(np.min(np.real_if_close(values).real))


def spectral_gap(H: np.ndarray) -> float:
    """Return the gap between the two lowest distinct eigenvalues."""

    values = np.sort(np.real_if_close(eigenvalues(H)).real)
    if values.size < 2:
        return 0.0
    unique = [values[0]]
    for value in values[1:]:
        if not np.isclose(value, unique[-1]):
            unique.append(value)
        if len(unique) == 2:
            return float(unique[1] - unique[0])
    return 0.0


def density_of_states(H: np.ndarray, bins: int = 50) -> tuple[np.ndarray, np.ndarray]:
    """Return histogram counts and bin edges for the Hamiltonian spectrum."""

    if bins < 1:
        raise ValueError("bins must be positive.")
    values = np.real_if_close(eigenvalues(H)).real
    counts, edges = np.histogram(values, bins=bins)
    return counts, edges


def _is_hermitian(matrix: np.ndarray) -> bool:
    return (
        matrix.ndim == 2
        and matrix.shape[0] == matrix.shape[1]
        and np.allclose(matrix, matrix.conj().T)
    )
