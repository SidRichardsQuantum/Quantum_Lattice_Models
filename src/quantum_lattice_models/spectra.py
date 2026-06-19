"""Spectral utilities for dense Hamiltonians."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from quantum_lattice_models._model_utils import as_dense


def eigenvalues(H: np.ndarray) -> np.ndarray:
    """Return sorted eigenvalues, using Hermitian routines when applicable."""

    matrix = as_dense(H)
    if _is_hermitian(matrix):
        return np.linalg.eigvalsh(matrix)
    return np.linalg.eigvals(matrix)


def eigensystem(H: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return eigenvalues and eigenvectors as columns."""

    matrix = as_dense(H)
    if _is_hermitian(matrix):
        return np.linalg.eigh(matrix)
    return np.linalg.eig(matrix)


def ground_energy(H: np.ndarray) -> float:
    """Return the lowest real part of the spectrum."""

    values = eigenvalues(H)
    return float(np.min(np.real_if_close(values).real))


def lowest_eigenvalues(H: np.ndarray, k: int = 2) -> np.ndarray:
    """Return the lowest ``k`` eigenvalues, using sparse solvers for sparse inputs."""

    if k < 1:
        raise ValueError("k must be positive.")
    if sp.issparse(H):
        matrix = H.asfptype()
        n = matrix.shape[0]
        if k >= n:
            return np.sort(np.linalg.eigvalsh(matrix.toarray()))
        values = spla.eigsh(matrix, k=k, which="SA", return_eigenvectors=False)
        return np.sort(np.real_if_close(values))
    return np.sort(np.real_if_close(eigenvalues(H)).real)[:k]


def ground_state(H: np.ndarray) -> tuple[float, np.ndarray]:
    """Return the ground-state energy and eigenvector."""

    if sp.issparse(H):
        matrix = H.asfptype()
        if matrix.shape[0] <= 2:
            values, vectors = np.linalg.eigh(matrix.toarray())
        else:
            values, vectors = spla.eigsh(matrix, k=1, which="SA")
    else:
        values, vectors = eigensystem(H)
    index = int(np.argmin(np.real_if_close(values).real))
    return float(np.real_if_close(values[index]).real), vectors[:, index]


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
