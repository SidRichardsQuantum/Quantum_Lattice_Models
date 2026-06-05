"""Model builders for small dense lattice Hamiltonians."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from quantum_lattice_models.operators import (
    PAULI_X,
    PAULI_Y,
    PAULI_Z,
    local_operator,
    two_site_operator,
)
from quantum_lattice_models.types import DenseHamiltonian, PauliTerm


def transverse_field_ising(
    n_sites: int,
    j: float = 1.0,
    h: float = 0.5,
    periodic: bool = False,
) -> DenseHamiltonian:
    """Return the dense transverse-field Ising chain Hamiltonian.

    Convention:
        H = -J sum_i Z_i Z_{i+1} - h sum_i X_i
    """

    _validate_positive_int(n_sites, "n_sites")
    dim = 2**n_sites
    matrix = np.zeros((dim, dim), dtype=complex)
    terms: list[PauliTerm] = []

    for i, k in _nearest_neighbor_bonds(n_sites, periodic):
        matrix += -j * two_site_operator(PAULI_Z, PAULI_Z, i, k, n_sites)
        terms.append(PauliTerm(-j, _pauli_labels(n_sites, {i: "Z", k: "Z"})))

    for i in range(n_sites):
        matrix += -h * local_operator(PAULI_X, i, n_sites)
        terms.append(PauliTerm(-h, _pauli_labels(n_sites, {i: "X"})))

    return DenseHamiltonian(
        matrix, model_name="transverse_field_ising", n_sites=n_sites, terms=tuple(terms)
    )


def heisenberg_chain(
    n_sites: int,
    jx: float = 1.0,
    jy: float = 1.0,
    jz: float = 1.0,
    field: float = 0.0,
    periodic: bool = False,
) -> DenseHamiltonian:
    """Return the dense anisotropic Heisenberg spin-chain Hamiltonian.

    Convention:
        H = sum_i (Jx X_i X_{i+1} + Jy Y_i Y_{i+1} + Jz Z_i Z_{i+1})
            + field sum_i Z_i
    """

    _validate_positive_int(n_sites, "n_sites")
    dim = 2**n_sites
    matrix = np.zeros((dim, dim), dtype=complex)
    terms: list[PauliTerm] = []

    couplings = ((jx, PAULI_X, "X"), (jy, PAULI_Y, "Y"), (jz, PAULI_Z, "Z"))
    for i, k in _nearest_neighbor_bonds(n_sites, periodic):
        for coefficient, operator, label in couplings:
            matrix += coefficient * two_site_operator(operator, operator, i, k, n_sites)
            terms.append(PauliTerm(coefficient, _pauli_labels(n_sites, {i: label, k: label})))

    for i in range(n_sites):
        matrix += field * local_operator(PAULI_Z, i, n_sites)
        terms.append(PauliTerm(field, _pauli_labels(n_sites, {i: "Z"})))

    return DenseHamiltonian(
        matrix, model_name="heisenberg_chain", n_sites=n_sites, terms=tuple(terms)
    )


def ssh_model(
    n_cells: int,
    t1: float = 0.5,
    t2: float = 1.0,
    periodic: bool = False,
) -> np.ndarray:
    """Return the single-particle Su-Schrieffer-Heeger tight-binding matrix."""

    _validate_positive_int(n_cells, "n_cells")
    dim = 2 * n_cells
    matrix = np.zeros((dim, dim), dtype=complex)

    for cell in range(n_cells):
        a = 2 * cell
        b = a + 1
        _add_symmetric_hopping(matrix, a, b, -t1)

        if cell < n_cells - 1:
            next_a = 2 * (cell + 1)
            _add_symmetric_hopping(matrix, b, next_a, -t2)
        elif periodic and n_cells > 1:
            _add_symmetric_hopping(matrix, b, 0, -t2)

    return matrix


def tight_binding_chain(
    n_sites: int,
    hopping: float = 1.0,
    onsite: float | Iterable[float] = 0.0,
    periodic: bool = False,
) -> np.ndarray:
    """Return a generic one-dimensional single-particle tight-binding matrix."""

    _validate_positive_int(n_sites, "n_sites")
    matrix = np.zeros((n_sites, n_sites), dtype=complex)
    onsite_values = _onsite_values(onsite, n_sites)
    np.fill_diagonal(matrix, onsite_values)

    for i, k in _nearest_neighbor_bonds(n_sites, periodic):
        _add_symmetric_hopping(matrix, i, k, -hopping)

    return matrix


def ssh_edge_state_localization(
    eigenvector: np.ndarray, n_cells: int, edge_cells: int = 1
) -> float:
    """Return probability weight near the two ends of an SSH chain."""

    _validate_positive_int(n_cells, "n_cells")
    if edge_cells < 1 or edge_cells > n_cells:
        raise ValueError("edge_cells must satisfy 1 <= edge_cells <= n_cells.")

    vector = np.asarray(eigenvector, dtype=complex).reshape(-1)
    if vector.size != 2 * n_cells:
        raise ValueError("eigenvector length must equal 2 * n_cells.")

    probabilities = np.abs(vector) ** 2
    norm = probabilities.sum()
    if norm == 0:
        raise ValueError("eigenvector must have nonzero norm.")

    left = slice(0, 2 * edge_cells)
    right = slice(2 * (n_cells - edge_cells), 2 * n_cells)
    return float((probabilities[left].sum() + probabilities[right].sum()) / norm)


def ssh_edge_state_localizations(
    eigenvectors: np.ndarray, n_cells: int, edge_cells: int = 1
) -> np.ndarray:
    """Return edge-localization weights for eigenvectors stored as columns."""

    vectors = np.asarray(eigenvectors, dtype=complex)
    if vectors.ndim != 2 or vectors.shape[0] != 2 * n_cells:
        raise ValueError("eigenvectors must have shape (2 * n_cells, n_vectors).")
    return np.array(
        [
            ssh_edge_state_localization(vectors[:, col], n_cells, edge_cells)
            for col in range(vectors.shape[1])
        ]
    )


def _nearest_neighbor_bonds(n_sites: int, periodic: bool) -> list[tuple[int, int]]:
    bonds = [(i, i + 1) for i in range(n_sites - 1)]
    if periodic and n_sites > 2:
        bonds.append((n_sites - 1, 0))
    elif periodic and n_sites == 2:
        bonds.append((1, 0))
    return bonds


def _add_symmetric_hopping(matrix: np.ndarray, i: int, j: int, value: float) -> None:
    matrix[i, j] += value
    matrix[j, i] += np.conjugate(value)


def _validate_positive_int(value: int, name: str) -> None:
    if not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer.")


def _pauli_labels(n_sites: int, labels: dict[int, str]) -> tuple[str, ...]:
    return tuple(labels.get(i, "I") for i in range(n_sites))


def _onsite_values(onsite: float | Iterable[float], n_sites: int) -> np.ndarray:
    if np.isscalar(onsite):
        return np.full(n_sites, onsite, dtype=float)
    values = np.asarray(list(onsite), dtype=float)
    if values.shape != (n_sites,):
        raise ValueError("onsite must be a scalar or a length-n_sites iterable.")
    return values
