"""Single-particle tight-binding lattice model builders."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models._model_utils import (
    add_symmetric_hopping,
    kagome_index_raw,
    nearest_neighbor_bonds,
    onsite_values,
    square_lattice_index_raw,
    validate_positive_int,
    validate_sublattice_index,
    wrapped_cell,
)
from quantum_lattice_models.types import LatticeHamiltonian


def ssh_model(n_cells: int, t1: float = 0.5, t2: float = 1.0, periodic: bool = False) -> np.ndarray:
    """Return the single-particle Su-Schrieffer-Heeger tight-binding matrix."""

    validate_positive_int(n_cells, "n_cells")
    matrix = np.zeros((2 * n_cells, 2 * n_cells), dtype=complex)
    for cell in range(n_cells):
        a = 2 * cell
        b = a + 1
        add_symmetric_hopping(matrix, a, b, -t1)
        if cell < n_cells - 1:
            add_symmetric_hopping(matrix, b, 2 * (cell + 1), -t2)
        elif periodic and n_cells > 1:
            add_symmetric_hopping(matrix, b, 0, -t2)
    return LatticeHamiltonian(
        matrix,
        model_name="ssh_model",
        basis="single_particle",
        lattice_shape=(n_cells, 2),
        metadata={"periodic": periodic, "t1": t1, "t2": t2},
    )


def rice_mele_model(
    n_cells: int,
    hopping: float = 1.0,
    dimerization: float = 0.25,
    staggering: float = 0.5,
    periodic: bool = False,
) -> np.ndarray:
    """Return the single-particle Rice-Mele chain matrix."""

    validate_positive_int(n_cells, "n_cells")
    matrix = np.zeros((2 * n_cells, 2 * n_cells), dtype=complex)
    for cell in range(n_cells):
        a = 2 * cell
        b = a + 1
        matrix[a, a] = staggering
        matrix[b, b] = -staggering
        add_symmetric_hopping(matrix, a, b, -(hopping + dimerization))
        if cell < n_cells - 1:
            add_symmetric_hopping(matrix, b, 2 * (cell + 1), -(hopping - dimerization))
        elif periodic and n_cells > 1:
            add_symmetric_hopping(matrix, b, 0, -(hopping - dimerization))
    return LatticeHamiltonian(
        matrix,
        model_name="rice_mele_model",
        basis="single_particle",
        lattice_shape=(n_cells, 2),
        metadata={
            "periodic": periodic,
            "hopping": hopping,
            "dimerization": dimerization,
            "staggering": staggering,
        },
    )


def tight_binding_chain(
    n_sites: int,
    hopping: float = 1.0,
    onsite: float | Iterable[float] = 0.0,
    periodic: bool = False,
) -> np.ndarray:
    """Return a generic one-dimensional single-particle tight-binding matrix."""

    validate_positive_int(n_sites, "n_sites")
    matrix = np.zeros((n_sites, n_sites), dtype=complex)
    np.fill_diagonal(matrix, onsite_values(onsite, n_sites))
    for i, k in nearest_neighbor_bonds(n_sites, periodic):
        add_symmetric_hopping(matrix, i, k, -hopping)
    return LatticeHamiltonian(
        matrix,
        model_name="tight_binding_chain",
        basis="single_particle",
        lattice_shape=(n_sites,),
        metadata={"periodic": periodic, "hopping": hopping},
    )


def tight_binding_chain_sparse(
    n_sites: int,
    hopping: float = 1.0,
    onsite: float | Iterable[float] = 0.0,
    periodic: bool = False,
) -> sp.csr_matrix:
    """Return a sparse 1D single-particle tight-binding matrix."""

    validate_positive_int(n_sites, "n_sites")
    matrix = sp.diags(onsite_values(onsite, n_sites), format="lil", dtype=complex)
    for i, k in nearest_neighbor_bonds(n_sites, periodic):
        matrix[i, k] += -hopping
        matrix[k, i] += -hopping
    return matrix.tocsr()


def square_lattice_tight_binding(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
    onsite: float | Iterable[float] = 0.0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> np.ndarray:
    """Return a single-particle tight-binding matrix on a rectangular square lattice."""

    validate_positive_int(n_rows, "n_rows")
    validate_positive_int(n_cols, "n_cols")
    n_sites = n_rows * n_cols
    matrix = np.zeros((n_sites, n_sites), dtype=complex)
    np.fill_diagonal(matrix, onsite_values(onsite, n_sites))
    for row in range(n_rows):
        for col in range(n_cols):
            site = square_lattice_index_raw(row, col, n_cols)
            for n_row, n_col in ((row, col + 1), (row + 1, col)):
                wrapped_row, wrapped_col = wrapped_cell(
                    n_row, n_col, n_rows, n_cols, periodic_y, periodic_x
                )
                if wrapped_row is not None and wrapped_col is not None:
                    add_symmetric_hopping(
                        matrix,
                        site,
                        square_lattice_index_raw(wrapped_row, wrapped_col, n_cols),
                        -hopping,
                    )
    return LatticeHamiltonian(
        matrix,
        model_name="square_lattice_tight_binding",
        basis="single_particle",
        lattice_shape=(n_rows, n_cols),
        metadata={"periodic_x": periodic_x, "periodic_y": periodic_y, "hopping": hopping},
    )


def square_lattice_tight_binding_sparse(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
    onsite: float | Iterable[float] = 0.0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> sp.csr_matrix:
    """Return a sparse rectangular square-lattice tight-binding matrix."""

    dense = square_lattice_tight_binding(
        n_rows, n_cols, hopping=hopping, onsite=onsite, periodic_x=periodic_x, periodic_y=periodic_y
    )
    return sp.csr_matrix(np.asarray(dense))


def aubry_andre_harper_chain(
    n_sites: int,
    hopping: float = 1.0,
    potential: float = 1.0,
    beta: float = (np.sqrt(5.0) - 1.0) / 2.0,
    phase: float = 0.0,
    periodic: bool = False,
) -> np.ndarray:
    """Return the single-particle Aubry-Andre-Harper tight-binding matrix."""

    sites = np.arange(n_sites, dtype=float)
    onsite = potential * np.cos(2.0 * np.pi * beta * sites + phase)
    matrix = tight_binding_chain(n_sites, hopping=hopping, onsite=onsite, periodic=periodic)
    matrix.model_name = "aubry_andre_harper_chain"
    matrix.metadata.update({"potential": potential, "beta": beta, "phase": phase})
    return matrix


def triangular_lattice_tight_binding(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
    onsite: float | Iterable[float] = 0.0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> np.ndarray:
    """Return a single-particle tight-binding matrix on a triangular lattice."""

    validate_positive_int(n_rows, "n_rows")
    validate_positive_int(n_cols, "n_cols")
    n_sites = n_rows * n_cols
    matrix = np.zeros((n_sites, n_sites), dtype=complex)
    np.fill_diagonal(matrix, onsite_values(onsite, n_sites))
    for row in range(n_rows):
        for col in range(n_cols):
            site = square_lattice_index_raw(row, col, n_cols)
            for d_row, d_col in ((0, 1), (1, 0), (1, -1)):
                n_row, n_col = wrapped_cell(
                    row + d_row, col + d_col, n_rows, n_cols, periodic_y, periodic_x
                )
                if n_row is not None and n_col is not None:
                    add_symmetric_hopping(
                        matrix, site, square_lattice_index_raw(n_row, n_col, n_cols), -hopping
                    )
    return LatticeHamiltonian(
        matrix,
        model_name="triangular_lattice_tight_binding",
        basis="single_particle",
        lattice_shape=(n_rows, n_cols),
        metadata={"periodic_x": periodic_x, "periodic_y": periodic_y, "hopping": hopping},
    )


def triangular_lattice_tight_binding_sparse(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
    onsite: float | Iterable[float] = 0.0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> sp.csr_matrix:
    """Return a sparse triangular-lattice tight-binding matrix."""

    dense = triangular_lattice_tight_binding(
        n_rows,
        n_cols,
        hopping=hopping,
        onsite=onsite,
        periodic_x=periodic_x,
        periodic_y=periodic_y,
    )
    return sp.csr_matrix(np.asarray(dense))


def kagome_lattice_tight_binding(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
    onsite: float | Iterable[float] = 0.0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> np.ndarray:
    """Return a finite kagome-lattice single-particle tight-binding matrix."""

    validate_positive_int(n_rows, "n_rows")
    validate_positive_int(n_cols, "n_cols")
    n_sites = 3 * n_rows * n_cols
    matrix = np.zeros((n_sites, n_sites), dtype=complex)
    np.fill_diagonal(matrix, onsite_values(onsite, n_sites))
    for row in range(n_rows):
        for col in range(n_cols):
            a = kagome_index_raw(row, col, 0, n_cols)
            b = kagome_index_raw(row, col, 1, n_cols)
            c = kagome_index_raw(row, col, 2, n_cols)
            for i, k in ((a, b), (b, c), (c, a)):
                add_symmetric_hopping(matrix, i, k, -hopping)
            for sub_a, d_row, d_col, sub_b in ((1, 0, 1, 0), (2, 1, 0, 0), (2, 1, -1, 1)):
                n_row, n_col = wrapped_cell(
                    row + d_row, col + d_col, n_rows, n_cols, periodic_y, periodic_x
                )
                if n_row is not None and n_col is not None:
                    add_symmetric_hopping(
                        matrix,
                        kagome_index_raw(row, col, sub_a, n_cols),
                        kagome_index_raw(n_row, n_col, sub_b, n_cols),
                        -hopping,
                    )
    return LatticeHamiltonian(
        matrix,
        model_name="kagome_lattice_tight_binding",
        basis="single_particle",
        lattice_shape=(n_rows, n_cols, 3),
        metadata={"periodic_x": periodic_x, "periodic_y": periodic_y, "hopping": hopping},
    )


def kagome_lattice_tight_binding_sparse(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
    onsite: float | Iterable[float] = 0.0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> sp.csr_matrix:
    """Return a sparse kagome-lattice tight-binding matrix."""

    dense = kagome_lattice_tight_binding(
        n_rows,
        n_cols,
        hopping=hopping,
        onsite=onsite,
        periodic_x=periodic_x,
        periodic_y=periodic_y,
    )
    return sp.csr_matrix(np.asarray(dense))


def ssh_edge_state_localization(
    eigenvector: np.ndarray, n_cells: int, edge_cells: int = 1
) -> float:
    """Return probability weight near the two ends of an SSH chain."""

    validate_positive_int(n_cells, "n_cells")
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


def square_lattice_index(row: int, col: int, n_cols: int) -> int:
    """Return the row-major site index for a square/rectangular lattice."""

    validate_positive_int(n_cols, "n_cols")
    if row < 0 or col < 0 or col >= n_cols:
        raise ValueError("row and col must be nonnegative, and col must be less than n_cols.")
    return square_lattice_index_raw(row, col, n_cols)


def kagome_lattice_index(row: int, col: int, sublattice: int, n_cols: int) -> int:
    """Return the row-major index for a kagome unit-cell sublattice."""

    validate_sublattice_index(row, col, sublattice, n_cols, n_sublattices=3)
    return kagome_index_raw(row, col, sublattice, n_cols)
