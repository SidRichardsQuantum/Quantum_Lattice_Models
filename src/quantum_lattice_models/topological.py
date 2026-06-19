"""Topological lattice model builders."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models._model_utils import (
    add_symmetric_hopping,
    honeycomb_index_raw,
    nearest_neighbor_bonds,
    onsite_values,
    square_lattice_index_raw,
    validate_positive_int,
    validate_sublattice_index,
    wrapped_cell,
)
from quantum_lattice_models.types import LatticeHamiltonian


def harper_hofstadter_square_lattice(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
    flux: float = 0.25,
    onsite: float | Iterable[float] = 0.0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> np.ndarray:
    """Return a square-lattice Harper-Hofstadter single-particle Hamiltonian."""

    matrix = harper_hofstadter_square_lattice_sparse(
        n_rows,
        n_cols,
        hopping=hopping,
        flux=flux,
        onsite=onsite,
        periodic_x=periodic_x,
        periodic_y=periodic_y,
    ).toarray()
    return LatticeHamiltonian(
        matrix,
        model_name="harper_hofstadter_square_lattice",
        basis="single_particle",
        lattice_shape=(n_rows, n_cols),
        metadata={
            "periodic_x": periodic_x,
            "periodic_y": periodic_y,
            "hopping": hopping,
            "flux": flux,
        },
    )


def harper_hofstadter_square_lattice_sparse(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
    flux: float = 0.25,
    onsite: float | Iterable[float] = 0.0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> sp.csr_matrix:
    """Return a sparse Harper-Hofstadter square-lattice matrix."""

    validate_positive_int(n_rows, "n_rows")
    validate_positive_int(n_cols, "n_cols")
    n_sites = n_rows * n_cols
    matrix = sp.diags(onsite_values(onsite, n_sites), format="lil", dtype=complex)
    for row in range(n_rows):
        for col in range(n_cols):
            site = square_lattice_index_raw(row, col, n_cols)
            n_row, n_col = wrapped_cell(row, col + 1, n_rows, n_cols, periodic_y, periodic_x)
            if n_row is not None and n_col is not None:
                add_symmetric_hopping(
                    matrix,
                    site,
                    square_lattice_index_raw(n_row, n_col, n_cols),
                    -hopping,
                )
            phase = np.exp(2j * np.pi * flux * col)
            n_row, n_col = wrapped_cell(row + 1, col, n_rows, n_cols, periodic_y, periodic_x)
            if n_row is not None and n_col is not None:
                add_symmetric_hopping(
                    matrix,
                    site,
                    square_lattice_index_raw(n_row, n_col, n_cols),
                    -hopping * phase,
                )
    return matrix.tocsr()


def kitaev_chain_bdg(
    n_sites: int,
    hopping: float = 1.0,
    chemical_potential: float = 0.0,
    pairing: complex = 1.0,
    periodic: bool = False,
) -> np.ndarray:
    """Return a single-particle Bogoliubov-de Gennes matrix for the Kitaev chain."""

    validate_positive_int(n_sites, "n_sites")
    normal = np.zeros((n_sites, n_sites), dtype=complex)
    pairing_matrix = np.zeros((n_sites, n_sites), dtype=complex)
    np.fill_diagonal(normal, -chemical_potential)
    for i, k in nearest_neighbor_bonds(n_sites, periodic):
        add_symmetric_hopping(normal, i, k, -hopping)
        pairing_matrix[i, k] += pairing
        pairing_matrix[k, i] += -pairing
    matrix = np.block([[normal, pairing_matrix], [-pairing_matrix.conj(), -normal.conj()]])
    return LatticeHamiltonian(
        matrix,
        model_name="kitaev_chain_bdg",
        basis="nambu_single_particle",
        lattice_shape=(n_sites,),
        metadata={
            "periodic": periodic,
            "hopping": hopping,
            "chemical_potential": chemical_potential,
            "pairing": pairing,
        },
    )


def haldane_honeycomb_lattice(
    n_rows: int,
    n_cols: int,
    t1: float = 1.0,
    t2: float = 0.1,
    phi: float = np.pi / 2,
    sublattice_potential: float = 0.0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> np.ndarray:
    """Return a finite honeycomb-lattice Haldane-model matrix."""

    matrix = haldane_honeycomb_lattice_sparse(
        n_rows,
        n_cols,
        t1=t1,
        t2=t2,
        phi=phi,
        sublattice_potential=sublattice_potential,
        periodic_x=periodic_x,
        periodic_y=periodic_y,
    ).toarray()
    return LatticeHamiltonian(
        matrix,
        model_name="haldane_honeycomb_lattice",
        basis="single_particle",
        lattice_shape=(n_rows, n_cols, 2),
        metadata={
            "periodic_x": periodic_x,
            "periodic_y": periodic_y,
            "t1": t1,
            "t2": t2,
            "phi": phi,
            "sublattice_potential": sublattice_potential,
        },
    )


def haldane_honeycomb_lattice_sparse(
    n_rows: int,
    n_cols: int,
    t1: float = 1.0,
    t2: float = 0.1,
    phi: float = np.pi / 2,
    sublattice_potential: float = 0.0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> sp.csr_matrix:
    """Return a sparse finite Haldane honeycomb-lattice matrix."""

    validate_positive_int(n_rows, "n_rows")
    validate_positive_int(n_cols, "n_cols")
    n_sites = 2 * n_rows * n_cols
    matrix = sp.lil_matrix((n_sites, n_sites), dtype=complex)
    for row in range(n_rows):
        for col in range(n_cols):
            a = honeycomb_index_raw(row, col, 0, n_cols)
            b = honeycomb_index_raw(row, col, 1, n_cols)
            matrix[a, a] = sublattice_potential
            matrix[b, b] = -sublattice_potential
            add_symmetric_hopping(matrix, a, b, -t1)
            for n_row, n_col in (
                wrapped_cell(row, col - 1, n_rows, n_cols, periodic_y, periodic_x),
                wrapped_cell(row - 1, col, n_rows, n_cols, periodic_y, periodic_x),
            ):
                if n_row is not None and n_col is not None:
                    add_symmetric_hopping(
                        matrix,
                        a,
                        honeycomb_index_raw(n_row, n_col, 1, n_cols),
                        -t1,
                    )
            for d_row, d_col in ((0, 1), (1, 0), (1, -1)):
                n_row, n_col = wrapped_cell(
                    row + d_row, col + d_col, n_rows, n_cols, periodic_y, periodic_x
                )
                if n_row is None or n_col is None:
                    continue
                add_symmetric_hopping(
                    matrix,
                    a,
                    honeycomb_index_raw(n_row, n_col, 0, n_cols),
                    -t2 * np.exp(1j * phi),
                )
                add_symmetric_hopping(
                    matrix,
                    b,
                    honeycomb_index_raw(n_row, n_col, 1, n_cols),
                    -t2 * np.exp(-1j * phi),
                )
    return matrix.tocsr()


def honeycomb_lattice_index(row: int, col: int, sublattice: int, n_cols: int) -> int:
    """Return the row-major index for a honeycomb unit-cell sublattice."""

    validate_sublattice_index(row, col, sublattice, n_cols, n_sublattices=2)
    return honeycomb_index_raw(row, col, sublattice, n_cols)
