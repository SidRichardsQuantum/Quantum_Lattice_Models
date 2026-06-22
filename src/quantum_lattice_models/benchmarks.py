"""High-value benchmark lattice and spin model families."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models._model_utils import (
    add_symmetric_hopping,
    nearest_neighbor_bonds,
    validate_positive_int,
)
from quantum_lattice_models.geometry import honeycomb_lattice_positions
from quantum_lattice_models.spin import SpinField, SpinInteraction, graph_spin_hamiltonian
from quantum_lattice_models.topological import haldane_honeycomb_lattice_sparse
from quantum_lattice_models.types import LatticeHamiltonian


def anderson_chain(
    n_sites: int,
    hopping: float = 1.0,
    disorder: float = 1.0,
    seed: int = 0,
    periodic: bool = False,
) -> LatticeHamiltonian:
    """Return a reproducible one-dimensional Anderson model."""

    rng = np.random.default_rng(seed)
    onsite = rng.uniform(-disorder / 2, disorder / 2, n_sites)
    matrix = np.diag(onsite.astype(complex))
    for left, right in nearest_neighbor_bonds(n_sites, periodic):
        add_symmetric_hopping(matrix, left, right, -hopping)
    return LatticeHamiltonian(
        matrix,
        model_name="anderson_chain",
        basis="single_particle",
        lattice_shape=(n_sites,),
        metadata={"disorder": disorder, "seed": seed, "periodic": periodic},
    )


def graphene_lattice(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> LatticeHamiltonian:
    """Return a finite nearest-neighbor graphene/honeycomb Hamiltonian."""

    matrix = graphene_lattice_sparse(n_rows, n_cols, hopping, periodic_x, periodic_y).toarray()
    return LatticeHamiltonian(
        matrix,
        model_name="graphene_lattice",
        basis="single_particle",
        lattice_shape=(n_rows, n_cols, 2),
        metadata={
            "periodic_x": periodic_x,
            "periodic_y": periodic_y,
            "hopping": hopping,
            "positions": honeycomb_lattice_positions(n_rows, n_cols),
            "sublattice_labels": tuple(
                label for _ in range(n_rows * n_cols) for label in ("A", "B")
            ),
        },
    )


def graphene_lattice_sparse(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> sp.csr_matrix:
    """Return a sparse finite graphene/honeycomb Hamiltonian."""

    return haldane_honeycomb_lattice_sparse(
        n_rows,
        n_cols,
        t1=hopping,
        t2=0.0,
        phi=0.0,
        sublattice_potential=0.0,
        periodic_x=periodic_x,
        periodic_y=periodic_y,
    )


def anderson_square_lattice(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
    disorder: float = 1.0,
    seed: int = 0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> LatticeHamiltonian:
    """Return a reproducible two-dimensional Anderson model."""

    matrix = anderson_square_lattice_sparse(
        n_rows, n_cols, hopping, disorder, seed, periodic_x, periodic_y
    ).toarray()
    return LatticeHamiltonian(
        matrix,
        model_name="anderson_square_lattice",
        basis="single_particle",
        lattice_shape=(n_rows, n_cols),
        metadata={
            "disorder": disorder,
            "seed": seed,
            "periodic_x": periodic_x,
            "periodic_y": periodic_y,
            "positions": _square_positions(n_rows, n_cols),
        },
    )


def anderson_square_lattice_sparse(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
    disorder: float = 1.0,
    seed: int = 0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> sp.csr_matrix:
    """Return a sparse two-dimensional Anderson Hamiltonian."""

    validate_positive_int(n_rows, "n_rows")
    validate_positive_int(n_cols, "n_cols")
    rng = np.random.default_rng(seed)
    matrix = sp.diags(
        rng.uniform(-disorder / 2, disorder / 2, n_rows * n_cols),
        format="lil",
        dtype=complex,
    )

    def index(row: int, col: int) -> int:
        return row * n_cols + col

    for row in range(n_rows):
        for col in range(n_cols):
            for d_row, d_col, periodic in ((0, 1, periodic_x), (1, 0, periodic_y)):
                next_row, next_col = row + d_row, col + d_col
                if next_col >= n_cols:
                    if not periodic:
                        continue
                    next_col = 0
                if next_row >= n_rows:
                    if not periodic:
                        continue
                    next_row = 0
                if index(row, col) != index(next_row, next_col):
                    add_symmetric_hopping(
                        matrix, index(row, col), index(next_row, next_col), -hopping
                    )
    return matrix.tocsr()


def checkerboard_chern_insulator(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
    mass: float = 1.0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> LatticeHamiltonian:
    """Return a two-orbital square-lattice Chern-insulator benchmark."""

    matrix = checkerboard_chern_insulator_sparse(
        n_rows, n_cols, hopping, mass, periodic_x, periodic_y
    ).toarray()
    return LatticeHamiltonian(
        matrix,
        model_name="checkerboard_chern_insulator",
        basis="single_particle",
        lattice_shape=(n_rows, n_cols, 2),
        metadata={
            "mass": mass,
            "hopping": hopping,
            "convention": "Qi-Wu-Zhang-type two-orbital checkerboard representation",
            "positions": _two_orbital_square_positions(n_rows, n_cols),
            "orbital_labels": tuple(label for _ in range(n_rows * n_cols) for label in ("A", "B")),
        },
    )


def checkerboard_chern_insulator_sparse(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
    mass: float = 1.0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> sp.csr_matrix:
    """Return the sparse real-space two-band Chern-insulator matrix."""

    validate_positive_int(n_rows, "n_rows")
    validate_positive_int(n_cols, "n_cols")
    sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma_z = np.diag([1.0, -1.0]).astype(complex)
    matrix = sp.lil_matrix((2 * n_rows * n_cols, 2 * n_rows * n_cols), dtype=complex)

    def cell(row: int, col: int) -> int:
        return 2 * (row * n_cols + col)

    for row in range(n_rows):
        for col in range(n_cols):
            start = cell(row, col)
            matrix[start : start + 2, start : start + 2] += mass * sigma_z
            for d_row, d_col, periodic, pauli in (
                (0, 1, periodic_x, sigma_x),
                (1, 0, periodic_y, sigma_y),
            ):
                next_row, next_col = row + d_row, col + d_col
                if next_col >= n_cols:
                    if not periodic:
                        continue
                    next_col = 0
                if next_row >= n_rows:
                    if not periodic:
                        continue
                    next_row = 0
                target = cell(next_row, next_col)
                hopping_block = 0.5 * hopping * (sigma_z - 1j * pauli)
                matrix[start : start + 2, target : target + 2] += hopping_block
                matrix[target : target + 2, start : start + 2] += hopping_block.conj().T
    return matrix.tocsr()


def dice_lattice(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> LatticeHamiltonian:
    """Return a finite dice/T3 lattice with one hub and two rim sites per cell."""

    matrix = dice_lattice_sparse(n_rows, n_cols, hopping, periodic_x, periodic_y).toarray()
    return LatticeHamiltonian(
        matrix,
        model_name="dice_lattice",
        basis="single_particle",
        lattice_shape=(n_rows, n_cols, 3),
        metadata={
            "hopping": hopping,
            "periodic_x": periodic_x,
            "periodic_y": periodic_y,
            "positions": _dice_positions(n_rows, n_cols),
            "sublattice_labels": tuple(
                label for _ in range(n_rows * n_cols) for label in ("hub", "rim-A", "rim-B")
            ),
        },
    )


def dice_lattice_sparse(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
    periodic_x: bool = False,
    periodic_y: bool = False,
) -> sp.csr_matrix:
    """Return a sparse finite dice/T3 lattice Hamiltonian."""

    validate_positive_int(n_rows, "n_rows")
    validate_positive_int(n_cols, "n_cols")
    matrix = sp.lil_matrix((3 * n_rows * n_cols, 3 * n_rows * n_cols), dtype=complex)

    def index(row: int, col: int, orbital: int) -> int:
        return 3 * (row * n_cols + col) + orbital

    def wrapped(row: int, col: int) -> tuple[int, int] | None:
        if col < 0 or col >= n_cols:
            if not periodic_x:
                return None
            col %= n_cols
        if row < 0 or row >= n_rows:
            if not periodic_y:
                return None
            row %= n_rows
        return row, col

    for row in range(n_rows):
        for col in range(n_cols):
            hub = index(row, col, 0)
            for orbital, offsets in (
                (1, ((0, 0), (0, -1), (-1, 0))),
                (2, ((0, 0), (0, 1), (1, 0))),
            ):
                for d_row, d_col in offsets:
                    target_cell = wrapped(row + d_row, col + d_col)
                    if target_cell is not None:
                        add_symmetric_hopping(
                            matrix,
                            hub,
                            index(*target_cell, orbital),
                            -hopping,
                        )
    return matrix.tocsr()


def _square_positions(n_rows: int, n_cols: int) -> np.ndarray:
    return np.asarray([(float(col), float(-row)) for row in range(n_rows) for col in range(n_cols)])


def _two_orbital_square_positions(n_rows: int, n_cols: int) -> np.ndarray:
    return np.asarray(
        [
            (float(col) + offset, float(-row))
            for row in range(n_rows)
            for col in range(n_cols)
            for offset in (-0.12, 0.12)
        ]
    )


def _dice_positions(n_rows: int, n_cols: int) -> np.ndarray:
    return np.asarray(
        [
            (float(col) + x_offset, float(-row) + y_offset)
            for row in range(n_rows)
            for col in range(n_cols)
            for x_offset, y_offset in ((0.0, 0.0), (-0.25, 0.2), (0.25, -0.2))
        ]
    )


def long_range_tight_binding_chain(
    n_sites: int,
    hopping: float = 1.0,
    power: float = 3.0,
    onsite: float | Iterable[float] = 0.0,
    periodic: bool = False,
) -> LatticeHamiltonian:
    """Return a chain with hopping decaying as inverse distance to ``power``."""

    diagonal = np.full(n_sites, onsite) if np.isscalar(onsite) else np.asarray(tuple(onsite))
    matrix = np.diag(diagonal.astype(complex))
    for left in range(n_sites):
        for right in range(left + 1, n_sites):
            distance = right - left
            if periodic:
                distance = min(distance, n_sites - distance)
            add_symmetric_hopping(matrix, left, right, -hopping / distance**power)
    return LatticeHamiltonian(
        matrix,
        model_name="long_range_tight_binding_chain",
        basis="single_particle",
        lattice_shape=(n_sites,),
        metadata={"power": power, "periodic": periodic},
    )


def creutz_ladder(
    n_cells: int,
    hopping: float = 1.0,
    diagonal_hopping: float = 1.0,
    flux: float = np.pi,
    mass: float = 0.0,
    periodic: bool = False,
) -> LatticeHamiltonian:
    """Return a two-leg Creutz ladder with complex leg hoppings."""

    matrix = np.zeros((2 * n_cells, 2 * n_cells), dtype=complex)
    for cell in range(n_cells):
        upper, lower = 2 * cell, 2 * cell + 1
        matrix[upper, upper] = mass
        matrix[lower, lower] = -mass
        next_cell = (cell + 1) % n_cells
        if cell == n_cells - 1 and not periodic:
            continue
        next_upper, next_lower = 2 * next_cell, 2 * next_cell + 1
        add_symmetric_hopping(matrix, upper, next_upper, -hopping * np.exp(1j * flux / 2))
        add_symmetric_hopping(matrix, lower, next_lower, -hopping * np.exp(-1j * flux / 2))
        add_symmetric_hopping(matrix, upper, next_lower, -diagonal_hopping)
        add_symmetric_hopping(matrix, lower, next_upper, -diagonal_hopping)
    return LatticeHamiltonian(
        matrix,
        model_name="creutz_ladder",
        basis="single_particle",
        lattice_shape=(n_cells, 2),
        metadata={"flux": flux, "mass": mass, "periodic": periodic},
    )


def sawtooth_chain(
    n_cells: int,
    base_hopping: float = 1.0,
    tooth_hopping: float = np.sqrt(2.0),
    periodic: bool = False,
) -> LatticeHamiltonian:
    """Return a sawtooth/delta chain with two orbitals per cell."""

    matrix = np.zeros((2 * n_cells, 2 * n_cells), dtype=complex)
    for cell in range(n_cells):
        base, tooth = 2 * cell, 2 * cell + 1
        add_symmetric_hopping(matrix, base, tooth, -tooth_hopping)
        next_cell = (cell + 1) % n_cells
        if cell == n_cells - 1 and not periodic:
            continue
        next_base = 2 * next_cell
        add_symmetric_hopping(matrix, base, next_base, -base_hopping)
        add_symmetric_hopping(matrix, tooth, next_base, -tooth_hopping)
    return LatticeHamiltonian(
        matrix,
        model_name="sawtooth_chain",
        basis="single_particle",
        lattice_shape=(n_cells, 2),
        metadata={"periodic": periodic},
    )


def lieb_lattice(
    n_rows: int,
    n_cols: int,
    hopping: float = 1.0,
) -> LatticeHamiltonian:
    """Return an open finite Lieb lattice with three orbitals per unit cell."""

    n_sites = 3 * n_rows * n_cols
    matrix = np.zeros((n_sites, n_sites), dtype=complex)

    def index(row: int, col: int, orbital: int) -> int:
        return 3 * (row * n_cols + col) + orbital

    for row in range(n_rows):
        for col in range(n_cols):
            center = index(row, col, 0)
            horizontal = index(row, col, 1)
            vertical = index(row, col, 2)
            add_symmetric_hopping(matrix, center, horizontal, -hopping)
            add_symmetric_hopping(matrix, center, vertical, -hopping)
            if col + 1 < n_cols:
                add_symmetric_hopping(matrix, horizontal, index(row, col + 1, 0), -hopping)
            if row + 1 < n_rows:
                add_symmetric_hopping(matrix, vertical, index(row + 1, col, 0), -hopping)
    return LatticeHamiltonian(
        matrix,
        model_name="lieb_lattice",
        basis="single_particle",
        lattice_shape=(n_rows, n_cols, 3),
        metadata={"hopping": hopping},
    )


def xyz_chain(
    n_sites: int,
    jx: float = 1.0,
    jy: float = 0.8,
    jz: float = 0.6,
    field: float = 0.0,
    periodic: bool = False,
):
    """Return an XYZ spin chain."""

    interactions = tuple(
        SpinInteraction(left, right, axis, axis, coupling)
        for left, right in nearest_neighbor_bonds(n_sites, periodic)
        for axis, coupling in (("X", jx), ("Y", jy), ("Z", jz))
    )
    fields = tuple(SpinField(site, "Z", field) for site in range(n_sites))
    return graph_spin_hamiltonian(
        n_sites,
        interactions,
        fields,
        model_name="xyz_chain",
    )


def random_field_heisenberg_chain(
    n_sites: int,
    coupling: float = 1.0,
    disorder: float = 1.0,
    seed: int = 0,
    periodic: bool = False,
):
    """Return a Heisenberg chain with reproducible random longitudinal fields."""

    rng = np.random.default_rng(seed)
    interactions = tuple(
        SpinInteraction(left, right, axis, axis, coupling)
        for left, right in nearest_neighbor_bonds(n_sites, periodic)
        for axis in ("X", "Y", "Z")
    )
    fields = tuple(
        SpinField(site, "Z", value)
        for site, value in enumerate(rng.uniform(-disorder, disorder, n_sites))
    )
    return graph_spin_hamiltonian(
        n_sites,
        interactions,
        fields,
        model_name="random_field_heisenberg_chain",
    )
