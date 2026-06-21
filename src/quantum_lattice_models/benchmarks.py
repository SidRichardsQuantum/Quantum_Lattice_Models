"""High-value benchmark lattice and spin model families."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np

from quantum_lattice_models._model_utils import add_symmetric_hopping, nearest_neighbor_bonds
from quantum_lattice_models.spin import SpinField, SpinInteraction, graph_spin_hamiltonian
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
