"""Fixed-particle-number sectors for truncated Bose-Hubbard chains."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import product

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models._model_utils import nearest_neighbor_bonds


@dataclass(frozen=True)
class FixedParticleNumberBasis:
    n_sites: int
    n_particles: int
    max_occupancy: int
    states: tuple[tuple[int, ...], ...]

    @property
    def dimension(self) -> int:
        return len(self.states)

    @property
    def state_to_index(self) -> dict[tuple[int, ...], int]:
        return {state: index for index, state in enumerate(self.states)}

    def embed(self, state: np.ndarray) -> np.ndarray:
        local_dim = self.max_occupancy + 1
        vector = np.asarray(state, dtype=complex).reshape(-1)
        if vector.size != self.dimension:
            raise ValueError("Reduced state length must match sector dimension.")
        full = np.zeros(local_dim**self.n_sites, dtype=complex)
        for amplitude, occupations in zip(vector, self.states, strict=True):
            index = 0
            for occupation in occupations:
                index = index * local_dim + occupation
            full[index] = amplitude
        return full


@dataclass(frozen=True)
class BoseHubbardSector:
    matrix: sp.csr_matrix
    basis: FixedParticleNumberBasis
    parameters: dict[str, object]


def fixed_particle_number_basis(
    n_sites: int,
    n_particles: int,
    max_occupancy: int = 2,
) -> FixedParticleNumberBasis:
    if n_sites < 1 or max_occupancy < 1 or n_particles < 0:
        raise ValueError("Require positive sites/occupancy and nonnegative particle number.")
    states = tuple(
        state
        for state in product(range(max_occupancy + 1), repeat=n_sites)
        if sum(state) == n_particles
    )
    if not states:
        raise ValueError("The requested particle-number sector is empty.")
    return FixedParticleNumberBasis(n_sites, n_particles, max_occupancy, states)


def bose_hubbard_chain_sector(
    n_sites: int,
    n_particles: int,
    hopping: float = 1.0,
    interaction: float = 1.0,
    chemical_potential: float = 0.0,
    max_occupancy: int = 2,
    periodic: bool = False,
) -> BoseHubbardSector:
    basis = fixed_particle_number_basis(n_sites, n_particles, max_occupancy)
    mapping = basis.state_to_index
    matrix = sp.lil_matrix((basis.dimension, basis.dimension), dtype=complex)
    for column, occupations_tuple in enumerate(basis.states):
        occupations = list(occupations_tuple)
        matrix[column, column] = sum(
            0.5 * interaction * n * (n - 1) - chemical_potential * n for n in occupations
        )
        for left, right in nearest_neighbor_bonds(n_sites, periodic):
            for source, target in ((left, right), (right, left)):
                if occupations[source] == 0 or occupations[target] >= max_occupancy:
                    continue
                updated = occupations.copy()
                amplitude = -hopping * np.sqrt(occupations[source] * (occupations[target] + 1))
                updated[source] -= 1
                updated[target] += 1
                matrix[mapping[tuple(updated)], column] += amplitude
    return BoseHubbardSector(
        matrix.tocsr(),
        basis,
        {
            "n_sites": n_sites,
            "n_particles": n_particles,
            "hopping": hopping,
            "interaction": interaction,
            "chemical_potential": chemical_potential,
            "max_occupancy": max_occupancy,
            "periodic": periodic,
        },
    )


def boson_site_occupations(
    state: np.ndarray,
    basis: FixedParticleNumberBasis,
) -> np.ndarray:
    vector = np.asarray(state, dtype=complex).reshape(-1)
    if vector.size != basis.dimension:
        raise ValueError("State length must match sector dimension.")
    probabilities = np.abs(vector) ** 2
    probabilities /= probabilities.sum()
    return probabilities @ np.asarray(basis.states)
