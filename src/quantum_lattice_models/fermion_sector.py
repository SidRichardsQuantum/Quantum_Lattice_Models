"""Fixed-particle-number sectors for spinful Fermi-Hubbard chains."""

from __future__ import annotations

from dataclasses import dataclass
from math import comb

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models._model_utils import (
    bit_is_set,
    fermion_hop_sign,
    nearest_neighbor_bonds,
    validate_positive_int,
)
from quantum_lattice_models.reduced import ReducedBasisMapping


@dataclass(frozen=True)
class FermiHubbardBasis:
    """Occupation basis with fixed up- and down-spin particle numbers."""

    n_sites: int
    n_up: int
    n_down: int
    states: tuple[int, ...]

    @property
    def mapping(self) -> ReducedBasisMapping:
        return ReducedBasisMapping(
            kind="fixed_spin_particle_numbers",
            full_dimension=2 ** (2 * self.n_sites),
            states=self.states,
            quantum_numbers={"n_up": self.n_up, "n_down": self.n_down},
            metadata={
                "n_sites": self.n_sites,
                "orbital_order": "site0_up, site0_down, site1_up, site1_down, ...",
            },
        )

    @property
    def dimension(self) -> int:
        return len(self.states)

    @property
    def state_to_index(self) -> dict[int, int]:
        return {state: index for index, state in enumerate(self.states)}

    def embed(self, state: np.ndarray) -> np.ndarray:
        return self.mapping.embed(state)

    def project(self, state: np.ndarray) -> np.ndarray:
        return self.mapping.project(state)

    def to_metadata(self) -> dict[str, object]:
        return {
            "kind": "fixed_spin_particle_numbers",
            "n_sites": self.n_sites,
            "n_up": self.n_up,
            "n_down": self.n_down,
            "dimension": self.dimension,
            "basis_states": list(self.states),
            "orbital_order": "site0_up, site0_down, site1_up, site1_down, ...",
            "mapping": self.mapping.to_dict(),
        }


@dataclass(frozen=True)
class FermiHubbardSector:
    """Sparse Fermi-Hubbard Hamiltonian and its explicit reduced basis."""

    matrix: sp.csr_matrix
    basis: FermiHubbardBasis
    parameters: dict[str, object]
    model_name: str = "fermi_hubbard_chain_sector"

    @property
    def shape(self) -> tuple[int, int]:
        return int(self.matrix.shape[0]), int(self.matrix.shape[1])

    def to_metadata(self) -> dict[str, object]:
        return {
            "model_name": self.model_name,
            "sector": self.basis.to_metadata(),
            "parameters": dict(self.parameters),
        }


def fermi_hubbard_basis(n_sites: int, n_up: int, n_down: int) -> FermiHubbardBasis:
    """Return spinful occupation states with fixed ``N_up`` and ``N_down``."""

    validate_positive_int(n_sites, "n_sites")
    for value, name in ((n_up, "n_up"), (n_down, "n_down")):
        if not isinstance(value, int) or isinstance(value, bool):
            raise ValueError(f"{name} must be an integer.")
        if not 0 <= value <= n_sites:
            raise ValueError(f"{name} must satisfy 0 <= {name} <= n_sites.")
    states = tuple(
        state
        for state in range(2 ** (2 * n_sites))
        if sum(bit_is_set(state, 2 * site) for site in range(n_sites)) == n_up
        and sum(bit_is_set(state, 2 * site + 1) for site in range(n_sites)) == n_down
    )
    expected = comb(n_sites, n_up) * comb(n_sites, n_down)
    if len(states) != expected:
        raise RuntimeError("Fermi-Hubbard sector basis construction failed.")
    return FermiHubbardBasis(n_sites, n_up, n_down, states)


def fermi_hubbard_chain_sector(
    n_sites: int,
    n_up: int,
    n_down: int,
    hopping: float = 1.0,
    interaction: float = 1.0,
    chemical_potential: float = 0.0,
    periodic: bool = False,
) -> FermiHubbardSector:
    """Build a spinful Fermi-Hubbard chain in a fixed ``(N_up, N_down)`` sector."""

    basis = fermi_hubbard_basis(n_sites, n_up, n_down)
    mapping = basis.state_to_index
    matrix = sp.lil_matrix((basis.dimension, basis.dimension), dtype=complex)
    for column, state in enumerate(basis.states):
        diagonal = 0.0
        for site in range(n_sites):
            up = bit_is_set(state, 2 * site)
            down = bit_is_set(state, 2 * site + 1)
            diagonal += interaction * up * down - chemical_potential * (up + down)
        matrix[column, column] = diagonal
        for left, right in nearest_neighbor_bonds(n_sites, periodic):
            for spin in (0, 1):
                for source, target in (
                    (2 * left + spin, 2 * right + spin),
                    (2 * right + spin, 2 * left + spin),
                ):
                    if not bit_is_set(state, source) or bit_is_set(state, target):
                        continue
                    new_state = state ^ (1 << source) ^ (1 << target)
                    matrix[mapping[new_state], column] += -hopping * fermion_hop_sign(
                        state, source, target
                    )
    return FermiHubbardSector(
        matrix.tocsr(),
        basis,
        {
            "n_sites": n_sites,
            "n_up": n_up,
            "n_down": n_down,
            "hopping": hopping,
            "interaction": interaction,
            "chemical_potential": chemical_potential,
            "periodic": periodic,
        },
    )


def fermi_hubbard_chain_sector_sparse(
    n_sites: int,
    n_up: int,
    n_down: int,
    hopping: float = 1.0,
    interaction: float = 1.0,
    chemical_potential: float = 0.0,
    periodic: bool = False,
) -> FermiHubbardSector:
    """Registered sparse alias for :func:`fermi_hubbard_chain_sector`."""

    return fermi_hubbard_chain_sector(
        n_sites, n_up, n_down, hopping, interaction, chemical_potential, periodic
    )


def fermi_hubbard_observables(
    state: np.ndarray,
    basis: FermiHubbardBasis,
) -> dict[str, np.ndarray]:
    """Return site occupations, spin density, and double occupancy."""

    vector = np.asarray(state, dtype=complex).reshape(-1)
    if vector.size != basis.dimension:
        raise ValueError("State length must match sector dimension.")
    probabilities = np.abs(vector) ** 2
    norm = probabilities.sum()
    if norm == 0:
        raise ValueError("State must have nonzero norm.")
    probabilities /= norm
    up = np.zeros(basis.n_sites)
    down = np.zeros(basis.n_sites)
    double = np.zeros(basis.n_sites)
    for probability, occupation in zip(probabilities, basis.states, strict=True):
        for site in range(basis.n_sites):
            has_up = bit_is_set(occupation, 2 * site)
            has_down = bit_is_set(occupation, 2 * site + 1)
            up[site] += probability * has_up
            down[site] += probability * has_down
            double[site] += probability * has_up * has_down
    return {
        "up_occupation": up,
        "down_occupation": down,
        "total_occupation": up + down,
        "spin_density": 0.5 * (up - down),
        "double_occupancy": double,
    }
