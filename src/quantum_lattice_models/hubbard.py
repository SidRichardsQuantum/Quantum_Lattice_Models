"""Bose-Hubbard and Fermi-Hubbard Hamiltonian builders."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models._model_utils import (
    bit_is_set,
    fermion_hop_sign,
    nearest_neighbor_bonds,
    occupation_digits,
    occupations_to_index,
    validate_positive_int,
)
from quantum_lattice_models.types import LatticeHamiltonian


def bose_hubbard_chain(
    n_sites: int,
    hopping: float = 1.0,
    interaction: float = 1.0,
    chemical_potential: float = 0.0,
    max_occupancy: int = 2,
    periodic: bool = False,
) -> np.ndarray:
    """Return a dense truncated Bose-Hubbard chain Hamiltonian."""

    sparse = bose_hubbard_chain_sparse(
        n_sites=n_sites,
        hopping=hopping,
        interaction=interaction,
        chemical_potential=chemical_potential,
        max_occupancy=max_occupancy,
        periodic=periodic,
    )
    return LatticeHamiltonian(
        sparse.toarray(),
        model_name="bose_hubbard_chain",
        basis="truncated_boson_occupation",
        lattice_shape=(n_sites,),
        metadata={
            "periodic": periodic,
            "hopping": hopping,
            "interaction": interaction,
            "chemical_potential": chemical_potential,
            "max_occupancy": max_occupancy,
        },
    )


def bose_hubbard_chain_sparse(
    n_sites: int,
    hopping: float = 1.0,
    interaction: float = 1.0,
    chemical_potential: float = 0.0,
    max_occupancy: int = 2,
    periodic: bool = False,
) -> sp.csr_matrix:
    """Return a sparse truncated Bose-Hubbard chain Hamiltonian."""

    validate_positive_int(n_sites, "n_sites")
    validate_positive_int(max_occupancy, "max_occupancy")
    local_dim = max_occupancy + 1
    dim = local_dim**n_sites
    matrix = sp.lil_matrix((dim, dim), dtype=complex)

    for state in range(dim):
        occupations = occupation_digits(state, n_sites, local_dim)
        diagonal = sum(
            0.5 * interaction * n * (n - 1) - chemical_potential * n for n in occupations
        )
        matrix[state, state] = diagonal
        for i, k in nearest_neighbor_bonds(n_sites, periodic):
            for source, target in ((k, i), (i, k)):
                if occupations[source] == 0 or occupations[target] >= max_occupancy:
                    continue
                new_occupations = occupations.copy()
                amplitude = -hopping * np.sqrt(occupations[source] * (occupations[target] + 1))
                new_occupations[source] -= 1
                new_occupations[target] += 1
                new_state = occupations_to_index(new_occupations, local_dim)
                matrix[new_state, state] += amplitude

    return matrix.tocsr()


def fermi_hubbard_chain(
    n_sites: int,
    hopping: float = 1.0,
    interaction: float = 1.0,
    chemical_potential: float = 0.0,
    periodic: bool = False,
) -> np.ndarray:
    """Return a dense spinful Fermi-Hubbard chain in occupation-number basis."""

    sparse = fermi_hubbard_chain_sparse(
        n_sites=n_sites,
        hopping=hopping,
        interaction=interaction,
        chemical_potential=chemical_potential,
        periodic=periodic,
    )
    return LatticeHamiltonian(
        sparse.toarray(),
        model_name="fermi_hubbard_chain",
        basis="spinful_fermion_occupation",
        lattice_shape=(n_sites,),
        metadata={
            "periodic": periodic,
            "hopping": hopping,
            "interaction": interaction,
            "chemical_potential": chemical_potential,
            "orbital_order": "site0_up, site0_down, site1_up, site1_down, ...",
        },
    )


def fermi_hubbard_chain_sparse(
    n_sites: int,
    hopping: float = 1.0,
    interaction: float = 1.0,
    chemical_potential: float = 0.0,
    periodic: bool = False,
) -> sp.csr_matrix:
    """Return a sparse spinful Fermi-Hubbard chain Hamiltonian."""

    validate_positive_int(n_sites, "n_sites")
    n_orbitals = 2 * n_sites
    dim = 2**n_orbitals
    matrix = sp.lil_matrix((dim, dim), dtype=complex)

    for state in range(dim):
        diagonal = 0.0
        for site in range(n_sites):
            up = bit_is_set(state, 2 * site)
            down = bit_is_set(state, 2 * site + 1)
            diagonal += interaction * up * down
            diagonal += -chemical_potential * (up + down)
        matrix[state, state] = diagonal

        for i, k in nearest_neighbor_bonds(n_sites, periodic):
            for spin in (0, 1):
                for source, target in ((2 * i + spin, 2 * k + spin), (2 * k + spin, 2 * i + spin)):
                    if not bit_is_set(state, source) or bit_is_set(state, target):
                        continue
                    sign = fermion_hop_sign(state, source, target)
                    new_state = state ^ (1 << source) ^ (1 << target)
                    matrix[new_state, state] += -hopping * sign

    return matrix.tocsr()
