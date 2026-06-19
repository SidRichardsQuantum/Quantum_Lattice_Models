"""Graph-based dense and sparse spin-1/2 Hamiltonian builders."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from math import comb

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models._model_utils import (
    nearest_neighbor_bonds,
    next_nearest_neighbor_bonds,
    pauli_labels,
    validate_positive_int,
)
from quantum_lattice_models.operators import PAULI_MATRICES
from quantum_lattice_models.types import DenseHamiltonian, PauliTerm


@dataclass(frozen=True)
class SpinInteraction:
    """One two-site Pauli interaction on an arbitrary spin graph."""

    source: int
    target: int
    source_axis: str
    target_axis: str
    coefficient: complex


@dataclass(frozen=True)
class SpinField:
    """One site-resolved Pauli field."""

    site: int
    axis: str
    coefficient: complex


@dataclass(frozen=True)
class FixedMagnetizationBasis:
    """Computational basis restricted to a total Pauli-Z eigenvalue."""

    n_sites: int
    magnetization: int
    states: tuple[int, ...]

    @property
    def dimension(self) -> int:
        """Return the reduced Hilbert-space dimension."""

        return len(self.states)

    @property
    def state_to_index(self) -> dict[int, int]:
        """Map full computational-basis integer states to sector indices."""

        return {state: index for index, state in enumerate(self.states)}

    def embed(self, state: np.ndarray) -> np.ndarray:
        """Embed a reduced state vector into the full computational basis."""

        vector = np.asarray(state, dtype=complex).reshape(-1)
        if vector.size != self.dimension:
            raise ValueError("Reduced state length must match the sector dimension.")
        full = np.zeros(2**self.n_sites, dtype=complex)
        full[np.asarray(self.states, dtype=int)] = vector
        return full

    def project(self, state: np.ndarray) -> np.ndarray:
        """Project a full computational-basis state vector into this sector."""

        vector = np.asarray(state, dtype=complex).reshape(-1)
        if vector.size != 2**self.n_sites:
            raise ValueError("Full state length must equal 2**n_sites.")
        return vector[np.asarray(self.states, dtype=int)]

    def to_metadata(self) -> dict[str, object]:
        """Return portable sector metadata."""

        return {
            "kind": "fixed_magnetization",
            "n_sites": self.n_sites,
            "magnetization": self.magnetization,
            "dimension": self.dimension,
            "basis_states": list(self.states),
        }


@dataclass(frozen=True)
class SpinSectorHamiltonian:
    """Sparse Hamiltonian and explicit basis for a conserved spin sector."""

    matrix: sp.csr_matrix
    basis: FixedMagnetizationBasis
    model_name: str
    parameters: dict[str, object]

    @property
    def shape(self) -> tuple[int, int]:
        """Return the reduced matrix shape."""

        return self.matrix.shape

    def to_metadata(self) -> dict[str, object]:
        """Return portable construction and basis metadata."""

        return {
            "model_name": self.model_name,
            "sector": self.basis.to_metadata(),
            "parameters": dict(self.parameters),
        }


def fixed_magnetization_basis(
    n_sites: int,
    magnetization: int,
) -> FixedMagnetizationBasis:
    """Return basis states with ``sum_i Z_i = magnetization``."""

    validate_positive_int(n_sites, "n_sites")
    if not isinstance(magnetization, int) or isinstance(magnetization, bool):
        raise ValueError("magnetization must be an integer.")
    if abs(magnetization) > n_sites:
        raise ValueError("magnetization must satisfy abs(magnetization) <= n_sites.")
    if (n_sites - magnetization) % 2:
        raise ValueError("n_sites and magnetization must have the same parity.")
    n_ones = (n_sites - magnetization) // 2
    states = tuple(state for state in range(2**n_sites) if state.bit_count() == n_ones)
    if len(states) != comb(n_sites, n_ones):
        raise RuntimeError("Fixed-magnetization basis construction failed.")
    return FixedMagnetizationBasis(n_sites, magnetization, states)


def heisenberg_chain_sector(
    n_sites: int,
    magnetization: int,
    jx: float = 1.0,
    jy: float = 1.0,
    jz: float = 1.0,
    field: float = 0.0,
    periodic: bool = False,
) -> SpinSectorHamiltonian:
    """Build an anisotropic Heisenberg chain in a fixed-magnetization sector.

    Fixed magnetization is conserved only when ``jx == jy``.
    """

    if not np.isclose(jx, jy):
        raise ValueError("Fixed-magnetization Heisenberg sectors require jx == jy.")
    basis = fixed_magnetization_basis(n_sites, magnetization)
    matrix = _xxz_sector_matrix(
        basis,
        nearest_neighbor_bonds(n_sites, periodic),
        transverse_coupling=(jx + jy) / 2.0,
        longitudinal_coupling=jz,
        field=field,
    )
    return SpinSectorHamiltonian(
        matrix,
        basis,
        "heisenberg_chain_sector",
        {
            "n_sites": n_sites,
            "magnetization": magnetization,
            "jx": jx,
            "jy": jy,
            "jz": jz,
            "field": field,
            "periodic": periodic,
        },
    )


def xxz_chain_sector(
    n_sites: int,
    magnetization: int,
    coupling: float = 1.0,
    anisotropy: float = 1.0,
    field: float = 0.0,
    periodic: bool = False,
) -> SpinSectorHamiltonian:
    """Build an XXZ chain in a fixed total Pauli-Z sector."""

    basis = fixed_magnetization_basis(n_sites, magnetization)
    matrix = _xxz_sector_matrix(
        basis,
        nearest_neighbor_bonds(n_sites, periodic),
        transverse_coupling=coupling,
        longitudinal_coupling=coupling * anisotropy,
        field=field,
    )
    return SpinSectorHamiltonian(
        matrix,
        basis,
        "xxz_chain_sector",
        {
            "n_sites": n_sites,
            "magnetization": magnetization,
            "coupling": coupling,
            "anisotropy": anisotropy,
            "field": field,
            "periodic": periodic,
        },
    )


def heisenberg_chain_sector_sparse(
    n_sites: int,
    magnetization: int,
    jx: float = 1.0,
    jy: float = 1.0,
    jz: float = 1.0,
    field: float = 0.0,
    periodic: bool = False,
) -> SpinSectorHamiltonian:
    """Registered sparse alias for :func:`heisenberg_chain_sector`."""

    return heisenberg_chain_sector(n_sites, magnetization, jx, jy, jz, field, periodic)


def xxz_chain_sector_sparse(
    n_sites: int,
    magnetization: int,
    coupling: float = 1.0,
    anisotropy: float = 1.0,
    field: float = 0.0,
    periodic: bool = False,
) -> SpinSectorHamiltonian:
    """Registered sparse alias for :func:`xxz_chain_sector`."""

    return xxz_chain_sector(n_sites, magnetization, coupling, anisotropy, field, periodic)


def graph_spin_hamiltonian_sparse(
    n_sites: int,
    interactions: Iterable[SpinInteraction] = (),
    fields: Iterable[SpinField] = (),
) -> sp.csr_matrix:
    """Build a sparse spin-1/2 Hamiltonian directly from graph interactions."""

    terms = _graph_terms(n_sites, interactions, fields)
    dimension = 2**n_sites
    matrix = sp.csr_matrix((dimension, dimension), dtype=complex)
    for term in terms:
        matrix += term.coefficient * _sparse_pauli_string(term.operators)
    matrix.eliminate_zeros()
    return matrix


def graph_spin_hamiltonian(
    n_sites: int,
    interactions: Iterable[SpinInteraction] = (),
    fields: Iterable[SpinField] = (),
    *,
    model_name: str = "graph_spin_hamiltonian",
) -> DenseHamiltonian:
    """Build a dense graph spin Hamiltonian through the shared sparse backend."""

    terms = _graph_terms(n_sites, interactions, fields)
    matrix = _sparse_from_terms(n_sites, terms).toarray()
    return DenseHamiltonian(matrix, model_name=model_name, n_sites=n_sites, terms=terms)


def transverse_field_ising(
    n_sites: int, j: float = 1.0, h: float = 0.5, periodic: bool = False
) -> DenseHamiltonian:
    """Return ``H = -J sum ZZ - h sum X``."""

    interactions = _axis_interactions(nearest_neighbor_bonds(n_sites, periodic), "Z", -j)
    fields = tuple(SpinField(i, "X", -h) for i in range(n_sites))
    return graph_spin_hamiltonian(
        n_sites, interactions, fields, model_name="transverse_field_ising"
    )


def transverse_field_ising_sparse(
    n_sites: int, j: float = 1.0, h: float = 0.5, periodic: bool = False
) -> sp.csr_matrix:
    """Return the sparse transverse-field Ising Hamiltonian."""

    interactions = _axis_interactions(nearest_neighbor_bonds(n_sites, periodic), "Z", -j)
    fields = tuple(SpinField(i, "X", -h) for i in range(n_sites))
    return graph_spin_hamiltonian_sparse(n_sites, interactions, fields)


def longitudinal_field_ising(
    n_sites: int,
    j: float = 1.0,
    h_x: float = 0.5,
    h_z: float = 0.0,
    periodic: bool = False,
) -> DenseHamiltonian:
    """Return ``H = -J sum ZZ - h_x sum X - h_z sum Z``."""

    interactions = _axis_interactions(nearest_neighbor_bonds(n_sites, periodic), "Z", -j)
    fields = tuple(
        field
        for i in range(n_sites)
        for field in (SpinField(i, "X", -h_x), SpinField(i, "Z", -h_z))
    )
    return graph_spin_hamiltonian(
        n_sites, interactions, fields, model_name="longitudinal_field_ising"
    )


def longitudinal_field_ising_sparse(
    n_sites: int,
    j: float = 1.0,
    h_x: float = 0.5,
    h_z: float = 0.0,
    periodic: bool = False,
) -> sp.csr_matrix:
    """Return the sparse Ising Hamiltonian with transverse and longitudinal fields."""

    interactions = _axis_interactions(nearest_neighbor_bonds(n_sites, periodic), "Z", -j)
    fields = tuple(
        field
        for i in range(n_sites)
        for field in (SpinField(i, "X", -h_x), SpinField(i, "Z", -h_z))
    )
    return graph_spin_hamiltonian_sparse(n_sites, interactions, fields)


def next_nearest_neighbor_ising(
    n_sites: int,
    j1: float = 1.0,
    j2: float = 0.5,
    h: float = 0.5,
    periodic: bool = False,
) -> DenseHamiltonian:
    """Return ``H = -J1 sum ZiZi+1 - J2 sum ZiZi+2 - h sum X``."""

    interactions = _axis_interactions(nearest_neighbor_bonds(n_sites, periodic), "Z", -j1)
    interactions += _axis_interactions(next_nearest_neighbor_bonds(n_sites, periodic), "Z", -j2)
    fields = tuple(SpinField(i, "X", -h) for i in range(n_sites))
    return graph_spin_hamiltonian(
        n_sites, interactions, fields, model_name="next_nearest_neighbor_ising"
    )


def next_nearest_neighbor_ising_sparse(
    n_sites: int,
    j1: float = 1.0,
    j2: float = 0.5,
    h: float = 0.5,
    periodic: bool = False,
) -> sp.csr_matrix:
    """Return the sparse next-nearest-neighbor Ising Hamiltonian."""

    interactions = _axis_interactions(nearest_neighbor_bonds(n_sites, periodic), "Z", -j1)
    interactions += _axis_interactions(next_nearest_neighbor_bonds(n_sites, periodic), "Z", -j2)
    fields = tuple(SpinField(i, "X", -h) for i in range(n_sites))
    return graph_spin_hamiltonian_sparse(n_sites, interactions, fields)


def heisenberg_chain(
    n_sites: int,
    jx: float = 1.0,
    jy: float = 1.0,
    jz: float = 1.0,
    field: float = 0.0,
    periodic: bool = False,
) -> DenseHamiltonian:
    """Return the dense anisotropic Heisenberg spin-chain Hamiltonian."""

    interactions = _xyz_interactions(nearest_neighbor_bonds(n_sites, periodic), jx, jy, jz)
    fields = tuple(SpinField(i, "Z", field) for i in range(n_sites))
    return graph_spin_hamiltonian(n_sites, interactions, fields, model_name="heisenberg_chain")


def heisenberg_chain_sparse(
    n_sites: int,
    jx: float = 1.0,
    jy: float = 1.0,
    jz: float = 1.0,
    field: float = 0.0,
    periodic: bool = False,
) -> sp.csr_matrix:
    """Return the sparse anisotropic Heisenberg spin-chain Hamiltonian."""

    interactions = _xyz_interactions(nearest_neighbor_bonds(n_sites, periodic), jx, jy, jz)
    fields = tuple(SpinField(i, "Z", field) for i in range(n_sites))
    return graph_spin_hamiltonian_sparse(n_sites, interactions, fields)


def xy_chain(
    n_sites: int,
    coupling: float = 1.0,
    anisotropy: float = 0.0,
    field: float = 0.0,
    periodic: bool = False,
) -> DenseHamiltonian:
    """Return the dense anisotropic XY spin-chain Hamiltonian."""

    xx = -coupling * (1.0 + anisotropy) / 2.0
    yy = -coupling * (1.0 - anisotropy) / 2.0
    bonds = nearest_neighbor_bonds(n_sites, periodic)
    interactions = _axis_interactions(bonds, "X", xx)
    interactions += _axis_interactions(bonds, "Y", yy)
    fields = tuple(SpinField(i, "Z", -field) for i in range(n_sites))
    return graph_spin_hamiltonian(n_sites, interactions, fields, model_name="xy_chain")


def xy_chain_sparse(
    n_sites: int,
    coupling: float = 1.0,
    anisotropy: float = 0.0,
    field: float = 0.0,
    periodic: bool = False,
) -> sp.csr_matrix:
    """Return the sparse anisotropic XY spin-chain Hamiltonian."""

    xx = -coupling * (1.0 + anisotropy) / 2.0
    yy = -coupling * (1.0 - anisotropy) / 2.0
    bonds = nearest_neighbor_bonds(n_sites, periodic)
    interactions = _axis_interactions(bonds, "X", xx)
    interactions += _axis_interactions(bonds, "Y", yy)
    fields = tuple(SpinField(i, "Z", -field) for i in range(n_sites))
    return graph_spin_hamiltonian_sparse(n_sites, interactions, fields)


def xxz_chain(
    n_sites: int,
    coupling: float = 1.0,
    anisotropy: float = 1.0,
    field: float = 0.0,
    periodic: bool = False,
) -> DenseHamiltonian:
    """Return the dense XXZ spin-chain Hamiltonian."""

    matrix = heisenberg_chain(n_sites, coupling, coupling, coupling * anisotropy, field, periodic)
    matrix.model_name = "xxz_chain"
    return matrix


def xxz_chain_sparse(
    n_sites: int,
    coupling: float = 1.0,
    anisotropy: float = 1.0,
    field: float = 0.0,
    periodic: bool = False,
) -> sp.csr_matrix:
    """Return the sparse XXZ spin-chain Hamiltonian."""

    return heisenberg_chain_sparse(
        n_sites, coupling, coupling, coupling * anisotropy, field, periodic
    )


def j1_j2_heisenberg_chain(
    n_sites: int,
    j1: float = 1.0,
    j2: float = 0.5,
    field: float = 0.0,
    periodic: bool = False,
) -> DenseHamiltonian:
    """Return a dense frustrated J1-J2 Heisenberg spin-chain Hamiltonian."""

    interactions = _xyz_interactions(nearest_neighbor_bonds(n_sites, periodic), j1, j1, j1)
    interactions += _xyz_interactions(next_nearest_neighbor_bonds(n_sites, periodic), j2, j2, j2)
    fields = tuple(SpinField(i, "Z", field) for i in range(n_sites))
    return graph_spin_hamiltonian(
        n_sites, interactions, fields, model_name="j1_j2_heisenberg_chain"
    )


def j1_j2_heisenberg_chain_sparse(
    n_sites: int,
    j1: float = 1.0,
    j2: float = 0.5,
    field: float = 0.0,
    periodic: bool = False,
) -> sp.csr_matrix:
    """Return the sparse frustrated J1-J2 Heisenberg spin-chain Hamiltonian."""

    interactions = _xyz_interactions(nearest_neighbor_bonds(n_sites, periodic), j1, j1, j1)
    interactions += _xyz_interactions(next_nearest_neighbor_bonds(n_sites, periodic), j2, j2, j2)
    fields = tuple(SpinField(i, "Z", field) for i in range(n_sites))
    return graph_spin_hamiltonian_sparse(n_sites, interactions, fields)


def heisenberg_ladder(
    n_rungs: int,
    leg_coupling: float = 1.0,
    rung_coupling: float = 1.0,
    field: float = 0.0,
    periodic: bool = False,
) -> DenseHamiltonian:
    """Return a dense two-leg Heisenberg ladder with ``2 * n_rungs`` spins."""

    n_sites, interactions, fields = _ladder_graph(
        n_rungs, leg_coupling, rung_coupling, field, periodic
    )
    return graph_spin_hamiltonian(n_sites, interactions, fields, model_name="heisenberg_ladder")


def heisenberg_ladder_sparse(
    n_rungs: int,
    leg_coupling: float = 1.0,
    rung_coupling: float = 1.0,
    field: float = 0.0,
    periodic: bool = False,
) -> sp.csr_matrix:
    """Return the sparse two-leg Heisenberg ladder Hamiltonian."""

    n_sites, interactions, fields = _ladder_graph(
        n_rungs, leg_coupling, rung_coupling, field, periodic
    )
    return graph_spin_hamiltonian_sparse(n_sites, interactions, fields)


def _graph_terms(
    n_sites: int,
    interactions: Iterable[SpinInteraction],
    fields: Iterable[SpinField],
) -> tuple[PauliTerm, ...]:
    validate_positive_int(n_sites, "n_sites")
    terms: list[PauliTerm] = []
    for interaction in interactions:
        _validate_axis(interaction.source_axis)
        _validate_axis(interaction.target_axis)
        _validate_site(interaction.source, n_sites)
        _validate_site(interaction.target, n_sites)
        if interaction.source == interaction.target:
            raise ValueError("Spin interactions require distinct source and target sites.")
        terms.append(
            PauliTerm(
                complex(interaction.coefficient),
                pauli_labels(
                    n_sites,
                    {
                        interaction.source: interaction.source_axis,
                        interaction.target: interaction.target_axis,
                    },
                ),
            )
        )
    for field in fields:
        _validate_axis(field.axis)
        _validate_site(field.site, n_sites)
        terms.append(
            PauliTerm(
                complex(field.coefficient),
                pauli_labels(n_sites, {field.site: field.axis}),
            )
        )
    return tuple(terms)


def _sparse_from_terms(n_sites: int, terms: tuple[PauliTerm, ...]) -> sp.csr_matrix:
    dimension = 2**n_sites
    matrix = sp.csr_matrix((dimension, dimension), dtype=complex)
    for term in terms:
        matrix += term.coefficient * _sparse_pauli_string(term.operators)
    matrix.eliminate_zeros()
    return matrix


def _sparse_pauli_string(operators: tuple[str, ...]) -> sp.csr_matrix:
    result = sp.csr_matrix([[1.0 + 0.0j]])
    for label in operators:
        result = sp.kron(result, sp.csr_matrix(PAULI_MATRICES[label]), format="csr")
    return result


def _xxz_sector_matrix(
    basis: FixedMagnetizationBasis,
    bonds: Iterable[tuple[int, int]],
    *,
    transverse_coupling: complex,
    longitudinal_coupling: complex,
    field: complex,
) -> sp.csr_matrix:
    state_to_index = basis.state_to_index
    bond_list = tuple(bonds)
    rows: list[int] = []
    columns: list[int] = []
    values: list[complex] = []
    for column, state in enumerate(basis.states):
        diagonal = complex(field * basis.magnetization)
        for source, target in bond_list:
            source_mask = 1 << (basis.n_sites - 1 - source)
            target_mask = 1 << (basis.n_sites - 1 - target)
            source_bit = bool(state & source_mask)
            target_bit = bool(state & target_mask)
            diagonal += longitudinal_coupling * (1 if source_bit == target_bit else -1)
            if source_bit != target_bit:
                flipped = state ^ source_mask ^ target_mask
                rows.append(state_to_index[flipped])
                columns.append(column)
                values.append(2.0 * transverse_coupling)
        rows.append(column)
        columns.append(column)
        values.append(diagonal)
    matrix = sp.coo_matrix(
        (values, (rows, columns)),
        shape=(basis.dimension, basis.dimension),
        dtype=complex,
    ).tocsr()
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    return matrix


def _axis_interactions(
    bonds: Iterable[tuple[int, int]], axis: str, coefficient: complex
) -> tuple[SpinInteraction, ...]:
    return tuple(SpinInteraction(i, j, axis, axis, coefficient) for i, j in bonds)


def _xyz_interactions(
    bonds: Iterable[tuple[int, int]], jx: complex, jy: complex, jz: complex
) -> tuple[SpinInteraction, ...]:
    return tuple(
        interaction
        for i, j in bonds
        for interaction in (
            SpinInteraction(i, j, "X", "X", jx),
            SpinInteraction(i, j, "Y", "Y", jy),
            SpinInteraction(i, j, "Z", "Z", jz),
        )
    )


def _ladder_graph(
    n_rungs: int,
    leg_coupling: float,
    rung_coupling: float,
    field: float,
    periodic: bool,
) -> tuple[int, tuple[SpinInteraction, ...], tuple[SpinField, ...]]:
    validate_positive_int(n_rungs, "n_rungs")
    n_sites = 2 * n_rungs
    leg_bonds = nearest_neighbor_bonds(n_rungs, periodic)
    bonds = [(i, j, leg_coupling) for i, j in leg_bonds]
    bonds.extend((n_rungs + i, n_rungs + j, leg_coupling) for i, j in leg_bonds)
    bonds.extend((rung, n_rungs + rung, rung_coupling) for rung in range(n_rungs))
    interactions = tuple(
        interaction
        for i, j, coupling in bonds
        for interaction in (
            SpinInteraction(i, j, "X", "X", coupling),
            SpinInteraction(i, j, "Y", "Y", coupling),
            SpinInteraction(i, j, "Z", "Z", coupling),
        )
    )
    fields = tuple(SpinField(i, "Z", field) for i in range(n_sites))
    return n_sites, interactions, fields


def _validate_axis(axis: str) -> None:
    if axis not in {"X", "Y", "Z"}:
        raise ValueError("Spin axes must be 'X', 'Y', or 'Z'.")


def _validate_site(site: int, n_sites: int) -> None:
    if not isinstance(site, int) or not 0 <= site < n_sites:
        raise ValueError(f"Spin site must satisfy 0 <= site < {n_sites}.")
