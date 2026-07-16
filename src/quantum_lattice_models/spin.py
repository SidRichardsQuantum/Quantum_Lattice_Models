"""Graph-based dense and sparse spin-1/2 Hamiltonian builders."""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from math import comb
from typing import TYPE_CHECKING, Any, cast

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models._model_utils import (
    nearest_neighbor_bonds,
    next_nearest_neighbor_bonds,
    pauli_labels,
    validate_positive_int,
)
from quantum_lattice_models.operators import PAULI_MATRICES
from quantum_lattice_models.reduced import ReducedBasisMapping, reduced_operator
from quantum_lattice_models.types import DenseHamiltonian, PauliTerm

if TYPE_CHECKING:
    from quantum_lattice_models.specs import ModelSpec


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
    def mapping(self) -> ReducedBasisMapping:
        return ReducedBasisMapping(
            kind="fixed_magnetization",
            full_dimension=2**self.n_sites,
            states=self.states,
            quantum_numbers={"magnetization": self.magnetization},
            metadata={"n_sites": self.n_sites},
        )

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

        return self.mapping.embed(state)

    def project(self, state: np.ndarray) -> np.ndarray:
        """Project a full computational-basis state vector into this sector."""

        return self.mapping.project(state)

    def to_metadata(self) -> dict[str, object]:
        """Return portable sector metadata."""

        return {
            "kind": "fixed_magnetization",
            "n_sites": self.n_sites,
            "magnetization": self.magnetization,
            "dimension": self.dimension,
            "basis_states": list(self.states),
            "mapping": self.mapping.to_dict(),
        }


@dataclass(frozen=True)
class SpinParityBasis:
    """Global spin-flip parity basis made from complementary bitstrings."""

    n_sites: int
    parity: int
    representatives: tuple[int, ...]

    @property
    def states(self) -> tuple[int, ...]:
        """Return canonical representative states for compatibility."""

        return self.representatives

    @property
    def mapping(self) -> ReducedBasisMapping:
        normalization = 1.0 / np.sqrt(2.0)
        mask = 2**self.n_sites - 1
        return ReducedBasisMapping(
            kind="spin_flip_parity",
            full_dimension=2**self.n_sites,
            states=self.representatives,
            quantum_numbers={"parity": self.parity},
            labels=tuple(
                f"(|{state:0{self.n_sites}b}> "
                f"{'+' if self.parity == 1 else '-'} "
                f"|{state ^ mask:0{self.n_sites}b}>)/sqrt(2)"
                for state in self.representatives
            ),
            components=tuple(
                (
                    (state, normalization),
                    (state ^ mask, self.parity * normalization),
                )
                for state in self.representatives
            ),
            metadata={
                "n_sites": self.n_sites,
                "symmetry": "global_spin_flip",
                "operator": "X^tensor_n",
            },
        )

    @property
    def dimension(self) -> int:
        """Return the parity-sector dimension."""

        return len(self.representatives)

    def embed(self, state: np.ndarray) -> np.ndarray:
        """Embed a parity-sector vector into the full computational basis."""

        return self.mapping.embed(state)

    def project(self, state: np.ndarray) -> np.ndarray:
        """Project a full state onto this parity sector."""

        return self.mapping.project(state)

    def to_metadata(self) -> dict[str, object]:
        """Return portable parity-sector metadata."""

        return {
            "kind": "spin_flip_parity",
            "n_sites": self.n_sites,
            "parity": self.parity,
            "dimension": self.dimension,
            "basis_states": list(self.representatives),
            "mapping": self.mapping.to_dict(),
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

        return int(self.matrix.shape[0]), int(self.matrix.shape[1])

    def to_metadata(self) -> dict[str, object]:
        """Return portable construction and basis metadata."""

        return {
            "model_name": self.model_name,
            "sector": self.basis.to_metadata(),
            "parameters": dict(self.parameters),
        }


@dataclass(frozen=True)
class SpinParitySectorHamiltonian:
    """Sparse Hamiltonian and explicit basis for global spin-flip parity."""

    matrix: sp.csr_matrix
    basis: SpinParityBasis
    model_name: str
    parameters: dict[str, object]

    @property
    def shape(self) -> tuple[int, int]:
        """Return the reduced matrix shape."""

        return int(self.matrix.shape[0]), int(self.matrix.shape[1])

    def to_metadata(self) -> dict[str, object]:
        """Return portable construction and parity-basis metadata."""

        return {
            "model_name": self.model_name,
            "sector": self.basis.to_metadata(),
            "parameters": dict(self.parameters),
        }


def graph_spin_sector(
    n_sites: int,
    magnetization: int,
    interactions: Iterable[SpinInteraction] = (),
    fields: Iterable[SpinField] = (),
    *,
    tolerance: float = 1e-12,
) -> SpinSectorHamiltonian:
    """Project a magnetization-conserving arbitrary spin graph without densifying."""

    basis = fixed_magnetization_basis(n_sites, magnetization)
    terms = _graph_terms(n_sites, interactions, fields)
    state_to_index = basis.state_to_index
    entries: dict[tuple[int, int], complex] = {}
    leakage: dict[tuple[int, int], complex] = {}
    for column, state in enumerate(basis.states):
        for term in terms:
            target, amplitude = _apply_pauli_labels(state, term.operators)
            value = term.coefficient * amplitude
            if abs(value) <= tolerance:
                continue
            row = state_to_index.get(target)
            destination = entries if row is not None else leakage
            key = (row if row is not None else target, column)
            destination[key] = destination.get(key, 0.0j) + value
    broken = {key: value for key, value in leakage.items() if abs(value) > tolerance}
    if broken:
        raise ValueError(
            "Graph interactions do not conserve the requested total magnetization; "
            f"found {len(broken)} nonzero transitions outside the sector."
        )
    rows = [row for row, _ in entries]
    columns = [column for _, column in entries]
    values = list(entries.values())
    matrix = sp.coo_matrix(
        (values, (rows, columns)),
        shape=(basis.dimension, basis.dimension),
        dtype=complex,
    ).tocsr()
    matrix.sum_duplicates()
    matrix.eliminate_zeros()
    return SpinSectorHamiltonian(
        matrix,
        basis,
        "graph_spin_sector",
        {
            "n_sites": n_sites,
            "magnetization": magnetization,
            "interaction_count": len(terms),
        },
    )


def graph_spin_model_sector(model: ModelSpec, magnetization: int) -> SpinSectorHamiltonian:
    """Build a fixed-magnetization sector directly from a portable graph-spin model."""

    from quantum_lattice_models.specs import GRAPH_SPIN_FAMILY

    if model.family != GRAPH_SPIN_FAMILY:
        raise ValueError("graph_spin_model_sector requires a graph_spin ModelSpec.")
    interactions = []
    fields = []
    for term in model.interactions:
        if len(term.degrees) == 1:
            fields.append(SpinField(term.degrees[0], term.operators[0], term.coefficient))
        else:
            interactions.append(
                SpinInteraction(
                    term.degrees[0],
                    term.degrees[1],
                    term.operators[0],
                    term.operators[1],
                    term.coefficient,
                )
            )
    return graph_spin_sector(
        int(cast(Any, model.parameters["n_sites"])),
        magnetization,
        interactions,
        fields,
    )


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


def spin_flip_parity_basis(n_sites: int, parity: int) -> SpinParityBasis:
    """Return eigenstates of the global spin-flip operator ``prod_i X_i``."""

    validate_positive_int(n_sites, "n_sites")
    if parity not in {-1, 1}:
        raise ValueError("parity must be +1 or -1.")
    mask = 2**n_sites - 1
    representatives = tuple(state for state in range(2**n_sites) if state < (state ^ mask))
    if len(representatives) != 2 ** (n_sites - 1):
        raise RuntimeError("Spin-flip parity basis construction failed.")
    return SpinParityBasis(n_sites, parity, representatives)


def transverse_field_ising_parity_sector(
    n_sites: int,
    parity: int,
    j: float = 1.0,
    h: float = 0.5,
    periodic: bool = False,
) -> SpinParitySectorHamiltonian:
    """Build the transverse-field Ising chain in a spin-flip parity sector."""

    basis = spin_flip_parity_basis(n_sites, parity)
    full = transverse_field_ising_sparse(n_sites, j=j, h=h, periodic=periodic)
    matrix = reduced_operator(full, basis.mapping)
    if not sp.issparse(matrix):
        raise RuntimeError("Parity reduction unexpectedly returned a dense matrix.")
    return SpinParitySectorHamiltonian(
        sp.csr_matrix(matrix),
        basis,
        "transverse_field_ising_parity_sector",
        {
            "n_sites": n_sites,
            "parity": parity,
            "j": j,
            "h": h,
            "periodic": periodic,
        },
    )


def transverse_field_ising_parity_sector_sparse(
    n_sites: int,
    parity: int,
    j: float = 1.0,
    h: float = 0.5,
    periodic: bool = False,
) -> SpinParitySectorHamiltonian:
    """Registered sparse alias for the transverse-field Ising parity sector."""

    return transverse_field_ising_parity_sector(n_sites, parity, j, h, periodic)


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


def heisenberg_ladder_sector(
    n_rungs: int,
    magnetization: int,
    leg_coupling: float = 1.0,
    rung_coupling: float = 1.0,
    field: float = 0.0,
    periodic: bool = False,
) -> SpinSectorHamiltonian:
    """Build a two-leg Heisenberg ladder in a fixed-magnetization sector."""

    validate_positive_int(n_rungs, "n_rungs")
    n_sites = 2 * n_rungs
    basis = fixed_magnetization_basis(n_sites, magnetization)
    leg_bonds = nearest_neighbor_bonds(n_rungs, periodic)
    weighted_bonds = [(i, j, leg_coupling) for i, j in leg_bonds]
    weighted_bonds.extend((n_rungs + i, n_rungs + j, leg_coupling) for i, j in leg_bonds)
    weighted_bonds.extend((rung, n_rungs + rung, rung_coupling) for rung in range(n_rungs))
    matrix = sp.csr_matrix((basis.dimension, basis.dimension), dtype=complex)
    for source, target, coupling in weighted_bonds:
        matrix += _xxz_sector_matrix(
            basis,
            ((source, target),),
            transverse_coupling=coupling,
            longitudinal_coupling=coupling,
            field=0.0,
        )
    if field:
        matrix += sp.eye(basis.dimension, format="csr") * field * magnetization
    return SpinSectorHamiltonian(
        matrix,
        basis,
        "heisenberg_ladder_sector",
        {
            "n_rungs": n_rungs,
            "magnetization": magnetization,
            "leg_coupling": leg_coupling,
            "rung_coupling": rung_coupling,
            "field": field,
            "periodic": periodic,
        },
    )


def heisenberg_ladder_sector_sparse(
    n_rungs: int,
    magnetization: int,
    leg_coupling: float = 1.0,
    rung_coupling: float = 1.0,
    field: float = 0.0,
    periodic: bool = False,
) -> SpinSectorHamiltonian:
    """Registered sparse alias for :func:`heisenberg_ladder_sector`."""

    return heisenberg_ladder_sector(
        n_rungs, magnetization, leg_coupling, rung_coupling, field, periodic
    )


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


def _apply_pauli_labels(state: int, operators: tuple[str, ...]) -> tuple[int, complex]:
    target = state
    amplitude = 1.0 + 0.0j
    n_sites = len(operators)
    for site, label in enumerate(operators):
        mask = 1 << (n_sites - 1 - site)
        bit = bool(target & mask)
        if label == "X":
            target ^= mask
        elif label == "Y":
            amplitude *= -1j if bit else 1j
            target ^= mask
        elif label == "Z":
            amplitude *= -1.0 if bit else 1.0
    return target, amplitude


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
