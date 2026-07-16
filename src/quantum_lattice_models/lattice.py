"""User-facing lattice containers and generic tight-binding builders."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from typing import Any, cast

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models._model_utils import onsite_values, validate_positive_int
from quantum_lattice_models.types import LatticeHamiltonian


@dataclass(frozen=True)
class Bond:
    """A directed edge between two lattice sites with an optional matrix element."""

    source: int
    target: int
    value: complex | None = None


@dataclass(frozen=True)
class Lattice:
    """Finite lattice geometry for user-defined single-particle models."""

    n_sites: int
    bonds: tuple[Bond, ...] = ()
    positions: np.ndarray | None = None
    metadata: dict[str, object] | None = None

    def __init__(
        self,
        n_sites: int | None = None,
        *,
        bonds: Iterable[Bond | Sequence[object]] = (),
        positions: Iterable[Sequence[float]] | None = None,
        metadata: dict[str, object] | None = None,
    ) -> None:
        parsed_bonds = tuple(_parse_bond(bond) for bond in bonds)
        parsed_positions = _parse_positions(positions)
        inferred_n_sites = _infer_n_sites(n_sites, parsed_bonds, parsed_positions)
        _validate_site_indices(parsed_bonds, inferred_n_sites)

        object.__setattr__(self, "n_sites", inferred_n_sites)
        object.__setattr__(self, "bonds", parsed_bonds)
        object.__setattr__(self, "positions", parsed_positions)
        object.__setattr__(self, "metadata", dict(metadata or {}))

    @classmethod
    def from_edges(
        cls,
        n_sites: int,
        edges: Iterable[Bond | Sequence[object]],
        *,
        positions: Iterable[Sequence[float]] | None = None,
        metadata: dict[str, object] | None = None,
    ) -> Lattice:
        """Create a lattice from site count and edge-like bond records."""

        return cls(n_sites=n_sites, bonds=edges, positions=positions, metadata=metadata)


class TightBindingModel:
    """Build single-particle tight-binding Hamiltonians on a finite lattice."""

    def __init__(self, lattice: Lattice) -> None:
        self.lattice = lattice

    def hamiltonian(
        self,
        *,
        hopping: complex = 1.0,
        onsite: float | Iterable[float] = 0.0,
        hermitian: bool = True,
        model_name: str = "custom_tight_binding",
    ) -> LatticeHamiltonian:
        """Return a dense Hamiltonian for this lattice."""

        return custom_tight_binding(
            lattice=self.lattice,
            hopping=hopping,
            onsite=onsite,
            hermitian=hermitian,
            model_name=model_name,
        )

    def sparse_hamiltonian(
        self,
        *,
        hopping: complex = 1.0,
        onsite: float | Iterable[float] = 0.0,
        hermitian: bool = True,
    ) -> sp.csr_matrix:
        """Return a sparse Hamiltonian for this lattice."""

        return custom_tight_binding_sparse(
            lattice=self.lattice,
            hopping=hopping,
            onsite=onsite,
            hermitian=hermitian,
        )


def custom_tight_binding(
    n_sites: int | None = None,
    bonds: Iterable[Bond | Sequence[object]] = (),
    *,
    lattice: Lattice | None = None,
    positions: Iterable[Sequence[float]] | None = None,
    hopping: complex = 1.0,
    onsite: float | Iterable[float] = 0.0,
    hermitian: bool = True,
    model_name: str = "custom_tight_binding",
) -> LatticeHamiltonian:
    """Return a dense single-particle Hamiltonian on a user-defined graph.

    Two-item bond records use ``-hopping`` as their matrix element. Three-item
    bond records use the third item directly as the matrix element.
    """

    lattice = _coerce_lattice(n_sites, bonds, lattice, positions)
    matrix = custom_tight_binding_sparse(
        lattice=lattice,
        hopping=hopping,
        onsite=onsite,
        hermitian=hermitian,
    ).toarray()

    metadata = dict(lattice.metadata or {})
    metadata.update({"hermitian": hermitian, "hopping": hopping})
    if lattice.positions is not None:
        metadata["positions"] = lattice.positions.copy()

    return LatticeHamiltonian(
        matrix,
        model_name=model_name,
        basis="single_particle",
        lattice_shape=(lattice.n_sites,),
        metadata=metadata,
    )


def custom_tight_binding_sparse(
    n_sites: int | None = None,
    bonds: Iterable[Bond | Sequence[object]] = (),
    *,
    lattice: Lattice | None = None,
    positions: Iterable[Sequence[float]] | None = None,
    hopping: complex = 1.0,
    onsite: float | Iterable[float] = 0.0,
    hermitian: bool = True,
) -> sp.csr_matrix:
    """Return a sparse single-particle Hamiltonian on a user-defined graph."""

    lattice = _coerce_lattice(n_sites, bonds, lattice, positions)
    matrix = sp.diags(
        onsite_values(onsite, lattice.n_sites),
        format="lil",
        dtype=complex,
    )
    for bond in lattice.bonds:
        value = _bond_value(bond, hopping)
        matrix[bond.source, bond.target] += value
        if hermitian and bond.source != bond.target:
            matrix[bond.target, bond.source] += np.conjugate(value)
    return matrix.tocsr()


def _coerce_lattice(
    n_sites: int | None,
    bonds: Iterable[Bond | Sequence[object]],
    lattice: Lattice | None,
    positions: Iterable[Sequence[float]] | None,
) -> Lattice:
    if lattice is not None:
        if n_sites is not None or positions is not None:
            raise ValueError("Pass either lattice or n_sites/bonds/positions, not both.")
        return lattice
    return Lattice(n_sites=n_sites, bonds=bonds, positions=positions)


def _parse_bond(record: Bond | Sequence[object]) -> Bond:
    if isinstance(record, Bond):
        bond = record
    else:
        values = tuple(record)
        if len(values) not in (2, 3):
            raise ValueError("Each bond must contain source, target, and optional value.")
        value = values[2] if len(values) == 3 else None
        bond = Bond(
            int(cast(Any, values[0])),
            int(cast(Any, values[1])),
            None if value is None else complex(cast(Any, value)),
        )
    if bond.source < 0 or bond.target < 0:
        raise ValueError("Bond site indices must be nonnegative.")
    return bond


def _parse_positions(positions: Iterable[Sequence[float]] | None) -> np.ndarray | None:
    if positions is None:
        return None
    parsed = np.asarray(list(positions), dtype=float)
    if parsed.ndim != 2 or parsed.shape[1] not in (2, 3):
        raise ValueError("positions must be an iterable of 2D or 3D coordinate records.")
    return parsed


def _infer_n_sites(
    n_sites: int | None, bonds: tuple[Bond, ...], positions: np.ndarray | None
) -> int:
    if n_sites is not None:
        validate_positive_int(n_sites, "n_sites")
        return n_sites
    if positions is not None:
        validate_positive_int(len(positions), "n_sites")
        return len(positions)
    if bonds:
        return max(max(bond.source, bond.target) for bond in bonds) + 1
    raise ValueError("n_sites is required when positions and bonds do not define any sites.")


def _validate_site_indices(bonds: tuple[Bond, ...], n_sites: int) -> None:
    for bond in bonds:
        if bond.source >= n_sites or bond.target >= n_sites:
            raise ValueError("Bond site indices must be less than n_sites.")


def _bond_value(bond: Bond, hopping: complex) -> complex:
    if bond.value is not None:
        return complex(bond.value)
    return -complex(hopping)
