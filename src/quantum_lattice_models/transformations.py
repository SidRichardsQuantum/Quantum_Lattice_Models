"""Immutable transformations for portable finite-lattice specifications."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import replace

import numpy as np

from quantum_lattice_models.lattice import Bond
from quantum_lattice_models.specs import LatticeSpec


def relabel_lattice(lattice: LatticeSpec, order: Sequence[int]) -> LatticeSpec:
    """Return a lattice whose new site order is given by old site indices."""

    lattice.validate()
    old_sites = tuple(order)
    if sorted(old_sites) != list(range(lattice.n_sites)):
        raise ValueError("order must be a permutation of all lattice site indices.")
    old_to_new = {old: new for new, old in enumerate(old_sites)}
    bonds = tuple(
        Bond(old_to_new[bond.source], old_to_new[bond.target], bond.value) for bond in lattice.bonds
    )
    return _transformed(
        lattice,
        "relabel",
        {"order": list(old_sites)},
        positions=_select(lattice.positions, old_sites),
        bonds=bonds,
        site_labels=_select(lattice.site_labels, old_sites),
        orbital_labels=_select(lattice.orbital_labels, old_sites),
        sublattice_labels=_select(lattice.sublattice_labels, old_sites),
        unit_cells=_select(lattice.unit_cells, old_sites),
    )


def lattice_subgraph(lattice: LatticeSpec, sites: Iterable[int]) -> LatticeSpec:
    """Return the induced subgraph on selected sites, preserving their order."""

    lattice.validate()
    selected = tuple(sites)
    if not selected:
        raise ValueError("sites must contain at least one site.")
    if len(set(selected)) != len(selected):
        raise ValueError("sites must be unique.")
    if any(not isinstance(site, int) or not 0 <= site < lattice.n_sites for site in selected):
        raise ValueError("sites must contain valid lattice indices.")
    old_to_new = {old: new for new, old in enumerate(selected)}
    bonds = tuple(
        Bond(old_to_new[bond.source], old_to_new[bond.target], bond.value)
        for bond in lattice.bonds
        if bond.source in old_to_new and bond.target in old_to_new
    )
    return _transformed(
        lattice,
        "subgraph",
        {"sites": list(selected)},
        n_sites=len(selected),
        positions=_select(lattice.positions, selected),
        bonds=bonds,
        site_labels=_select(lattice.site_labels, selected),
        orbital_labels=_select(lattice.orbital_labels, selected),
        sublattice_labels=_select(lattice.sublattice_labels, selected),
        unit_cells=_select(lattice.unit_cells, selected),
    )


def remove_lattice_sites(lattice: LatticeSpec, sites: Iterable[int]) -> LatticeSpec:
    """Remove vacancy sites and all incident bonds."""

    removed = tuple(sorted(set(sites)))
    if any(not isinstance(site, int) or not 0 <= site < lattice.n_sites for site in removed):
        raise ValueError("sites must contain valid lattice indices.")
    retained = tuple(site for site in range(lattice.n_sites) if site not in removed)
    if not retained:
        raise ValueError("A transformed lattice must retain at least one site.")
    transformed = lattice_subgraph(lattice, retained)
    provenance = lattice.provenance + (
        {"operation": "remove_sites", "parameters": {"sites": list(removed)}},
    )
    return replace(transformed, provenance=provenance)


def with_boundary_conditions(
    lattice: LatticeSpec,
    boundary_conditions: dict[str, str],
) -> LatticeSpec:
    """Return a lattice with replaced boundary-condition metadata."""

    return _transformed(
        lattice,
        "set_boundary_conditions",
        {"boundary_conditions": dict(boundary_conditions)},
        boundary_conditions=dict(boundary_conditions),
    )


def add_lattice_bonds(
    lattice: LatticeSpec,
    bonds: Iterable[Bond | Sequence[object]],
) -> LatticeSpec:
    """Append validated bonds to a lattice."""

    additions = tuple(_bond(record) for record in bonds)
    return _transformed(
        lattice,
        "add_bonds",
        {"bonds": [_bond_record(bond) for bond in additions]},
        bonds=lattice.bonds + additions,
    )


def remove_lattice_bonds(
    lattice: LatticeSpec,
    bonds: Iterable[tuple[int, int]],
) -> LatticeSpec:
    """Remove every bond matching a source-target pair."""

    removed = {tuple(pair) for pair in bonds}
    if any(len(pair) != 2 for pair in removed):
        raise ValueError("bonds must contain source-target pairs.")
    retained = tuple(bond for bond in lattice.bonds if (bond.source, bond.target) not in removed)
    return _transformed(
        lattice,
        "remove_bonds",
        {"bonds": [list(pair) for pair in sorted(removed)]},
        bonds=retained,
    )


def apply_onsite_disorder(
    lattice: LatticeSpec,
    strength: float,
    *,
    seed: int,
    distribution: str = "uniform",
) -> LatticeSpec:
    """Attach reproducible onsite disorder values to lattice metadata."""

    values = _disorder_values(lattice.n_sites, strength, seed, distribution)
    metadata = dict(lattice.metadata)
    metadata["onsite_potential"] = values.tolist()
    return _transformed(
        lattice,
        "onsite_disorder",
        {"strength": strength, "seed": seed, "distribution": distribution},
        metadata=metadata,
    )


def apply_bond_disorder(
    lattice: LatticeSpec,
    strength: float,
    *,
    seed: int,
    distribution: str = "uniform",
    default_value: complex = -1.0,
) -> LatticeSpec:
    """Perturb explicit bond matrix elements reproducibly."""

    offsets = _disorder_values(len(lattice.bonds), strength, seed, distribution)
    bonds = tuple(
        Bond(
            bond.source,
            bond.target,
            complex(default_value if bond.value is None else bond.value) + offset,
        )
        for bond, offset in zip(lattice.bonds, offsets, strict=True)
    )
    return _transformed(
        lattice,
        "bond_disorder",
        {
            "strength": strength,
            "seed": seed,
            "distribution": distribution,
            "default_value": _portable_complex(default_value),
        },
        bonds=bonds,
    )


def _transformed(
    lattice: LatticeSpec,
    operation: str,
    parameters: dict[str, object],
    **changes: object,
) -> LatticeSpec:
    provenance = lattice.provenance + ({"operation": operation, "parameters": parameters},)
    transformed = replace(lattice, provenance=provenance, **changes)
    transformed.validate()
    return transformed


def _select(values: tuple, indices: Sequence[int]) -> tuple:
    return tuple(values[index] for index in indices) if values else ()


def _bond(record: Bond | Sequence[object]) -> Bond:
    if isinstance(record, Bond):
        return record
    values = tuple(record)
    if len(values) not in {2, 3}:
        raise ValueError("Bond records require source, target, and optional value.")
    return Bond(
        int(values[0]),
        int(values[1]),
        None if len(values) == 2 else complex(values[2]),
    )


def _bond_record(bond: Bond) -> dict[str, object]:
    return {
        "source": bond.source,
        "target": bond.target,
        "value": None if bond.value is None else _portable_complex(bond.value),
    }


def _portable_complex(value: complex) -> dict[str, list[float]]:
    number = complex(value)
    return {"__complex__": [number.real, number.imag]}


def _disorder_values(
    count: int,
    strength: float,
    seed: int,
    distribution: str,
) -> np.ndarray:
    if not isinstance(strength, (int, float)) or strength < 0:
        raise ValueError("strength must be a nonnegative number.")
    if not isinstance(seed, int):
        raise ValueError("seed must be an integer.")
    rng = np.random.default_rng(seed)
    if distribution == "uniform":
        return rng.uniform(-strength, strength, size=count)
    if distribution == "normal":
        return rng.normal(0.0, strength, size=count)
    raise ValueError("distribution must be 'uniform' or 'normal'.")
