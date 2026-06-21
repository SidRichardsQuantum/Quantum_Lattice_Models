"""Immutable transformations for portable finite-lattice specifications."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from dataclasses import replace

import numpy as np

from quantum_lattice_models.lattice import Bond
from quantum_lattice_models.periodic import PeriodicLatticeSpec
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


def repeat_unit_cell(
    periodic: PeriodicLatticeSpec,
    repeats: int | Sequence[int],
    *,
    periodic_boundaries: Sequence[bool] | None = None,
) -> LatticeSpec:
    """Expand a periodic unit cell into a finite portable lattice."""

    periodic.validate()
    shape = (
        (repeats,) * periodic.dimension
        if isinstance(repeats, int)
        else tuple(int(value) for value in repeats)
    )
    if len(shape) != periodic.dimension or any(value < 1 for value in shape):
        raise ValueError("repeats must contain one positive count per lattice dimension.")
    boundaries = (
        (False,) * periodic.dimension
        if periodic_boundaries is None
        else tuple(bool(value) for value in periodic_boundaries)
    )
    if len(boundaries) != periodic.dimension:
        raise ValueError("periodic_boundaries must match the lattice dimension.")

    vectors = np.asarray(periodic.primitive_vectors)
    offsets = np.asarray(periodic.orbital_positions)
    cells = tuple(np.ndindex(shape))
    cell_index = {cell: index for index, cell in enumerate(cells)}
    positions = tuple(
        tuple(np.asarray(cell) @ vectors + offsets[orbital])
        for cell in cells
        for orbital in range(periodic.n_orbitals)
    )
    bonds: list[Bond] = []
    for cell in cells:
        for bond in periodic.bonds:
            target_cell = tuple(
                cell[axis] + bond.displacement[axis] for axis in range(periodic.dimension)
            )
            valid = True
            wrapped = list(target_cell)
            for axis, count in enumerate(shape):
                if 0 <= wrapped[axis] < count:
                    continue
                if boundaries[axis]:
                    wrapped[axis] %= count
                else:
                    valid = False
                    break
            if not valid:
                continue
            source = cell_index[cell] * periodic.n_orbitals + bond.source
            target = cell_index[tuple(wrapped)] * periodic.n_orbitals + bond.target
            bonds.append(Bond(source, target, bond.value))
    labels = periodic.orbital_labels or tuple(
        f"orbital {index}" for index in range(periodic.n_orbitals)
    )
    return LatticeSpec(
        n_sites=len(cells) * periodic.n_orbitals,
        positions=positions,
        bonds=tuple(bonds),
        site_labels=tuple(
            f"{labels[orbital]}@{','.join(map(str, cell))}"
            for cell in cells
            for orbital in range(periodic.n_orbitals)
        ),
        orbital_labels=tuple(labels[orbital] for _ in cells for orbital in range(len(labels))),
        sublattice_labels=(
            tuple(
                periodic.sublattice_labels[orbital]
                for _ in cells
                for orbital in range(periodic.n_orbitals)
            )
            if periodic.sublattice_labels
            else ()
        ),
        unit_cells=tuple(cell_index[cell] for cell in cells for _ in range(periodic.n_orbitals)),
        boundary_conditions={
            f"axis_{axis}": "periodic" if boundaries[axis] else "open"
            for axis in range(periodic.dimension)
        },
        units=dict(periodic.units),
        conventions={
            **periodic.conventions,
            "source": "repeated periodic unit cell",
        },
        references=periodic.references,
        provenance=periodic.provenance
        + (
            {
                "operation": "repeat_unit_cell",
                "parameters": {
                    "repeats": list(shape),
                    "periodic_boundaries": list(boundaries),
                },
            },
        ),
        metadata={**periodic.metadata, "onsite_potential": list(periodic.onsite)},
    )


def apply_impurity_potential(
    lattice: LatticeSpec,
    potentials: dict[int, complex],
) -> LatticeSpec:
    """Add site-selective impurity potentials to portable onsite metadata."""

    lattice.validate()
    if any(not isinstance(site, int) or not 0 <= site < lattice.n_sites for site in potentials):
        raise ValueError("Impurity sites must be valid lattice indices.")
    onsite = np.asarray(lattice.metadata.get("onsite_potential", [0.0] * lattice.n_sites), complex)
    if onsite.shape != (lattice.n_sites,):
        raise ValueError("Existing onsite_potential metadata must contain one value per site.")
    for site, value in potentials.items():
        onsite[site] += complex(value)
    metadata = dict(lattice.metadata)
    metadata["onsite_potential"] = [_portable_number(value) for value in onsite]
    return _transformed(
        lattice,
        "impurity_potential",
        {
            "potentials": {
                str(site): _portable_complex(value) for site, value in sorted(potentials.items())
            }
        },
        metadata=metadata,
    )


def apply_domain_wall(
    lattice: LatticeSpec,
    *,
    axis: int = 0,
    position: float = 0.0,
    left: complex = -1.0,
    right: complex = 1.0,
) -> LatticeSpec:
    """Add a step-like onsite potential across a spatial domain wall."""

    lattice.validate()
    if not lattice.positions:
        raise ValueError("Domain walls require lattice positions.")
    if not 0 <= axis < len(lattice.positions[0]):
        raise ValueError("axis is out of range for lattice positions.")
    values = {
        site: left if coordinate[axis] < position else right
        for site, coordinate in enumerate(lattice.positions)
    }
    transformed = apply_impurity_potential(lattice, values)
    provenance = lattice.provenance + (
        {
            "operation": "domain_wall",
            "parameters": {
                "axis": axis,
                "position": position,
                "left": _portable_complex(left),
                "right": _portable_complex(right),
            },
        },
    )
    return replace(transformed, provenance=provenance)


def apply_twisted_boundary(
    lattice: LatticeSpec,
    angle: float,
    *,
    axis: int = 0,
    crossing_bonds: Iterable[tuple[int, int]] | None = None,
    default_value: complex = -1.0,
) -> LatticeSpec:
    """Apply a phase twist to bonds crossing a selected finite-lattice boundary."""

    lattice.validate()
    selected = set(crossing_bonds or ())
    if crossing_bonds is None:
        if not lattice.positions:
            raise ValueError("Automatic twisted-boundary detection requires positions.")
        coordinates = np.asarray(lattice.positions)
        if not 0 <= axis < coordinates.shape[1]:
            raise ValueError("axis is out of range for lattice positions.")
        span = float(np.ptp(coordinates[:, axis]))
        selected = {
            (bond.source, bond.target)
            for bond in lattice.bonds
            if span > 0
            and abs(coordinates[bond.target, axis] - coordinates[bond.source, axis]) > span / 2
        }
    phase = np.exp(1j * float(angle))
    bonds = tuple(
        (
            Bond(
                bond.source,
                bond.target,
                complex(default_value if bond.value is None else bond.value) * phase,
            )
            if (bond.source, bond.target) in selected
            else bond
        )
        for bond in lattice.bonds
    )
    return _transformed(
        lattice,
        "twisted_boundary",
        {"angle": float(angle), "axis": axis, "bonds": [list(pair) for pair in sorted(selected)]},
        bonds=bonds,
    )


def add_long_range_bonds(
    lattice: LatticeSpec,
    *,
    max_distance: float,
    strength: complex = -1.0,
    power: float = 0.0,
    min_distance: float = 0.0,
    include_existing: bool = False,
) -> LatticeSpec:
    """Generate directed bonds by Euclidean distance with optional power-law decay."""

    lattice.validate()
    if not lattice.positions:
        raise ValueError("Long-range bond generation requires lattice positions.")
    if max_distance <= 0 or min_distance < 0 or min_distance >= max_distance:
        raise ValueError("Require 0 <= min_distance < max_distance.")
    if power < 0:
        raise ValueError("power must be nonnegative.")
    coordinates = np.asarray(lattice.positions)
    existing = {(bond.source, bond.target) for bond in lattice.bonds}
    additions = []
    for source in range(lattice.n_sites):
        for target in range(source + 1, lattice.n_sites):
            distance = float(np.linalg.norm(coordinates[target] - coordinates[source]))
            if not min_distance < distance <= max_distance:
                continue
            if not include_existing and (
                (source, target) in existing or (target, source) in existing
            ):
                continue
            value = complex(strength) / (distance**power if power else 1.0)
            additions.append(Bond(source, target, value))
    return _transformed(
        lattice,
        "long_range_bonds",
        {
            "max_distance": max_distance,
            "min_distance": min_distance,
            "strength": _portable_complex(strength),
            "power": power,
        },
        bonds=lattice.bonds + tuple(additions),
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


def _portable_number(value: complex) -> object:
    number = complex(value)
    if abs(number.imag) < 1e-15:
        return float(number.real)
    return _portable_complex(number)


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
