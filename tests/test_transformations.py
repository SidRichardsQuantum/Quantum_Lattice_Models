from __future__ import annotations

import numpy as np
import pytest

from quantum_lattice_models import (
    LatticeSpec,
    add_lattice_bonds,
    add_long_range_bonds,
    apply_bond_disorder,
    apply_domain_wall,
    apply_impurity_potential,
    apply_onsite_disorder,
    apply_twisted_boundary,
    lattice_subgraph,
    relabel_lattice,
    remove_lattice_bonds,
    remove_lattice_sites,
    with_boundary_conditions,
)
from quantum_lattice_models.lattice import Bond


def _chain() -> LatticeSpec:
    return LatticeSpec(
        n_sites=4,
        positions=((0.0, 0.0), (1.0, 0.0), (2.0, 0.0), (3.0, 0.0)),
        bonds=(Bond(0, 1, -1.0), Bond(1, 2, -1.0), Bond(2, 3, -1.0)),
        site_labels=("a", "b", "c", "d"),
        unit_cells=(0, 0, 1, 1),
        metadata={"name": "chain"},
    )


def test_relabel_subgraph_and_vacancy_preserve_site_metadata() -> None:
    lattice = _chain()
    relabeled = relabel_lattice(lattice, [3, 2, 1, 0])
    subgraph = lattice_subgraph(lattice, [0, 2, 3])
    vacancy = remove_lattice_sites(lattice, [1])

    assert relabeled.site_labels == ("d", "c", "b", "a")
    assert relabeled.bonds[0] == Bond(3, 2, -1.0)
    assert subgraph.site_labels == ("a", "c", "d")
    assert subgraph.bonds == (Bond(1, 2, -1.0),)
    assert vacancy.n_sites == subgraph.n_sites
    assert vacancy.positions == subgraph.positions
    assert vacancy.bonds == subgraph.bonds
    assert vacancy.site_labels == subgraph.site_labels
    assert vacancy.provenance[-1]["operation"] == "remove_sites"
    assert lattice.provenance == ()


def test_bond_and_boundary_transformations_are_immutable() -> None:
    lattice = _chain()
    added = add_lattice_bonds(lattice, [(3, 0, 0.5j)])
    removed = remove_lattice_bonds(added, [(1, 2)])
    periodic = with_boundary_conditions(removed, {"x": "periodic"})

    assert len(lattice.bonds) == 3
    assert added.bonds[-1] == Bond(3, 0, 0.5j)
    assert (1, 2) not in {(bond.source, bond.target) for bond in removed.bonds}
    assert periodic.boundary_conditions == {"x": "periodic"}
    assert [record["operation"] for record in periodic.provenance] == [
        "add_bonds",
        "remove_bonds",
        "set_boundary_conditions",
    ]


def test_disorder_is_reproducible_and_records_provenance() -> None:
    lattice = _chain()
    onsite_a = apply_onsite_disorder(lattice, 0.4, seed=12)
    onsite_b = apply_onsite_disorder(lattice, 0.4, seed=12)
    bonds_a = apply_bond_disorder(lattice, 0.2, seed=9, distribution="normal")
    bonds_b = apply_bond_disorder(lattice, 0.2, seed=9, distribution="normal")

    assert onsite_a == onsite_b
    assert np.max(np.abs(onsite_a.metadata["onsite_potential"])) <= 0.4
    assert bonds_a == bonds_b
    assert bonds_a.bonds != lattice.bonds
    assert onsite_a.provenance[-1]["parameters"]["seed"] == 12


def test_transformation_validation_failures() -> None:
    lattice = _chain()
    with pytest.raises(ValueError, match="permutation"):
        relabel_lattice(lattice, [0, 1])
    with pytest.raises(ValueError, match="retain at least one"):
        remove_lattice_sites(lattice, range(4))
    with pytest.raises(ValueError, match="distribution"):
        apply_onsite_disorder(lattice, 1.0, seed=1, distribution="invalid")


def test_impurity_domain_wall_twist_and_long_range_transformations() -> None:
    lattice = _chain()
    impurity = apply_impurity_potential(lattice, {1: 0.7, 3: -0.2})
    wall = apply_domain_wall(lattice, position=1.5, left=-0.4, right=0.4)
    ring = add_lattice_bonds(lattice, [(3, 0, -1.0)])
    twisted = apply_twisted_boundary(ring, np.pi / 2)
    long_range = add_long_range_bonds(lattice, min_distance=1.1, max_distance=2.1)

    assert impurity.metadata["onsite_potential"] == [0.0, 0.7, 0.0, -0.2]
    assert wall.metadata["onsite_potential"] == [-0.4, -0.4, 0.4, 0.4]
    assert twisted.bonds[-1].value == pytest.approx(-1j)
    assert {(bond.source, bond.target) for bond in long_range.bonds} >= {(0, 2), (1, 3)}
