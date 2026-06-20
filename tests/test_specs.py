from __future__ import annotations

import json

import numpy as np
import pytest
import scipy.sparse as sp

from quantum_lattice_models import (
    BasisIndexMapping,
    HamiltonianResult,
    InteractionTerm,
    LatticeSpec,
    LocalDegreeOfFreedom,
    ModelSpec,
    create_graph_spin_spec,
    create_model_spec,
    load_model,
    migrate_spec_data,
)
from quantum_lattice_models.cli import main
from quantum_lattice_models.lattice import Bond
from quantum_lattice_models.spin import SpinField, SpinInteraction


def test_lattice_spec_round_trip_preserves_complex_bonds() -> None:
    original = LatticeSpec(
        n_sites=3,
        positions=((0.0, 0.0), (1.0, 0.0), (0.5, 0.8)),
        bonds=(Bond(0, 1), Bond(1, 2, 0.25j)),
        sublattice_labels=("A", "B", "A"),
        unit_cells=(0, 0, 1),
        boundary_conditions={"x": "open"},
        units={"position": "lattice_constant"},
        conventions={"bond_sign": "explicit matrix element"},
        references=("doi:10.0000/example",),
        provenance=({"operation": "created", "parameters": {}},),
        metadata={"label": "triangle fragment"},
    )

    restored = LatticeSpec.from_dict(original.to_dict())

    assert restored == original
    assert restored.to_lattice().bonds[1].value == 0.25j


def test_model_spec_json_round_trip_and_dense_sparse_construction(tmp_path) -> None:
    spec = create_model_spec(
        "tight_binding_chain",
        parameters={"n_sites": 5, "hopping": 0.7, "periodic": True},
        units={"hopping": "eV"},
        conventions={"hopping_sign": "H_ij=-t"},
        references=("doi:10.0000/example",),
        provenance={"source": "test"},
        metadata={"purpose": "round-trip"},
    )
    path = spec.save(tmp_path / "chain.json")
    restored = load_model(path)

    dense = restored.hamiltonian()
    sparse = restored.hamiltonian(sparse=True)

    assert restored == spec
    assert dense.shape == (5, 5)
    assert sp.issparse(sparse)
    assert np.allclose(sparse.toarray(), dense)
    assert json.loads(path.read_text())["schema_version"] == "1.0"
    assert restored.parameters["onsite"] == 0.0
    assert restored.units["hopping"] == "eV"
    assert restored.conventions["hopping_sign"] == "H_ij=-t"

    onsite_spec = create_model_spec(
        "tight_binding_chain",
        parameters={"n_sites": 3, "onsite": (0.1, 0.2, 0.3)},
    )
    assert np.allclose(np.diag(onsite_spec.hamiltonian()), [0.1, 0.2, 0.3])

    kitaev = create_model_spec(
        "kitaev_chain_bdg",
        parameters={"n_sites": 3, "pairing": 0.25j},
    )
    kitaev_path = kitaev.save(tmp_path / "kitaev.json")
    assert load_model(kitaev_path).parameters["pairing"] == 0.25j


def test_custom_lattice_model_spec_builds_dense_and_sparse() -> None:
    lattice = LatticeSpec(
        n_sites=3,
        bonds=(Bond(0, 1), Bond(1, 2, 0.5j)),
    )
    spec = create_model_spec(
        "custom_tight_binding",
        parameters={"hopping": 0.8, "hermitian": True},
        lattice=lattice,
    )

    dense = spec.hamiltonian()
    sparse = spec.hamiltonian(sparse=True)

    assert dense[0, 1] == -0.8
    assert dense[1, 2] == 0.5j
    assert np.allclose(sparse.toarray(), dense)


def test_model_spec_validation_rejects_invalid_schema_and_representation() -> None:
    with pytest.raises(ValueError, match="schema_version"):
        ModelSpec(family="ssh_model", schema_version="2.0").validate()
    with pytest.raises(ValueError, match="representation"):
        ModelSpec(family="ssh_model", representation="gpu").validate()
    ModelSpec(
        family="transverse_field_ising",
        parameters={"n_sites": 2},
        representation="sparse",
    ).validate()
    with pytest.raises(ValueError, match="must have type int"):
        ModelSpec(
            family="ssh_model",
            parameters={"n_cells": "four"},
        ).validate()
    with pytest.raises(ValueError, match="registered basis"):
        ModelSpec(
            family="ssh_model",
            parameters={"n_cells": 2},
            basis="qubit",
        ).validate()


def test_create_inspect_and_validate_cli(tmp_path, capsys) -> None:
    path = tmp_path / "ssh.json"
    assert (
        main(
            [
                "create",
                "ssh_model",
                "--n-cells",
                "4",
                "--t1",
                "0.3",
                "--representation",
                "dense",
                "--output",
                str(path),
            ]
        )
        == 0
    )
    assert path.exists()
    capsys.readouterr()

    assert main(["inspect", str(path)]) == 0
    inspection = json.loads(capsys.readouterr().out)
    assert inspection["family"] == "ssh_model"
    assert inspection["dimension"] == 8

    assert main(["validate", str(path)]) == 0
    assert capsys.readouterr().out.startswith("valid\tssh_model\t1.0")


def test_schema_migration_and_strict_file_validation() -> None:
    legacy = {
        "family": "ssh_model",
        "parameters": {"n_cells": 2},
        "basis": "single particle",
    }
    migrated = migrate_spec_data(legacy, kind="model")
    assert migrated["schema_version"] == "1.0"
    assert ModelSpec.from_dict(legacy).schema_version == "1.0"

    with pytest.raises(ValueError, match="no migration"):
        ModelSpec.from_dict({**legacy, "schema_version": "2.0"})
    with pytest.raises(ValueError, match="unknown fields"):
        ModelSpec.from_dict({**legacy, "unexpected": True})
    with pytest.raises(ValueError, match="entries require 'target'"):
        LatticeSpec.from_dict({"n_sites": 2, "bonds": [{"source": 0}]})
    with pytest.raises(ValueError, match="finite coordinates"):
        LatticeSpec(n_sites=1, positions=((float("nan"), 0.0),)).validate()
    with pytest.raises(ValueError, match="conventions"):
        ModelSpec(family="ssh_model", conventions={"gauge": 2}).validate()
    with pytest.raises(ValueError, match="provenance"):
        LatticeSpec(n_sites=1, provenance=("invalid",)).validate()


def test_model_spec_build_result_preserves_model_and_representation() -> None:
    spec = create_model_spec(
        "tight_binding_chain",
        parameters={"n_sites": 4, "periodic": True},
        representation="sparse",
    )
    result = spec.build_result()

    assert isinstance(result, HamiltonianResult)
    assert sp.issparse(result.matrix)
    assert result.model == spec
    assert result.basis == "single particle"
    assert result.representation == "sparse"


@pytest.mark.parametrize(
    ("family", "parameters", "degree_kind", "mapping_role"),
    [
        ("transverse_field_ising", {"n_sites": 3}, "spin", "tensor_factor"),
        ("heisenberg_chain", {"n_sites": 3}, "spin", "tensor_factor"),
        ("xxz_chain", {"n_sites": 3}, "spin", "tensor_factor"),
        ("ssh_model", {"n_cells": 2}, "orbital", "single_particle_state"),
    ],
)
def test_supported_models_include_portable_physical_system(
    family, parameters, degree_kind, mapping_role
) -> None:
    spec = create_model_spec(family, parameters=parameters)
    restored = ModelSpec.from_dict(spec.to_dict())

    assert restored == spec
    assert spec.lattice is not None
    assert all(degree.kind == degree_kind for degree in spec.local_degrees)
    assert all(mapping.role == mapping_role for mapping in spec.basis_mappings)
    assert spec.interactions
    assert spec.summary()["local_degree_count"] == len(spec.local_degrees)


def test_spin_physical_terms_preserve_axes_coefficients_and_geometry() -> None:
    spec = create_model_spec(
        "xxz_chain",
        parameters={
            "n_sites": 3,
            "coupling": 2.0,
            "anisotropy": 0.25,
            "field": 0.3,
        },
    )

    pair_terms = [term for term in spec.interactions if len(term.degrees) == 2]
    onsite_terms = [term for term in spec.interactions if len(term.degrees) == 1]

    assert spec.lattice.positions == ((0.0, 0.0), (1.0, 0.0), (2.0, 0.0))
    assert {term.operators for term in pair_terms} == {("X", "X"), ("Y", "Y"), ("Z", "Z")}
    assert {term.coefficient for term in pair_terms} == {2.0, 0.5}
    assert all(term.operators == ("Z",) and term.coefficient == 0.3 for term in onsite_terms)


def test_ssh_and_custom_tight_binding_physical_terms() -> None:
    ssh = create_model_spec(
        "ssh_model",
        parameters={"n_cells": 2, "t1": 0.4, "t2": 1.2},
    )
    assert ssh.lattice.sublattice_labels == ("A", "B", "A", "B")
    assert [term.coefficient for term in ssh.interactions] == [-0.4, -1.2, -0.4]
    assert all(term.metadata["hermitian_conjugate"] is True for term in ssh.interactions)

    lattice = LatticeSpec(
        n_sites=2,
        positions=((0.0, 0.0), (1.0, 0.0)),
        bonds=(Bond(0, 1, 0.25j),),
        site_labels=("left", "right"),
        orbital_labels=("s", "p"),
    )
    custom = create_model_spec(
        "custom_tight_binding",
        lattice=lattice,
        parameters={"onsite": (0.2, -0.1), "hermitian": True},
    )
    assert [degree.orbital for degree in custom.local_degrees] == ["s", "p"]
    assert custom.interactions[0].coefficient == 0.25j
    assert {term.kind for term in custom.interactions[1:]} == {"onsite"}


def test_physical_system_validation_rejects_invalid_indices_and_operators() -> None:
    degrees = (
        LocalDegreeOfFreedom(0, 0, "spin", 2, "s0"),
        LocalDegreeOfFreedom(1, 1, "spin", 2, "s1"),
    )
    with pytest.raises(ValueError, match="incompatible"):
        ModelSpec(
            family="transverse_field_ising",
            parameters={"n_sites": 2},
            basis="qubit",
            local_degrees=degrees,
            basis_mappings=(
                BasisIndexMapping(0, 0, "tensor_factor"),
                BasisIndexMapping(1, 1, "tensor_factor"),
            ),
            interactions=(InteractionTerm((0,), ("number",), 1.0, "invalid"),),
        ).validate()


@pytest.mark.parametrize(
    ("family", "parameters", "kinds", "mapping_role"),
    [
        (
            "bose_hubbard_chain",
            {"n_sites": 2, "max_occupancy": 3},
            {"boson"},
            "tensor_factor",
        ),
        (
            "fermi_hubbard_chain",
            {"n_sites": 2},
            {"fermion"},
            "mode",
        ),
        (
            "kitaev_chain_bdg",
            {"n_sites": 2, "pairing": 0.4},
            {"nambu"},
            "single_particle_state",
        ),
    ],
)
def test_hubbard_and_bdg_specs_include_physical_records(
    family, parameters, kinds, mapping_role
) -> None:
    spec = create_model_spec(family, parameters=parameters)
    restored = ModelSpec.from_dict(spec.to_dict())

    assert restored == spec
    assert {degree.kind for degree in spec.local_degrees} == kinds
    assert all(mapping.role == mapping_role for mapping in spec.basis_mappings)
    assert spec.lattice is not None
    assert spec.interactions


def test_hubbard_physical_terms_match_basis_conventions() -> None:
    bose = create_model_spec(
        "bose_hubbard_chain",
        parameters={
            "n_sites": 2,
            "hopping": 0.7,
            "interaction": 2.0,
            "chemical_potential": 0.25,
            "max_occupancy": 3,
        },
    )
    assert [degree.local_dimension for degree in bose.local_degrees] == [4, 4]
    assert {term.operators for term in bose.interactions} >= {
        ("create", "annihilate"),
        ("number_pair",),
        ("number",),
    }
    assert any(
        term.kind == "onsite_interaction" and term.coefficient == 1.0 for term in bose.interactions
    )

    fermi = create_model_spec(
        "fermi_hubbard_chain",
        parameters={
            "n_sites": 2,
            "hopping": 0.6,
            "interaction": 3.0,
            "chemical_potential": 0.2,
        },
    )
    assert [degree.component for degree in fermi.local_degrees] == [
        "up",
        "down",
        "up",
        "down",
    ]
    onsite = [term for term in fermi.interactions if term.kind == "onsite_interaction"]
    assert [term.degrees for term in onsite] == [(0, 1), (2, 3)]
    assert all(term.operators == ("number", "number") for term in onsite)


def test_kitaev_bdg_records_particle_hole_order_and_pairing() -> None:
    spec = create_model_spec(
        "kitaev_chain_bdg",
        parameters={
            "n_sites": 3,
            "hopping": 0.8,
            "chemical_potential": 0.2,
            "pairing": 0.5j,
        },
    )

    assert [degree.component for degree in spec.local_degrees] == [
        "particle",
        "particle",
        "particle",
        "hole",
        "hole",
        "hole",
    ]
    assert [mapping.basis_index for mapping in spec.basis_mappings] == list(range(6))
    pairing = [term for term in spec.interactions if term.kind == "bdg_pairing"]
    assert {term.coefficient for term in pairing} == {0.5j, -0.5j}
    assert all(term.operators == ("particle", "hole") for term in pairing)


def test_portable_graph_spin_spec_round_trip_and_construction() -> None:
    interactions = (
        SpinInteraction(0, 1, "X", "Y", 0.25),
        SpinInteraction(1, 2, "Z", "Z", -0.7),
    )
    fields = (SpinField(2, "X", 0.4),)
    spec = create_graph_spin_spec(
        3,
        interactions=interactions,
        fields=fields,
        positions=((0.0, 0.0), (1.0, 0.0), (0.5, 0.8)),
        site_labels=("left", "right", "top"),
    )
    restored = ModelSpec.from_dict(spec.to_dict())

    assert restored == spec
    assert restored.family == "graph_spin"
    assert restored.summary()["dimension"] == 8
    assert np.allclose(restored.hamiltonian(), spec.hamiltonian())
    assert np.allclose(restored.hamiltonian(sparse=True).toarray(), spec.hamiltonian())
    assert [degree.label for degree in spec.local_degrees] == ["left", "right", "top"]
    with pytest.raises(ValueError, match="contiguous"):
        ModelSpec(
            family="transverse_field_ising",
            parameters={"n_sites": 2},
            basis="qubit",
            local_degrees=(LocalDegreeOfFreedom(1, 0, "spin", 2),),
        ).validate()
