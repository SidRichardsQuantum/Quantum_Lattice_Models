from __future__ import annotations

from math import comb

import numpy as np
import pytest

from quantum_lattice_models import (
    ConservedQuantity,
    LatticeSpec,
    apply_interface,
    apply_spatial_potential,
    berry_curvature_convergence,
    commutator_diagnostic,
    create_graph_spin_spec,
    create_model_spec,
    estimate_model_dimension,
    evolve_time_dependent,
    export_lattice_dot,
    export_lattice_xyz,
    export_model_yaml,
    fermi_hubbard_basis,
    fermi_hubbard_chain,
    fermi_hubbard_chain_sector,
    fermi_hubbard_observables,
    gap_closing_points,
    graph_spin_subgraph,
    haldane_unit_cell,
    heisenberg_ladder_sector,
    heisenberg_ladder_sparse,
    import_lattice_xyz,
    import_model_yaml,
    load_model_plugins,
    sector_compatibility,
)
from quantum_lattice_models.lattice import Bond
from quantum_lattice_models.spin import SpinField, SpinInteraction


@pytest.mark.parametrize("periodic", [False, True])
def test_fermi_hubbard_sector_matches_full_space_block(periodic: bool) -> None:
    sector = fermi_hubbard_chain_sector(
        3, 2, 1, hopping=0.7, interaction=1.3, chemical_potential=0.2, periodic=periodic
    )
    full = fermi_hubbard_chain(
        3, hopping=0.7, interaction=1.3, chemical_potential=0.2, periodic=periodic
    )
    states = np.asarray(sector.basis.states)
    assert sector.basis.dimension == comb(3, 2) * comb(3, 1)
    assert np.allclose(sector.matrix.toarray(), full[np.ix_(states, states)])

    vector = np.ones(sector.basis.dimension) / np.sqrt(sector.basis.dimension)
    values = fermi_hubbard_observables(vector, sector.basis)
    assert values["up_occupation"].sum() == pytest.approx(2)
    assert values["down_occupation"].sum() == pytest.approx(1)


def test_fermi_sector_registry_spec_and_mapping() -> None:
    basis = fermi_hubbard_basis(3, 1, 2)
    reduced = np.arange(basis.dimension, dtype=complex)
    assert np.allclose(basis.project(basis.embed(reduced)), reduced)
    assert (
        estimate_model_dimension("fermi_hubbard_chain_sector_sparse", n_sites=4, n_up=2, n_down=1)
        == 24
    )
    spec = create_model_spec(
        "fermi_hubbard_chain_sector_sparse",
        parameters={"n_sites": 3, "n_up": 1, "n_down": 1},
    )
    result = spec.build_result()
    assert result.shape == (9, 9)
    assert result.metadata["sector"]["n_up"] == 1


def test_ladder_sector_and_commutator_diagnostics() -> None:
    sector = heisenberg_ladder_sector(3, 0, leg_coupling=0.8, rung_coupling=1.2)
    full = heisenberg_ladder_sparse(3, leg_coupling=0.8, rung_coupling=1.2)
    states = np.asarray(sector.basis.states)
    assert np.allclose(sector.matrix.toarray(), full[states][:, states].toarray())

    hamiltonian = np.diag([0.0, 1.0])
    conserved = ConservedQuantity("Z", np.diag([1.0, -1.0]))
    result = commutator_diagnostic(hamiltonian, conserved)
    assert result.metadata["conserved"]
    assert sector_compatibility(hamiltonian, conserved)
    assert not sector_compatibility(np.array([[0.0, 1.0], [1.0, 0.0]]), conserved)


def test_time_dependent_gap_and_topology_extensions() -> None:
    times = np.linspace(0, 1, 6)
    initial = np.array([1.0, 0.0])
    result = evolve_time_dependent(
        lambda time: time * np.array([[0.0, 1.0], [1.0, 0.0]]),
        initial,
        times,
        method="rk4",
        substeps=4,
    )
    assert np.allclose(result.values["state_norm"], 1.0)
    closings = gap_closing_points([-1, 0, 1], [1, 0, 1], threshold=0.1)
    assert closings.metadata["threshold_indices"] == [1]
    convergence = berry_curvature_convergence(haldane_unit_cell(t2=0.2, phi=np.pi / 2), [7, 11])
    assert convergence.values["chern_number"].shape == (2,)


def test_spatial_transformations_and_text_adapters(tmp_path) -> None:
    lattice = LatticeSpec(
        n_sites=3,
        positions=((-1.0, 0.0), (0.0, 0.0), (1.0, 0.0)),
        bonds=(Bond(0, 1), Bond(1, 2)),
        site_labels=("L", "C", "R"),
    )
    varied = apply_spatial_potential(lattice, lambda position: position[0] ** 2)
    interface = apply_interface(lattice, left=-1, right=1, width=2)
    assert varied.metadata["onsite_potential"] == [1.0, 0.0, 1.0]
    assert interface.metadata["onsite_potential"] == [-1.0, 0.0, 1.0]
    assert "graph lattice" in export_lattice_dot(lattice, tmp_path / "lattice.dot").read_text()
    xyz = export_lattice_xyz(lattice, tmp_path / "lattice.xyz")
    restored = import_lattice_xyz(xyz)
    assert restored.positions == ((-1.0, 0.0, 0.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.0))

    graph = create_graph_spin_spec(
        3,
        interactions=(
            SpinInteraction(0, 1, "Z", "Z", 1.0),
            SpinInteraction(1, 2, "X", "X", 0.5),
        ),
        fields=(SpinField(2, "Z", 0.2),),
        positions=lattice.positions,
        site_labels=lattice.site_labels,
    )
    reduced = graph_spin_subgraph(graph, (1, 2))
    assert reduced.parameters["n_sites"] == 2
    assert len(reduced.interactions) == 2


def test_yaml_model_adapter(tmp_path) -> None:
    pytest.importorskip("yaml")
    model = create_model_spec("ssh_model", parameters={"n_cells": 3})
    path = export_model_yaml(model, tmp_path / "model.yaml")
    assert import_model_yaml(path) == model


def test_plugin_discovery_reports_versioned_registration(monkeypatch) -> None:
    from quantum_lattice_models import registry

    calls = []

    class EntryPoint:
        name = "test-plugin"

        @staticmethod
        def load():
            def callback(register, *, api_version):
                calls.append((register, api_version))

            return callback

    monkeypatch.setattr(registry, "entry_points", lambda **kwargs: [EntryPoint()])
    assert load_model_plugins() == {"test-plugin": "loaded"}
    assert calls[0][1] == "1.0"
