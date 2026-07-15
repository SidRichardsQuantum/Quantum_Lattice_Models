from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pytest

from quantum_lattice_models import (
    MultiPanelPlotSpec,
    ReducedBasisMapping,
    anderson_square_lattice,
    anderson_square_lattice_sparse,
    checkerboard_chern_insulator,
    checkerboard_chern_insulator_sparse,
    create_graph_spin_spec,
    create_model_spec,
    dice_lattice,
    dice_lattice_sparse,
    fermi_hubbard_basis,
    fixed_magnetization_basis,
    fixed_particle_number_basis,
    graph_spin_model_sector,
    graph_spin_sector,
    graphene_lattice,
    graphene_lattice_sparse,
    lattice_annotation_data,
    matrix_block_structure,
    particle_model_subgraph,
    particle_model_vacancies,
    reduced_expectation,
    reduced_operator,
    spin_texture_plot_data,
    substitute_particle_bonds,
    to_styled_networkx,
)
from quantum_lattice_models.plotting import plot_spin_texture
from quantum_lattice_models.spin import SpinField, SpinInteraction, graph_spin_hamiltonian_sparse


def test_generic_reduced_mapping_shared_by_all_sector_bases() -> None:
    spin = fixed_magnetization_basis(4, 0)
    boson = fixed_particle_number_basis(3, 2)
    fermion = fermi_hubbard_basis(3, 1, 1)
    for mapping in (spin.mapping, boson.mapping, fermion.mapping):
        restored = ReducedBasisMapping.from_dict(mapping.to_dict())
        vector = np.arange(mapping.dimension, dtype=complex)
        assert restored == mapping
        assert np.allclose(mapping.project(mapping.embed(vector)), vector)

    full_operator = np.diag(np.arange(spin.mapping.full_dimension))
    projected = reduced_operator(full_operator, spin.mapping)
    state = np.ones(spin.dimension) / np.sqrt(spin.dimension)
    assert reduced_expectation(state, full_operator, spin.mapping) == pytest.approx(
        np.vdot(state, projected @ state)
    )


def test_graph_spin_sector_matches_full_block_and_rejects_leakage() -> None:
    interactions = (
        SpinInteraction(0, 1, "X", "X", 0.7),
        SpinInteraction(0, 1, "Y", "Y", 0.7),
        SpinInteraction(1, 2, "Z", "Z", 1.2),
    )
    fields = (SpinField(2, "Z", 0.3),)
    sector = graph_spin_sector(3, 1, interactions, fields)
    full = graph_spin_hamiltonian_sparse(3, interactions, fields)
    states = np.asarray(sector.basis.states)
    assert np.allclose(sector.matrix.toarray(), full[states][:, states].toarray())

    model = create_graph_spin_spec(3, interactions=interactions, fields=fields)
    assert np.allclose(graph_spin_model_sector(model, 1).matrix.toarray(), sector.matrix.toarray())
    with pytest.raises(ValueError, match="do not conserve"):
        graph_spin_sector(3, 1, fields=(SpinField(0, "X", 1.0),))


def test_particle_model_transformations_preserve_selected_matrix_blocks() -> None:
    model = create_model_spec("ssh_model", parameters={"n_cells": 3, "t1": 0.4, "t2": 1.1})
    full = model.hamiltonian()
    selected = (0, 1, 4, 5)
    transformed = particle_model_subgraph(model, selected)
    assert np.allclose(
        transformed.hamiltonian(),
        full[np.ix_(selected, selected)],
    )
    vacancy = particle_model_vacancies(model, (2, 3))
    assert np.allclose(vacancy.hamiltonian(), full[np.ix_(selected, selected)])

    bond = transformed.lattice.bonds[0]
    substituted = substitute_particle_bonds(
        transformed, {(bond.source, bond.target): 2 * complex(bond.value)}
    )
    assert substituted.lattice.bonds[0].value == 2 * bond.value


def test_visual_metadata_spin_texture_blocks_and_styled_graph() -> None:
    pytest.importorskip("networkx")
    model = create_graph_spin_spec(
        2,
        interactions=(SpinInteraction(0, 1, "Z", "Z", 1.0),),
        fields=(SpinField(0, "X", 0.2),),
        positions=((0.0, 0.0), (1.0, 0.0)),
    )
    texture = spin_texture_plot_data(model, np.array([[1, 0, 0.2], [0, 1, -0.3]]))
    assert texture["sites"][0]["vector"] == [1.0, 0.0, 0.2]
    ax = plot_spin_texture(model, np.array([[1, 0, 0.2], [0, 1, -0.3]]))
    assert ax.get_title() == "Local spin texture"
    plt.close(ax.figure)

    graph = to_styled_networkx(model)
    assert graph.edges[0, 1, 0]["color"] == "#D55E00"
    blocks = matrix_block_structure(np.eye(4), (("left", 0, 2), ("right", 2, 4)))
    assert blocks["blocks"][0]["frobenius_norm"] == pytest.approx(np.sqrt(2))
    panels = MultiPanelPlotSpec(({"kind": "spectrum"}, {"kind": "heatmap"}))
    panels.validate()
    assert panels.to_dict()["columns"] == 2
    assert lattice_annotation_data(model.lattice)["operations"] == []


def test_new_benchmark_models_dense_sparse_and_reference_properties() -> None:
    graphene = graphene_lattice(2, 3)
    assert np.allclose(graphene, graphene_lattice_sparse(2, 3).toarray())
    assert graphene.metadata["positions"].shape == (12, 2)
    graphene_values = np.linalg.eigvalsh(graphene)
    assert np.allclose(graphene_values, -graphene_values[::-1])

    anderson = anderson_square_lattice(3, 4, seed=7)
    assert np.allclose(anderson, anderson_square_lattice(3, 4, seed=7))
    assert np.allclose(anderson, anderson_square_lattice_sparse(3, 4, seed=7).toarray())

    checkerboard = checkerboard_chern_insulator(2, 3, mass=0.5)
    assert np.allclose(checkerboard, checkerboard_chern_insulator_sparse(2, 3, mass=0.5).toarray())
    assert np.allclose(checkerboard, checkerboard.conj().T)
    assert checkerboard.metadata["positions"].shape == (12, 2)

    dice = dice_lattice(2, 3)
    assert np.allclose(dice, dice_lattice_sparse(2, 3).toarray())
    assert np.count_nonzero(np.isclose(np.linalg.eigvalsh(dice), 0.0, atol=1e-10)) >= 6
    assert dice.metadata["positions"].shape == (18, 2)
