from __future__ import annotations

import matplotlib
import numpy as np

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from quantum_lattice_models import create_model_spec
from quantum_lattice_models.models import (
    Lattice,
    custom_tight_binding,
    harper_hofstadter_square_lattice,
    tight_binding_chain,
)
from quantum_lattice_models.plotting import (
    plot_density,
    plot_hamiltonian_matrix,
    plot_hofstadter_butterfly,
    plot_interaction_graph,
    plot_lattice_graph,
    plot_lattice_state,
    plot_parameter_sweep,
    plot_site_probabilities,
    plot_spectrum,
)


def test_new_plotting_helpers_return_axes() -> None:
    H = tight_binding_chain(n_sites=4)
    ax = plot_spectrum(H, highlight_gap=True, zero_line=True)
    assert ax.get_xlabel() == "Eigenvalue index"
    assert len(ax.lines) >= 2
    plt.close(ax.figure)

    ax = plot_density(np.arange(5), bins=2)
    assert ax.get_ylabel() == "Count"
    plt.close(ax.figure)

    ax = plot_site_probabilities(np.ones(4) / 2)
    assert ax.get_ylabel() == "Probability"
    plt.close(ax.figure)

    positions = np.array([[0, 0], [1, 0], [2, 0], [3, 0]], dtype=float)
    ax = plot_lattice_graph(
        H,
        positions,
        sublattices=np.array([0, 1, 0, 1]),
        unit_cells=np.array([0, 0, 1, 1]),
        show_unit_cells=True,
    )
    assert ax.get_title() == "Lattice graph"
    assert len(ax.patches) == 2
    plt.close(ax.figure)


def test_metadata_aware_lattice_graph_and_state_plotters() -> None:
    lattice = Lattice(
        positions=[(0, 0), (1, 0), (0.5, 0.75)],
        bonds=[(0, 1), (1, 2, 0.5j), (2, 0)],
    )
    H = custom_tight_binding(lattice=lattice)

    ax = plot_lattice_graph(H, show_colorbar=True, labels=True)
    assert ax.get_title() == "Lattice graph"
    plt.close(ax.figure)

    ax = plot_lattice_state(H, np.ones(3, dtype=complex))
    assert ax.get_title() == "Lattice state"
    assert [tick.get_text() for tick in ax.figure.axes[1].get_yticklabels()] == [
        r"$-\pi$",
        "0",
        r"$\pi$",
    ]
    plt.close(ax.figure)


def test_hamiltonian_matrix_plotter_modes() -> None:
    H = tight_binding_chain(n_sites=4)
    for mode in ("real", "imag", "magnitude", "phase"):
        ax = plot_hamiltonian_matrix(H, mode=mode)
        assert ax.get_xlabel() == "Ket index"
        plt.close(ax.figure)


def test_hofstadter_butterfly_plotter() -> None:
    ax = plot_hofstadter_butterfly(
        lambda flux: harper_hofstadter_square_lattice(n_rows=2, n_cols=2, flux=flux),
        np.linspace(0.0, 0.5, 3),
    )
    assert ax.get_xlabel() == "Flux per plaquette"
    plt.close(ax.figure)


def test_parameter_sweep_plotter() -> None:
    ax = plot_parameter_sweep(
        lambda hopping: tight_binding_chain(n_sites=3, hopping=hopping),
        np.linspace(0.5, 1.0, 3),
        parameter_name="Hopping",
    )
    assert ax.get_xlabel() == "Hopping"
    plt.close(ax.figure)


def test_interaction_graph_consumes_portable_spin_and_hopping_terms() -> None:
    spin = create_model_spec(
        "xxz_chain",
        parameters={"n_sites": 3, "coupling": 1.0, "anisotropy": 0.5, "field": 0.2},
    )
    ax = plot_interaction_graph(spin)
    assert ax.get_title() == "Physical interaction graph"
    assert len(ax.lines) == 2
    assert any("XX=" in text.get_text() for text in ax.texts)
    assert any("Z=0.2" in text.get_text() for text in ax.texts)
    plt.close(ax.figure)

    ssh = create_model_spec("ssh_model", parameters={"n_cells": 2})
    ax = plot_interaction_graph(ssh, show_coefficients=False)
    assert len(ax.lines) == 3
    assert {text.get_text() for text in ax.texts} >= {"A0", "B0", "A1", "B1"}
    plt.close(ax.figure)


def test_interaction_graph_offsets_multiple_modes_on_one_site() -> None:
    fermi = create_model_spec("fermi_hubbard_chain", parameters={"n_sites": 2})
    ax = plot_interaction_graph(fermi, show_coefficients=False)
    offsets = ax.collections[0].get_offsets()

    assert offsets.shape == (4, 2)
    assert offsets[0, 1] != offsets[1, 1]
    assert offsets[2, 1] != offsets[3, 1]
    plt.close(ax.figure)

    bdg = create_model_spec("kitaev_chain_bdg", parameters={"n_sites": 2})
    ax = plot_interaction_graph(bdg, show_coefficients=False)
    assert ax.collections[0].get_offsets().shape == (4, 2)
    plt.close(ax.figure)
