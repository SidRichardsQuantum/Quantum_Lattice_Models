from __future__ import annotations

import matplotlib
import numpy as np

matplotlib.use("Agg")

import matplotlib.pyplot as plt

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
    plot_lattice_graph,
    plot_lattice_spectrum,
    plot_lattice_state,
    plot_parameter_sweep,
    plot_site_probabilities,
)


def test_new_plotting_helpers_return_axes() -> None:
    H = tight_binding_chain(n_sites=4)
    ax = plot_lattice_spectrum(H)
    assert ax.get_xlabel() == "Eigenvalue index"
    plt.close(ax.figure)

    ax = plot_density(np.arange(5), bins=2)
    assert ax.get_ylabel() == "Count"
    plt.close(ax.figure)

    ax = plot_site_probabilities(np.ones(4) / 2)
    assert ax.get_ylabel() == "Probability"
    plt.close(ax.figure)

    positions = np.array([[0, 0], [1, 0], [2, 0], [3, 0]], dtype=float)
    ax = plot_lattice_graph(H, positions)
    assert ax.get_title() == "Lattice graph"
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
