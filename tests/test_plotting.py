from __future__ import annotations

import matplotlib
import numpy as np

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from quantum_lattice_models.models import harper_hofstadter_square_lattice, tight_binding_chain
from quantum_lattice_models.plotting import (
    plot_density,
    plot_hofstadter_butterfly,
    plot_lattice_graph,
    plot_lattice_spectrum,
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


def test_hofstadter_butterfly_plotter() -> None:
    ax = plot_hofstadter_butterfly(
        lambda flux: harper_hofstadter_square_lattice(n_rows=2, n_cols=2, flux=flux),
        np.linspace(0.0, 0.5, 3),
    )
    assert ax.get_xlabel() == "Flux per plaquette"
    plt.close(ax.figure)
