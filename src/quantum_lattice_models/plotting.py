"""Matplotlib helpers for spectra and SSH edge states."""

from __future__ import annotations

import numpy as np

from quantum_lattice_models.models import ssh_edge_state_localization
from quantum_lattice_models.spectra import density_of_states, eigenvalues


def plot_spectrum(H: np.ndarray, ax=None, **scatter_kwargs):
    """Plot sorted eigenvalues and return the Matplotlib axes."""

    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    values = np.real_if_close(eigenvalues(H)).real
    kwargs = {"s": 24, "color": "tab:blue"}
    kwargs.update(scatter_kwargs)
    ax.scatter(np.arange(values.size), np.sort(values), **kwargs)
    ax.set_xlabel("Eigenvalue index")
    ax.set_ylabel("Energy")
    ax.set_title("Spectrum")
    return ax


def plot_density_of_states(H: np.ndarray, bins: int = 50, ax=None, **bar_kwargs):
    """Plot a histogram approximation to the density of states."""

    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    counts, edges = density_of_states(H, bins=bins)
    widths = np.diff(edges)
    kwargs = {
        "align": "edge",
        "color": "tab:green",
        "edgecolor": "black",
        "alpha": 0.75,
    }
    kwargs.update(bar_kwargs)
    ax.bar(edges[:-1], counts, width=widths, **kwargs)
    ax.set_xlabel("Energy")
    ax.set_ylabel("Count")
    ax.set_title("Density of states")
    return ax


def plot_ssh_edge_state(
    eigenvector: np.ndarray,
    n_cells: int,
    ax=None,
    *,
    edge_cells: int = 1,
    **line_kwargs,
):
    """Plot site probability for an SSH eigenvector."""

    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    vector = np.asarray(eigenvector, dtype=complex).reshape(-1)
    if vector.size != 2 * n_cells:
        raise ValueError("eigenvector length must equal 2 * n_cells.")

    probability = np.abs(vector) ** 2
    probability = probability / probability.sum()
    localization = ssh_edge_state_localization(vector, n_cells, edge_cells=edge_cells)

    kwargs = {"marker": "o", "color": "tab:red"}
    kwargs.update(line_kwargs)
    ax.plot(np.arange(vector.size), probability, **kwargs)
    ax.set_xlabel("SSH site")
    ax.set_ylabel("Probability")
    ax.set_title(f"SSH edge-state weight: {localization:.3f}")
    return ax
