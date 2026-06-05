"""Matplotlib helpers for spectra, densities, probabilities, and lattice graphs."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models.models import ssh_edge_state_localization
from quantum_lattice_models.spectra import density_of_states, eigenvalues


def plot_spectrum(H: np.ndarray, ax=None, **scatter_kwargs):
    """Plot sorted eigenvalues and return the Matplotlib axes."""

    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    values = np.real_if_close(eigenvalues(_as_dense(H))).real
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
    counts, edges = density_of_states(_as_dense(H), bins=bins)
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


def plot_lattice_spectrum(H: np.ndarray, ax=None, **scatter_kwargs):
    """Plot a lattice Hamiltonian spectrum and return the Matplotlib axes."""

    return plot_spectrum(H, ax=ax, **scatter_kwargs)


def plot_density(values: np.ndarray, ax=None, *, bins: int = 50, **bar_kwargs):
    """Plot a histogram for a one-dimensional array of values."""

    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    data = np.asarray(values, dtype=float).reshape(-1)
    counts, edges = np.histogram(data, bins=bins)
    kwargs = {"align": "edge", "color": "tab:purple", "edgecolor": "black", "alpha": 0.75}
    kwargs.update(bar_kwargs)
    ax.bar(edges[:-1], counts, width=np.diff(edges), **kwargs)
    ax.set_xlabel("Value")
    ax.set_ylabel("Count")
    ax.set_title("Density")
    return ax


def plot_site_probabilities(
    state: np.ndarray,
    ax=None,
    *,
    normalize: bool = True,
    title: str = "Site probabilities",
    **line_kwargs,
):
    """Plot site probabilities for a state vector or single-particle eigenvector."""

    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    vector = np.asarray(state, dtype=complex).reshape(-1)
    probabilities = np.abs(vector) ** 2
    norm = probabilities.sum()
    if normalize:
        if norm == 0:
            raise ValueError("state must have nonzero norm when normalize=True.")
        probabilities = probabilities / norm
    kwargs = {"marker": "o", "color": "tab:red"}
    kwargs.update(line_kwargs)
    ax.plot(np.arange(probabilities.size), probabilities, **kwargs)
    ax.set_xlabel("Site")
    ax.set_ylabel("Probability")
    ax.set_title(title)
    return ax


def plot_hofstadter_butterfly(
    builder,
    flux_values: np.ndarray,
    ax=None,
    **scatter_kwargs,
):
    """Plot eigenvalues over flux values for a Harper-Hofstadter-style builder.

    ``builder`` must accept a single flux value and return a Hamiltonian matrix.
    """

    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    fluxes = np.asarray(flux_values, dtype=float).reshape(-1)
    kwargs = {"s": 8, "color": "tab:blue", "alpha": 0.65}
    kwargs.update(scatter_kwargs)
    for flux in fluxes:
        values = np.real_if_close(eigenvalues(_as_dense(builder(float(flux))))).real
        ax.scatter(np.full(values.size, flux), values, **kwargs)
    ax.set_xlabel("Flux per plaquette")
    ax.set_ylabel("Energy")
    ax.set_title("Hofstadter spectrum")
    return ax


def plot_lattice_graph(
    H: np.ndarray,
    positions: np.ndarray,
    ax=None,
    *,
    threshold: float = 1e-12,
    node_size: float = 34,
    **line_kwargs,
):
    """Plot graph connectivity implied by nonzero off-diagonal matrix entries."""

    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    matrix = _as_dense(H)
    coords = np.asarray(positions, dtype=float)
    if coords.shape != (matrix.shape[0], 2):
        raise ValueError("positions must have shape (n_sites, 2).")

    kwargs = {"color": "0.55", "linewidth": 1.0, "alpha": 0.75}
    kwargs.update(line_kwargs)
    for i in range(matrix.shape[0]):
        for j in range(i + 1, matrix.shape[1]):
            if abs(matrix[i, j]) > threshold:
                ax.plot([coords[i, 0], coords[j, 0]], [coords[i, 1], coords[j, 1]], **kwargs)
    ax.scatter(coords[:, 0], coords[:, 1], s=node_size, color="tab:blue", zorder=3)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Lattice graph")
    return ax


def _as_dense(H: np.ndarray) -> np.ndarray:
    if sp.issparse(H):
        return H.toarray()
    return np.asarray(H, dtype=complex)
