"""Matplotlib helpers for spectra, densities, probabilities, and lattice graphs."""

from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models.models import ssh_edge_state_localization
from quantum_lattice_models.spectra import density_of_states, eigenvalues

COLORBLIND_PALETTE = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "red": "#D55E00",
    "purple": "#CC79A7",
    "sky": "#56B4E9",
    "yellow": "#F0E442",
    "black": "#000000",
}


def plot_spectrum(
    H: np.ndarray,
    ax=None,
    *,
    highlight_gap: bool = False,
    zero_line: bool = False,
    **scatter_kwargs,
):
    """Plot sorted eigenvalues and return the Matplotlib axes."""

    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    values = np.real_if_close(eigenvalues(_as_dense(H))).real
    values = np.sort(values)
    kwargs = {"s": 24, "color": COLORBLIND_PALETTE["blue"]}
    kwargs.update(scatter_kwargs)
    ax.scatter(np.arange(values.size), values, **kwargs)
    if highlight_gap and values.size:
        count = min(2, values.size)
        ax.scatter(
            np.arange(count),
            values[:count],
            s=max(float(kwargs.get("s", 24)) * 2.2, 44),
            color=[COLORBLIND_PALETTE["red"], COLORBLIND_PALETTE["orange"]][:count],
            edgecolor="black",
            linewidth=0.55,
            zorder=4,
        )
        if values.size > 1:
            ax.annotate(
                "",
                xy=(0.5, values[1]),
                xytext=(0.5, values[0]),
                arrowprops={"arrowstyle": "<->", "color": "0.35", "linewidth": 0.9},
            )
            ax.annotate(
                f"E₀ = {values[0]:.3g}\ngap = {values[1] - values[0]:.3g}",
                xy=(0.5, 0.5 * (values[0] + values[1])),
                xytext=(8, 0),
                textcoords="offset points",
                va="center",
                fontsize=9,
                color="0.25",
            )
            ax.plot([0, 1], [values[0], values[0]], color="0.45", linewidth=0.8, alpha=0.7)
    if zero_line:
        ax.axhline(0.0, color="0.35", linewidth=0.9, linestyle="--", alpha=0.8)
    ax.set_xlabel("Eigenvalue index")
    ax.set_ylabel("Energy")
    ax.set_title("Spectrum")
    apply_plot_style(ax)
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
        "color": COLORBLIND_PALETTE["green"],
        "edgecolor": "black",
        "alpha": 0.75,
    }
    kwargs.update(bar_kwargs)
    ax.bar(edges[:-1], counts, width=widths, **kwargs)
    ax.set_xlabel("Energy")
    ax.set_ylabel("Count")
    ax.set_title("Density of states")
    apply_plot_style(ax)
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

    kwargs = {"marker": "o", "color": COLORBLIND_PALETTE["red"]}
    kwargs.update(line_kwargs)
    ax.plot(np.arange(vector.size), probability, **kwargs)
    ax.set_xlabel("SSH site")
    ax.set_ylabel("Probability")
    ax.set_title(f"SSH edge-state weight: {localization:.3f}")
    apply_plot_style(ax)
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
    kwargs = {
        "align": "edge",
        "color": COLORBLIND_PALETTE["purple"],
        "edgecolor": "black",
        "alpha": 0.75,
    }
    kwargs.update(bar_kwargs)
    ax.bar(edges[:-1], counts, width=np.diff(edges), **kwargs)
    ax.set_xlabel("Value")
    ax.set_ylabel("Count")
    ax.set_title("Density")
    apply_plot_style(ax)
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
    kwargs = {"marker": "o", "color": COLORBLIND_PALETTE["red"]}
    kwargs.update(line_kwargs)
    ax.plot(np.arange(probabilities.size), probabilities, **kwargs)
    ax.set_xlabel("Site")
    ax.set_ylabel("Probability")
    ax.set_title(title)
    apply_plot_style(ax)
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

    return plot_parameter_sweep(
        builder,
        flux_values,
        parameter_name="Flux per plaquette",
        title="Hofstadter spectrum",
        ax=ax,
        **scatter_kwargs,
    )


def plot_parameter_sweep(
    builder,
    values: np.ndarray,
    ax=None,
    *,
    parameter_name: str = "Parameter",
    title: str = "Parameter sweep spectrum",
    **scatter_kwargs,
):
    """Plot eigenvalues over a one-dimensional model parameter sweep.

    ``builder`` must accept one parameter value and return a Hamiltonian matrix.
    """

    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    parameters = np.asarray(values, dtype=float).reshape(-1)
    kwargs = {"s": 8, "color": COLORBLIND_PALETTE["blue"], "alpha": 0.65}
    kwargs.update(scatter_kwargs)
    for parameter in parameters:
        eigs = np.real_if_close(eigenvalues(_as_dense(builder(float(parameter))))).real
        ax.scatter(np.full(eigs.size, parameter), eigs, **kwargs)
    ax.set_xlabel(parameter_name)
    ax.set_ylabel("Energy")
    ax.set_title(title)
    apply_plot_style(ax)
    return ax


def plot_lattice_graph(
    H: np.ndarray,
    positions: np.ndarray | None = None,
    ax=None,
    *,
    threshold: float = 1e-12,
    node_size: float = 34,
    node_color: str | list[str] | np.ndarray = COLORBLIND_PALETTE["blue"],
    sublattices: np.ndarray | None = None,
    sublattice_colors: tuple[str, ...] = (
        COLORBLIND_PALETTE["blue"],
        COLORBLIND_PALETTE["orange"],
        COLORBLIND_PALETTE["green"],
    ),
    unit_cells: np.ndarray | None = None,
    show_unit_cells: bool = False,
    show_sublattice_legend: bool = False,
    edge_color_by: str = "phase",
    scale_edges: bool = True,
    show_colorbar: bool = False,
    labels: bool = False,
    arrows: bool = False,
    **line_kwargs,
):
    """Plot graph connectivity implied by nonzero off-diagonal matrix entries.

    When ``positions`` is omitted, the function uses ``H.metadata["positions"]``
    if present. Edge widths can encode hopping magnitude, and edge colors can
    encode hopping phase for complex-valued models. Optional per-site
    ``sublattices`` and ``unit_cells`` arrays control node colors and dashed
    unit-cell outlines.
    """

    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection

    if ax is None:
        _, ax = plt.subplots()
    matrix = _as_dense(H)
    coords = _positions_for_plot(H, positions, matrix.shape[0])
    explicit_color = "color" in line_kwargs
    kwargs = {"linewidth": 1.0, "alpha": 0.78}
    if edge_color_by not in ("phase", "magnitude", "none"):
        raise ValueError("edge_color_by must be 'phase', 'magnitude', or 'none'.")
    if explicit_color:
        kwargs["color"] = line_kwargs.pop("color")
    elif edge_color_by == "none":
        kwargs["color"] = "0.55"
    kwargs.update(line_kwargs)

    segments = []
    phases = []
    magnitudes = []
    arrow_edges = []
    for i in range(matrix.shape[0]):
        for j in range(i + 1, matrix.shape[1]):
            value = matrix[i, j]
            if abs(value) > threshold:
                segments.append([(coords[i, 0], coords[i, 1]), (coords[j, 0], coords[j, 1])])
                phases.append(np.angle(value))
                magnitudes.append(abs(value))
                arrow_edges.append((i, j, value))

    if segments:
        linewidths = _scaled_linewidths(
            np.asarray(magnitudes), kwargs.pop("linewidth"), scale_edges
        )
        collection = LineCollection(segments, linewidths=linewidths, **kwargs)
        if not explicit_color and edge_color_by in ("phase", "magnitude"):
            if edge_color_by == "phase":
                collection.set_array(np.asarray(phases))
                collection.set_cmap("twilight")
                collection.set_clim(-np.pi, np.pi)
            else:
                collection.set_array(np.asarray(magnitudes))
                collection.set_cmap("viridis")
        ax.add_collection(collection)
        if show_colorbar and not explicit_color and edge_color_by in ("phase", "magnitude"):
            cbar = ax.figure.colorbar(collection, ax=ax, fraction=0.046, pad=0.04)
            cbar.set_label("Hopping phase" if edge_color_by == "phase" else "Hopping magnitude")
            if edge_color_by == "phase":
                _set_phase_colorbar_ticks(cbar)
        if arrows:
            _draw_arrows(ax, coords, arrow_edges, threshold)

    if sublattices is not None:
        labels_array = np.asarray(sublattices).reshape(-1)
        if labels_array.size != matrix.shape[0]:
            raise ValueError("sublattices must contain one label per lattice site.")
        unique_labels = list(dict.fromkeys(labels_array.tolist()))
        color_map = {
            label: sublattice_colors[index % len(sublattice_colors)]
            for index, label in enumerate(unique_labels)
        }
        node_color = [color_map[label] for label in labels_array]
        if show_sublattice_legend:
            from matplotlib.lines import Line2D

            handles = [
                Line2D(
                    [],
                    [],
                    marker="o",
                    linestyle="none",
                    markerfacecolor=color_map[label],
                    markeredgecolor="none",
                    label=f"sublattice {label}",
                )
                for label in unique_labels
            ]
            ax.legend(handles=handles, loc="best", frameon=False, fontsize=8)
    if show_unit_cells:
        _draw_unit_cell_outlines(ax, coords, unit_cells)
    ax.scatter(coords[:, 0], coords[:, 1], s=node_size, color=node_color, zorder=3)
    if labels:
        for index, (x_coord, y_coord) in enumerate(coords):
            ax.text(x_coord, y_coord, str(index), ha="center", va="center", fontsize=8, zorder=4)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Lattice graph")
    ax.autoscale_view()
    return ax


def plot_lattice_state(
    H: np.ndarray,
    state: np.ndarray,
    positions: np.ndarray | None = None,
    ax=None,
    *,
    normalize: bool = True,
    node_size: float = 220,
    min_node_size: float = 24,
    show_colorbar: bool = True,
    draw_edges: bool = True,
):
    """Plot a single-particle eigenstate or state vector on lattice positions."""

    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    matrix = _as_dense(H)
    coords = _positions_for_plot(H, positions, matrix.shape[0])
    vector = np.asarray(state, dtype=complex).reshape(-1)
    if vector.size != matrix.shape[0]:
        raise ValueError("state length must match the Hamiltonian dimension.")

    probabilities = np.abs(vector) ** 2
    norm = probabilities.sum()
    if normalize:
        if norm == 0:
            raise ValueError("state must have nonzero norm when normalize=True.")
        probabilities = probabilities / norm
    phases = np.angle(vector)

    if draw_edges:
        plot_lattice_graph(H, coords, ax=ax, node_size=0, color="0.82")

    max_probability = probabilities.max()
    if max_probability == 0:
        sizes = np.full(probabilities.shape, min_node_size)
    else:
        sizes = min_node_size + node_size * probabilities / max_probability
    scatter = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        s=sizes,
        c=phases,
        cmap="twilight",
        vmin=-np.pi,
        vmax=np.pi,
        edgecolors="black",
        linewidths=0.45,
        zorder=4,
    )
    if show_colorbar:
        cbar = ax.figure.colorbar(scatter, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("State phase")
        _set_phase_colorbar_ticks(cbar)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title("Lattice state")
    ax.autoscale_view()
    return ax


def plot_hamiltonian_matrix(
    H: np.ndarray,
    ax=None,
    *,
    mode: str = "magnitude",
    show_colorbar: bool = True,
    **imshow_kwargs,
):
    """Visualize Hamiltonian matrix entries."""

    import matplotlib.pyplot as plt

    if ax is None:
        _, ax = plt.subplots()
    matrix = _as_dense(H)
    data, label, cmap = _matrix_image_data(matrix, mode)
    kwargs = {"origin": "lower", "interpolation": "nearest", "cmap": cmap}
    kwargs.update(imshow_kwargs)
    image = ax.imshow(data, **kwargs)
    if show_colorbar:
        cbar = ax.figure.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label(label)
        if mode == "phase":
            image.set_clim(-np.pi, np.pi)
            _set_phase_colorbar_ticks(cbar)
    ax.set_xlabel("Ket index")
    ax.set_ylabel("Bra index")
    ax.set_title(f"Hamiltonian {mode}")
    return ax


def apply_plot_style(ax, *, grid: bool = True):
    """Apply a compact, consistent style to Matplotlib axes."""

    if grid:
        ax.grid(True, color="0.88", linewidth=0.8, alpha=0.8)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.tick_params(direction="out", length=3, width=0.8)
    return ax


def _as_dense(H: np.ndarray) -> np.ndarray:
    if sp.issparse(H):
        return H.toarray()
    return np.asarray(H, dtype=complex)


def _positions_for_plot(H: np.ndarray, positions: np.ndarray | None, n_sites: int) -> np.ndarray:
    if positions is None:
        positions = getattr(H, "metadata", {}).get("positions")
    if positions is None:
        raise ValueError("positions must be provided or stored in H.metadata['positions'].")
    coords = np.asarray(positions, dtype=float)
    if coords.ndim != 2 or coords.shape[0] != n_sites or coords.shape[1] not in (2, 3):
        raise ValueError("positions must have shape (n_sites, 2) or (n_sites, 3).")
    return coords[:, :2]


def _scaled_linewidths(values: np.ndarray, base: float, scale_edges: bool) -> np.ndarray:
    if not scale_edges or values.size == 0:
        return np.full(values.shape, base)
    maximum = values.max()
    if maximum == 0:
        return np.full(values.shape, base)
    return base * (0.75 + 1.75 * values / maximum)


def _draw_arrows(
    ax, coords: np.ndarray, arrow_edges: list[tuple[int, int, complex]], threshold: float
) -> None:
    from matplotlib.patches import FancyArrowPatch

    for source, target, value in arrow_edges:
        reverse = value.conjugate()
        if abs(reverse) <= threshold:
            patch = FancyArrowPatch(
                coords[source],
                coords[target],
                arrowstyle="-|>",
                mutation_scale=8,
                color="0.35",
                linewidth=0.8,
                alpha=0.65,
                zorder=2,
            )
            ax.add_patch(patch)


def _set_phase_colorbar_ticks(colorbar) -> None:
    colorbar.set_ticks([-np.pi, 0.0, np.pi])
    colorbar.set_ticklabels([r"$-\pi$", "0", r"$\pi$"])


def _draw_unit_cell_outlines(ax, coords: np.ndarray, unit_cells: np.ndarray | None) -> None:
    from matplotlib.patches import FancyBboxPatch

    if unit_cells is None:
        raise ValueError("unit_cells must be provided when show_unit_cells=True.")
    cell_labels = np.asarray(unit_cells).reshape(-1)
    if cell_labels.size != coords.shape[0]:
        raise ValueError("unit_cells must contain one label per lattice site.")
    for label in dict.fromkeys(cell_labels.tolist()):
        cell_coords = coords[cell_labels == label]
        x_min, y_min = cell_coords.min(axis=0)
        x_max, y_max = cell_coords.max(axis=0)
        padding = 0.16
        outline = FancyBboxPatch(
            (x_min - padding, y_min - padding),
            max(x_max - x_min + 2 * padding, 2 * padding),
            max(y_max - y_min + 2 * padding, 2 * padding),
            boxstyle="round,pad=0.02",
            fill=False,
            edgecolor="0.45",
            linewidth=0.8,
            linestyle="--",
            alpha=0.65,
            zorder=1,
        )
        ax.add_patch(outline)


def _matrix_image_data(matrix: np.ndarray, mode: str) -> tuple[np.ndarray, str, object]:
    if mode == "real":
        return matrix.real, "Real part", "coolwarm"
    if mode == "imag":
        return matrix.imag, "Imaginary part", "coolwarm"
    if mode == "magnitude":
        return np.abs(matrix), "Magnitude", "viridis"
    if mode == "phase":
        from matplotlib import colormaps

        cmap = colormaps["twilight"].copy()
        cmap.set_bad("#f5f5f5")
        phases = np.ma.masked_where(np.abs(matrix) <= 1e-14, np.angle(matrix))
        return phases, "Phase", cmap
    raise ValueError("mode must be one of 'real', 'imag', 'magnitude', or 'phase'.")
