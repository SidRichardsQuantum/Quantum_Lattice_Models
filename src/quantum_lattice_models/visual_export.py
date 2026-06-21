"""Deterministic SVG and plot-data exports for portable lattice workflows."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from quantum_lattice_models.periodic import BandStructure, PeriodicLatticeSpec
from quantum_lattice_models.specs import LatticeSpec, ModelSpec


def export_lattice_plot_data(
    lattice: LatticeSpec,
    path: str | Path,
    *,
    format: str | None = None,
) -> Path:
    """Export lattice site and bond coordinates as JSON or a single CSV table."""

    lattice.validate()
    if not lattice.positions:
        raise ValueError("Lattice plot-data export requires positions.")
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output_format = (format or output.suffix.lstrip(".")).lower()
    records = _lattice_records(lattice)
    if output_format == "json":
        output.write_text(json.dumps(records, indent=2, sort_keys=True) + "\n")
    elif output_format == "csv":
        fields = ("kind", "index", "source", "target", "x", "y", "target_x", "target_y")
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            writer.writerows(records["sites"] + records["bonds"])
    else:
        raise ValueError("Plot-data format must be 'json' or 'csv'.")
    return output


def export_band_data(
    bands: BandStructure,
    path: str | Path,
    *,
    format: str | None = None,
) -> Path:
    """Export band coordinates and energies as deterministic JSON or CSV."""

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output_format = (format or output.suffix.lstrip(".")).lower()
    if output_format == "json":
        output.write_text(json.dumps(bands.to_dict(), indent=2, sort_keys=True) + "\n")
    elif output_format == "csv":
        fields = ["point", "distance", *[f"k{axis}" for axis in range(bands.momenta.shape[1])]]
        fields.extend(f"band_{band}" for band in range(bands.energies.shape[1]))
        with output.open("w", newline="", encoding="utf-8") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields)
            writer.writeheader()
            for point, (momentum, distance, energies) in enumerate(
                zip(bands.momenta, bands.distances, bands.energies, strict=True)
            ):
                row = {"point": point, "distance": distance}
                row.update({f"k{axis}": value for axis, value in enumerate(momentum)})
                row.update({f"band_{band}": value for band, value in enumerate(energies)})
                writer.writerow(row)
    else:
        raise ValueError("Band-data format must be 'json' or 'csv'.")
    return output


def export_analysis_svg(result, path: str | Path, *, width: int = 720, height: int = 480) -> Path:
    """Render any supported portable analysis plot to SVG."""

    import matplotlib.pyplot as plt

    from quantum_lattice_models.plotting import plot_analysis_result

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    figure, axis = plt.subplots(figsize=(width / 100, height / 100))
    plot_analysis_result(result, ax=axis)
    figure.tight_layout()
    figure.savefig(output, format="svg")
    plt.close(figure)
    return output


def export_analysis_pdf(result, path: str | Path) -> Path:
    """Render any supported portable analysis plot to PDF."""

    import matplotlib.pyplot as plt

    from quantum_lattice_models.plotting import plot_analysis_result

    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    axis = plot_analysis_result(result)
    axis.figure.tight_layout()
    axis.figure.savefig(output, format="pdf")
    plt.close(axis.figure)
    return output


def export_matrix_plot_data(matrix: np.ndarray, path: str | Path) -> Path:
    """Export matrix real, imaginary, magnitude, phase, and sparsity arrays as JSON."""

    values = np.asarray(matrix, dtype=complex)
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(
            {
                "shape": list(values.shape),
                "real": values.real.tolist(),
                "imaginary": values.imag.tolist(),
                "magnitude": np.abs(values).tolist(),
                "phase": np.angle(values).tolist(),
                "nonzero": (np.abs(values) > 0).astype(int).tolist(),
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    return output


def export_lattice_svg(
    lattice: LatticeSpec,
    path: str | Path,
    *,
    width: int = 720,
    height: int = 480,
    labels: bool = True,
) -> Path:
    """Export a finite lattice diagram as dependency-free deterministic SVG."""

    lattice.validate()
    if not lattice.positions:
        raise ValueError("Lattice SVG export requires positions.")
    coordinates = np.asarray(lattice.positions, dtype=float)[:, :2]
    projected = _project(coordinates, width, height)
    lines = _svg_header(width, height)
    for bond in lattice.bonds:
        start, end = projected[bond.source], projected[bond.target]
        lines.append(
            f'<line x1="{start[0]:.3f}" y1="{start[1]:.3f}" '
            f'x2="{end[0]:.3f}" y2="{end[1]:.3f}" stroke="#666" stroke-width="2"/>'
        )
    for index, (x_coord, y_coord) in enumerate(projected):
        lines.append(
            f'<circle cx="{x_coord:.3f}" cy="{y_coord:.3f}" r="7" '
            'fill="#0072B2" stroke="#fff" stroke-width="1"/>'
        )
        if labels:
            label = lattice.site_labels[index] if lattice.site_labels else str(index)
            lines.append(
                f'<text x="{x_coord + 10:.3f}" y="{y_coord - 10:.3f}" '
                f'font-family="sans-serif" font-size="12">{_escape(label)}</text>'
            )
    lines.append("</svg>")
    return _write_svg(path, lines)


def export_periodic_svg(
    lattice: PeriodicLatticeSpec,
    path: str | Path,
    *,
    repeats: tuple[int, int] = (3, 3),
    width: int = 720,
    height: int = 480,
) -> Path:
    """Export a repeated two-dimensional unit-cell diagram as SVG."""

    from quantum_lattice_models.transformations import repeat_unit_cell

    if lattice.dimension != 2:
        raise ValueError("Periodic SVG export currently requires a two-dimensional lattice.")
    finite = repeat_unit_cell(lattice, repeats)
    return export_lattice_svg(finite, path, width=width, height=height, labels=False)


def export_interaction_plot_data(
    model: ModelSpec,
    path: str | Path,
) -> Path:
    """Export physical degrees and interaction edges as deterministic JSON."""

    model.validate()
    if not model.local_degrees:
        raise ValueError("Interaction plot-data export requires local degrees.")
    coordinates = _model_degree_coordinates(model)
    data = {
        "degrees": [
            {
                "index": degree.index,
                "site": degree.site,
                "kind": degree.kind,
                "label": degree.label,
                "x": coordinates[degree.index, 0],
                "y": coordinates[degree.index, 1],
            }
            for degree in model.local_degrees
        ],
        "interactions": [
            {
                "degrees": list(interaction.degrees),
                "operators": list(interaction.operators),
                "coefficient": {
                    "real": complex(interaction.coefficient).real,
                    "imag": complex(interaction.coefficient).imag,
                },
                "kind": interaction.kind,
                "label": interaction.label,
            }
            for interaction in model.interactions
        ],
    }
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")
    return output


def export_interaction_svg(
    model: ModelSpec,
    path: str | Path,
    *,
    width: int = 720,
    height: int = 480,
) -> Path:
    """Export a portable physical interaction graph as deterministic SVG."""

    model.validate()
    coordinates = _model_degree_coordinates(model)
    projected = _project(coordinates, width, height)
    lines = _svg_header(width, height)
    for interaction in model.interactions:
        if len(interaction.degrees) != 2:
            continue
        source, target = interaction.degrees
        start, end = projected[source], projected[target]
        lines.append(
            f'<line x1="{start[0]:.3f}" y1="{start[1]:.3f}" '
            f'x2="{end[0]:.3f}" y2="{end[1]:.3f}" stroke="#666" stroke-width="2"/>'
        )
    for degree in model.local_degrees:
        x_coord, y_coord = projected[degree.index]
        lines.append(
            f'<circle cx="{x_coord:.3f}" cy="{y_coord:.3f}" r="7" '
            'fill="#0072B2" stroke="#fff" stroke-width="1"/>'
        )
        lines.append(
            f'<text x="{x_coord + 10:.3f}" y="{y_coord - 10:.3f}" '
            f'font-family="sans-serif" font-size="12">{_escape(degree.label)}</text>'
        )
    lines.append("</svg>")
    return _write_svg(path, lines)


def _lattice_records(lattice: LatticeSpec) -> dict[str, list[dict[str, object]]]:
    positions = np.asarray(lattice.positions)
    sites = [
        {
            "kind": "site",
            "index": index,
            "source": "",
            "target": "",
            "x": coordinate[0],
            "y": coordinate[1],
            "target_x": "",
            "target_y": "",
        }
        for index, coordinate in enumerate(positions)
    ]
    bonds = [
        {
            "kind": "bond",
            "index": index,
            "source": bond.source,
            "target": bond.target,
            "x": positions[bond.source, 0],
            "y": positions[bond.source, 1],
            "target_x": positions[bond.target, 0],
            "target_y": positions[bond.target, 1],
        }
        for index, bond in enumerate(lattice.bonds)
    ]
    return {"sites": sites, "bonds": bonds}


def _model_degree_coordinates(model: ModelSpec) -> np.ndarray:
    if not model.local_degrees:
        raise ValueError("Interaction SVG export requires local degrees.")
    if model.lattice is not None and model.lattice.positions:
        site_coordinates = np.asarray(model.lattice.positions, dtype=float)[:, :2]
    else:
        site_count = max(degree.site for degree in model.local_degrees) + 1
        site_coordinates = np.column_stack(
            (np.arange(site_count, dtype=float), np.zeros(site_count))
        )
    grouped: dict[int, list[int]] = {}
    for degree in model.local_degrees:
        grouped.setdefault(degree.site, []).append(degree.index)
    coordinates = np.zeros((len(model.local_degrees), 2), dtype=float)
    for site, indices in grouped.items():
        offsets = np.linspace(-0.12, 0.12, len(indices)) if len(indices) > 1 else np.zeros(1)
        for index, offset in zip(indices, offsets, strict=True):
            coordinates[index] = site_coordinates[site] + (0.0, offset)
    return coordinates


def _project(coordinates: np.ndarray, width: int, height: int) -> np.ndarray:
    minimum = coordinates.min(axis=0)
    span = np.ptp(coordinates, axis=0)
    span[span == 0] = 1.0
    normalized = (coordinates - minimum) / span
    margin = 40.0
    result = np.empty_like(normalized)
    result[:, 0] = margin + normalized[:, 0] * (width - 2 * margin)
    result[:, 1] = height - margin - normalized[:, 1] * (height - 2 * margin)
    return result


def _svg_header(width: int, height: int) -> list[str]:
    return [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
    ]


def _write_svg(path: str | Path, lines: list[str]) -> Path:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines) + "\n")
    return output


def _escape(value: str) -> str:
    return (
        value.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")
    )
