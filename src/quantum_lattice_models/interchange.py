"""CSV, NetworkX, and GraphML interchange for portable lattice specifications."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

from quantum_lattice_models.lattice import Bond
from quantum_lattice_models.specs import LatticeSpec

SITE_COLUMNS = (
    "site",
    "x",
    "y",
    "z",
    "site_label",
    "orbital_label",
    "sublattice_label",
    "unit_cell",
)
BOND_COLUMNS = ("source", "target", "has_value", "value_real", "value_imag")


def export_lattice_csv(
    lattice: LatticeSpec,
    sites_path: str | Path,
    bonds_path: str | Path,
    *,
    metadata_path: str | Path | None = None,
) -> tuple[Path, Path, Path]:
    """Export site and bond tables plus a metadata sidecar."""

    lattice.validate()
    sites = Path(sites_path)
    bonds = Path(bonds_path)
    metadata = Path(metadata_path) if metadata_path is not None else _metadata_path(sites)
    for output in (sites, bonds, metadata):
        output.parent.mkdir(parents=True, exist_ok=True)

    with sites.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=SITE_COLUMNS)
        writer.writeheader()
        for site in range(lattice.n_sites):
            position = lattice.positions[site] if lattice.positions else ()
            writer.writerow(
                {
                    "site": site,
                    "x": position[0] if position else "",
                    "y": position[1] if position else "",
                    "z": position[2] if len(position) == 3 else "",
                    "site_label": _item(lattice.site_labels, site),
                    "orbital_label": _item(lattice.orbital_labels, site),
                    "sublattice_label": _item(lattice.sublattice_labels, site),
                    "unit_cell": _item(lattice.unit_cells, site),
                }
            )

    with bonds.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=BOND_COLUMNS)
        writer.writeheader()
        for bond in lattice.bonds:
            value = 0.0j if bond.value is None else complex(bond.value)
            writer.writerow(
                {
                    "source": bond.source,
                    "target": bond.target,
                    "has_value": "0" if bond.value is None else "1",
                    "value_real": value.real,
                    "value_imag": value.imag,
                }
            )

    sidecar = lattice.to_dict()
    for field in (
        "n_sites",
        "positions",
        "bonds",
        "site_labels",
        "orbital_labels",
        "sublattice_labels",
        "unit_cells",
    ):
        sidecar.pop(field, None)
    metadata.write_text(json.dumps(sidecar, indent=2, sort_keys=True) + "\n")
    return sites, bonds, metadata


def import_lattice_csv(
    sites_path: str | Path,
    bonds_path: str | Path,
    *,
    metadata_path: str | Path | None = None,
) -> LatticeSpec:
    """Import paired site and bond tables with an optional metadata sidecar."""

    sites = Path(sites_path)
    bonds = Path(bonds_path)
    metadata = Path(metadata_path) if metadata_path is not None else _metadata_path(sites)
    site_rows = _read_csv(sites, SITE_COLUMNS)
    if not site_rows:
        raise ValueError("Site CSV must contain at least one site.")
    site_rows.sort(key=lambda row: _csv_int(row, "site", "site CSV"))
    site_indices = [_csv_int(row, "site", "site CSV") for row in site_rows]
    if site_indices != list(range(len(site_rows))):
        raise ValueError("Site CSV indices must be contiguous and start at zero.")

    dimensions = _coordinate_dimensions(site_rows)
    positions = (
        tuple(tuple(float(row[axis]) for axis in ("x", "y", "z")[:dimensions]) for row in site_rows)
        if dimensions
        else ()
    )
    parsed_bonds = tuple(
        Bond(
            _csv_int(row, "source", "bond CSV"),
            _csv_int(row, "target", "bond CSV"),
            _csv_bond_value(row),
        )
        for row in _read_csv(bonds, BOND_COLUMNS)
    )
    sidecar: dict[str, object] = {}
    if metadata.exists():
        decoded = json.loads(metadata.read_text())
        if not isinstance(decoded, dict):
            raise ValueError("Lattice metadata sidecar must contain a JSON object.")
        sidecar = decoded

    data = {
        **sidecar,
        "n_sites": len(site_rows),
        "positions": positions,
        "bonds": [
            {"source": bond.source, "target": bond.target, "value": _complex_json(bond.value)}
            for bond in parsed_bonds
        ],
        "site_labels": _optional_column(site_rows, "site_label", str),
        "orbital_labels": _optional_column(site_rows, "orbital_label", str),
        "sublattice_labels": _optional_column(site_rows, "sublattice_label", str),
        "unit_cells": _optional_column(site_rows, "unit_cell", int),
    }
    return LatticeSpec.from_dict(data)


def to_networkx(lattice: LatticeSpec):
    """Return a NetworkX ``MultiDiGraph`` preserving lattice attributes."""

    nx = _networkx()
    lattice.validate()
    graph = nx.MultiDiGraph()
    graph.graph["qlm_schema_version"] = lattice.schema_version
    graph.graph["qlm_lattice_json"] = json.dumps(lattice.to_dict(), sort_keys=True)
    for site in range(lattice.n_sites):
        attributes: dict[str, object] = {"qlm_index": site}
        if lattice.positions:
            for axis, value in zip(("x", "y", "z"), lattice.positions[site], strict=False):
                attributes[axis] = value
        for name, values in (
            ("site_label", lattice.site_labels),
            ("orbital_label", lattice.orbital_labels),
            ("sublattice_label", lattice.sublattice_labels),
            ("unit_cell", lattice.unit_cells),
        ):
            if values:
                attributes[name] = values[site]
        graph.add_node(site, **attributes)
    for bond in lattice.bonds:
        value = 0.0j if bond.value is None else complex(bond.value)
        graph.add_edge(
            bond.source,
            bond.target,
            has_value=bond.value is not None,
            value_real=value.real,
            value_imag=value.imag,
        )
    return graph


def from_networkx(graph) -> LatticeSpec:
    """Construct a portable lattice specification from a NetworkX graph."""

    _networkx()
    nodes = list(graph.nodes(data=True))
    if not nodes:
        raise ValueError("NetworkX graph must contain at least one node.")
    nodes.sort(key=lambda item: int(item[1].get("qlm_index", 0)))
    mapping = {node: index for index, (node, _) in enumerate(nodes)}
    attributes = [data for _, data in nodes]
    dimensions = _networkx_coordinate_dimensions(attributes)
    positions = (
        tuple(
            tuple(float(data[axis]) for axis in ("x", "y", "z")[:dimensions]) for data in attributes
        )
        if dimensions
        else ()
    )
    edges = graph.edges(data=True, keys=True) if graph.is_multigraph() else graph.edges(data=True)
    bonds = []
    for edge in edges:
        source, target, data = (edge[0], edge[1], edge[-1])
        has_value = _as_bool(data.get("has_value", False))
        value = (
            complex(float(data.get("value_real", 0.0)), float(data.get("value_imag", 0.0)))
            if has_value
            else None
        )
        bonds.append(Bond(mapping[source], mapping[target], value))

    base: dict[str, object] = {}
    encoded = graph.graph.get("qlm_lattice_json")
    if isinstance(encoded, str):
        decoded = json.loads(encoded)
        if isinstance(decoded, dict):
            base = decoded
    base.update(
        {
            "n_sites": len(nodes),
            "positions": positions,
            "bonds": [
                {"source": bond.source, "target": bond.target, "value": _complex_json(bond.value)}
                for bond in bonds
            ],
            "site_labels": _networkx_optional(attributes, "site_label", str),
            "orbital_labels": _networkx_optional(attributes, "orbital_label", str),
            "sublattice_labels": _networkx_optional(attributes, "sublattice_label", str),
            "unit_cells": _networkx_optional(attributes, "unit_cell", int),
        }
    )
    if not graph.is_directed():
        conventions = dict(base.get("conventions", {}))
        conventions["edge_interpretation"] = "undirected"
        base["conventions"] = conventions
    return LatticeSpec.from_dict(base)


def export_graphml(lattice: LatticeSpec, path: str | Path) -> Path:
    """Export a lattice as GraphML through the optional NetworkX dependency."""

    nx = _networkx()
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    graph = to_networkx(lattice)
    for _, _, data in graph.edges(data=True):
        data["has_value"] = int(bool(data["has_value"]))
    nx.write_graphml(graph, output)
    return output


def import_graphml(path: str | Path) -> LatticeSpec:
    """Import a GraphML lattice through the optional NetworkX dependency."""

    nx = _networkx()
    return from_networkx(nx.read_graphml(Path(path), node_type=int, force_multigraph=True))


def _networkx():
    try:
        import networkx as nx
    except ImportError as exc:
        raise ImportError(
            "NetworkX/GraphML support requires the optional 'networkx' extra."
        ) from exc
    return nx


def _read_csv(path: Path, expected: tuple[str, ...]) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = sorted(set(expected) - set(reader.fieldnames or ()))
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(missing)}.")
        return [dict(row) for row in reader]


def _coordinate_dimensions(rows: list[dict[str, str]]) -> int:
    x_present = [bool(row["x"].strip()) for row in rows]
    y_present = [bool(row["y"].strip()) for row in rows]
    z_present = [bool(row["z"].strip()) for row in rows]
    if not any(x_present + y_present + z_present):
        return 0
    if not all(x_present) or not all(y_present):
        raise ValueError("Site CSV coordinates require x and y for every site.")
    if any(z_present) and not all(z_present):
        raise ValueError("Site CSV z coordinates must be present for every site or none.")
    return 3 if all(z_present) else 2


def _networkx_coordinate_dimensions(attributes: list[dict[str, Any]]) -> int:
    present = [{axis: axis in data for axis in ("x", "y", "z")} for data in attributes]
    if not any(item["x"] or item["y"] or item["z"] for item in present):
        return 0
    if not all(item["x"] and item["y"] for item in present):
        raise ValueError("NetworkX nodes require x and y coordinates consistently.")
    if any(item["z"] for item in present) and not all(item["z"] for item in present):
        raise ValueError("NetworkX z coordinates must be present for every node or none.")
    return 3 if all(item["z"] for item in present) else 2


def _optional_column(rows: list[dict[str, str]], name: str, converter) -> list[object]:
    values = [row[name].strip() for row in rows]
    if not any(values):
        return []
    if not all(values):
        raise ValueError(f"CSV column {name!r} must be complete or entirely empty.")
    return [converter(value) for value in values]


def _networkx_optional(attributes: list[dict[str, Any]], name: str, converter) -> list[object]:
    values = [data.get(name) for data in attributes]
    if not any(value is not None for value in values):
        return []
    if any(value is None for value in values):
        raise ValueError(f"NetworkX node attribute {name!r} must be present on every node.")
    return [converter(value) for value in values]


def _csv_int(row: dict[str, str], name: str, source: str) -> int:
    try:
        return int(row[name])
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{source} field {name!r} must be an integer.") from exc


def _csv_bond_value(row: dict[str, str]) -> complex | None:
    if not _as_bool(row["has_value"]):
        return None
    try:
        return complex(float(row["value_real"]), float(row["value_imag"]))
    except ValueError as exc:
        raise ValueError("Bond CSV values must contain numeric real and imaginary parts.") from exc


def _as_bool(value: object) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes"}


def _item(values: tuple[object, ...], index: int) -> object:
    return values[index] if values else ""


def _metadata_path(sites_path: Path) -> Path:
    return sites_path.with_suffix(sites_path.suffix + ".json")


def _complex_json(value: complex | None) -> object:
    if value is None:
        return None
    number = complex(value)
    return {"__complex__": [number.real, number.imag]}
