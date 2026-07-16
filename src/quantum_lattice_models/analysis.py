"""Portable numerical analysis records and result-producing helpers."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from quantum_lattice_models._model_utils import as_dense
from quantum_lattice_models.observables import (
    site_magnetization_z,
    spin_correlation_matrix,
)
from quantum_lattice_models.spectra import eigenvalues

if TYPE_CHECKING:
    from quantum_lattice_models.periodic import BandStructure, PeriodicLatticeSpec
    from quantum_lattice_models.specs import ModelSpec
    from quantum_lattice_models.spin import FixedMagnetizationBasis, SpinParityBasis

ANALYSIS_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class AnalysisResult:
    """Portable numerical result with coordinates, provenance, and plot metadata."""

    analysis: str
    values: dict[str, np.ndarray]
    coordinates: dict[str, np.ndarray] = field(default_factory=dict)
    parameters: dict[str, object] = field(default_factory=dict)
    source: dict[str, object] = field(default_factory=dict)
    solver: dict[str, object] = field(default_factory=dict)
    warnings: tuple[str, ...] = ()
    plot: dict[str, object] = field(default_factory=dict)
    provenance: dict[str, object] = field(default_factory=dict)
    metadata: dict[str, object] = field(default_factory=dict)
    package_version: str = field(default_factory=lambda: _package_version())
    created_at: str | None = None
    schema_version: str = ANALYSIS_SCHEMA_VERSION

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "values",
            {str(name): np.asarray(value) for name, value in self.values.items()},
        )
        object.__setattr__(
            self,
            "coordinates",
            {str(name): np.asarray(value) for name, value in self.coordinates.items()},
        )
        self.validate()

    def validate(self) -> None:
        if self.schema_version != ANALYSIS_SCHEMA_VERSION:
            raise ValueError(f"Unsupported analysis schema_version {self.schema_version!r}.")
        if not isinstance(self.analysis, str) or not self.analysis:
            raise ValueError("analysis must be a nonempty string.")
        if not self.values:
            raise ValueError("Analysis results require at least one values array.")
        for group_name, arrays in (("values", self.values), ("coordinates", self.coordinates)):
            if not all(isinstance(name, str) and name for name in arrays):
                raise ValueError(f"{group_name} keys must be nonempty strings.")
            for name, array in arrays.items():
                if array.dtype == object:
                    raise ValueError(f"{group_name}.{name} must not use object dtype.")
                if np.issubdtype(array.dtype, np.number) and not np.all(np.isfinite(array)):
                    raise ValueError(f"{group_name}.{name} must contain only finite values.")
        if not all(isinstance(item, str) for item in self.warnings):
            raise ValueError("warnings must contain strings.")
        for name, value in (
            ("parameters", self.parameters),
            ("source", self.source),
            ("solver", self.solver),
            ("plot", self.plot),
            ("provenance", self.provenance),
            ("metadata", self.metadata),
        ):
            if not isinstance(value, dict):
                raise ValueError(f"{name} must be a dictionary.")
            _json_value(value)

    def summary(self) -> dict[str, object]:
        return {
            "schema_version": self.schema_version,
            "analysis": self.analysis,
            "value_shapes": {name: list(value.shape) for name, value in self.values.items()},
            "coordinate_shapes": {
                name: list(value.shape) for name, value in self.coordinates.items()
            },
            "parameters": _json_value(self.parameters),
            "source_kind": self.source.get("kind"),
            "solver": _json_value(self.solver),
            "warnings": list(self.warnings),
            "plot_kind": self.plot.get("kind"),
            "package_version": self.package_version,
            "created_at": self.created_at,
            "metadata": _json_value(self.metadata),
        }

    def to_dict(self) -> dict[str, object]:
        self.validate()
        return {
            "schema_version": self.schema_version,
            "analysis": self.analysis,
            "values": {name: _array_dict(value) for name, value in self.values.items()},
            "coordinates": {name: _array_dict(value) for name, value in self.coordinates.items()},
            "parameters": _json_value(self.parameters),
            "source": _json_value(self.source),
            "solver": _json_value(self.solver),
            "warnings": list(self.warnings),
            "plot": _json_value(self.plot),
            "provenance": _json_value(self.provenance),
            "metadata": _json_value(self.metadata),
            "package_version": self.package_version,
            "created_at": self.created_at,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> AnalysisResult:
        allowed = {
            "schema_version",
            "analysis",
            "values",
            "coordinates",
            "parameters",
            "source",
            "solver",
            "warnings",
            "plot",
            "provenance",
            "metadata",
            "package_version",
            "created_at",
        }
        unknown = sorted(set(data) - allowed)
        if unknown:
            raise ValueError(f"Analysis result contains unknown fields: {', '.join(unknown)}.")
        values = _array_mapping(data.get("values"), "values")
        coordinates = _array_mapping(data.get("coordinates", {}), "coordinates")
        return cls(
            schema_version=str(data.get("schema_version", "")),
            analysis=str(data.get("analysis", "")),
            values=values,
            coordinates=coordinates,
            parameters=_mapping(data.get("parameters", {}), "parameters"),
            source=_mapping(data.get("source", {}), "source"),
            solver=_mapping(data.get("solver", {}), "solver"),
            warnings=tuple(str(item) for item in data.get("warnings", [])),
            plot=_mapping(data.get("plot", {}), "plot"),
            provenance=_mapping(data.get("provenance", {}), "provenance"),
            metadata=_mapping(data.get("metadata", {}), "metadata"),
            package_version=str(data.get("package_version", "")),
            created_at=(None if data.get("created_at") is None else str(data.get("created_at"))),
        )


def spectrum_result(
    hamiltonian: np.ndarray,
    *,
    model: ModelSpec | None = None,
    parameters: dict[str, object] | None = None,
) -> AnalysisResult:
    """Return a portable complete-spectrum result."""

    matrix = as_dense(hamiltonian)
    values = np.sort(np.real_if_close(eigenvalues(matrix)).real)
    return AnalysisResult(
        analysis="spectrum",
        coordinates={"index": np.arange(values.size)},
        values={"eigenvalues": values},
        parameters=dict(parameters or {}),
        source=_model_source(model),
        solver={
            "method": (
                "numpy.linalg.eigvalsh"
                if np.allclose(matrix, matrix.conj().T)
                else "numpy.linalg.eigvals"
            ),
            "exact": True,
        },
        plot={"kind": "spectrum", "x": "index", "y": "eigenvalues"},
        metadata={"dimension": matrix.shape[0]},
    )


def band_structure_result(
    bands: BandStructure,
    *,
    periodic: PeriodicLatticeSpec | None = None,
    parameters: dict[str, object] | None = None,
) -> AnalysisResult:
    """Convert a band structure into a portable analysis result."""

    values = {"energies": bands.energies}
    if bands.eigenvectors is not None:
        values["eigenvectors"] = bands.eigenvectors
    return AnalysisResult(
        analysis="band_structure",
        coordinates={"momenta": bands.momenta, "distance": bands.distances},
        values=values,
        parameters=dict(parameters or {}),
        source=(
            {"kind": "periodic_lattice", "periodic_lattice": periodic.to_dict()}
            if periodic is not None
            else {}
        ),
        solver={"method": "numpy.linalg.eigh", "exact": True},
        plot={
            "kind": "bands",
            "x": "distance",
            "y": "energies",
            "labels": list(bands.labels),
            "label_indices": list(bands.label_indices),
        },
        metadata=dict(bands.metadata),
    )


def topology_result(
    invariant: str,
    value: float | int,
    *,
    periodic: PeriodicLatticeSpec | None = None,
    parameters: dict[str, object] | None = None,
    solver: dict[str, object] | None = None,
) -> AnalysisResult:
    """Return a portable scalar topological-invariant result."""

    if invariant not in {"zak", "winding", "chern"}:
        raise ValueError("invariant must be 'zak', 'winding', or 'chern'.")
    return AnalysisResult(
        analysis=f"{invariant}_phase" if invariant == "zak" else f"{invariant}_number",
        values={"value": np.asarray(value)},
        parameters=dict(parameters or {}),
        source=(
            {"kind": "periodic_lattice", "periodic_lattice": periodic.to_dict()}
            if periodic is not None
            else {}
        ),
        solver=dict(solver or {}),
        plot={"kind": "scalar", "label": invariant},
    )


def spin_observables_result(
    state: np.ndarray,
    n_sites: int,
    *,
    axis: str = "Z",
    connected: bool = False,
    basis: FixedMagnetizationBasis | SpinParityBasis | None = None,
    model: ModelSpec | None = None,
) -> AnalysisResult:
    """Return site magnetization and same-axis correlation records."""

    magnetization = site_magnetization_z(state, n_sites, basis=basis)
    correlations = spin_correlation_matrix(
        state,
        n_sites,
        axis=axis,
        basis=basis,
        connected=connected,
    )
    return AnalysisResult(
        analysis="spin_observables",
        coordinates={"site": np.arange(n_sites)},
        values={"site_magnetization_z": magnetization, "correlations": correlations},
        parameters={"axis": axis.upper(), "connected": connected, "n_sites": n_sites},
        source=_model_source(model),
        solver={"method": "direct expectation values", "exact": True},
        plot={
            "kind": "observable",
            "primary": "site_magnetization_z",
            "matrix": "correlations",
        },
        metadata={
            "basis": "full" if basis is None else "fixed_magnetization",
        },
    )


def save_analysis_result(
    result: AnalysisResult,
    path: str | Path,
    *,
    format: str | None = None,
) -> Path:
    """Save an analysis result as readable JSON or self-contained NPZ."""

    output_format = (format or Path(path).suffix.lstrip(".") or "json").lower()
    if output_format not in {"json", "npz"}:
        raise ValueError("Analysis result format must be 'json' or 'npz'.")
    output = _with_suffix(Path(path), output_format)
    output.parent.mkdir(parents=True, exist_ok=True)
    if output_format == "json":
        output.write_text(json.dumps(result.to_dict(), indent=2, sort_keys=True) + "\n")
        return output

    metadata = result.to_dict()
    metadata["values"] = {name: {"array_key": f"value__{name}"} for name in result.values}
    metadata["coordinates"] = {
        name: {"array_key": f"coordinate__{name}"} for name in result.coordinates
    }
    arrays = {f"value__{name}": value for name, value in result.values.items()}
    arrays.update({f"coordinate__{name}": value for name, value in result.coordinates.items()})
    np.savez_compressed(
        output,
        metadata=np.array(json.dumps(metadata, sort_keys=True)),
        **arrays,
    )
    return output


def load_analysis_result(path: str | Path) -> AnalysisResult:
    """Load a JSON or NPZ analysis result."""

    source = Path(path)
    if source.suffix.lower() == ".json":
        data = json.loads(source.read_text())
        if not isinstance(data, dict):
            raise ValueError("Analysis result JSON must contain an object.")
        return AnalysisResult.from_dict(data)
    if source.suffix.lower() != ".npz":
        raise ValueError("Analysis result path must end in .json or .npz.")
    with np.load(source, allow_pickle=False) as archive:
        metadata = json.loads(str(archive["metadata"].item()))
        metadata["values"] = {
            name: _array_dict(archive[record["array_key"]])
            for name, record in metadata["values"].items()
        }
        metadata["coordinates"] = {
            name: _array_dict(archive[record["array_key"]])
            for name, record in metadata["coordinates"].items()
        }
    return AnalysisResult.from_dict(metadata)


def _model_source(model: ModelSpec | None) -> dict[str, object]:
    return {} if model is None else {"kind": "model", "model": model.to_dict()}


def _array_dict(value: np.ndarray) -> dict[str, object]:
    array = np.asarray(value)
    return {"dtype": str(array.dtype), "shape": list(array.shape), "data": _json_value(array)}


def _array_mapping(value: object, name: str) -> dict[str, np.ndarray]:
    mapping = _mapping(value, name)
    return {key: _array_from_dict(record, f"{name}.{key}") for key, record in mapping.items()}


def _array_from_dict(value: object, name: str) -> np.ndarray:
    record = _mapping(value, name)
    if set(record) != {"dtype", "shape", "data"}:
        raise ValueError(f"{name} requires dtype, shape, and data.")
    array = np.asarray(_decode_value(record["data"]), dtype=str(record["dtype"]))
    shape = tuple(int(item) for item in record["shape"])
    if array.shape != shape:
        raise ValueError(f"{name} data shape does not match its declared shape.")
    return array


def _mapping(value: object, name: str) -> dict[str, object]:
    if not isinstance(value, dict) or not all(isinstance(key, str) for key in value):
        raise ValueError(f"{name} must be an object with string keys.")
    return dict(value)


def _json_value(value: object) -> object:
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, complex):
        return {"__complex__": [value.real, value.imag]}
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise ValueError(f"Value of type {type(value).__name__} is not JSON-compatible.")


def _decode_value(value: object) -> object:
    if isinstance(value, dict):
        if set(value) == {"__complex__"}:
            parts = value["__complex__"]
            if not isinstance(parts, list) or len(parts) != 2:
                raise ValueError("Complex values require [real, imaginary].")
            return complex(float(parts[0]), float(parts[1]))
        return {str(key): _decode_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return [_decode_value(item) for item in value]
    return value


def _with_suffix(path: Path, output_format: str) -> Path:
    suffix = f".{output_format}"
    return path if path.suffix.lower() == suffix else path.with_suffix(suffix)


def _package_version() -> str:
    try:
        return version("quantum-lattice-models")
    except PackageNotFoundError:
        return "0+unknown"
