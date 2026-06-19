"""Versioned portable specifications for lattices and registered models."""

from __future__ import annotations

import json
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Any

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models.lattice import Bond, Lattice

SPEC_SCHEMA_VERSION = "1.0"
_MODEL_FIELDS = {
    "schema_version",
    "family",
    "parameters",
    "lattice",
    "basis",
    "representation",
    "units",
    "conventions",
    "references",
    "provenance",
    "metadata",
}
_LATTICE_FIELDS = {
    "schema_version",
    "n_sites",
    "positions",
    "bonds",
    "site_labels",
    "orbital_labels",
    "sublattice_labels",
    "unit_cells",
    "boundary_conditions",
    "units",
    "conventions",
    "references",
    "provenance",
    "metadata",
}


@dataclass(frozen=True)
class LatticeSpec:
    """Portable finite-lattice geometry and metadata."""

    n_sites: int
    positions: tuple[tuple[float, ...], ...] = ()
    bonds: tuple[Bond, ...] = ()
    site_labels: tuple[str, ...] = ()
    orbital_labels: tuple[str, ...] = ()
    sublattice_labels: tuple[str, ...] = ()
    unit_cells: tuple[int, ...] = ()
    boundary_conditions: dict[str, str] = field(default_factory=dict)
    units: dict[str, str] = field(default_factory=dict)
    conventions: dict[str, str] = field(default_factory=dict)
    references: tuple[str, ...] = ()
    provenance: tuple[dict[str, object], ...] = ()
    metadata: dict[str, object] = field(default_factory=dict)
    schema_version: str = SPEC_SCHEMA_VERSION

    def validate(self) -> None:
        """Validate geometry, labels, bonds, and schema version."""

        if self.schema_version != SPEC_SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported lattice schema_version {self.schema_version!r}; "
                f"expected {SPEC_SCHEMA_VERSION!r}."
            )
        if not isinstance(self.n_sites, int) or self.n_sites < 1:
            raise ValueError("lattice.n_sites must be a positive integer.")
        if self.positions:
            if len(self.positions) != self.n_sites:
                raise ValueError("lattice.positions must contain one coordinate per site.")
            dimensions = {len(position) for position in self.positions}
            if dimensions not in ({2}, {3}):
                raise ValueError(
                    "lattice.positions must contain consistently 2D or 3D coordinates."
                )
            if not all(np.isfinite(value) for position in self.positions for value in position):
                raise ValueError("lattice.positions must contain only finite coordinates.")
        for bond in self.bonds:
            if not isinstance(bond.source, int) or not isinstance(bond.target, int):
                raise ValueError("lattice bond indices must be integers.")
            if not 0 <= bond.source < self.n_sites or not 0 <= bond.target < self.n_sites:
                raise ValueError("lattice bond indices must satisfy 0 <= index < n_sites.")
            if bond.value is not None and not isinstance(
                bond.value, (int, float, complex, np.number)
            ):
                raise ValueError("lattice bond values must be numeric or null.")
        for name, labels in (
            ("site_labels", self.site_labels),
            ("orbital_labels", self.orbital_labels),
            ("sublattice_labels", self.sublattice_labels),
            ("unit_cells", self.unit_cells),
        ):
            if labels and len(labels) != self.n_sites:
                raise ValueError(f"lattice.{name} must contain one value per site.")
        invalid_boundaries = {
            axis: value
            for axis, value in self.boundary_conditions.items()
            if value not in {"open", "periodic"}
        }
        if invalid_boundaries:
            raise ValueError("lattice boundary conditions must be 'open' or 'periodic'.")
        if not all(isinstance(axis, str) and axis for axis in self.boundary_conditions):
            raise ValueError("lattice boundary-condition axes must be nonempty strings.")
        _validate_string_mapping(self.units, "lattice.units")
        _validate_string_mapping(self.conventions, "lattice.conventions")
        if not all(isinstance(reference, str) for reference in self.references):
            raise ValueError("lattice.references must contain strings.")
        if not all(isinstance(record, dict) for record in self.provenance):
            raise ValueError("lattice.provenance must contain objects.")
        _validate_json_value(self.provenance, "lattice.provenance")
        _validate_json_value(self.metadata, "lattice.metadata")

    def to_lattice(self) -> Lattice:
        """Return the runtime ``Lattice`` represented by this specification."""

        self.validate()
        metadata = dict(self.metadata)
        if self.site_labels:
            metadata["site_labels"] = self.site_labels
        if self.orbital_labels:
            metadata["orbital_labels"] = self.orbital_labels
        if self.sublattice_labels:
            metadata["sublattice_labels"] = self.sublattice_labels
        if self.unit_cells:
            metadata["unit_cells"] = self.unit_cells
        if self.boundary_conditions:
            metadata["boundary_conditions"] = dict(self.boundary_conditions)
        if self.units:
            metadata["units"] = dict(self.units)
        if self.conventions:
            metadata["conventions"] = dict(self.conventions)
        if self.references:
            metadata["references"] = self.references
        if self.provenance:
            metadata["provenance"] = self.provenance
        return Lattice(
            n_sites=self.n_sites,
            positions=self.positions or None,
            bonds=self.bonds,
            metadata=metadata,
        )

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-compatible dictionary."""

        self.validate()
        return {
            "schema_version": self.schema_version,
            "n_sites": self.n_sites,
            "positions": [list(position) for position in self.positions],
            "bonds": [
                {
                    "source": bond.source,
                    "target": bond.target,
                    "value": _encode_value(bond.value),
                }
                for bond in self.bonds
            ],
            "site_labels": list(self.site_labels),
            "orbital_labels": list(self.orbital_labels),
            "sublattice_labels": list(self.sublattice_labels),
            "unit_cells": list(self.unit_cells),
            "boundary_conditions": dict(self.boundary_conditions),
            "units": dict(self.units),
            "conventions": dict(self.conventions),
            "references": list(self.references),
            "provenance": _encode_value(self.provenance),
            "metadata": _encode_value(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> LatticeSpec:
        """Construct a lattice specification from decoded JSON data."""

        data = migrate_spec_data(data, kind="lattice")
        _validate_fields(data, _LATTICE_FIELDS, "lattice")
        _require_fields(data, {"n_sites"}, "lattice")
        _require_type(data["n_sites"], int, "lattice.n_sites")
        bonds = tuple(
            _bond_from_record(record)
            for record in _record_list(data.get("bonds", []), "lattice.bonds")
        )
        spec = cls(
            schema_version=str(data.get("schema_version", "")),
            n_sites=data["n_sites"],
            positions=tuple(
                tuple(float(value) for value in position) for position in data.get("positions", [])
            ),
            bonds=bonds,
            site_labels=tuple(str(value) for value in data.get("site_labels", [])),
            orbital_labels=tuple(str(value) for value in data.get("orbital_labels", [])),
            sublattice_labels=tuple(str(value) for value in data.get("sublattice_labels", [])),
            unit_cells=tuple(int(value) for value in data.get("unit_cells", [])),
            boundary_conditions={
                str(key): str(value)
                for key, value in _mapping(data.get("boundary_conditions", {})).items()
            },
            units={str(key): str(value) for key, value in _mapping(data.get("units", {})).items()},
            conventions={
                str(key): str(value) for key, value in _mapping(data.get("conventions", {})).items()
            },
            references=tuple(str(value) for value in data.get("references", [])),
            provenance=tuple(
                _mapping(_decode_value(record)) for record in data.get("provenance", [])
            ),
            metadata=_mapping(_decode_value(data.get("metadata", {}))),
        )
        spec.validate()
        return spec


@dataclass(frozen=True)
class ModelSpec:
    """Portable registered-model parameters and optional lattice geometry."""

    family: str
    parameters: dict[str, object] = field(default_factory=dict)
    lattice: LatticeSpec | None = None
    basis: str | None = None
    representation: str = "dense"
    units: dict[str, str] = field(default_factory=dict)
    conventions: dict[str, str] = field(default_factory=dict)
    references: tuple[str, ...] = ()
    provenance: dict[str, object] = field(default_factory=dict)
    metadata: dict[str, object] = field(default_factory=dict)
    schema_version: str = SPEC_SCHEMA_VERSION

    def validate(self) -> None:
        """Validate schema, model family, parameters, and lattice requirements."""

        from quantum_lattice_models.registry import get_model_info, validate_parameter

        if self.schema_version != SPEC_SCHEMA_VERSION:
            raise ValueError(
                f"Unsupported model schema_version {self.schema_version!r}; "
                f"expected {SPEC_SCHEMA_VERSION!r}."
            )
        if self.representation not in {"dense", "sparse"}:
            raise ValueError("model.representation must be 'dense' or 'sparse'.")
        if not isinstance(self.family, str) or not self.family:
            raise ValueError("model.family must be a nonempty string.")
        if not isinstance(self.parameters, dict) or not all(
            isinstance(name, str) for name in self.parameters
        ):
            raise ValueError("model.parameters must be an object with string keys.")
        try:
            info = get_model_info(self.family)
        except KeyError as exc:
            raise ValueError(f"Unknown model family {self.family!r}.") from exc
        accepted = {parameter.name: parameter for parameter in info.parameters}
        if self.basis is not None and self.basis != info.basis:
            raise ValueError(
                f"model.basis must match registered basis {info.basis!r} for {self.family!r}."
            )
        unsupported = sorted(set(self.parameters) - set(accepted))
        if unsupported:
            raise ValueError(
                f"Model {self.family!r} does not accept parameters: {', '.join(unsupported)}."
            )
        values = dict(info.defaults)
        values.update(self.parameters)
        for name, parameter in accepted.items():
            validate_parameter(parameter, values.get(name))
        if self.lattice is not None:
            self.lattice.validate()
        if _base_family(self.family) == "custom_tight_binding" and self.lattice is None:
            if "n_sites" not in values and "bonds" not in values:
                raise ValueError("Custom tight-binding specifications require a lattice.")
        _resolve_builder_name(self.family, self.representation)
        _validate_string_mapping(self.units, "model.units")
        _validate_string_mapping(self.conventions, "model.conventions")
        if not all(isinstance(reference, str) for reference in self.references):
            raise ValueError("model.references must contain strings.")
        _validate_json_value(self.parameters, "model.parameters")
        _validate_json_value(self.provenance, "model.provenance")
        _validate_json_value(self.metadata, "model.metadata")

    def hamiltonian(self, *, sparse: bool | None = None) -> np.ndarray | sp.csr_matrix:
        """Build the represented Hamiltonian in dense or sparse form."""

        representation = (
            self.representation if sparse is None else ("sparse" if sparse else "dense")
        )
        built = self._build_value(representation)
        return getattr(built, "matrix", built)

    def build_result(self, *, sparse: bool | None = None):
        """Build a Hamiltonian with its model, basis, and construction metadata."""

        from quantum_lattice_models.types import HamiltonianResult

        representation = (
            self.representation if sparse is None else ("sparse" if sparse else "dense")
        )
        built = self._build_value(representation)
        matrix = getattr(built, "matrix", built)
        result_model = (
            self
            if representation == self.representation
            else replace(self, representation=representation)
        )
        matrix_metadata = dict(getattr(matrix, "metadata", {}))
        if built is not matrix and hasattr(built, "to_metadata"):
            matrix_metadata.update(built.to_metadata())
        for name in ("model_name", "n_sites", "lattice_shape"):
            value = getattr(matrix, name, None)
            if value is not None:
                matrix_metadata[name] = value
        terms = getattr(matrix, "terms", ())
        if terms:
            matrix_metadata["pauli_terms"] = [
                {
                    "coefficient": _encode_value(term.coefficient),
                    "operators": list(term.operators),
                }
                for term in terms
            ]
        return HamiltonianResult(
            matrix=matrix,
            model=result_model,
            basis=self.basis or "",
            representation=representation,
            metadata=matrix_metadata,
        )

    def _build_value(self, representation: str) -> object:
        from quantum_lattice_models.registry import get_model_info

        self.validate()
        builder_name = _resolve_builder_name(self.family, representation)
        info = get_model_info(builder_name)
        if info.builder is None:
            raise ValueError(f"Registered model {builder_name!r} does not define a builder.")
        kwargs = dict(info.defaults)
        kwargs.update(self.parameters)
        if self.lattice is not None and _base_family(builder_name) == "custom_tight_binding":
            kwargs.pop("n_sites", None)
            kwargs.pop("bonds", None)
            kwargs["lattice"] = self.lattice.to_lattice()
        return info.builder(**kwargs)

    def to_dict(self) -> dict[str, object]:
        """Return a JSON-compatible dictionary."""

        self.validate()
        return {
            "schema_version": self.schema_version,
            "family": self.family,
            "parameters": _encode_value(self.parameters),
            "lattice": None if self.lattice is None else self.lattice.to_dict(),
            "basis": self.basis,
            "representation": self.representation,
            "units": dict(self.units),
            "conventions": dict(self.conventions),
            "references": list(self.references),
            "provenance": _encode_value(self.provenance),
            "metadata": _encode_value(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> ModelSpec:
        """Construct a model specification from decoded JSON data."""

        data = migrate_spec_data(data, kind="model")
        _validate_fields(data, _MODEL_FIELDS, "model")
        _require_fields(data, {"family"}, "model")
        _require_type(data["family"], str, "model.family")
        _require_type(data.get("parameters", {}), dict, "model.parameters")
        lattice_data = data.get("lattice")
        if lattice_data is not None:
            _require_type(lattice_data, dict, "model.lattice")
        spec = cls(
            schema_version=str(data.get("schema_version", "")),
            family=str(data["family"]),
            parameters=_mapping(_decode_value(data.get("parameters", {}))),
            lattice=(
                None if lattice_data is None else LatticeSpec.from_dict(_mapping(lattice_data))
            ),
            basis=None if data.get("basis") is None else str(data["basis"]),
            representation=str(data.get("representation", "dense")),
            units={str(key): str(value) for key, value in _mapping(data.get("units", {})).items()},
            conventions={
                str(key): str(value) for key, value in _mapping(data.get("conventions", {})).items()
            },
            references=tuple(str(value) for value in data.get("references", [])),
            provenance=_mapping(_decode_value(data.get("provenance", {}))),
            metadata=_mapping(_decode_value(data.get("metadata", {}))),
        )
        spec.validate()
        return spec

    def save(self, path: str | Path) -> Path:
        """Save this specification as canonical indented JSON."""

        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n")
        return output

    def summary(self) -> dict[str, object]:
        """Return a deterministic, human-readable inspection summary."""

        from quantum_lattice_models.diagnostics import (
            estimate_dense_memory,
            estimate_model_dimension,
        )

        summary: dict[str, object] = {
            "schema_version": self.schema_version,
            "family": self.family,
            "basis": self.basis,
            "representation": self.representation,
            "parameters": dict(sorted(self.parameters.items())),
            "has_lattice": self.lattice is not None,
            "units": dict(sorted(self.units.items())),
            "conventions": dict(sorted(self.conventions.items())),
            "references": list(self.references),
            "provenance": dict(self.provenance),
            "metadata": dict(self.metadata),
        }
        if self.lattice is not None:
            summary["lattice"] = {
                "n_sites": self.lattice.n_sites,
                "n_bonds": len(self.lattice.bonds),
                "coordinate_dimension": (
                    len(self.lattice.positions[0]) if self.lattice.positions else None
                ),
                "units": dict(sorted(self.lattice.units.items())),
                "conventions": dict(sorted(self.lattice.conventions.items())),
                "references": list(self.lattice.references),
                "provenance": list(self.lattice.provenance),
            }
        dimension_parameters = dict(self.parameters)
        if self.lattice is not None:
            dimension_parameters.setdefault("n_sites", self.lattice.n_sites)
        try:
            dimension = estimate_model_dimension(self.family, **dimension_parameters)
        except ValueError:
            dimension = None
        if dimension is not None:
            summary["dimension"] = dimension
            summary["dense_memory_bytes"] = estimate_dense_memory(dimension)
        return summary


def create_model_spec(
    family: str,
    *,
    parameters: dict[str, object] | None = None,
    lattice: LatticeSpec | None = None,
    representation: str | None = None,
    units: dict[str, str] | None = None,
    conventions: dict[str, str] | None = None,
    references: tuple[str, ...] = (),
    provenance: dict[str, object] | None = None,
    metadata: dict[str, object] | None = None,
) -> ModelSpec:
    """Create a validated specification using registered defaults."""

    from quantum_lattice_models.registry import get_model_info

    info = get_model_info(family)
    resolved = {
        parameter.name: parameter.default for parameter in info.parameters if not parameter.required
    }
    resolved.update(info.defaults)
    resolved.update(parameters or {})
    if lattice is not None and _base_family(family) == "custom_tight_binding":
        resolved.pop("n_sites", None)
        resolved.pop("bonds", None)
    inferred_representation = representation or (
        "sparse" if family.endswith("_sparse") else "dense"
    )
    spec = ModelSpec(
        family=family,
        parameters=resolved,
        lattice=lattice,
        basis=info.basis,
        representation=inferred_representation,
        units=dict(units or {}),
        conventions=dict(conventions or {}),
        references=tuple(references),
        provenance=dict(provenance or {}),
        metadata=dict(metadata or {}),
    )
    spec.validate()
    return spec


def create_model_from_preset(
    preset_name: str,
    *,
    parameters: dict[str, object] | None = None,
    representation: str | None = None,
) -> ModelSpec:
    """Create a model specification from a named preset with optional overrides."""

    from quantum_lattice_models.registry import get_preset

    preset = get_preset(preset_name)
    resolved = dict(preset.parameters)
    resolved.update(parameters or {})
    return create_model_spec(
        preset.model,
        parameters=resolved,
        representation=representation,
        provenance={"preset": preset.name},
        metadata={"preset_description": preset.description},
    )


def load_model(path: str | Path) -> ModelSpec:
    """Load and validate a model specification from JSON."""

    data = json.loads(Path(path).read_text())
    if not isinstance(data, dict):
        raise ValueError("Model specification JSON must contain an object.")
    return ModelSpec.from_dict(data)


def migrate_spec_data(
    data: dict[str, object],
    *,
    kind: str,
) -> dict[str, object]:
    """Migrate decoded model or lattice data to the current schema.

    Unversioned files from the pre-schema prototype are treated as legacy
    version ``0`` and receive the current version marker. Unknown past and all
    future versions are rejected until an explicit migration is implemented.
    """

    if kind not in {"model", "lattice"}:
        raise ValueError("kind must be 'model' or 'lattice'.")
    if not isinstance(data, dict):
        raise ValueError(f"{kind} specification must contain an object.")
    migrated = dict(data)
    version = migrated.get("schema_version")
    if version is None:
        migrated["schema_version"] = SPEC_SCHEMA_VERSION
        return migrated
    if not isinstance(version, str):
        raise ValueError(f"{kind}.schema_version must be a string.")
    if version == SPEC_SCHEMA_VERSION:
        return migrated
    raise ValueError(
        f"Unsupported {kind} schema_version {version!r}; no migration to "
        f"{SPEC_SCHEMA_VERSION!r} is registered."
    )


def _resolve_builder_name(family: str, representation: str) -> str:
    from quantum_lattice_models.registry import MODEL_REGISTRY

    base = _base_family(family)
    candidate = f"{base}_sparse" if representation == "sparse" else base
    if candidate not in MODEL_REGISTRY:
        raise ValueError(f"Model {base!r} does not provide a {representation} builder.")
    return candidate


def _base_family(family: str) -> str:
    return family.removesuffix("_sparse")


def _encode_value(value: object) -> object:
    if isinstance(value, complex):
        return {"__complex__": [value.real, value.imag]}
    if isinstance(value, np.generic):
        return _encode_value(value.item())
    if isinstance(value, dict):
        return {str(key): _encode_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_encode_value(item) for item in value]
    return value


def _decode_value(value: object) -> Any:
    if isinstance(value, dict):
        if set(value) == {"__complex__"}:
            parts = value["__complex__"]
            if not isinstance(parts, list) or len(parts) != 2:
                raise ValueError("Complex values must contain [real, imaginary].")
            return complex(float(parts[0]), float(parts[1]))
        return {str(key): _decode_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return tuple(_decode_value(item) for item in value)
    return value


def _mapping(value: object) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError("Expected a JSON object.")
    return dict(value)


def _record_list(value: object, field_name: str) -> list[dict[str, Any]]:
    if not isinstance(value, list) or not all(isinstance(item, dict) for item in value):
        raise ValueError(f"{field_name} must be a list of objects.")
    return value


def _validate_fields(data: dict[str, object], allowed: set[str], field_name: str) -> None:
    unknown = sorted(set(data) - allowed)
    if unknown:
        raise ValueError(f"{field_name} contains unknown fields: {', '.join(unknown)}.")


def _require_fields(data: dict[str, object], required: set[str], field_name: str) -> None:
    missing = sorted(required - set(data))
    if missing:
        raise ValueError(f"{field_name} is missing required fields: {', '.join(missing)}.")


def _require_type(value: object, expected: type, field_name: str) -> None:
    if expected is int:
        valid = isinstance(value, int) and not isinstance(value, bool)
    else:
        valid = isinstance(value, expected)
    if not valid:
        raise ValueError(f"{field_name} must have type {expected.__name__}.")


def _required_int(record: dict[str, Any], name: str, field_name: str) -> int:
    if name not in record:
        raise ValueError(f"{field_name} entries require {name!r}.")
    _require_type(record[name], int, f"{field_name}.{name}")
    return record[name]


def _bond_from_record(record: dict[str, Any]) -> Bond:
    _validate_fields(record, {"source", "target", "value"}, "lattice.bonds entry")
    return Bond(
        source=_required_int(record, "source", "lattice.bonds"),
        target=_required_int(record, "target", "lattice.bonds"),
        value=_decode_value(record.get("value")),
    )


def _validate_json_value(value: object, field_name: str) -> None:
    try:
        json.dumps(_encode_value(value))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must contain portable JSON values.") from exc


def _validate_string_mapping(value: object, field_name: str) -> None:
    if not isinstance(value, dict) or not all(
        isinstance(key, str) and isinstance(item, str) for key, item in value.items()
    ):
        raise ValueError(f"{field_name} must map strings to strings.")
