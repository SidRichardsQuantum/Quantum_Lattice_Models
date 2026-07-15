"""Versioned portable specifications for lattices and registered models."""

from __future__ import annotations

import json
from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models._schema_codec import (
    bond_from_record as _bond_from_record,
)
from quantum_lattice_models._schema_codec import (
    decode_value as _decode_value,
)
from quantum_lattice_models._schema_codec import (
    encode_value as _encode_value,
)
from quantum_lattice_models._schema_codec import (
    mapping as _mapping,
)
from quantum_lattice_models._schema_codec import (
    record_list as _record_list,
)
from quantum_lattice_models._schema_codec import (
    require_fields as _require_fields,
)
from quantum_lattice_models._schema_codec import (
    require_type as _require_type,
)
from quantum_lattice_models._schema_codec import (
    validate_fields as _validate_fields,
)
from quantum_lattice_models._schema_codec import (
    validate_json_value as _validate_json_value,
)
from quantum_lattice_models._schema_codec import (
    validate_string_mapping as _validate_string_mapping,
)
from quantum_lattice_models.lattice import Bond, Lattice
from quantum_lattice_models.physical import (
    BasisIndexMapping,
    InteractionTerm,
    LocalDegreeOfFreedom,
    validate_physical_system,
)

SPEC_SCHEMA_VERSION = "1.0"
EXTERNAL_MATRIX_FAMILY = "external_matrix"
GRAPH_SPIN_FAMILY = "graph_spin"
_MODEL_FIELDS = {
    "schema_version",
    "family",
    "parameters",
    "lattice",
    "basis",
    "local_degrees",
    "basis_mappings",
    "interactions",
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
    local_degrees: tuple[LocalDegreeOfFreedom, ...] = ()
    basis_mappings: tuple[BasisIndexMapping, ...] = ()
    interactions: tuple[InteractionTerm, ...] = ()
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
        validate_physical_system(
            self.local_degrees,
            self.basis_mappings,
            self.interactions,
            n_sites=None if self.lattice is None else self.lattice.n_sites,
        )
        if self.family == EXTERNAL_MATRIX_FAMILY:
            if self.parameters:
                raise ValueError("External matrix specifications do not accept parameters.")
            if not isinstance(self.basis, str) or not self.basis.strip():
                raise ValueError("External matrix specifications require a nonempty basis.")
            if self.lattice is not None:
                self.lattice.validate()
            _validate_string_mapping(self.units, "model.units")
            _validate_string_mapping(self.conventions, "model.conventions")
            if not all(isinstance(reference, str) for reference in self.references):
                raise ValueError("model.references must contain strings.")
            _validate_json_value(self.provenance, "model.provenance")
            _validate_json_value(self.metadata, "model.metadata")
            return
        if self.family == GRAPH_SPIN_FAMILY:
            if set(self.parameters) != {"n_sites"}:
                raise ValueError("Graph-spin specifications require only the n_sites parameter.")
            n_sites = self.parameters["n_sites"]
            if not isinstance(n_sites, int) or isinstance(n_sites, bool) or n_sites < 1:
                raise ValueError("Graph-spin n_sites must be a positive integer.")
            if self.basis not in (None, "qubit"):
                raise ValueError("Graph-spin specifications require the 'qubit' basis.")
            if self.lattice is not None:
                self.lattice.validate()
                if self.lattice.n_sites != n_sites:
                    raise ValueError("Graph-spin lattice.n_sites must equal parameters.n_sites.")
            if len(self.local_degrees) != n_sites:
                raise ValueError("Graph-spin specifications require one local spin per site.")
            if any(degree.kind != "spin" for degree in self.local_degrees):
                raise ValueError("Graph-spin local degrees must have kind 'spin'.")
            _validate_string_mapping(self.units, "model.units")
            _validate_string_mapping(self.conventions, "model.conventions")
            if not all(isinstance(reference, str) for reference in self.references):
                raise ValueError("model.references must contain strings.")
            _validate_json_value(self.provenance, "model.provenance")
            _validate_json_value(self.metadata, "model.metadata")
            return
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
        if self.family == EXTERNAL_MATRIX_FAMILY:
            raise ValueError(
                "External matrix specifications cannot reconstruct a matrix; "
                "load the persisted Hamiltonian file instead."
            )
        if self.family == GRAPH_SPIN_FAMILY:
            from quantum_lattice_models.spin import (
                SpinField,
                SpinInteraction,
                graph_spin_hamiltonian,
                graph_spin_hamiltonian_sparse,
            )

            graph_interactions = []
            graph_fields = []
            for term in self.interactions:
                if len(term.degrees) == 1:
                    graph_fields.append(
                        SpinField(term.degrees[0], term.operators[0], term.coefficient)
                    )
                else:
                    graph_interactions.append(
                        SpinInteraction(
                            term.degrees[0],
                            term.degrees[1],
                            term.operators[0],
                            term.operators[1],
                            term.coefficient,
                        )
                    )
            builder = (
                graph_spin_hamiltonian_sparse
                if representation == "sparse"
                else graph_spin_hamiltonian
            )
            return builder(
                int(self.parameters["n_sites"]),
                graph_interactions,
                graph_fields,
            )
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
            "local_degrees": [degree.to_dict() for degree in self.local_degrees],
            "basis_mappings": [mapping.to_dict() for mapping in self.basis_mappings],
            "interactions": [interaction.to_dict() for interaction in self.interactions],
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
            local_degrees=tuple(
                LocalDegreeOfFreedom.from_dict(record)
                for record in _record_list(data.get("local_degrees", []), "model.local_degrees")
            ),
            basis_mappings=tuple(
                BasisIndexMapping.from_dict(record)
                for record in _record_list(data.get("basis_mappings", []), "model.basis_mappings")
            ),
            interactions=tuple(
                InteractionTerm.from_dict(record)
                for record in _record_list(data.get("interactions", []), "model.interactions")
            ),
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
            "local_degree_count": len(self.local_degrees),
            "interaction_count": len(self.interactions),
            "basis_mapping_count": len(self.basis_mappings),
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
        if self.family == EXTERNAL_MATRIX_FAMILY:
            dimension = self.metadata.get("matrix_dimension")
            if not isinstance(dimension, int) or dimension < 1:
                dimension = None
        elif self.family == GRAPH_SPIN_FAMILY:
            dimension = 2 ** int(self.parameters["n_sites"])
        else:
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
    local_degrees: tuple[LocalDegreeOfFreedom, ...] | None = None,
    basis_mappings: tuple[BasisIndexMapping, ...] | None = None,
    interactions: tuple[InteractionTerm, ...] | None = None,
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
    inferred_lattice, inferred_degrees, inferred_mappings, inferred_interactions = (
        _infer_physical_system(family, resolved, lattice)
    )
    spec = ModelSpec(
        family=family,
        parameters=resolved,
        lattice=inferred_lattice,
        basis=info.basis,
        local_degrees=(inferred_degrees if local_degrees is None else tuple(local_degrees)),
        basis_mappings=(inferred_mappings if basis_mappings is None else tuple(basis_mappings)),
        interactions=(inferred_interactions if interactions is None else tuple(interactions)),
        representation=inferred_representation,
        units=dict(units or {}),
        conventions=dict(conventions or {}),
        references=tuple(references),
        provenance=dict(provenance or {}),
        metadata=dict(metadata or {}),
    )
    spec.validate()
    return spec


def create_graph_spin_spec(
    n_sites: int,
    *,
    interactions: tuple[object, ...] = (),
    fields: tuple[object, ...] = (),
    lattice: LatticeSpec | None = None,
    positions: tuple[tuple[float, ...], ...] | None = None,
    site_labels: tuple[str, ...] | None = None,
    representation: str = "dense",
    units: dict[str, str] | None = None,
    conventions: dict[str, str] | None = None,
    references: tuple[str, ...] = (),
    provenance: dict[str, object] | None = None,
    metadata: dict[str, object] | None = None,
) -> ModelSpec:
    """Create a portable arbitrary spin-1/2 interaction-graph specification."""

    from quantum_lattice_models.spin import SpinField, SpinInteraction

    if not isinstance(n_sites, int) or isinstance(n_sites, bool) or n_sites < 1:
        raise ValueError("n_sites must be a positive integer.")
    if representation not in {"dense", "sparse"}:
        raise ValueError("representation must be 'dense' or 'sparse'.")
    if lattice is not None and (positions is not None or site_labels is not None):
        raise ValueError("Use either lattice or positions/site_labels, not both.")
    portable_terms: list[InteractionTerm] = []
    bonds: list[Bond] = []
    for interaction in interactions:
        if not isinstance(interaction, SpinInteraction):
            raise TypeError("interactions must contain SpinInteraction records.")
        portable_terms.append(
            InteractionTerm(
                (interaction.source, interaction.target),
                (interaction.source_axis, interaction.target_axis),
                interaction.coefficient,
                "spin_exchange",
                label=(
                    f"{interaction.source_axis}{interaction.target_axis} "
                    f"{interaction.source}-{interaction.target}"
                ),
            )
        )
        bond = Bond(interaction.source, interaction.target)
        if bond not in bonds:
            bonds.append(bond)
    for spin_field in fields:
        if not isinstance(spin_field, SpinField):
            raise TypeError("fields must contain SpinField records.")
        portable_terms.append(
            InteractionTerm(
                (spin_field.site,),
                (spin_field.axis,),
                spin_field.coefficient,
                "field",
                label=f"{spin_field.axis} field {spin_field.site}",
            )
        )
    if lattice is None:
        resolved_positions = (
            tuple(positions)
            if positions is not None
            else tuple((float(site), 0.0) for site in range(n_sites))
        )
        resolved_labels = (
            tuple(site_labels)
            if site_labels is not None
            else tuple(f"spin {site}" for site in range(n_sites))
        )
        lattice = LatticeSpec(
            n_sites=n_sites,
            positions=resolved_positions,
            bonds=tuple(bonds),
            site_labels=resolved_labels,
            conventions={"edge_interpretation": "spin interaction graph"},
        )
    degrees = tuple(
        LocalDegreeOfFreedom(
            site,
            site,
            "spin",
            2,
            lattice.site_labels[site] if lattice.site_labels else f"spin {site}",
            component="spin-1/2",
        )
        for site in range(n_sites)
    )
    spec = ModelSpec(
        family=GRAPH_SPIN_FAMILY,
        parameters={"n_sites": n_sites},
        lattice=lattice,
        basis="qubit",
        local_degrees=degrees,
        basis_mappings=tuple(
            BasisIndexMapping(site, site, "tensor_factor", f"qubit {site}")
            for site in range(n_sites)
        ),
        interactions=tuple(portable_terms),
        representation=representation,
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


def _infer_physical_system(
    family: str,
    parameters: dict[str, object],
    lattice: LatticeSpec | None,
) -> tuple[
    LatticeSpec | None,
    tuple[LocalDegreeOfFreedom, ...],
    tuple[BasisIndexMapping, ...],
    tuple[InteractionTerm, ...],
]:
    base = _base_family(family)
    if base in {
        "transverse_field_ising",
        "longitudinal_field_ising",
        "next_nearest_neighbor_ising",
        "heisenberg_chain",
        "heisenberg_chain_sector",
        "xxz_chain",
        "xxz_chain_sector",
        "xy_chain",
        "j1_j2_heisenberg_chain",
        "xyz_chain",
        "random_field_heisenberg_chain",
    }:
        return _spin_chain_physical_system(base, parameters, lattice)
    if base in {"heisenberg_ladder", "heisenberg_ladder_sector"}:
        return _spin_ladder_physical_system(parameters, lattice)
    if base == "ssh_model":
        return _ssh_physical_system(parameters, lattice)
    if base in {
        "anderson_chain",
        "aubry_andre_harper_chain",
        "long_range_tight_binding_chain",
        "rice_mele_model",
        "sawtooth_chain",
        "creutz_ladder",
        "tight_binding_chain",
    }:
        return _chain_particle_physical_system(base, parameters, lattice)
    if base in {
        "anderson_square_lattice",
        "checkerboard_chern_insulator",
        "dice_lattice",
        "graphene_lattice",
        "haldane_honeycomb_lattice",
        "harper_hofstadter_square_lattice",
        "kagome_lattice_tight_binding",
        "lieb_lattice",
        "square_lattice_tight_binding",
        "triangular_lattice_tight_binding",
    }:
        return _planar_particle_physical_system(base, parameters, lattice)
    if base == "bose_hubbard_chain":
        return _bose_hubbard_physical_system(parameters, lattice)
    if base in {"fermi_hubbard_chain", "fermi_hubbard_chain_sector"}:
        return _fermi_hubbard_physical_system(parameters, lattice)
    if base == "kitaev_chain_bdg":
        return _kitaev_bdg_physical_system(parameters, lattice)
    if base == "custom_tight_binding":
        return _custom_tight_binding_physical_system(
            parameters,
            lattice or _custom_lattice_from_parameters(parameters),
        )
    return lattice, (), (), ()


def _spin_chain_physical_system(
    family: str,
    parameters: dict[str, object],
    lattice: LatticeSpec | None,
) -> tuple[
    LatticeSpec,
    tuple[LocalDegreeOfFreedom, ...],
    tuple[BasisIndexMapping, ...],
    tuple[InteractionTerm, ...],
]:
    n_sites = int(parameters["n_sites"])
    periodic = bool(parameters.get("periodic", False))
    nearest = _chain_bonds(n_sites, periodic)
    if lattice is None:
        lattice = LatticeSpec(
            n_sites=n_sites,
            positions=tuple((float(site), 0.0) for site in range(n_sites)),
            bonds=tuple(Bond(source, target) for source, target in nearest),
            site_labels=tuple(f"spin {site}" for site in range(n_sites)),
            boundary_conditions={"x": "periodic" if periodic else "open"},
            conventions={"site_ordering": "left-to-right chain order"},
        )
    degrees = tuple(
        LocalDegreeOfFreedom(
            index=site,
            site=site,
            kind="spin",
            local_dimension=2,
            label=f"spin {site}",
            component="spin-1/2",
        )
        for site in range(n_sites)
    )
    mappings = tuple(
        BasisIndexMapping(site, site, "tensor_factor", f"qubit {site}") for site in range(n_sites)
    )
    interactions: list[InteractionTerm] = []

    def pair_terms(
        bonds: tuple[tuple[int, int], ...],
        couplings: tuple[tuple[str, complex], ...],
        kind: str,
    ) -> None:
        for source, target in bonds:
            for axis, coefficient in couplings:
                interactions.append(
                    InteractionTerm(
                        (source, target),
                        (axis, axis),
                        coefficient,
                        kind,
                        label=f"{axis}{axis} {source}-{target}",
                    )
                )

    def fields(axis: str, coefficient: complex) -> None:
        for site in range(n_sites):
            interactions.append(
                InteractionTerm(
                    (site,),
                    (axis,),
                    coefficient,
                    "field",
                    label=f"{axis} field {site}",
                )
            )

    if family == "transverse_field_ising":
        pair_terms(nearest, (("Z", -complex(parameters["j"])),), "ising_exchange")
        fields("X", -complex(parameters["h"]))
    elif family == "longitudinal_field_ising":
        pair_terms(nearest, (("Z", -complex(parameters["j"])),), "ising_exchange")
        fields("X", -complex(parameters["h_x"]))
        fields("Z", -complex(parameters["h_z"]))
    elif family == "next_nearest_neighbor_ising":
        pair_terms(nearest, (("Z", -complex(parameters["j1"])),), "nearest_exchange")
        pair_terms(
            _next_nearest_chain_bonds(n_sites, periodic),
            (("Z", -complex(parameters["j2"])),),
            "next_nearest_exchange",
        )
        fields("X", -complex(parameters["h"]))
    elif family in {"heisenberg_chain", "heisenberg_chain_sector", "xyz_chain"}:
        pair_terms(
            nearest,
            (
                ("X", complex(parameters["jx"])),
                ("Y", complex(parameters["jy"])),
                ("Z", complex(parameters["jz"])),
            ),
            "heisenberg_exchange",
        )
        fields("Z", complex(parameters.get("field", 0.0)))
    elif family == "xy_chain":
        coupling = complex(parameters["coupling"])
        anisotropy = complex(parameters["anisotropy"])
        pair_terms(
            nearest,
            (
                ("X", -coupling * (1.0 + anisotropy) / 2.0),
                ("Y", -coupling * (1.0 - anisotropy) / 2.0),
            ),
            "xy_exchange",
        )
        fields("Z", -complex(parameters.get("field", 0.0)))
    elif family == "j1_j2_heisenberg_chain":
        pair_terms(
            nearest,
            tuple((axis, complex(parameters["j1"])) for axis in ("X", "Y", "Z")),
            "nearest_exchange",
        )
        pair_terms(
            _next_nearest_chain_bonds(n_sites, periodic),
            tuple((axis, complex(parameters["j2"])) for axis in ("X", "Y", "Z")),
            "next_nearest_exchange",
        )
        fields("Z", complex(parameters.get("field", 0.0)))
    elif family == "random_field_heisenberg_chain":
        pair_terms(
            nearest,
            tuple((axis, complex(parameters["coupling"])) for axis in ("X", "Y", "Z")),
            "heisenberg_exchange",
        )
        rng = np.random.default_rng(int(parameters["seed"]))
        for site, coefficient in enumerate(
            rng.uniform(
                -float(parameters["disorder"]),
                float(parameters["disorder"]),
                n_sites,
            )
        ):
            interactions.append(
                InteractionTerm(
                    (site,),
                    ("Z",),
                    complex(coefficient),
                    "random_field",
                    label=f"Z field {site}",
                    metadata={"seed": int(parameters["seed"])},
                )
            )
    else:
        coupling = complex(parameters["coupling"])
        pair_terms(
            nearest,
            (
                ("X", coupling),
                ("Y", coupling),
                ("Z", coupling * complex(parameters["anisotropy"])),
            ),
            "xxz_exchange",
        )
        fields("Z", complex(parameters.get("field", 0.0)))
    return lattice, degrees, mappings, tuple(interactions)


def _spin_ladder_physical_system(
    parameters: dict[str, object],
    lattice: LatticeSpec | None,
) -> tuple[
    LatticeSpec,
    tuple[LocalDegreeOfFreedom, ...],
    tuple[BasisIndexMapping, ...],
    tuple[InteractionTerm, ...],
]:
    n_rungs = int(parameters["n_rungs"])
    n_sites = 2 * n_rungs
    periodic = bool(parameters.get("periodic", False))
    rung_bonds = tuple((2 * rung, 2 * rung + 1) for rung in range(n_rungs))
    leg_bonds = tuple(
        (2 * rung + leg, 2 * next_rung + leg)
        for leg in range(2)
        for rung, next_rung in _chain_bonds(n_rungs, periodic)
    )
    all_bonds = rung_bonds + leg_bonds
    if lattice is None:
        lattice = LatticeSpec(
            n_sites=n_sites,
            positions=tuple(
                (float(rung), float(leg)) for rung in range(n_rungs) for leg in range(2)
            ),
            bonds=tuple(Bond(source, target) for source, target in all_bonds),
            site_labels=tuple(
                f"rung {rung} leg {leg}" for rung in range(n_rungs) for leg in range(2)
            ),
            unit_cells=tuple(rung for rung in range(n_rungs) for _ in range(2)),
            boundary_conditions={"x": "periodic" if periodic else "open", "y": "open"},
            conventions={"site_ordering": "rung-major, leg-minor"},
        )
    degrees = tuple(
        LocalDegreeOfFreedom(site, site, "spin", 2, lattice.site_labels[site], component="spin-1/2")
        for site in range(n_sites)
    )
    mappings = tuple(
        BasisIndexMapping(site, site, "tensor_factor", f"qubit {site}") for site in range(n_sites)
    )
    interactions = tuple(
        InteractionTerm(
            (source, target),
            (axis, axis),
            complex(coupling),
            kind,
            label=f"{axis}{axis} {source}-{target}",
        )
        for bonds, coupling, kind in (
            (leg_bonds, parameters["leg_coupling"], "leg_exchange"),
            (rung_bonds, parameters["rung_coupling"], "rung_exchange"),
        )
        for source, target in bonds
        for axis in ("X", "Y", "Z")
    ) + tuple(
        InteractionTerm(
            (site,),
            ("Z",),
            complex(parameters.get("field", 0.0)),
            "field",
            label=f"Z field {site}",
        )
        for site in range(n_sites)
    )
    return lattice, degrees, mappings, interactions


def _ssh_physical_system(
    parameters: dict[str, object],
    lattice: LatticeSpec | None,
) -> tuple[
    LatticeSpec,
    tuple[LocalDegreeOfFreedom, ...],
    tuple[BasisIndexMapping, ...],
    tuple[InteractionTerm, ...],
]:
    n_cells = int(parameters["n_cells"])
    periodic = bool(parameters.get("periodic", False))
    t1 = complex(parameters["t1"])
    t2 = complex(parameters["t2"])
    bonds: list[Bond] = []
    interactions: list[InteractionTerm] = []
    for cell in range(n_cells):
        a = 2 * cell
        b = a + 1
        bonds.append(Bond(a, b, -t1))
        interactions.append(
            _hopping_interaction(a, b, -t1, "intracell_hopping", f"cell {cell} A-B")
        )
        if cell < n_cells - 1:
            target = 2 * (cell + 1)
            bonds.append(Bond(b, target, -t2))
            interactions.append(
                _hopping_interaction(b, target, -t2, "intercell_hopping", f"cell {cell} B-A")
            )
        elif periodic and n_cells > 1:
            bonds.append(Bond(b, 0, -t2))
            interactions.append(
                _hopping_interaction(b, 0, -t2, "intercell_hopping", "periodic B-A")
            )
    if lattice is None:
        lattice = LatticeSpec(
            n_sites=2 * n_cells,
            positions=tuple(
                (float(cell) + (0.0 if sublattice == 0 else 0.45), 0.0)
                for cell in range(n_cells)
                for sublattice in range(2)
            ),
            bonds=tuple(bonds),
            site_labels=tuple(
                f"{'A' if sublattice == 0 else 'B'}{cell}"
                for cell in range(n_cells)
                for sublattice in range(2)
            ),
            orbital_labels=tuple("orbital" for _ in range(2 * n_cells)),
            sublattice_labels=tuple("A" if index % 2 == 0 else "B" for index in range(2 * n_cells)),
            unit_cells=tuple(cell for cell in range(n_cells) for _ in range(2)),
            boundary_conditions={"x": "periodic" if periodic else "open"},
            conventions={"basis_ordering": "A, B within each unit cell"},
        )
    degrees = _orbital_degrees(lattice)
    mappings = _single_particle_mappings(degrees)
    return lattice, degrees, mappings, tuple(interactions)


def _chain_particle_physical_system(
    family: str,
    parameters: dict[str, object],
    lattice: LatticeSpec | None,
) -> tuple[
    LatticeSpec,
    tuple[LocalDegreeOfFreedom, ...],
    tuple[BasisIndexMapping, ...],
    tuple[InteractionTerm, ...],
]:
    periodic = bool(parameters.get("periodic", False))
    hoppings: dict[tuple[int, int], complex] = {}
    onsite: dict[int, complex] = {}

    if family in {"rice_mele_model", "sawtooth_chain", "creutz_ladder"}:
        n_cells = int(parameters["n_cells"])
        n_sites = 2 * n_cells
        positions = tuple(
            (float(cell), 0.35 if orbital else 0.0)
            for cell in range(n_cells)
            for orbital in range(2)
        )
        unit_cells = tuple(cell for cell in range(n_cells) for _ in range(2))
        if family == "rice_mele_model":
            sublattices = tuple(label for _ in range(n_cells) for label in ("A", "B"))
            hopping = complex(parameters["hopping"])
            dimerization = complex(parameters["dimerization"])
            staggering = complex(parameters["staggering"])
            for cell in range(n_cells):
                a, b = 2 * cell, 2 * cell + 1
                onsite[a], onsite[b] = staggering, -staggering
                _accumulate_hopping(hoppings, onsite, a, b, -(hopping + dimerization))
                if cell < n_cells - 1:
                    _accumulate_hopping(
                        hoppings, onsite, b, 2 * (cell + 1), -(hopping - dimerization)
                    )
                elif periodic and n_cells > 1:
                    _accumulate_hopping(hoppings, onsite, b, 0, -(hopping - dimerization))
        elif family == "sawtooth_chain":
            sublattices = tuple(label for _ in range(n_cells) for label in ("base", "tooth"))
            for cell in range(n_cells):
                base, tooth = 2 * cell, 2 * cell + 1
                _accumulate_hopping(
                    hoppings, onsite, base, tooth, -complex(parameters["tooth_hopping"])
                )
                if cell < n_cells - 1 or periodic:
                    next_base = 2 * ((cell + 1) % n_cells)
                    _accumulate_hopping(
                        hoppings, onsite, base, next_base, -complex(parameters["base_hopping"])
                    )
                    _accumulate_hopping(
                        hoppings, onsite, tooth, next_base, -complex(parameters["tooth_hopping"])
                    )
        else:
            sublattices = tuple(label for _ in range(n_cells) for label in ("upper", "lower"))
            phase = float(parameters["flux"]) / 2.0
            for cell in range(n_cells):
                upper, lower = 2 * cell, 2 * cell + 1
                onsite[upper] = complex(parameters.get("mass", 0.0))
                onsite[lower] = -complex(parameters.get("mass", 0.0))
                if cell < n_cells - 1 or periodic:
                    next_upper = 2 * ((cell + 1) % n_cells)
                    next_lower = next_upper + 1
                    _accumulate_hopping(
                        hoppings,
                        onsite,
                        upper,
                        next_upper,
                        -complex(parameters["hopping"]) * np.exp(1j * phase),
                    )
                    _accumulate_hopping(
                        hoppings,
                        onsite,
                        lower,
                        next_lower,
                        -complex(parameters["hopping"]) * np.exp(-1j * phase),
                    )
                    _accumulate_hopping(
                        hoppings,
                        onsite,
                        upper,
                        next_lower,
                        -complex(parameters["diagonal_hopping"]),
                    )
                    _accumulate_hopping(
                        hoppings,
                        onsite,
                        lower,
                        next_upper,
                        -complex(parameters["diagonal_hopping"]),
                    )
        labels = tuple(f"{sublattices[index]} {index // 2}" for index in range(n_sites))
        orbital_labels = sublattices
    else:
        n_sites = int(parameters["n_sites"])
        positions = tuple((float(site), 0.0) for site in range(n_sites))
        unit_cells = tuple(range(n_sites))
        sublattices = ()
        labels = tuple(f"site {site}" for site in range(n_sites))
        orbital_labels = tuple("orbital" for _ in range(n_sites))
        if family == "long_range_tight_binding_chain":
            for left in range(n_sites):
                for right in range(left + 1, n_sites):
                    distance = right - left
                    if periodic:
                        distance = min(distance, n_sites - distance)
                    _accumulate_hopping(
                        hoppings,
                        onsite,
                        left,
                        right,
                        -complex(parameters["hopping"]) / distance ** float(parameters["power"]),
                    )
        else:
            for source, target in _chain_bonds(n_sites, periodic):
                _accumulate_hopping(
                    hoppings, onsite, source, target, -complex(parameters["hopping"])
                )
        if family == "anderson_chain":
            values = np.random.default_rng(int(parameters["seed"])).uniform(
                -float(parameters["disorder"]) / 2.0,
                float(parameters["disorder"]) / 2.0,
                n_sites,
            )
        elif family == "aubry_andre_harper_chain":
            sites = np.arange(n_sites, dtype=float)
            values = float(parameters["potential"]) * np.cos(
                2.0 * np.pi * float(parameters["beta"]) * sites + float(parameters["phase"])
            )
        else:
            values = _site_parameter_values(parameters.get("onsite", 0.0), n_sites)
        onsite.update({site: complex(value) for site, value in enumerate(values)})

    if lattice is None:
        lattice = LatticeSpec(
            n_sites=n_sites,
            positions=positions,
            bonds=tuple(
                Bond(source, target, value) for (source, target), value in hoppings.items()
            ),
            site_labels=labels,
            orbital_labels=orbital_labels,
            sublattice_labels=sublattices,
            unit_cells=unit_cells,
            boundary_conditions={"x": "periodic" if periodic else "open"},
            conventions={"basis_ordering": "site/orbital index order"},
        )
    return _particle_records(lattice, hoppings, onsite)


def _planar_particle_physical_system(
    family: str,
    parameters: dict[str, object],
    lattice: LatticeSpec | None,
) -> tuple[
    LatticeSpec,
    tuple[LocalDegreeOfFreedom, ...],
    tuple[BasisIndexMapping, ...],
    tuple[InteractionTerm, ...],
]:
    n_rows = int(parameters["n_rows"])
    n_cols = int(parameters["n_cols"])
    periodic_x = bool(parameters.get("periodic_x", False))
    periodic_y = bool(parameters.get("periodic_y", False))
    orbitals = 1
    if family in {"graphene_lattice", "haldane_honeycomb_lattice", "checkerboard_chern_insulator"}:
        orbitals = 2
    elif family in {"kagome_lattice_tight_binding", "lieb_lattice", "dice_lattice"}:
        orbitals = 3
    n_sites = orbitals * n_rows * n_cols
    hoppings: dict[tuple[int, int], complex] = {}
    onsite: dict[int, complex] = {}

    def index(row: int, col: int, orbital: int = 0) -> int:
        return orbitals * (row * n_cols + col) + orbital

    if family in {
        "square_lattice_tight_binding",
        "anderson_square_lattice",
        "harper_hofstadter_square_lattice",
    }:
        for row in range(n_rows):
            for col in range(n_cols):
                source = index(row, col)
                for d_row, d_col in ((0, 1), (1, 0)):
                    target_cell = _wrapped_cell(
                        row + d_row,
                        col + d_col,
                        n_rows,
                        n_cols,
                        periodic_y,
                        periodic_x,
                    )
                    if target_cell is None:
                        continue
                    coefficient = -complex(parameters["hopping"])
                    if family == "harper_hofstadter_square_lattice" and d_row == 1:
                        coefficient *= np.exp(2j * np.pi * float(parameters["flux"]) * col)
                    _accumulate_hopping(hoppings, onsite, source, index(*target_cell), coefficient)
        if family == "anderson_square_lattice":
            values = np.random.default_rng(int(parameters["seed"])).uniform(
                -float(parameters["disorder"]) / 2.0,
                float(parameters["disorder"]) / 2.0,
                n_sites,
            )
        else:
            values = _site_parameter_values(parameters.get("onsite", 0.0), n_sites)
        onsite.update({site: complex(value) for site, value in enumerate(values)})
    elif family == "triangular_lattice_tight_binding":
        for row in range(n_rows):
            for col in range(n_cols):
                for d_row, d_col in ((0, 1), (1, 0), (1, -1)):
                    target = _wrapped_cell(
                        row + d_row,
                        col + d_col,
                        n_rows,
                        n_cols,
                        periodic_y,
                        periodic_x,
                    )
                    if target is not None:
                        _accumulate_hopping(
                            hoppings,
                            onsite,
                            index(row, col),
                            index(*target),
                            -complex(parameters["hopping"]),
                        )
        values = _site_parameter_values(parameters.get("onsite", 0.0), n_sites)
        onsite.update({site: complex(value) for site, value in enumerate(values)})
    elif family == "kagome_lattice_tight_binding":
        for row in range(n_rows):
            for col in range(n_cols):
                for source_orbital, target_orbital in ((0, 1), (1, 2), (2, 0)):
                    _accumulate_hopping(
                        hoppings,
                        onsite,
                        index(row, col, source_orbital),
                        index(row, col, target_orbital),
                        -complex(parameters["hopping"]),
                    )
                for source_orbital, d_row, d_col, target_orbital in (
                    (1, 0, 1, 0),
                    (2, 1, 0, 0),
                    (2, 1, -1, 1),
                ):
                    target = _wrapped_cell(
                        row + d_row,
                        col + d_col,
                        n_rows,
                        n_cols,
                        periodic_y,
                        periodic_x,
                    )
                    if target is not None:
                        _accumulate_hopping(
                            hoppings,
                            onsite,
                            index(row, col, source_orbital),
                            index(*target, target_orbital),
                            -complex(parameters["hopping"]),
                        )
        values = _site_parameter_values(parameters.get("onsite", 0.0), n_sites)
        onsite.update({site: complex(value) for site, value in enumerate(values)})
    elif family in {"graphene_lattice", "haldane_honeycomb_lattice"}:
        t1 = complex(parameters.get("t1", parameters.get("hopping", 1.0)))
        t2 = complex(parameters.get("t2", 0.0))
        phi = float(parameters.get("phi", 0.0))
        mass = complex(parameters.get("sublattice_potential", 0.0))
        for row in range(n_rows):
            for col in range(n_cols):
                a, b = index(row, col, 0), index(row, col, 1)
                onsite[a], onsite[b] = mass, -mass
                _accumulate_hopping(hoppings, onsite, a, b, -t1)
                for target in (
                    _wrapped_cell(row, col - 1, n_rows, n_cols, periodic_y, periodic_x),
                    _wrapped_cell(row - 1, col, n_rows, n_cols, periodic_y, periodic_x),
                ):
                    if target is not None:
                        _accumulate_hopping(hoppings, onsite, a, index(*target, 1), -t1)
                for d_row, d_col in ((0, 1), (1, 0), (1, -1)):
                    target = _wrapped_cell(
                        row + d_row,
                        col + d_col,
                        n_rows,
                        n_cols,
                        periodic_y,
                        periodic_x,
                    )
                    if target is not None and t2 != 0:
                        _accumulate_hopping(
                            hoppings, onsite, a, index(*target, 0), -t2 * np.exp(1j * phi)
                        )
                        _accumulate_hopping(
                            hoppings, onsite, b, index(*target, 1), -t2 * np.exp(-1j * phi)
                        )
    elif family == "checkerboard_chern_insulator":
        sigma_x = np.array([[0, 1], [1, 0]], dtype=complex)
        sigma_y = np.array([[0, -1j], [1j, 0]], dtype=complex)
        sigma_z = np.diag([1.0, -1.0]).astype(complex)
        for row in range(n_rows):
            for col in range(n_cols):
                onsite[index(row, col, 0)] = complex(parameters["mass"])
                onsite[index(row, col, 1)] = -complex(parameters["mass"])
                for d_row, d_col, pauli in ((0, 1, sigma_x), (1, 0, sigma_y)):
                    target = _wrapped_cell(
                        row + d_row,
                        col + d_col,
                        n_rows,
                        n_cols,
                        periodic_y,
                        periodic_x,
                    )
                    if target is None:
                        continue
                    block = 0.5 * complex(parameters["hopping"]) * (sigma_z - 1j * pauli)
                    for source_orbital in range(2):
                        for target_orbital in range(2):
                            if block[source_orbital, target_orbital] != 0:
                                _accumulate_hopping(
                                    hoppings,
                                    onsite,
                                    index(row, col, source_orbital),
                                    index(*target, target_orbital),
                                    complex(block[source_orbital, target_orbital]),
                                )
    elif family == "lieb_lattice":
        for row in range(n_rows):
            for col in range(n_cols):
                center = index(row, col, 0)
                for target in (index(row, col, 1), index(row, col, 2)):
                    _accumulate_hopping(
                        hoppings, onsite, center, target, -complex(parameters["hopping"])
                    )
                if col + 1 < n_cols:
                    _accumulate_hopping(
                        hoppings,
                        onsite,
                        index(row, col, 1),
                        index(row, col + 1, 0),
                        -complex(parameters["hopping"]),
                    )
                if row + 1 < n_rows:
                    _accumulate_hopping(
                        hoppings,
                        onsite,
                        index(row, col, 2),
                        index(row + 1, col, 0),
                        -complex(parameters["hopping"]),
                    )
    else:  # dice_lattice
        for row in range(n_rows):
            for col in range(n_cols):
                hub = index(row, col, 0)
                for orbital, offsets in (
                    (1, ((0, 0), (0, -1), (-1, 0))),
                    (2, ((0, 0), (0, 1), (1, 0))),
                ):
                    for d_row, d_col in offsets:
                        target = _wrapped_cell(
                            row + d_row,
                            col + d_col,
                            n_rows,
                            n_cols,
                            periodic_y,
                            periodic_x,
                        )
                        if target is not None:
                            _accumulate_hopping(
                                hoppings,
                                onsite,
                                hub,
                                index(*target, orbital),
                                -complex(parameters["hopping"]),
                            )

    if orbitals == 1:
        positions = tuple(
            (float(col), float(-row)) for row in range(n_rows) for col in range(n_cols)
        )
        sublattices: tuple[str, ...] = ()
        orbital_labels = tuple("orbital" for _ in range(n_sites))
    else:
        names = {
            "checkerboard_chern_insulator": ("A", "B"),
            "graphene_lattice": ("A", "B"),
            "haldane_honeycomb_lattice": ("A", "B"),
            "kagome_lattice_tight_binding": ("A", "B", "C"),
            "lieb_lattice": ("center", "horizontal", "vertical"),
            "dice_lattice": ("hub", "rim-A", "rim-B"),
        }[family]
        offsets = (
            ((0.0, 0.0), (0.0, 0.35)) if orbitals == 2 else ((0.0, 0.0), (-0.25, 0.2), (0.25, -0.2))
        )
        positions = tuple(
            (float(col) + offsets[orbital][0], float(-row) + offsets[orbital][1])
            for row in range(n_rows)
            for col in range(n_cols)
            for orbital in range(orbitals)
        )
        sublattices = tuple(name for _ in range(n_rows * n_cols) for name in names)
        orbital_labels = sublattices
    unit_cells = tuple(cell for cell in range(n_rows * n_cols) for _ in range(orbitals))
    labels = tuple(
        f"{sublattices[site]} cell {unit_cells[site]}" if sublattices else f"site {site}"
        for site in range(n_sites)
    )
    if lattice is None:
        lattice = LatticeSpec(
            n_sites=n_sites,
            positions=positions,
            bonds=tuple(
                Bond(source, target, value) for (source, target), value in hoppings.items()
            ),
            site_labels=labels,
            orbital_labels=orbital_labels,
            sublattice_labels=sublattices,
            unit_cells=unit_cells,
            boundary_conditions={
                "x": "periodic" if periodic_x else "open",
                "y": "periodic" if periodic_y else "open",
            },
            conventions={"basis_ordering": "row-major cell order, then orbital"},
        )
    return _particle_records(lattice, hoppings, onsite)


def _particle_records(
    lattice: LatticeSpec,
    hoppings: dict[tuple[int, int], complex],
    onsite: dict[int, complex],
) -> tuple[
    LatticeSpec,
    tuple[LocalDegreeOfFreedom, ...],
    tuple[BasisIndexMapping, ...],
    tuple[InteractionTerm, ...],
]:
    degrees = _orbital_degrees(lattice)
    interactions = tuple(
        InteractionTerm(
            (source, target),
            ("create", "annihilate"),
            coefficient,
            "hopping",
            label=f"{source}-{target}",
            metadata={"hermitian_conjugate": True},
        )
        for (source, target), coefficient in sorted(hoppings.items())
        if coefficient != 0
    ) + tuple(
        InteractionTerm(
            (site,),
            ("number",),
            coefficient,
            "onsite",
            label=f"onsite {site}",
        )
        for site, coefficient in sorted(onsite.items())
        if coefficient != 0
    )
    return lattice, degrees, _single_particle_mappings(degrees), interactions


def _accumulate_hopping(
    hoppings: dict[tuple[int, int], complex],
    onsite: dict[int, complex],
    source: int,
    target: int,
    coefficient: complex,
) -> None:
    coefficient = complex(coefficient)
    if source == target:
        onsite[source] = onsite.get(source, 0.0j) + coefficient + coefficient.conjugate()
        return
    if source > target:
        source, target = target, source
        coefficient = coefficient.conjugate()
    hoppings[(source, target)] = hoppings.get((source, target), 0.0j) + coefficient


def _wrapped_cell(
    row: int,
    col: int,
    n_rows: int,
    n_cols: int,
    periodic_y: bool,
    periodic_x: bool,
) -> tuple[int, int] | None:
    if col < 0 or col >= n_cols:
        if not periodic_x:
            return None
        col %= n_cols
    if row < 0 or row >= n_rows:
        if not periodic_y:
            return None
        row %= n_rows
    return row, col


def _site_parameter_values(value: object, n_sites: int) -> tuple[complex, ...]:
    if isinstance(value, (tuple, list, np.ndarray)):
        values = tuple(complex(item) for item in value)
        if len(values) != n_sites:
            raise ValueError("onsite must contain one value per site.")
        return values
    return tuple(complex(value) for _ in range(n_sites))


def _custom_lattice_from_parameters(parameters: dict[str, object]) -> LatticeSpec:
    n_sites = int(parameters["n_sites"])
    bonds = tuple(
        Bond(
            int(record[0]),
            int(record[1]),
            None if len(record) == 2 else complex(record[2]),
        )
        for record in parameters.get("bonds", ())
    )
    return LatticeSpec(
        n_sites=n_sites,
        positions=tuple((float(site), 0.0) for site in range(n_sites)),
        bonds=bonds,
        site_labels=tuple(f"site {site}" for site in range(n_sites)),
        orbital_labels=tuple("orbital" for _ in range(n_sites)),
        conventions={"basis_ordering": "site index order"},
    )


def _custom_tight_binding_physical_system(
    parameters: dict[str, object],
    lattice: LatticeSpec,
) -> tuple[
    LatticeSpec,
    tuple[LocalDegreeOfFreedom, ...],
    tuple[BasisIndexMapping, ...],
    tuple[InteractionTerm, ...],
]:
    degrees = _orbital_degrees(lattice)
    interactions: list[InteractionTerm] = []
    hopping = complex(parameters.get("hopping", 1.0))
    hermitian = bool(parameters.get("hermitian", True))
    for bond in lattice.bonds:
        coefficient = -hopping if bond.value is None else complex(bond.value)
        interactions.append(
            InteractionTerm(
                (bond.source, bond.target),
                ("create", "annihilate"),
                coefficient,
                "hopping",
                label=f"{bond.source}-{bond.target}",
                metadata={"hermitian_conjugate": hermitian},
            )
        )
    onsite = parameters.get("onsite", 0.0)
    onsite_values = (
        tuple(complex(value) for value in onsite)
        if isinstance(onsite, (tuple, list))
        else tuple(complex(onsite) for _ in range(lattice.n_sites))
    )
    for site, coefficient in enumerate(onsite_values):
        if coefficient != 0:
            interactions.append(
                InteractionTerm(
                    (site,),
                    ("number",),
                    coefficient,
                    "onsite",
                    label=f"onsite {site}",
                )
            )
    return lattice, degrees, _single_particle_mappings(degrees), tuple(interactions)


def _bose_hubbard_physical_system(
    parameters: dict[str, object],
    lattice: LatticeSpec | None,
) -> tuple[
    LatticeSpec,
    tuple[LocalDegreeOfFreedom, ...],
    tuple[BasisIndexMapping, ...],
    tuple[InteractionTerm, ...],
]:
    n_sites = int(parameters["n_sites"])
    periodic = bool(parameters.get("periodic", False))
    lattice = lattice or _chain_lattice(n_sites, periodic, "boson")
    max_occupancy = int(parameters["max_occupancy"])
    degrees = tuple(
        LocalDegreeOfFreedom(
            site,
            site,
            "boson",
            max_occupancy + 1,
            f"boson {site}",
            component="truncated occupation",
            metadata={"max_occupancy": max_occupancy},
        )
        for site in range(n_sites)
    )
    mappings = tuple(
        BasisIndexMapping(site, site, "tensor_factor", f"occupation digit {site}")
        for site in range(n_sites)
    )
    interactions: list[InteractionTerm] = []
    for source, target in _chain_bonds(n_sites, periodic):
        interactions.append(
            InteractionTerm(
                (source, target),
                ("create", "annihilate"),
                -complex(parameters["hopping"]),
                "boson_hopping",
                label=f"{source}-{target}",
                metadata={"hermitian_conjugate": True},
            )
        )
    for site in range(n_sites):
        interactions.extend(
            (
                InteractionTerm(
                    (site,),
                    ("number_pair",),
                    0.5 * complex(parameters["interaction"]),
                    "onsite_interaction",
                    label=f"U {site}",
                ),
                InteractionTerm(
                    (site,),
                    ("number",),
                    -complex(parameters.get("chemical_potential", 0.0)),
                    "chemical_potential",
                    label=f"mu {site}",
                ),
            )
        )
    return lattice, degrees, mappings, tuple(interactions)


def _fermi_hubbard_physical_system(
    parameters: dict[str, object],
    lattice: LatticeSpec | None,
) -> tuple[
    LatticeSpec,
    tuple[LocalDegreeOfFreedom, ...],
    tuple[BasisIndexMapping, ...],
    tuple[InteractionTerm, ...],
]:
    n_sites = int(parameters["n_sites"])
    periodic = bool(parameters.get("periodic", False))
    lattice = lattice or _chain_lattice(n_sites, periodic, "fermion")
    degrees = tuple(
        LocalDegreeOfFreedom(
            2 * site + spin_index,
            site,
            "fermion",
            2,
            f"{site} {spin}",
            component=spin,
            orbital=spin,
        )
        for site in range(n_sites)
        for spin_index, spin in enumerate(("up", "down"))
    )
    mappings = tuple(
        BasisIndexMapping(degree.index, degree.index, "mode", degree.label) for degree in degrees
    )
    interactions: list[InteractionTerm] = []
    for source, target in _chain_bonds(n_sites, periodic):
        for spin_index, spin in enumerate(("up", "down")):
            source_mode = 2 * source + spin_index
            target_mode = 2 * target + spin_index
            interactions.append(
                InteractionTerm(
                    (source_mode, target_mode),
                    ("create", "annihilate"),
                    -complex(parameters["hopping"]),
                    "fermion_hopping",
                    label=f"{spin} {source}-{target}",
                    metadata={"hermitian_conjugate": True},
                )
            )
    for site in range(n_sites):
        up = 2 * site
        down = up + 1
        interactions.append(
            InteractionTerm(
                (up, down),
                ("number", "number"),
                complex(parameters["interaction"]),
                "onsite_interaction",
                label=f"U {site}",
            )
        )
        for mode in (up, down):
            interactions.append(
                InteractionTerm(
                    (mode,),
                    ("number",),
                    -complex(parameters.get("chemical_potential", 0.0)),
                    "chemical_potential",
                    label=f"mu {degrees[mode].label}",
                )
            )
    return lattice, degrees, mappings, tuple(interactions)


def _kitaev_bdg_physical_system(
    parameters: dict[str, object],
    lattice: LatticeSpec | None,
) -> tuple[
    LatticeSpec,
    tuple[LocalDegreeOfFreedom, ...],
    tuple[BasisIndexMapping, ...],
    tuple[InteractionTerm, ...],
]:
    n_sites = int(parameters["n_sites"])
    periodic = bool(parameters.get("periodic", False))
    lattice = lattice or _chain_lattice(n_sites, periodic, "nambu")
    degrees = tuple(
        LocalDegreeOfFreedom(
            index=(0 if component == "particle" else n_sites) + site,
            site=site,
            kind="nambu",
            local_dimension=2,
            label=f"{component} {site}",
            component=component,
        )
        for component in ("particle", "hole")
        for site in range(n_sites)
    )
    mappings = tuple(
        BasisIndexMapping(degree.index, degree.index, "single_particle_state", degree.label)
        for degree in degrees
    )
    interactions: list[InteractionTerm] = []
    chemical_potential = complex(parameters.get("chemical_potential", 0.0))
    for site in range(n_sites):
        interactions.extend(
            (
                InteractionTerm(
                    (site,),
                    ("particle",),
                    -chemical_potential,
                    "bdg_onsite",
                    label=f"particle {site}",
                ),
                InteractionTerm(
                    (n_sites + site,),
                    ("hole",),
                    chemical_potential.conjugate(),
                    "bdg_onsite",
                    label=f"hole {site}",
                ),
            )
        )
    hopping = complex(parameters["hopping"])
    pairing = complex(parameters["pairing"])
    for source, target in _chain_bonds(n_sites, periodic):
        interactions.extend(
            (
                InteractionTerm(
                    (source, target),
                    ("particle", "particle"),
                    -hopping,
                    "bdg_hopping",
                    label=f"particle {source}-{target}",
                    metadata={"hermitian_conjugate": True},
                ),
                InteractionTerm(
                    (n_sites + source, n_sites + target),
                    ("hole", "hole"),
                    hopping.conjugate(),
                    "bdg_hopping",
                    label=f"hole {source}-{target}",
                    metadata={"hermitian_conjugate": True},
                ),
                InteractionTerm(
                    (source, n_sites + target),
                    ("particle", "hole"),
                    pairing,
                    "bdg_pairing",
                    label=f"pair {source}-{target}",
                    metadata={"antisymmetric_partner": [target, n_sites + source]},
                ),
                InteractionTerm(
                    (target, n_sites + source),
                    ("particle", "hole"),
                    -pairing,
                    "bdg_pairing",
                    label=f"pair {target}-{source}",
                    metadata={"antisymmetric_partner": [source, n_sites + target]},
                ),
            )
        )
    return lattice, degrees, mappings, tuple(interactions)


def _chain_lattice(n_sites: int, periodic: bool, label: str) -> LatticeSpec:
    return LatticeSpec(
        n_sites=n_sites,
        positions=tuple((float(site), 0.0) for site in range(n_sites)),
        bonds=tuple(Bond(source, target) for source, target in _chain_bonds(n_sites, periodic)),
        site_labels=tuple(f"{label} site {site}" for site in range(n_sites)),
        boundary_conditions={"x": "periodic" if periodic else "open"},
        conventions={"site_ordering": "left-to-right chain order"},
    )


def _orbital_degrees(lattice: LatticeSpec) -> tuple[LocalDegreeOfFreedom, ...]:
    return tuple(
        LocalDegreeOfFreedom(
            index=site,
            site=site,
            kind="orbital",
            local_dimension=2,
            label=lattice.site_labels[site] if lattice.site_labels else f"site {site}",
            orbital=lattice.orbital_labels[site] if lattice.orbital_labels else None,
            component=(lattice.sublattice_labels[site] if lattice.sublattice_labels else None),
        )
        for site in range(lattice.n_sites)
    )


def _single_particle_mappings(
    degrees: tuple[LocalDegreeOfFreedom, ...],
) -> tuple[BasisIndexMapping, ...]:
    return tuple(
        BasisIndexMapping(degree.index, degree.index, "single_particle_state", degree.label)
        for degree in degrees
    )


def _hopping_interaction(
    source: int,
    target: int,
    coefficient: complex,
    kind: str,
    label: str,
) -> InteractionTerm:
    return InteractionTerm(
        (source, target),
        ("create", "annihilate"),
        coefficient,
        kind,
        label=label,
        metadata={"hermitian_conjugate": True},
    )


def _chain_bonds(n_sites: int, periodic: bool) -> tuple[tuple[int, int], ...]:
    bonds = [(site, site + 1) for site in range(n_sites - 1)]
    if periodic and n_sites > 2:
        bonds.append((n_sites - 1, 0))
    return tuple(bonds)


def _next_nearest_chain_bonds(
    n_sites: int,
    periodic: bool,
) -> tuple[tuple[int, int], ...]:
    bonds = {(site, site + 2) for site in range(n_sites - 2)}
    if periodic and n_sites > 3:
        bonds.update(tuple(sorted((site, (site + 2) % n_sites))) for site in range(n_sites))
    return tuple(sorted(bonds))


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
