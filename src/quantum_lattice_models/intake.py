"""Model-intake summaries, linting, and adapter capability reports."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models.diagnostics import diagnose_matrix, estimate_dense_memory
from quantum_lattice_models.specs import EXTERNAL_MATRIX_FAMILY, ModelSpec
from quantum_lattice_models.types import HamiltonianResult


@dataclass(frozen=True)
class ModelSummary:
    """User-facing summary of a portable model or Hamiltonian result."""

    family: str
    basis: str | None
    representation: str
    dimension: int | None
    lattice_sites: int | None
    lattice_bonds: int | None
    local_degrees: int
    interactions: int
    symmetry_actions: int
    basis_mappings: int
    boundary_conditions: dict[str, str]
    units: dict[str, str]
    conventions: dict[str, str]
    reconstructable: bool
    matrix_shape: tuple[int, int] | None = None
    warnings: tuple[str, ...] = ()

    def to_dict(self) -> dict[str, object]:
        """Return a deterministic JSON-compatible summary."""

        return {
            "family": self.family,
            "basis": self.basis,
            "representation": self.representation,
            "dimension": self.dimension,
            "lattice_sites": self.lattice_sites,
            "lattice_bonds": self.lattice_bonds,
            "local_degrees": self.local_degrees,
            "interactions": self.interactions,
            "symmetry_actions": self.symmetry_actions,
            "basis_mappings": self.basis_mappings,
            "boundary_conditions": dict(sorted(self.boundary_conditions.items())),
            "units": dict(sorted(self.units.items())),
            "conventions": dict(sorted(self.conventions.items())),
            "reconstructable": self.reconstructable,
            "matrix_shape": None if self.matrix_shape is None else list(self.matrix_shape),
            "warnings": list(self.warnings),
        }


@dataclass(frozen=True)
class ModelLintReport:
    """Validation and consistency report for model intake workflows."""

    valid: bool
    errors: tuple[str, ...]
    warnings: tuple[str, ...]
    suggestions: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        """Return a deterministic JSON-compatible lint report."""

        return {
            "valid": self.valid,
            "errors": list(self.errors),
            "warnings": list(self.warnings),
            "suggestions": list(self.suggestions),
        }


@dataclass(frozen=True)
class AdapterCapabilityReport:
    """Report what a model translation target can preserve."""

    target: str
    supported: bool
    preserved: tuple[str, ...]
    lost: tuple[str, ...]
    warnings: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        """Return a deterministic JSON-compatible capability report."""

        return {
            "target": self.target,
            "supported": self.supported,
            "preserved": list(self.preserved),
            "lost": list(self.lost),
            "warnings": list(self.warnings),
        }


def describe_model(source: ModelSpec | HamiltonianResult) -> ModelSummary:
    """Return a compact user-facing description of a model or Hamiltonian result."""

    model, matrix = _model_and_matrix(source)
    raw_summary = model.summary()
    lattice = model.lattice
    warnings: list[str] = []
    if model.family == EXTERNAL_MATRIX_FAMILY:
        warnings.append("external matrices cannot be reconstructed from a registered builder")
    if not model.local_degrees:
        warnings.append("physical local-degree records are not available")
    if lattice is not None and not lattice.positions:
        warnings.append("lattice coordinates are not available")
    matrix_shape = None
    if matrix is not None:
        matrix_shape = (int(matrix.shape[0]), int(matrix.shape[1]))
    dimension = raw_summary.get("dimension")
    resolved_dimension = (
        dimension if isinstance(dimension, int) else matrix_shape[0] if matrix_shape else None
    )
    return ModelSummary(
        family=model.family,
        basis=model.basis,
        representation=model.representation,
        dimension=resolved_dimension,
        lattice_sites=None if lattice is None else lattice.n_sites,
        lattice_bonds=None if lattice is None else len(lattice.bonds),
        local_degrees=len(model.local_degrees),
        interactions=len(model.interactions),
        symmetry_actions=len(model.symmetry_actions),
        basis_mappings=len(model.basis_mappings),
        boundary_conditions={} if lattice is None else dict(lattice.boundary_conditions),
        units={**({} if lattice is None else lattice.units), **model.units},
        conventions={**({} if lattice is None else lattice.conventions), **model.conventions},
        reconstructable=model.family != EXTERNAL_MATRIX_FAMILY,
        matrix_shape=matrix_shape,
        warnings=tuple(warnings),
    )


def lint_model(source: ModelSpec | HamiltonianResult) -> ModelLintReport:
    """Return validation errors, warnings, and follow-up checks for a model."""

    model, matrix = _model_and_matrix(source)
    errors: list[str] = []
    warnings: list[str] = []
    suggestions: list[str] = []
    try:
        model.validate()
    except ValueError as exc:
        errors.append(str(exc))

    lattice = model.lattice
    if lattice is None:
        warnings.append("model has no portable lattice geometry")
    else:
        _extend_duplicate_warnings(warnings, "site labels", lattice.site_labels)
        _extend_duplicate_warnings(warnings, "orbital labels", lattice.orbital_labels)
        if not lattice.positions:
            warnings.append("lattice positions are missing")
        if not lattice.boundary_conditions:
            warnings.append("boundary conditions are unspecified")

    degree_labels = tuple(degree.label for degree in model.local_degrees if degree.label)
    _extend_duplicate_warnings(warnings, "local-degree labels", degree_labels)
    if model.local_degrees and not model.basis_mappings:
        warnings.append("local degrees are present without basis mappings")
    if model.interactions and not model.local_degrees:
        warnings.append("interactions are present without local-degree records")
    if model.symmetry_actions and matrix is None:
        suggestions.append("validate advertised symmetry actions after construction")
    if model.family == EXTERNAL_MATRIX_FAMILY:
        warnings.append("external matrix records cannot reconstruct a Hamiltonian builder")

    if matrix is not None:
        diagnostics = diagnose_matrix(matrix)
        if diagnostics.shape[0] != diagnostics.shape[1]:
            errors.append("matrix is not square")
        if not diagnostics.hermitian:
            warnings.append("matrix is not Hermitian")
        if diagnostics.storage_bytes >= 64 * 1024**2:
            warnings.append("matrix storage is at least 64 MiB")
    else:
        dimension = describe_model(model).dimension
        if dimension is not None and estimate_dense_memory(dimension) >= 64 * 1024**2:
            warnings.append("dense construction estimate is at least 64 MiB")

    suggestions.extend(_validation_suggestions(model, matrix is not None))
    return ModelLintReport(
        valid=not errors,
        errors=tuple(errors),
        warnings=tuple(dict.fromkeys(warnings)),
        suggestions=tuple(dict.fromkeys(suggestions)),
    )


def adapter_capability_report(
    source: ModelSpec | HamiltonianResult,
    target: str,
) -> AdapterCapabilityReport:
    """Report whether a target adapter can preserve key model semantics."""

    model, _matrix = _model_and_matrix(source)
    normalized = target.lower().replace("_", "-")
    features = _model_features(model)
    preserved: set[str] = set()
    lost: set[str] = set()
    warnings: list[str] = []
    supported = True

    if normalized in {"yaml", "model-yaml", "json", "model-json"}:
        preserved = set(features)
    elif normalized in {"csv", "graphml", "networkx"}:
        supported = model.lattice is not None
        preserved = features & {"geometry", "site_labels", "bonds", "units", "provenance"}
        lost = features - preserved
        if not supported:
            warnings.append("lattice export requires portable lattice data")
    elif normalized in {"xyz", "ase"}:
        supported = model.lattice is not None and bool(model.lattice.positions)
        preserved = features & {"geometry", "site_labels", "units", "provenance"}
        lost = features - preserved
        warnings.append("coordinate formats do not preserve bonds or Hamiltonian terms")
    elif normalized in {"dot", "svg", "plot-data-json"}:
        supported = model.lattice is not None
        preserved = features & {"geometry", "site_labels", "bonds"}
        lost = features - preserved
        warnings.append("visual exports are derived artifacts, not reconstructable models")
    elif normalized in {"qiskit", "quspin"}:
        supported = bool(model.local_degrees) and all(
            degree.kind == "spin" for degree in model.local_degrees
        )
        preserved = features & {"interactions", "basis", "local_degrees"}
        lost = features - preserved
        if not supported:
            warnings.append("target supports portable spin Pauli terms only")
    elif normalized == "openfermion":
        supported = bool(model.local_degrees) and all(
            degree.kind == "fermion" for degree in model.local_degrees
        )
        preserved = features & {"interactions", "basis", "local_degrees"}
        lost = features - preserved
        if not supported:
            warnings.append("target supports portable fermionic terms only")
    elif normalized == "pennylane":
        supported = model.family != EXTERNAL_MATRIX_FAMILY
        preserved = features & {"interactions", "basis"}
        lost = features - preserved
        warnings.append("PennyLane export uses operator terms and does not preserve geometry")
    elif normalized == "qutip":
        supported = model.family != EXTERNAL_MATRIX_FAMILY
        preserved = features & {"basis", "local_degrees"}
        lost = features - preserved
        warnings.append(
            "QuTiP export constructs a matrix object and does not preserve full metadata"
        )
    elif normalized == "netket":
        supported = model.lattice is not None and bool(model.lattice.bonds)
        preserved = features & {"geometry", "bonds"}
        lost = features - preserved
        if not supported:
            warnings.append("NetKet graph export requires lattice bonds")
    else:
        raise ValueError(f"Unknown adapter target {target!r}.")

    if not lost:
        lost = features - preserved
    full_spec_targets = {"yaml", "model-yaml", "json", "model-json"}
    if model.family == EXTERNAL_MATRIX_FAMILY and normalized not in full_spec_targets:
        warnings.append("external matrix imports may not be reconstructable after translation")
    return AdapterCapabilityReport(
        target=normalized,
        supported=supported,
        preserved=tuple(sorted(preserved)),
        lost=tuple(sorted(lost)),
        warnings=tuple(dict.fromkeys(warnings)),
    )


def _model_and_matrix(
    source: ModelSpec | HamiltonianResult,
) -> tuple[ModelSpec, np.ndarray | sp.spmatrix | None]:
    if isinstance(source, HamiltonianResult):
        return source.model, source.matrix
    if isinstance(source, ModelSpec):
        return source, None
    raise TypeError("source must be a ModelSpec or HamiltonianResult.")


def _extend_duplicate_warnings(warnings: list[str], name: str, values: tuple[Any, ...]) -> None:
    if values and len(set(values)) != len(values):
        warnings.append(f"duplicate {name} detected")


def _model_features(model: ModelSpec) -> set[str]:
    features = {"basis", "parameters", "schema_version"}
    if model.lattice is not None:
        features.add("geometry")
        if model.lattice.bonds:
            features.add("bonds")
        if model.lattice.site_labels:
            features.add("site_labels")
        if model.lattice.units or model.units:
            features.add("units")
    if model.local_degrees:
        features.add("local_degrees")
    if model.basis_mappings:
        features.add("basis_mappings")
    if model.interactions:
        features.add("interactions")
    if model.symmetry_actions:
        features.add("symmetry_actions")
    if model.provenance or (model.lattice is not None and model.lattice.provenance):
        features.add("provenance")
    if model.references or (model.lattice is not None and model.lattice.references):
        features.add("references")
    if model.metadata or (model.lattice is not None and model.lattice.metadata):
        features.add("metadata")
    return features


def _validation_suggestions(model: ModelSpec, has_matrix: bool) -> tuple[str, ...]:
    suggestions = ["run model.validate() before persistence or export"]
    if has_matrix:
        suggestions.append("check Hermiticity before spectral analysis")
    if model.interactions:
        suggestions.append("check conserved-quantity commutators for advertised sectors")
    if model.symmetry_actions:
        suggestions.append("check symmetry-action commutators and sector eigenvalues")
    if model.lattice is not None and model.lattice.bonds:
        suggestions.append("inspect connected components for imported lattice graphs")
    if model.family == EXTERNAL_MATRIX_FAMILY:
        suggestions.append("record basis conventions and provenance for external matrices")
    return tuple(suggestions)
