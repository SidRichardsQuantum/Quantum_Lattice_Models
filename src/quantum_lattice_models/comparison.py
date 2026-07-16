"""Deterministic summaries comparing portable model specifications."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models.diagnostics import estimate_model_dimension
from quantum_lattice_models.specs import ModelSpec
from quantum_lattice_models.spectra import eigenvalues, spectral_gap


@dataclass(frozen=True)
class ModelComparison:
    """Parameters, matrix, spectrum, and gap comparison for two models."""

    left_family: str
    right_family: str
    same_basis: bool
    parameter_differences: dict[str, dict[str, object]]
    left_dimension: int | None
    right_dimension: int | None
    same_shape: bool | None
    matrix_frobenius_difference: float | None
    matrix_maximum_difference: float | None
    spectrum_maximum_difference: float | None
    left_gap: float | None
    right_gap: float | None
    warnings: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        """Return a deterministic machine-readable summary."""

        return {
            "left_family": self.left_family,
            "right_family": self.right_family,
            "same_basis": self.same_basis,
            "parameter_differences": {
                key: value for key, value in sorted(self.parameter_differences.items())
            },
            "left_dimension": self.left_dimension,
            "right_dimension": self.right_dimension,
            "same_shape": self.same_shape,
            "matrix_frobenius_difference": self.matrix_frobenius_difference,
            "matrix_maximum_difference": self.matrix_maximum_difference,
            "spectrum_maximum_difference": self.spectrum_maximum_difference,
            "left_gap": self.left_gap,
            "right_gap": self.right_gap,
            "warnings": list(self.warnings),
        }


def compare_models(
    left: ModelSpec,
    right: ModelSpec,
    *,
    max_dimension: int = 2048,
) -> ModelComparison:
    """Compare model metadata and, when safe, matrices and spectra."""

    if max_dimension < 1:
        raise ValueError("max_dimension must be positive.")
    parameter_differences = _parameter_differences(left.parameters, right.parameters)
    left_dimension = _dimension(left)
    right_dimension = _dimension(right)
    warnings: list[str] = []
    if left.basis != right.basis:
        warnings.append("models use different basis conventions")
    if left_dimension is None or right_dimension is None:
        warnings.append("one or both model dimensions could not be estimated")
    elif max(left_dimension, right_dimension) > max_dimension:
        warnings.append(f"matrix and spectrum comparison skipped above dimension {max_dimension}")
        return ModelComparison(
            left.family,
            right.family,
            left.basis == right.basis,
            parameter_differences,
            left_dimension,
            right_dimension,
            left_dimension == right_dimension,
            None,
            None,
            None,
            None,
            None,
            tuple(warnings),
        )

    left_matrix = left.hamiltonian()
    right_matrix = right.hamiltonian()
    same_shape = left_matrix.shape == right_matrix.shape
    frobenius = maximum = spectrum_difference = None
    if same_shape:
        difference = _dense(left_matrix) - _dense(right_matrix)
        frobenius = float(np.linalg.norm(difference))
        maximum = float(np.max(np.abs(difference), initial=0.0))
        left_values = np.sort(np.real_if_close(eigenvalues(left_matrix)).real)
        right_values = np.sort(np.real_if_close(eigenvalues(right_matrix)).real)
        spectrum_difference = float(np.max(np.abs(left_values - right_values), initial=0.0))
    else:
        warnings.append("matrix and spectrum differences require matching shapes")
    return ModelComparison(
        left.family,
        right.family,
        left.basis == right.basis,
        parameter_differences,
        left_dimension,
        right_dimension,
        same_shape,
        frobenius,
        maximum,
        spectrum_difference,
        float(spectral_gap(left_matrix)),
        float(spectral_gap(right_matrix)),
        tuple(warnings),
    )


def _parameter_differences(
    left: dict[str, object], right: dict[str, object]
) -> dict[str, dict[str, object]]:
    differences = {}
    for name in sorted(set(left) | set(right)):
        left_value = left.get(name)
        right_value = right.get(name)
        if left_value != right_value:
            differences[name] = {"left": left_value, "right": right_value}
    return differences


def _dimension(spec: ModelSpec) -> int | None:
    parameters = dict(spec.parameters)
    if spec.lattice is not None:
        parameters.setdefault("n_sites", spec.lattice.n_sites)
    try:
        return estimate_model_dimension(spec.family, **parameters)
    except ValueError:
        return None


def _dense(matrix: np.ndarray | sp.spmatrix) -> np.ndarray:
    return sp.csr_matrix(matrix).toarray() if sp.issparse(matrix) else np.asarray(matrix)
