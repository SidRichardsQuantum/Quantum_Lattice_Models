"""Computational estimates and numerical symmetry diagnostics."""

from __future__ import annotations

import ast
from dataclasses import dataclass
from math import comb

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models.registry import get_model_info, supports_sparse


@dataclass(frozen=True)
class MatrixDiagnostics:
    """Summary of matrix structure, storage, and Hermiticity."""

    shape: tuple[int, int]
    sparse: bool
    nonzero_entries: int
    density: float
    storage_bytes: int
    hermitian: bool


@dataclass(frozen=True)
class ModelInspection:
    """Pre-build model dimension, storage, basis, and warning report."""

    model: str
    category: str
    basis: str
    representation: str
    dimension: int
    dense_memory_bytes: int
    supports_sparse: bool
    validation_status: str
    parameters: dict[str, object]
    warnings: tuple[str, ...]

    def to_dict(self) -> dict[str, object]:
        """Return a deterministic machine-readable report."""

        return {
            "model": self.model,
            "category": self.category,
            "basis": self.basis,
            "representation": self.representation,
            "dimension": self.dimension,
            "dense_memory_bytes": self.dense_memory_bytes,
            "supports_sparse": self.supports_sparse,
            "validation_status": self.validation_status,
            "parameters": dict(sorted(self.parameters.items())),
            "warnings": list(self.warnings),
        }


def inspect_model(
    model_name: str,
    *,
    representation: str | None = None,
    **parameters: object,
) -> ModelInspection:
    """Return a construction-free resource and risk report."""

    info = get_model_info(model_name)
    values = dict(info.defaults)
    values.update(parameters)
    resolved_representation = representation or (
        "sparse" if model_name.endswith("_sparse") else "dense"
    )
    if resolved_representation not in {"dense", "sparse"}:
        raise ValueError("representation must be 'dense' or 'sparse'.")
    sparse_available = supports_sparse(model_name)
    if resolved_representation == "sparse" and not sparse_available:
        raise ValueError(f"Model {model_name!r} does not support sparse construction.")
    dimension = estimate_model_dimension(model_name, **values)
    memory = estimate_dense_memory(dimension)
    warnings: list[str] = []
    if memory >= 512 * 1024**2:
        warnings.append("dense matrix estimate is at least 512 MiB")
    elif memory >= 64 * 1024**2:
        warnings.append("dense matrix estimate is at least 64 MiB")
    if "**n_sites" in info.dimension or "**(2*n_sites)" in info.dimension:
        warnings.append("Hilbert-space dimension grows exponentially with site count")
    if resolved_representation == "dense" and sparse_available and dimension >= 1024:
        warnings.append("a sparse builder is available and may use substantially less memory")
    if info.validation_status == "unvalidated":
        warnings.append("model has no registered package validation status")
    return ModelInspection(
        model=model_name,
        category=info.category,
        basis=info.basis,
        representation=resolved_representation,
        dimension=dimension,
        dense_memory_bytes=memory,
        supports_sparse=sparse_available,
        validation_status=info.validation_status,
        parameters=values,
        warnings=tuple(warnings),
    )


def estimate_model_dimension(model_name: str, **parameters: object) -> int:
    """Estimate a registered model's matrix dimension from its dimension expression."""

    info = get_model_info(model_name)
    values = dict(info.defaults)
    values.update(parameters)
    if "magnetization" in values:
        n_sites = values.get("n_sites")
        magnetization = values["magnetization"]
        if not isinstance(n_sites, int) or not isinstance(magnetization, int):
            raise ValueError("n_sites and magnetization must be integers.")
        if abs(magnetization) > n_sites or (n_sites - magnetization) % 2:
            raise ValueError(
                "magnetization must satisfy abs(magnetization) <= n_sites "
                "and have the same parity as n_sites."
            )
    if "parity" in values and values["parity"] not in {-1, 1}:
        raise ValueError("parity must be +1 or -1.")
    dimension = _evaluate_dimension(info.dimension, values)
    if dimension < 1:
        raise ValueError("Estimated model dimension must be positive.")
    return dimension


def estimate_dense_memory(
    dimension: int,
    *,
    dtype: np.dtype | type = np.complex128,
) -> int:
    """Return bytes required for a square dense matrix."""

    if not isinstance(dimension, int) or dimension < 1:
        raise ValueError("dimension must be a positive integer.")
    return dimension * dimension * np.dtype(dtype).itemsize


def matrix_density(matrix: np.ndarray | sp.spmatrix) -> float:
    """Return the fraction of nonzero matrix entries."""

    rows, cols = _matrix_shape(matrix)
    nonzero = matrix.nnz if sp.issparse(matrix) else np.count_nonzero(matrix)
    return float(nonzero / (rows * cols))


def matrix_storage_bytes(matrix: np.ndarray | sp.spmatrix) -> int:
    """Return bytes used by dense data or a SciPy sparse matrix's core arrays."""

    if sp.issparse(matrix):
        compressed = matrix.tocsr()
        return int(compressed.data.nbytes + compressed.indices.nbytes + compressed.indptr.nbytes)
    return int(np.asarray(matrix).nbytes)


def is_hermitian(
    matrix: np.ndarray | sp.spmatrix,
    *,
    atol: float = 1e-10,
    rtol: float = 1e-10,
) -> bool:
    """Return whether a square dense or sparse matrix is Hermitian."""

    rows, cols = _matrix_shape(matrix)
    if rows != cols:
        return False
    if sp.issparse(matrix):
        difference = (matrix - matrix.getH()).tocsr()
        if difference.nnz == 0:
            return True
        scale = max(float(np.max(np.abs(matrix.data), initial=0.0)), 1.0)
        return bool(np.max(np.abs(difference.data), initial=0.0) <= atol + rtol * scale)
    dense = np.asarray(matrix)
    return bool(np.allclose(dense, dense.conj().T, atol=atol, rtol=rtol))


def has_particle_hole_symmetric_spectrum(
    matrix: np.ndarray | sp.spmatrix,
    *,
    atol: float = 1e-9,
    rtol: float = 1e-9,
) -> bool:
    """Return whether a Hermitian spectrum is paired as ``E`` and ``-E``."""

    if not is_hermitian(matrix, atol=atol, rtol=rtol):
        return False
    dense = matrix.toarray() if sp.issparse(matrix) else np.asarray(matrix)
    values = np.linalg.eigvalsh(dense)
    return bool(np.allclose(values, -values[::-1], atol=atol, rtol=rtol))


def diagnose_matrix(matrix: np.ndarray | sp.spmatrix) -> MatrixDiagnostics:
    """Return a compact structural diagnostic summary."""

    rows, cols = _matrix_shape(matrix)
    nonzero = int(matrix.nnz if sp.issparse(matrix) else np.count_nonzero(matrix))
    return MatrixDiagnostics(
        shape=(rows, cols),
        sparse=bool(sp.issparse(matrix)),
        nonzero_entries=nonzero,
        density=float(nonzero / (rows * cols)),
        storage_bytes=matrix_storage_bytes(matrix),
        hermitian=is_hermitian(matrix),
    )


def _matrix_shape(matrix: np.ndarray | sp.spmatrix) -> tuple[int, int]:
    shape = matrix.shape
    if len(shape) != 2 or shape[0] < 1 or shape[1] < 1:
        raise ValueError("matrix must be a nonempty two-dimensional array.")
    return int(shape[0]), int(shape[1])


def _evaluate_dimension(expression: str, values: dict[str, object]) -> int:
    tree = ast.parse(expression, mode="eval")

    def evaluate(node: ast.AST) -> int:
        if isinstance(node, ast.Expression):
            return evaluate(node.body)
        if isinstance(node, ast.Constant) and isinstance(node.value, int):
            return node.value
        if isinstance(node, ast.Name):
            try:
                value = values[node.id]
            except KeyError as exc:
                raise ValueError(f"Missing dimension parameter {node.id!r}.") from exc
            if not isinstance(value, int):
                raise ValueError(f"Dimension parameter {node.id!r} must be an integer.")
            return value
        if isinstance(node, ast.BinOp) and isinstance(
            node.op, (ast.Add, ast.Sub, ast.Mult, ast.FloorDiv, ast.Pow)
        ):
            left = evaluate(node.left)
            right = evaluate(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Sub):
                return left - right
            if isinstance(node.op, ast.Mult):
                return left * right
            if isinstance(node.op, ast.FloorDiv):
                return left // right
            return left**right
        if (
            isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "comb"
            and len(node.args) == 2
            and not node.keywords
        ):
            n = evaluate(node.args[0])
            k = evaluate(node.args[1])
            try:
                return comb(n, k)
            except ValueError as exc:
                raise ValueError("Invalid fixed-sector dimension parameters.") from exc
        raise ValueError(f"Unsupported dimension expression {expression!r}.")

    return evaluate(tree)
