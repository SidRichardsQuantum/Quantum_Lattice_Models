"""Computational estimates and numerical symmetry diagnostics."""

from __future__ import annotations

import ast
from dataclasses import dataclass

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models.registry import get_model_info


@dataclass(frozen=True)
class MatrixDiagnostics:
    """Summary of matrix structure, storage, and Hermiticity."""

    shape: tuple[int, int]
    sparse: bool
    nonzero_entries: int
    density: float
    storage_bytes: int
    hermitian: bool


def estimate_model_dimension(model_name: str, **parameters: object) -> int:
    """Estimate a registered model's matrix dimension from its dimension expression."""

    info = get_model_info(model_name)
    values = dict(info.defaults)
    values.update(parameters)
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
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Add, ast.Mult, ast.Pow)):
            left = evaluate(node.left)
            right = evaluate(node.right)
            if isinstance(node.op, ast.Add):
                return left + right
            if isinstance(node.op, ast.Mult):
                return left * right
            return left**right
        raise ValueError(f"Unsupported dimension expression {expression!r}.")

    return evaluate(tree)
