"""Coordinate helpers for plotting finite lattice geometries."""

from __future__ import annotations

import numpy as np

from quantum_lattice_models._model_utils import validate_positive_int


def square_lattice_positions(n_rows: int, n_cols: int) -> np.ndarray:
    """Return row-major ``(x, y)`` coordinates for a square lattice."""

    validate_positive_int(n_rows, "n_rows")
    validate_positive_int(n_cols, "n_cols")
    return np.array([(col, row) for row in range(n_rows) for col in range(n_cols)], dtype=float)


def triangular_lattice_positions(n_rows: int, n_cols: int) -> np.ndarray:
    """Return row-major ``(x, y)`` coordinates for a triangular lattice."""

    validate_positive_int(n_rows, "n_rows")
    validate_positive_int(n_cols, "n_cols")
    return np.array(
        [
            (col + 0.5 * row, np.sqrt(3.0) * row / 2.0)
            for row in range(n_rows)
            for col in range(n_cols)
        ],
        dtype=float,
    )


def honeycomb_lattice_positions(n_rows: int, n_cols: int) -> np.ndarray:
    """Return row-major two-sublattice coordinates for a honeycomb lattice."""

    validate_positive_int(n_rows, "n_rows")
    validate_positive_int(n_cols, "n_cols")
    offsets = np.array([[0.0, 0.0], [0.0, 1.0 / np.sqrt(3.0)]])
    positions = []
    for row in range(n_rows):
        for col in range(n_cols):
            origin = np.array([1.5 * col, np.sqrt(3.0) * row + (np.sqrt(3.0) / 2.0) * col])
            positions.extend(origin + offsets)
    return np.asarray(positions, dtype=float)


def kagome_lattice_positions(n_rows: int, n_cols: int) -> np.ndarray:
    """Return row-major three-sublattice coordinates for a kagome lattice."""

    validate_positive_int(n_rows, "n_rows")
    validate_positive_int(n_cols, "n_cols")
    offsets = np.array([[0.0, 0.0], [0.5, 0.0], [0.25, np.sqrt(3.0) / 4.0]])
    positions = []
    for row in range(n_rows):
        for col in range(n_cols):
            origin = np.array([col + 0.5 * row, np.sqrt(3.0) * row / 2.0])
            positions.extend(origin + offsets)
    return np.asarray(positions, dtype=float)
