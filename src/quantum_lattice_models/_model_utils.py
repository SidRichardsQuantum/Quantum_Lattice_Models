"""Shared helpers for lattice model builders."""

from __future__ import annotations

from collections.abc import Iterable

import numpy as np


def nearest_neighbor_bonds(n_sites: int, periodic: bool) -> list[tuple[int, int]]:
    bonds = [(i, i + 1) for i in range(n_sites - 1)]
    if periodic and n_sites > 2:
        bonds.append((n_sites - 1, 0))
    elif periodic and n_sites == 2:
        bonds.append((1, 0))
    return bonds


def next_nearest_neighbor_bonds(n_sites: int, periodic: bool) -> list[tuple[int, int]]:
    bonds = [(i, i + 2) for i in range(n_sites - 2)]
    if periodic and n_sites > 3:
        bonds.extend([(n_sites - 2, 0), (n_sites - 1, 1)])
    return bonds


def add_symmetric_hopping(matrix: np.ndarray, i: int, j: int, value: complex) -> None:
    matrix[i, j] += value
    matrix[j, i] += np.conjugate(value)


def validate_positive_int(value: int, name: str) -> None:
    if not isinstance(value, int) or value < 1:
        raise ValueError(f"{name} must be a positive integer.")


def pauli_labels(n_sites: int, labels: dict[int, str]) -> tuple[str, ...]:
    return tuple(labels.get(i, "I") for i in range(n_sites))


def onsite_values(onsite: float | Iterable[float], n_sites: int) -> np.ndarray:
    if np.isscalar(onsite):
        return np.full(n_sites, onsite, dtype=float)
    values = np.asarray(list(onsite), dtype=float)
    if values.shape != (n_sites,):
        raise ValueError("onsite must be a scalar or a length-n_sites iterable.")
    return values


def square_lattice_index_raw(row: int, col: int, n_cols: int) -> int:
    return row * n_cols + col


def honeycomb_index_raw(row: int, col: int, sublattice: int, n_cols: int) -> int:
    return 2 * (row * n_cols + col) + sublattice


def kagome_index_raw(row: int, col: int, sublattice: int, n_cols: int) -> int:
    return 3 * (row * n_cols + col) + sublattice


def wrapped_cell(
    row: int,
    col: int,
    n_rows: int,
    n_cols: int,
    periodic_y: bool,
    periodic_x: bool,
) -> tuple[int | None, int | None]:
    if periodic_y:
        row %= n_rows
    if periodic_x:
        col %= n_cols
    if not (0 <= row < n_rows and 0 <= col < n_cols):
        return None, None
    return row, col


def validate_sublattice_index(
    row: int, col: int, sublattice: int, n_cols: int, n_sublattices: int
) -> None:
    validate_positive_int(n_cols, "n_cols")
    if row < 0 or col < 0 or col >= n_cols:
        raise ValueError("row and col must be nonnegative, and col must be less than n_cols.")
    if not 0 <= sublattice < n_sublattices:
        raise ValueError(f"sublattice must satisfy 0 <= sublattice < {n_sublattices}.")


def bit_is_set(state: int, bit: int) -> int:
    return (state >> bit) & 1


def fermion_hop_sign(state: int, source: int, target: int) -> int:
    lower = min(source, target) + 1
    upper = max(source, target)
    parity = sum(bit_is_set(state, orbital) for orbital in range(lower, upper))
    return -1 if parity % 2 else 1


def occupation_digits(state: int, n_sites: int, local_dim: int) -> list[int]:
    occupations = []
    for _ in range(n_sites):
        occupations.append(state % local_dim)
        state //= local_dim
    return occupations


def occupations_to_index(occupations: list[int], local_dim: int) -> int:
    index = 0
    multiplier = 1
    for occupation in occupations:
        index += occupation * multiplier
        multiplier *= local_dim
    return index
