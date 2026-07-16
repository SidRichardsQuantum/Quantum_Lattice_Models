"""Gauge-tolerant topological invariants for Bloch Hamiltonians."""

from __future__ import annotations

from collections.abc import Callable, Sequence
from typing import TYPE_CHECKING

import numpy as np

from quantum_lattice_models.periodic import PeriodicLatticeSpec, bloch_function

if TYPE_CHECKING:
    from quantum_lattice_models.analysis import AnalysisResult


def zak_phase(
    source: PeriodicLatticeSpec | Callable[[np.ndarray], np.ndarray],
    *,
    band: int = 0,
    n_points: int = 401,
    reciprocal_period: Sequence[float] | None = None,
) -> float:
    """Return the discretized Zak phase for one isolated one-dimensional band."""

    if n_points < 3:
        raise ValueError("n_points must be at least 3.")
    builder = bloch_function(source)
    if reciprocal_period is None:
        if not isinstance(source, PeriodicLatticeSpec) or source.dimension != 1:
            raise ValueError("reciprocal_period is required for a callable Bloch Hamiltonian.")
        period = source.reciprocal_vectors()[0]
    else:
        period = np.asarray(reciprocal_period, dtype=float)
    states = []
    for fraction in np.arange(n_points) / n_points:
        _, vectors = np.linalg.eigh(builder(fraction * period))
        if not 0 <= band < vectors.shape[1]:
            raise ValueError("band index is out of range.")
        states.append(vectors[:, band])
    product = 1.0 + 0.0j
    for index, state in enumerate(states):
        overlap = np.vdot(state, states[(index + 1) % len(states)])
        if abs(overlap) < 1e-12:
            raise ValueError("Zak phase is undefined because adjacent band states are orthogonal.")
        product *= overlap / abs(overlap)
    return float(np.angle(product))


def winding_number(
    source: PeriodicLatticeSpec | Callable[[np.ndarray], np.ndarray],
    *,
    n_points: int = 801,
    reciprocal_period: Sequence[float] | None = None,
    tolerance: float = 1e-10,
) -> int:
    """Return the winding of a two-band chiral Hamiltonian's off-diagonal term."""

    builder = bloch_function(source)
    if reciprocal_period is None:
        if not isinstance(source, PeriodicLatticeSpec) or source.dimension != 1:
            raise ValueError("reciprocal_period is required for a callable Bloch Hamiltonian.")
        period = source.reciprocal_vectors()[0]
    else:
        period = np.asarray(reciprocal_period, dtype=float)
    phases = []
    for fraction in np.linspace(0.0, 1.0, n_points):
        matrix = np.asarray(builder(fraction * period), dtype=complex)
        if matrix.shape != (2, 2):
            raise ValueError("Winding number currently requires a two-band Hamiltonian.")
        if abs(matrix[0, 0] - matrix[1, 1]) > tolerance:
            raise ValueError("Winding number requires equal diagonal terms (chiral symmetry).")
        value = matrix[0, 1]
        if abs(value) < tolerance:
            raise ValueError("Winding number is undefined at a bulk gap closing.")
        phases.append(np.angle(value))
    return int(np.rint((np.unwrap(phases)[-1] - np.unwrap(phases)[0]) / (2.0 * np.pi)))


def chern_number(
    source: PeriodicLatticeSpec | Callable[[np.ndarray], np.ndarray],
    *,
    occupied_bands: int = 1,
    mesh: tuple[int, int] = (31, 31),
    reciprocal_vectors: np.ndarray | None = None,
) -> float:
    """Return the Fukui-Hatsugai-Suzuki Chern number on a periodic momentum mesh."""

    if min(mesh) < 3:
        raise ValueError("Chern mesh dimensions must be at least 3.")
    builder = bloch_function(source)
    if reciprocal_vectors is None:
        if not isinstance(source, PeriodicLatticeSpec) or source.dimension != 2:
            raise ValueError("reciprocal_vectors are required for a callable Hamiltonian.")
        reciprocal = source.reciprocal_vectors()
    else:
        reciprocal = np.asarray(reciprocal_vectors, dtype=float)
    if reciprocal.shape != (2, 2):
        raise ValueError("reciprocal_vectors must have shape (2, 2).")
    states: list[list[np.ndarray]] = []
    for i in range(mesh[0]):
        row = []
        for j in range(mesh[1]):
            momentum = (i / mesh[0]) * reciprocal[0] + (j / mesh[1]) * reciprocal[1]
            _, vectors = np.linalg.eigh(builder(momentum))
            if not 1 <= occupied_bands < vectors.shape[1]:
                raise ValueError("occupied_bands must select a nonempty proper band subspace.")
            row.append(vectors[:, :occupied_bands])
        states.append(row)

    curvature = 0.0
    for i in range(mesh[0]):
        for j in range(mesh[1]):
            state = states[i][j]
            state_x = states[(i + 1) % mesh[0]][j]
            state_y = states[i][(j + 1) % mesh[1]]
            state_xy = states[(i + 1) % mesh[0]][(j + 1) % mesh[1]]
            ux = _link(state, state_x)
            uy = _link(state, state_y)
            uy_x = _link(state_x, state_xy)
            ux_y = _link(state_y, state_xy)
            curvature += np.angle(ux * uy_x / (ux_y * uy))
    return float(curvature / (2.0 * np.pi))


def zak_phase_result(
    source: PeriodicLatticeSpec,
    *,
    band: int = 0,
    n_points: int = 401,
) -> AnalysisResult:
    """Return the Zak phase as a portable analysis result."""

    from quantum_lattice_models.analysis import topology_result

    value = zak_phase(source, band=band, n_points=n_points)
    return topology_result(
        "zak",
        value,
        periodic=source,
        parameters={"band": band, "n_points": n_points},
        solver={"method": "discrete Wilson loop", "gauge_tolerant": True},
    )


def winding_number_result(
    source: PeriodicLatticeSpec,
    *,
    n_points: int = 801,
    tolerance: float = 1e-10,
) -> AnalysisResult:
    """Return the chiral winding number as a portable analysis result."""

    from quantum_lattice_models.analysis import topology_result

    value = winding_number(source, n_points=n_points, tolerance=tolerance)
    return topology_result(
        "winding",
        value,
        periodic=source,
        parameters={"n_points": n_points, "tolerance": tolerance},
        solver={"method": "unwrapped off-diagonal phase", "gauge_tolerant": True},
    )


def chern_number_result(
    source: PeriodicLatticeSpec,
    *,
    occupied_bands: int = 1,
    mesh: tuple[int, int] = (31, 31),
) -> AnalysisResult:
    """Return the occupied-subspace Chern number as a portable analysis result."""

    from quantum_lattice_models.analysis import topology_result

    value = chern_number(source, occupied_bands=occupied_bands, mesh=mesh)
    return topology_result(
        "chern",
        value,
        periodic=source,
        parameters={"occupied_bands": occupied_bands, "mesh": list(mesh)},
        solver={"method": "Fukui-Hatsugai-Suzuki", "gauge_tolerant": True},
    )


def _link(left: np.ndarray, right: np.ndarray) -> complex:
    determinant = np.linalg.det(left.conj().T @ right)
    if abs(determinant) < 1e-14:
        raise ValueError("Topological invariant encountered a singular occupied-band overlap.")
    return complex(determinant / abs(determinant))
