"""Structured solvers, dynamics, sweeps, thermal tools, and reciprocal analysis."""

from __future__ import annotations

from collections.abc import Callable, Sequence

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from quantum_lattice_models.analysis import AnalysisResult
from quantum_lattice_models.periodic import PeriodicLatticeSpec


def solve_eigenpairs(
    matrix: np.ndarray | sp.spmatrix,
    *,
    k: int | None = None,
    which: str = "SA",
    tolerance: float = 0.0,
    max_iterations: int | None = None,
    dense_threshold: int = 512,
    allow_dense: bool = False,
) -> AnalysisResult:
    """Solve eigenpairs with explicit sparse-to-dense safeguards and residuals."""

    n = matrix.shape[0]
    if matrix.shape != (n, n):
        raise ValueError("matrix must be square.")
    sparse = sp.issparse(matrix)
    full = k is None or k >= n
    warnings = []
    if sparse and full and n > dense_threshold and not allow_dense:
        raise ValueError(
            "Full diagonalization would densify a large sparse matrix; "
            "provide k or allow_dense=True."
        )
    if k is not None and not 1 <= k <= n:
        raise ValueError("k must satisfy 1 <= k <= matrix dimension.")
    method = "dense"
    exact = True
    if sparse and not full:
        method = "scipy.sparse.linalg.eigsh"
        exact = False
        values, vectors = spla.eigsh(
            matrix.asfptype(),
            k=k,
            which=which,
            tol=tolerance,
            maxiter=max_iterations,
        )
        order = np.argsort(values.real)
        values, vectors = values[order], vectors[:, order]
    else:
        dense = matrix.toarray() if sparse else np.asarray(matrix)
        if sparse:
            warnings.append("Sparse matrix converted to dense for complete diagonalization.")
        if np.allclose(dense, dense.conj().T):
            values, vectors = np.linalg.eigh(dense)
            method = "numpy.linalg.eigh"
        else:
            values, vectors = np.linalg.eig(dense)
            order = np.argsort(values.real)
            values, vectors = values[order], vectors[:, order]
            method = "numpy.linalg.eig"
    residuals = np.asarray(
        [
            np.linalg.norm(matrix @ vectors[:, index] - values[index] * vectors[:, index])
            for index in range(vectors.shape[1])
        ]
    )
    degeneracies = _degeneracy_groups(values)
    return AnalysisResult(
        analysis="eigensolver",
        coordinates={"index": np.arange(len(values))},
        values={
            "eigenvalues": values,
            "eigenvectors": vectors,
            "residuals": residuals,
        },
        parameters={
            "k": k,
            "which": which,
            "tolerance": tolerance,
            "max_iterations": max_iterations,
        },
        solver={
            "method": method,
            "exact": exact,
            "converged": True,
            "max_residual": float(residuals.max(initial=0.0)),
        },
        warnings=tuple(warnings),
        plot={"kind": "spectrum", "x": "index", "y": "eigenvalues"},
        metadata={
            "dimension": n,
            "sparse_input": sparse,
            "estimated_dense_bytes": n * n * 16,
            "degeneracy_groups": degeneracies,
        },
    )


def evolve_state(
    hamiltonian: np.ndarray | sp.spmatrix,
    initial_state: np.ndarray,
    times: Sequence[float],
    *,
    observables: dict[str, np.ndarray | sp.spmatrix] | None = None,
) -> AnalysisResult:
    """Evolve a state under a time-independent Hamiltonian."""

    time_values = np.asarray(times, dtype=float)
    if time_values.ndim != 1 or len(time_values) < 1 or np.any(np.diff(time_values) < 0):
        raise ValueError("times must be a nondecreasing one-dimensional sequence.")
    state = np.asarray(initial_state, dtype=complex).reshape(-1)
    if hamiltonian.shape != (state.size, state.size):
        raise ValueError("Hamiltonian shape must match initial-state dimension.")
    if not np.isclose(np.vdot(state, state).real, 1.0):
        norm = np.linalg.norm(state)
        if norm == 0:
            raise ValueError("initial_state must have nonzero norm.")
        state = state / norm
    if sp.issparse(hamiltonian):
        states = np.asarray(
            [spla.expm_multiply(-1j * time * hamiltonian, state) for time in time_values]
        )
        method = "scipy.sparse.linalg.expm_multiply"
    else:
        matrix = np.asarray(hamiltonian)
        states = np.asarray([la.expm(-1j * time * matrix) @ state for time in time_values])
        method = "scipy.linalg.expm"
    values: dict[str, np.ndarray] = {"states": states}
    for name, operator in (observables or {}).items():
        values[name] = np.asarray([np.vdot(vector, operator @ vector) for vector in states])
    overlaps = states @ state.conj()
    values["loschmidt_amplitude"] = overlaps
    values["loschmidt_echo"] = np.abs(overlaps) ** 2
    primary = next(iter(observables), "loschmidt_echo") if observables else "loschmidt_echo"
    return AnalysisResult(
        analysis="time_evolution",
        coordinates={"time": time_values},
        values=values,
        solver={"method": method, "exact": not sp.issparse(hamiltonian)},
        plot={"kind": "time_series", "x": "time", "y": primary},
        metadata={"dimension": state.size},
    )


def evolve_time_dependent(
    hamiltonian: Callable[[float], np.ndarray | sp.spmatrix],
    initial_state: np.ndarray,
    times: Sequence[float],
    *,
    observables: dict[str, np.ndarray | sp.spmatrix] | None = None,
    method: str = "midpoint",
    substeps: int = 1,
) -> AnalysisResult:
    """Evolve under ``H(t)`` using piecewise midpoint or fourth-order propagation."""

    time_values = np.asarray(times, dtype=float)
    if time_values.ndim != 1 or len(time_values) < 1 or np.any(np.diff(time_values) < 0):
        raise ValueError("times must be a nondecreasing one-dimensional sequence.")
    if method not in {"midpoint", "rk4"}:
        raise ValueError("method must be 'midpoint' or 'rk4'.")
    if not isinstance(substeps, int) or substeps < 1:
        raise ValueError("substeps must be a positive integer.")
    state = np.asarray(initial_state, dtype=complex).reshape(-1)
    norm = np.linalg.norm(state)
    if norm == 0:
        raise ValueError("initial_state must have nonzero norm.")
    state = state / norm
    states = [state.copy()]
    for start, stop in zip(time_values[:-1], time_values[1:], strict=True):
        step = (stop - start) / substeps
        for substep in range(substeps):
            left = start + substep * step
            midpoint = left + 0.5 * step
            if method == "midpoint":
                matrix = hamiltonian(midpoint)
                state = spla.expm_multiply(-1j * step * matrix, state)
            else:

                def derivative(time: float, vector: np.ndarray) -> np.ndarray:
                    return -1j * (hamiltonian(time) @ vector)

                k1 = derivative(left, state)
                k2 = derivative(midpoint, state + 0.5 * step * k1)
                k3 = derivative(midpoint, state + 0.5 * step * k2)
                k4 = derivative(left + step, state + step * k3)
                state = state + step * (k1 + 2 * k2 + 2 * k3 + k4) / 6
                state /= np.linalg.norm(state)
        states.append(state.copy())
    state_array = np.asarray(states)
    values: dict[str, np.ndarray] = {
        "states": state_array,
        "state_norm": np.linalg.norm(state_array, axis=1),
    }
    for name, operator in (observables or {}).items():
        values[name] = np.asarray([np.vdot(vector, operator @ vector) for vector in state_array])
    return AnalysisResult(
        analysis="time_dependent_evolution",
        coordinates={"time": time_values},
        values=values,
        parameters={"method": method, "substeps": substeps},
        solver={"method": method, "order": 4 if method == "rk4" else 2},
        plot={
            "kind": "time_series",
            "x": "time",
            "y": next(iter(observables), "state_norm") if observables else "state_norm",
        },
        metadata={"dimension": state.size},
    )


def quench_dynamics(
    initial_hamiltonian: np.ndarray | sp.spmatrix,
    final_hamiltonian: np.ndarray | sp.spmatrix,
    times: Sequence[float],
    *,
    observables: dict[str, np.ndarray | sp.spmatrix] | None = None,
) -> AnalysisResult:
    """Evolve the initial Hamiltonian's ground state under a final Hamiltonian."""

    initial = solve_eigenpairs(initial_hamiltonian, k=1)
    state = initial.values["eigenvectors"][:, 0]
    result = evolve_state(final_hamiltonian, state, times, observables=observables)
    return AnalysisResult(
        analysis="quench_dynamics",
        coordinates=result.coordinates,
        values=result.values,
        parameters={"initial_ground_energy": initial.values["eigenvalues"][0]},
        solver=result.solver,
        warnings=result.warnings,
        plot=result.plot,
        metadata=result.metadata,
    )


def parameter_sweep(
    builder: Callable[..., np.ndarray | sp.spmatrix],
    parameter: str,
    values: Sequence[float],
    analyzer: Callable[[np.ndarray | sp.spmatrix], float],
    *,
    builder_parameters: dict[str, object] | None = None,
    analysis_name: str = "parameter_sweep",
) -> AnalysisResult:
    """Evaluate a scalar analysis across one model parameter."""

    coordinates = np.asarray(values)
    base = dict(builder_parameters or {})
    measured = []
    for value in coordinates:
        measured.append(analyzer(builder(**base, **{parameter: value.item()})))
    measured_array = np.asarray(measured)
    extrema = {
        "minimum_index": int(np.argmin(measured_array)),
        "maximum_index": int(np.argmax(measured_array)),
    }
    return AnalysisResult(
        analysis=analysis_name,
        coordinates={parameter: coordinates},
        values={"measurement": measured_array},
        parameters={"parameter": parameter, "builder_parameters": base},
        plot={"kind": "sweep", "x": parameter, "y": "measurement"},
        metadata=extrema,
    )


def two_parameter_sweep(
    builder: Callable[..., np.ndarray | sp.spmatrix],
    x_parameter: str,
    x_values: Sequence[float],
    y_parameter: str,
    y_values: Sequence[float],
    analyzer: Callable[[np.ndarray | sp.spmatrix], float],
    *,
    builder_parameters: dict[str, object] | None = None,
) -> AnalysisResult:
    """Evaluate a scalar analysis on a two-parameter grid."""

    x = np.asarray(x_values, dtype=float)
    y = np.asarray(y_values, dtype=float)
    base = dict(builder_parameters or {})
    measured = np.empty((len(y), len(x)), dtype=float)
    for row, y_value in enumerate(y):
        for column, x_value in enumerate(x):
            matrix = builder(
                **base,
                **{x_parameter: float(x_value), y_parameter: float(y_value)},
            )
            measured[row, column] = analyzer(matrix)
    return AnalysisResult(
        analysis="two_parameter_sweep",
        coordinates={x_parameter: x, y_parameter: y},
        values={"measurement": measured},
        parameters={"x_parameter": x_parameter, "y_parameter": y_parameter},
        plot={
            "kind": "heatmap",
            "x": x_parameter,
            "y": y_parameter,
            "values": "measurement",
        },
    )


def finite_size_sweep(
    builder: Callable[..., np.ndarray | sp.spmatrix],
    sizes: Sequence[int],
    analyzer: Callable[[np.ndarray | sp.spmatrix], float],
    *,
    size_parameter: str = "n_sites",
    builder_parameters: dict[str, object] | None = None,
) -> AnalysisResult:
    """Evaluate a scalar quantity over finite system sizes."""

    sizes_array = np.asarray(sizes, dtype=int)
    result = parameter_sweep(
        builder,
        size_parameter,
        sizes_array,
        analyzer,
        builder_parameters=builder_parameters,
        analysis_name="finite_size_sweep",
    )
    return result


def gap_closing_points(
    coordinates: Sequence[float],
    gaps: Sequence[float],
    *,
    threshold: float = 1e-6,
) -> AnalysisResult:
    """Identify threshold crossings and local minima in sampled spectral gaps."""

    x = np.asarray(coordinates, dtype=float)
    values = np.asarray(gaps, dtype=float)
    if x.ndim != 1 or values.shape != x.shape or x.size < 2:
        raise ValueError("coordinates and gaps must be matching one-dimensional arrays.")
    minima = [
        index
        for index in range(1, len(values) - 1)
        if values[index] <= values[index - 1] and values[index] <= values[index + 1]
    ]
    threshold_indices = np.flatnonzero(values <= threshold)
    return AnalysisResult(
        analysis="gap_closing_detection",
        coordinates={"parameter": x},
        values={"gap": values},
        parameters={"threshold": threshold},
        plot={"kind": "sweep", "x": "parameter", "y": "gap"},
        metadata={
            "minimum_index": int(np.argmin(values)),
            "local_minimum_indices": minima,
            "threshold_indices": threshold_indices.tolist(),
        },
    )


def susceptibility(
    builder: Callable[..., np.ndarray | sp.spmatrix],
    field_parameter: str,
    field_value: float,
    observable: np.ndarray | sp.spmatrix,
    *,
    step: float = 1e-4,
    builder_parameters: dict[str, object] | None = None,
) -> AnalysisResult:
    """Estimate zero-temperature response to an explicit conjugate field."""

    if step <= 0:
        raise ValueError("step must be positive.")
    base = dict(builder_parameters or {})

    def expectation_at(value: float) -> float:
        matrix = builder(**base, **{field_parameter: value})
        dense = matrix.toarray() if sp.issparse(matrix) else np.asarray(matrix)
        _, vectors = np.linalg.eigh(dense)
        ground = vectors[:, 0]
        return float(np.real(np.vdot(ground, observable @ ground)))

    minus = expectation_at(field_value - step)
    center = expectation_at(field_value)
    plus = expectation_at(field_value + step)
    response = (plus - minus) / (2 * step)
    return AnalysisResult(
        analysis="susceptibility",
        coordinates={"field": np.asarray([field_value - step, field_value, field_value + step])},
        values={
            "expectation": np.asarray([minus, center, plus]),
            "susceptibility": np.asarray([response]),
        },
        parameters={
            "field_parameter": field_parameter,
            "field_value": field_value,
            "step": step,
            "builder_parameters": base,
        },
        solver={"method": "central finite difference", "exact_ground_state": True},
        plot={"kind": "sweep", "x": "field", "y": "expectation"},
    )


def thermal_observables(
    hamiltonian: np.ndarray | sp.spmatrix,
    temperatures: Sequence[float],
    *,
    observable: np.ndarray | sp.spmatrix | None = None,
) -> AnalysisResult:
    """Return canonical small-system thermodynamics with stable shifted weights."""

    temperatures_array = np.asarray(temperatures, dtype=float)
    if np.any(temperatures_array <= 0):
        raise ValueError("temperatures must be positive.")
    dense = hamiltonian.toarray() if sp.issparse(hamiltonian) else np.asarray(hamiltonian)
    energies, vectors = np.linalg.eigh(dense)
    shifted = energies - energies.min()
    partition = []
    internal = []
    free = []
    entropy = []
    heat_capacity = []
    expectation_values = []
    diagonal_observable = None
    if observable is not None:
        operator = observable.toarray() if sp.issparse(observable) else np.asarray(observable)
        diagonal_observable = np.einsum("ij,ji->i", vectors.conj().T, operator @ vectors).real
    for temperature in temperatures_array:
        beta = 1.0 / temperature
        weights = np.exp(-beta * shifted)
        normalized = weights / weights.sum()
        mean = float(np.dot(normalized, energies))
        variance = float(np.dot(normalized, energies**2) - mean**2)
        log_z = -beta * energies.min() + np.log(weights.sum())
        partition.append(np.exp(min(log_z, 700.0)))
        internal.append(mean)
        free.append(-temperature * log_z)
        entropy.append((mean + temperature * log_z) / temperature)
        heat_capacity.append(variance / temperature**2)
        if diagonal_observable is not None:
            expectation_values.append(float(np.dot(normalized, diagonal_observable)))
    values = {
        "partition_function": np.asarray(partition),
        "internal_energy": np.asarray(internal),
        "free_energy": np.asarray(free),
        "entropy": np.asarray(entropy),
        "heat_capacity": np.asarray(heat_capacity),
    }
    if expectation_values:
        values["observable_expectation"] = np.asarray(expectation_values)
    return AnalysisResult(
        analysis="thermal_observables",
        coordinates={"temperature": temperatures_array},
        values=values,
        solver={"method": "complete eigendecomposition", "exact": True},
        plot={"kind": "time_series", "x": "temperature", "y": "internal_energy"},
        warnings=(
            "Intended for small systems requiring a complete spectrum.",
            "Partition functions are capped at exp(700) to remain finite.",
        ),
    )


def berry_curvature(
    lattice: PeriodicLatticeSpec,
    *,
    occupied_bands: int = 1,
    mesh: tuple[int, int] = (31, 31),
) -> AnalysisResult:
    """Return gauge-invariant plaquette Berry curvature over a reciprocal mesh."""

    reciprocal = lattice.reciprocal_vectors()
    states: list[list[np.ndarray]] = []
    kx = np.arange(mesh[0]) / mesh[0]
    ky = np.arange(mesh[1]) / mesh[1]
    for x in kx:
        row = []
        for y in ky:
            momentum = x * reciprocal[0] + y * reciprocal[1]
            _, vectors = np.linalg.eigh(lattice.bloch_hamiltonian(momentum))
            row.append(vectors[:, :occupied_bands])
        states.append(row)
    curvature = np.empty(mesh)
    for i in range(mesh[0]):
        for j in range(mesh[1]):
            u = states[i][j]
            ux = states[(i + 1) % mesh[0]][j]
            uy = states[i][(j + 1) % mesh[1]]
            uxy = states[(i + 1) % mesh[0]][(j + 1) % mesh[1]]
            curvature[i, j] = np.angle(
                _link(u, ux) * _link(ux, uxy) / (_link(uy, uxy) * _link(u, uy))
            )
    return AnalysisResult(
        analysis="berry_curvature",
        coordinates={"k1_fraction": kx, "k2_fraction": ky},
        values={"curvature": curvature},
        parameters={"occupied_bands": occupied_bands, "mesh": list(mesh)},
        source={"kind": "periodic_lattice", "periodic_lattice": lattice.to_dict()},
        solver={"method": "Fukui-Hatsugai-Suzuki plaquettes", "gauge_tolerant": True},
        plot={
            "kind": "heatmap",
            "x": "k1_fraction",
            "y": "k2_fraction",
            "values": "curvature",
        },
        metadata={"chern_number": float(curvature.sum() / (2 * np.pi))},
    )


def berry_curvature_convergence(
    lattice: PeriodicLatticeSpec,
    meshes: Sequence[int | tuple[int, int]],
    *,
    occupied_bands: int = 1,
) -> AnalysisResult:
    """Report Chern-number convergence across reciprocal meshes."""

    normalized = [(mesh, mesh) if isinstance(mesh, int) else tuple(mesh) for mesh in meshes]
    if not normalized or any(len(mesh) != 2 or min(mesh) < 2 for mesh in normalized):
        raise ValueError("meshes must contain positive two-dimensional mesh sizes.")
    chern_values = np.asarray(
        [
            berry_curvature(lattice, occupied_bands=occupied_bands, mesh=mesh).metadata[
                "chern_number"
            ]
            for mesh in normalized
        ],
        dtype=float,
    )
    return AnalysisResult(
        analysis="berry_curvature_convergence",
        coordinates={"mesh_points": np.asarray([a * b for a, b in normalized])},
        values={"chern_number": chern_values},
        parameters={
            "occupied_bands": occupied_bands,
            "meshes": [list(mesh) for mesh in normalized],
        },
        source={"kind": "periodic_lattice", "periodic_lattice": lattice.to_dict()},
        solver={"method": "Fukui-Hatsugai-Suzuki mesh sequence"},
        plot={"kind": "sweep", "x": "mesh_points", "y": "chern_number"},
        metadata={
            "successive_changes": np.abs(np.diff(chern_values)).tolist(),
            "final_value": float(chern_values[-1]),
        },
    )


def wilson_loop(
    lattice: PeriodicLatticeSpec,
    transverse_values: Sequence[float],
    *,
    occupied_bands: int = 1,
    loop_points: int = 101,
) -> AnalysisResult:
    """Return occupied-subspace Wilson-loop eigenphases along reciprocal vector one."""

    reciprocal = lattice.reciprocal_vectors()
    transverse = np.asarray(transverse_values, dtype=float)
    phases = np.empty((len(transverse), occupied_bands))
    for row, transverse_value in enumerate(transverse):
        product = np.eye(occupied_bands, dtype=complex)
        previous = None
        first = None
        for fraction in np.arange(loop_points) / loop_points:
            momentum = fraction * reciprocal[0] + transverse_value * reciprocal[1]
            _, vectors = np.linalg.eigh(lattice.bloch_hamiltonian(momentum))
            occupied = vectors[:, :occupied_bands]
            first = occupied if first is None else first
            if previous is not None:
                overlap = previous.conj().T @ occupied
                u, _, vh = np.linalg.svd(overlap)
                product = product @ (u @ vh)
            previous = occupied
        overlap = previous.conj().T @ first
        u, _, vh = np.linalg.svd(overlap)
        product = product @ (u @ vh)
        phases[row] = np.sort(np.angle(np.linalg.eigvals(product)))
    return AnalysisResult(
        analysis="wilson_loop",
        coordinates={"transverse_fraction": transverse},
        values={"phases": phases},
        parameters={"occupied_bands": occupied_bands, "loop_points": loop_points},
        source={"kind": "periodic_lattice", "periodic_lattice": lattice.to_dict()},
        solver={"method": "unitarized overlap Wilson loop", "gauge_tolerant": True},
        plot={"kind": "time_series", "x": "transverse_fraction", "y": "phases"},
    )


def wilson_loop_winding(
    lattice: PeriodicLatticeSpec,
    transverse_values: Sequence[float],
    *,
    occupied_bands: int = 1,
    loop_points: int = 101,
) -> AnalysisResult:
    """Estimate occupied-subspace Wilson-loop determinant winding."""

    loops = wilson_loop(
        lattice,
        transverse_values,
        occupied_bands=occupied_bands,
        loop_points=loop_points,
    )
    determinant_phase = np.unwrap(np.sum(loops.values["phases"], axis=1))
    winding = float((determinant_phase[-1] - determinant_phase[0]) / (2 * np.pi))
    return AnalysisResult(
        analysis="wilson_loop_winding",
        coordinates=loops.coordinates,
        values={"determinant_phase": determinant_phase, "winding": np.asarray([winding])},
        parameters=loops.parameters,
        source=loops.source,
        solver={"method": "occupied-subspace Wilson-loop determinant"},
        plot={
            "kind": "time_series",
            "x": "transverse_fraction",
            "y": "determinant_phase",
        },
        metadata={"nearest_integer": int(np.rint(winding))},
    )


def reciprocal_space_data(lattice: PeriodicLatticeSpec) -> AnalysisResult:
    """Return primitive reciprocal vectors and a first-zone polygon."""

    reciprocal = lattice.reciprocal_vectors()
    if lattice.dimension != 2:
        raise ValueError("Reciprocal-space diagrams currently require two dimensions.")
    candidates = np.asarray(
        [
            i * reciprocal[0] + j * reciprocal[1]
            for i in range(-1, 2)
            for j in range(-1, 2)
            if (i, j) != (0, 0)
        ]
    )
    angles = np.arctan2(candidates[:, 1], candidates[:, 0])
    polygon = candidates[np.argsort(angles)]
    return AnalysisResult(
        analysis="reciprocal_space",
        values={"reciprocal_vectors": reciprocal, "zone_points": polygon},
        source={"kind": "periodic_lattice", "periodic_lattice": lattice.to_dict()},
        plot={"kind": "reciprocal", "vectors": "reciprocal_vectors", "zone": "zone_points"},
    )


def _link(left: np.ndarray, right: np.ndarray) -> complex:
    determinant = np.linalg.det(left.conj().T @ right)
    return determinant / abs(determinant)


def _degeneracy_groups(values: np.ndarray) -> list[list[int]]:
    groups: list[list[int]] = []
    for index, value in enumerate(values):
        if not groups or not np.isclose(value, values[groups[-1][0]]):
            groups.append([index])
        else:
            groups[-1].append(index)
    return [group for group in groups if len(group) > 1]
