"""Spin observables, correlations, and entanglement for full and reduced bases."""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np

from quantum_lattice_models.spin import FixedMagnetizationBasis


def expectation(state: np.ndarray, operator: np.ndarray) -> complex:
    """Return ``<state|operator|state>`` for a state vector."""

    vector = np.asarray(state, dtype=complex).reshape(-1)
    matrix = np.asarray(operator, dtype=complex)
    if matrix.shape != (vector.size, vector.size):
        raise ValueError("operator shape must match the state-vector dimension.")
    return complex(np.vdot(vector, matrix @ vector))


def site_magnetization_z(
    state: np.ndarray,
    n_sites: int,
    *,
    basis: FixedMagnetizationBasis | None = None,
) -> np.ndarray:
    """Return the site-resolved expectations ``<Z_i>``."""

    return np.asarray(
        [
            _pauli_expectation(state, n_sites, {site: "Z"}, basis=basis).real
            for site in range(n_sites)
        ]
    )


def total_magnetization_z(
    state: np.ndarray,
    n_sites: int,
    *,
    basis: FixedMagnetizationBasis | None = None,
) -> float:
    """Return the total Pauli-Z expectation ``sum_i <Z_i>``."""

    return float(np.sum(site_magnetization_z(state, n_sites, basis=basis)))


def magnetization_z(
    state: np.ndarray,
    n_sites: int,
    *,
    basis: FixedMagnetizationBasis | None = None,
) -> float:
    """Return average Pauli-Z magnetization per site."""

    return float(np.mean(site_magnetization_z(state, n_sites, basis=basis)))


def spin_correlation_matrix(
    state: np.ndarray,
    n_sites: int,
    *,
    axis: str = "Z",
    basis: FixedMagnetizationBasis | None = None,
    connected: bool = False,
) -> np.ndarray:
    """Return all ``<P_i P_j>`` correlations for one Pauli axis.

    With ``connected=True``, subtract ``<P_i><P_j>``.
    """

    axis = _validated_axis(axis)
    one_point = np.asarray(
        [_pauli_expectation(state, n_sites, {site: axis}, basis=basis) for site in range(n_sites)],
        dtype=complex,
    )
    correlations = np.empty((n_sites, n_sites), dtype=complex)
    for i in range(n_sites):
        correlations[i, i] = _state_norm_squared(state)
        for j in range(i + 1, n_sites):
            value = _pauli_expectation(
                state,
                n_sites,
                {i: axis, j: axis},
                basis=basis,
            )
            correlations[i, j] = value
            correlations[j, i] = np.conjugate(value)
    if connected:
        correlations -= np.outer(one_point, one_point)
    return np.real_if_close(correlations)


def connected_spin_correlation_matrix(
    state: np.ndarray,
    n_sites: int,
    *,
    axis: str = "Z",
    basis: FixedMagnetizationBasis | None = None,
) -> np.ndarray:
    """Return connected same-axis spin correlations."""

    return spin_correlation_matrix(
        state,
        n_sites,
        axis=axis,
        basis=basis,
        connected=True,
    )


def correlation_zz(
    state: np.ndarray,
    n_sites: int,
    i: int,
    j: int,
    *,
    basis: FixedMagnetizationBasis | None = None,
) -> float:
    """Return the two-point ``<Z_i Z_j>`` correlation."""

    _validate_site(i, n_sites)
    _validate_site(j, n_sites)
    if i == j:
        return _state_norm_squared(state)
    value = _pauli_expectation(state, n_sites, {i: "Z", j: "Z"}, basis=basis)
    return float(np.real_if_close(value).real)


def static_spin_structure_factor(
    state: np.ndarray,
    n_sites: int,
    momenta: float | Sequence[float] | np.ndarray,
    *,
    axis: str = "Z",
    basis: FixedMagnetizationBasis | None = None,
    connected: bool = True,
    positions: Sequence[float] | np.ndarray | None = None,
) -> np.ndarray:
    """Return ``S(q) = N^-1 sum_ij exp[iq(r_i-r_j)] C_ij``."""

    q_values = np.atleast_1d(np.asarray(momenta, dtype=float))
    coordinates = (
        np.arange(n_sites, dtype=float)
        if positions is None
        else np.asarray(positions, dtype=float).reshape(-1)
    )
    if coordinates.size != n_sites:
        raise ValueError("positions must contain one coordinate per site.")
    correlations = spin_correlation_matrix(
        state,
        n_sites,
        axis=axis,
        basis=basis,
        connected=connected,
    )
    displacements = coordinates[:, None] - coordinates[None, :]
    factors = np.asarray(
        [
            np.sum(np.exp(1j * momentum * displacements) * correlations) / n_sites
            for momentum in q_values
        ]
    )
    return np.real_if_close(factors)


def reduced_density_matrix(
    state: np.ndarray,
    n_sites: int,
    subsystem: Sequence[int],
    *,
    basis: FixedMagnetizationBasis | None = None,
) -> np.ndarray:
    """Return the normalized reduced density matrix for selected spin sites.

    The calculation groups amplitudes by subsystem and environment bit
    patterns, so reduced-basis states do not need to be expanded to length
    ``2**n_sites``.
    """

    sites = _validated_subsystem(subsystem, n_sites)
    environment = tuple(site for site in range(n_sites) if site not in sites)
    vector, states = _state_and_basis_states(state, n_sites, basis)
    norm = float(np.vdot(vector, vector).real)
    if norm <= 0.0:
        raise ValueError("state must have nonzero norm.")
    amplitudes = np.zeros((2 ** len(sites), 2 ** len(environment)), dtype=complex)
    for amplitude, full_state in zip(vector, states, strict=True):
        subsystem_index = _extract_bits(full_state, sites, n_sites)
        environment_index = _extract_bits(full_state, environment, n_sites)
        amplitudes[subsystem_index, environment_index] += amplitude
    density = amplitudes @ amplitudes.conj().T / norm
    return (density + density.conj().T) / 2.0


def bipartite_entanglement_entropy(
    state: np.ndarray,
    n_sites: int,
    subsystem: Sequence[int],
    *,
    basis: FixedMagnetizationBasis | None = None,
    base: float = 2.0,
    tolerance: float = 1e-12,
) -> float:
    """Return the von Neumann entropy of a spin subsystem."""

    if base <= 0.0 or np.isclose(base, 1.0):
        raise ValueError("base must be positive and different from one.")
    if tolerance < 0.0:
        raise ValueError("tolerance must be nonnegative.")
    density = reduced_density_matrix(state, n_sites, subsystem, basis=basis)
    eigenvalues = np.linalg.eigvalsh(density)
    probabilities = eigenvalues[eigenvalues > tolerance]
    return float(-np.sum(probabilities * np.log(probabilities)) / np.log(base))


def inverse_participation_ratio(vector: np.ndarray) -> float:
    """Return sum_i |v_i|^4 / (sum_i |v_i|^2)^2."""

    values = np.asarray(vector, dtype=complex).reshape(-1)
    probabilities = np.abs(values) ** 2
    norm = probabilities.sum()
    if norm == 0:
        raise ValueError("vector must have nonzero norm.")
    return float(np.sum(probabilities**2) / norm**2)


def single_particle_occupations(state: np.ndarray) -> np.ndarray:
    """Return normalized site/orbital occupations for a single-particle state."""

    vector = np.asarray(state, dtype=complex).reshape(-1)
    probabilities = np.abs(vector) ** 2
    norm = probabilities.sum()
    if norm == 0:
        raise ValueError("state must have nonzero norm.")
    return probabilities / norm


def bond_current(
    state: np.ndarray,
    source: int,
    target: int,
    hopping: complex,
) -> float:
    """Return the oriented single-particle bond current ``2 Im(t ψ_i* ψ_j)``."""

    vector = np.asarray(state, dtype=complex).reshape(-1)
    if not 0 <= source < vector.size or not 0 <= target < vector.size:
        raise ValueError("Bond-current indices are out of range.")
    return float(2.0 * np.imag(complex(hopping) * np.conjugate(vector[source]) * vector[target]))


def mixed_spin_correlation_matrix(
    state: np.ndarray,
    n_sites: int,
    source_axis: str,
    target_axis: str,
    *,
    basis: FixedMagnetizationBasis | None = None,
) -> np.ndarray:
    """Return mixed-axis ``<P_i Q_j>`` spin correlations."""

    source_axis = _validated_axis(source_axis)
    target_axis = _validated_axis(target_axis)
    result = np.empty((n_sites, n_sites), dtype=complex)
    for source in range(n_sites):
        for target in range(n_sites):
            operators = {source: source_axis}
            if target == source:
                if source_axis != target_axis:
                    # Pauli products on one site: XY=iZ, YZ=iX, ZX=iY.
                    products = {
                        ("X", "Y"): (1j, "Z"),
                        ("Y", "Z"): (1j, "X"),
                        ("Z", "X"): (1j, "Y"),
                        ("Y", "X"): (-1j, "Z"),
                        ("Z", "Y"): (-1j, "X"),
                        ("X", "Z"): (-1j, "Y"),
                    }
                    phase, axis = products[(source_axis, target_axis)]
                    result[source, target] = phase * _pauli_expectation(
                        state, n_sites, {source: axis}, basis=basis
                    )
                else:
                    result[source, target] = _state_norm_squared(state)
            else:
                operators[target] = target_axis
                result[source, target] = _pauli_expectation(state, n_sites, operators, basis=basis)
    return np.real_if_close(result)


def total_spin_squared(
    state: np.ndarray,
    n_sites: int,
    *,
    basis: FixedMagnetizationBasis | None = None,
) -> float:
    """Return ``<S^2>`` for spin-1/2 sites with ``S=Σσ/2``."""

    value = 0.75 * n_sites * _state_norm_squared(state)
    for axis in ("X", "Y", "Z"):
        correlations = spin_correlation_matrix(state, n_sites, axis=axis, basis=basis)
        value += 0.5 * sum(
            correlations[left, right].real
            for left in range(n_sites)
            for right in range(left + 1, n_sites)
        )
    return float(value)


def local_density_of_states(
    hamiltonian: np.ndarray,
    energies: Sequence[float],
    *,
    broadening: float = 0.05,
) -> np.ndarray:
    """Return Lorentzian-broadened local density of states by site and energy."""

    if broadening <= 0:
        raise ValueError("broadening must be positive.")
    values, vectors = np.linalg.eigh(np.asarray(hamiltonian))
    grid = np.asarray(energies, dtype=float)
    lorentzian = broadening / np.pi / ((grid[:, None] - values[None, :]) ** 2 + broadening**2)
    return lorentzian @ (np.abs(vectors) ** 2).T


def _pauli_expectation(
    state: np.ndarray,
    n_sites: int,
    operators: dict[int, str],
    *,
    basis: FixedMagnetizationBasis | None,
) -> complex:
    vector, states = _state_and_basis_states(state, n_sites, basis)
    state_to_index = {full_state: index for index, full_state in enumerate(states)}
    value = 0.0j
    for column, full_state in enumerate(states):
        target = full_state
        phase = 1.0 + 0.0j
        for site, axis in operators.items():
            _validate_site(site, n_sites)
            axis = _validated_axis(axis)
            mask = 1 << (n_sites - 1 - site)
            bit = bool(full_state & mask)
            if axis == "Z":
                phase *= -1.0 if bit else 1.0
            elif axis == "X":
                target ^= mask
            else:
                phase *= -1j if bit else 1j
                target ^= mask
        row = state_to_index.get(target)
        if row is not None:
            value += np.conjugate(vector[row]) * phase * vector[column]
    return complex(value)


def _state_and_basis_states(
    state: np.ndarray,
    n_sites: int,
    basis: FixedMagnetizationBasis | None,
) -> tuple[np.ndarray, tuple[int, ...]]:
    if not isinstance(n_sites, int) or n_sites < 1:
        raise ValueError("n_sites must be a positive integer.")
    vector = np.asarray(state, dtype=complex).reshape(-1)
    if basis is None:
        expected = 2**n_sites
        states = tuple(range(expected))
    else:
        if basis.n_sites != n_sites:
            raise ValueError("basis.n_sites must match n_sites.")
        expected = basis.dimension
        states = basis.states
    if vector.size != expected:
        raise ValueError("state length must match the selected spin basis.")
    return vector, states


def _state_norm_squared(state: np.ndarray) -> float:
    vector = np.asarray(state, dtype=complex).reshape(-1)
    return float(np.vdot(vector, vector).real)


def _validated_axis(axis: str) -> str:
    normalized = axis.upper()
    if normalized not in {"X", "Y", "Z"}:
        raise ValueError("axis must be 'X', 'Y', or 'Z'.")
    return normalized


def _validate_site(site: int, n_sites: int) -> None:
    if not isinstance(site, int) or not 0 <= site < n_sites:
        raise ValueError(f"site must satisfy 0 <= site < {n_sites}.")


def _validated_subsystem(subsystem: Sequence[int], n_sites: int) -> tuple[int, ...]:
    sites = tuple(subsystem)
    if not sites:
        raise ValueError("subsystem must contain at least one site.")
    if len(set(sites)) != len(sites):
        raise ValueError("subsystem sites must be unique.")
    for site in sites:
        _validate_site(site, n_sites)
    return sites


def _extract_bits(state: int, sites: Sequence[int], n_sites: int) -> int:
    value = 0
    for site in sites:
        value = (value << 1) | ((state >> (n_sites - 1 - site)) & 1)
    return value
