"""Dense spin-chain and spin-ladder Hamiltonian builders."""

from __future__ import annotations

import numpy as np

from quantum_lattice_models._model_utils import (
    nearest_neighbor_bonds,
    next_nearest_neighbor_bonds,
    pauli_labels,
    validate_positive_int,
)
from quantum_lattice_models.operators import (
    PAULI_X,
    PAULI_Y,
    PAULI_Z,
    local_operator,
    two_site_operator,
)
from quantum_lattice_models.types import DenseHamiltonian, PauliTerm


def transverse_field_ising(
    n_sites: int, j: float = 1.0, h: float = 0.5, periodic: bool = False
) -> DenseHamiltonian:
    """Return ``H = -J sum ZZ - h sum X``."""

    validate_positive_int(n_sites, "n_sites")
    matrix = np.zeros((2**n_sites, 2**n_sites), dtype=complex)
    terms: list[PauliTerm] = []
    for i, k in nearest_neighbor_bonds(n_sites, periodic):
        matrix += -j * two_site_operator(PAULI_Z, PAULI_Z, i, k, n_sites)
        terms.append(PauliTerm(-j, pauli_labels(n_sites, {i: "Z", k: "Z"})))
    for i in range(n_sites):
        matrix += -h * local_operator(PAULI_X, i, n_sites)
        terms.append(PauliTerm(-h, pauli_labels(n_sites, {i: "X"})))
    return DenseHamiltonian(
        matrix, model_name="transverse_field_ising", n_sites=n_sites, terms=tuple(terms)
    )


def longitudinal_field_ising(
    n_sites: int,
    j: float = 1.0,
    h_x: float = 0.5,
    h_z: float = 0.0,
    periodic: bool = False,
) -> DenseHamiltonian:
    """Return ``H = -J sum ZZ - h_x sum X - h_z sum Z``."""

    validate_positive_int(n_sites, "n_sites")
    matrix = np.zeros((2**n_sites, 2**n_sites), dtype=complex)
    terms: list[PauliTerm] = []
    for i, k in nearest_neighbor_bonds(n_sites, periodic):
        matrix += -j * two_site_operator(PAULI_Z, PAULI_Z, i, k, n_sites)
        terms.append(PauliTerm(-j, pauli_labels(n_sites, {i: "Z", k: "Z"})))
    for i in range(n_sites):
        matrix += -h_x * local_operator(PAULI_X, i, n_sites)
        matrix += -h_z * local_operator(PAULI_Z, i, n_sites)
        terms.append(PauliTerm(-h_x, pauli_labels(n_sites, {i: "X"})))
        terms.append(PauliTerm(-h_z, pauli_labels(n_sites, {i: "Z"})))
    return DenseHamiltonian(
        matrix, model_name="longitudinal_field_ising", n_sites=n_sites, terms=tuple(terms)
    )


def next_nearest_neighbor_ising(
    n_sites: int,
    j1: float = 1.0,
    j2: float = 0.5,
    h: float = 0.5,
    periodic: bool = False,
) -> DenseHamiltonian:
    """Return ``H = -J1 sum ZiZi+1 - J2 sum ZiZi+2 - h sum X``."""

    validate_positive_int(n_sites, "n_sites")
    matrix = np.zeros((2**n_sites, 2**n_sites), dtype=complex)
    terms: list[PauliTerm] = []
    for coefficient, bonds in (
        (-j1, nearest_neighbor_bonds(n_sites, periodic)),
        (-j2, next_nearest_neighbor_bonds(n_sites, periodic)),
    ):
        for i, k in bonds:
            matrix += coefficient * two_site_operator(PAULI_Z, PAULI_Z, i, k, n_sites)
            terms.append(PauliTerm(coefficient, pauli_labels(n_sites, {i: "Z", k: "Z"})))
    for i in range(n_sites):
        matrix += -h * local_operator(PAULI_X, i, n_sites)
        terms.append(PauliTerm(-h, pauli_labels(n_sites, {i: "X"})))
    return DenseHamiltonian(
        matrix, model_name="next_nearest_neighbor_ising", n_sites=n_sites, terms=tuple(terms)
    )


def heisenberg_chain(
    n_sites: int,
    jx: float = 1.0,
    jy: float = 1.0,
    jz: float = 1.0,
    field: float = 0.0,
    periodic: bool = False,
) -> DenseHamiltonian:
    """Return the dense anisotropic Heisenberg spin-chain Hamiltonian."""

    validate_positive_int(n_sites, "n_sites")
    matrix = np.zeros((2**n_sites, 2**n_sites), dtype=complex)
    terms: list[PauliTerm] = []
    for i, k in nearest_neighbor_bonds(n_sites, periodic):
        for coefficient, operator, label in (
            (jx, PAULI_X, "X"),
            (jy, PAULI_Y, "Y"),
            (jz, PAULI_Z, "Z"),
        ):
            matrix += coefficient * two_site_operator(operator, operator, i, k, n_sites)
            terms.append(PauliTerm(coefficient, pauli_labels(n_sites, {i: label, k: label})))
    for i in range(n_sites):
        matrix += field * local_operator(PAULI_Z, i, n_sites)
        terms.append(PauliTerm(field, pauli_labels(n_sites, {i: "Z"})))
    return DenseHamiltonian(
        matrix, model_name="heisenberg_chain", n_sites=n_sites, terms=tuple(terms)
    )


def xy_chain(
    n_sites: int,
    coupling: float = 1.0,
    anisotropy: float = 0.0,
    field: float = 0.0,
    periodic: bool = False,
) -> DenseHamiltonian:
    """Return the dense anisotropic XY spin-chain Hamiltonian."""

    validate_positive_int(n_sites, "n_sites")
    matrix = np.zeros((2**n_sites, 2**n_sites), dtype=complex)
    terms: list[PauliTerm] = []
    xx_coefficient = -coupling * (1.0 + anisotropy) / 2.0
    yy_coefficient = -coupling * (1.0 - anisotropy) / 2.0
    for i, k in nearest_neighbor_bonds(n_sites, periodic):
        matrix += xx_coefficient * two_site_operator(PAULI_X, PAULI_X, i, k, n_sites)
        matrix += yy_coefficient * two_site_operator(PAULI_Y, PAULI_Y, i, k, n_sites)
        terms.append(PauliTerm(xx_coefficient, pauli_labels(n_sites, {i: "X", k: "X"})))
        terms.append(PauliTerm(yy_coefficient, pauli_labels(n_sites, {i: "Y", k: "Y"})))
    for i in range(n_sites):
        matrix += -field * local_operator(PAULI_Z, i, n_sites)
        terms.append(PauliTerm(-field, pauli_labels(n_sites, {i: "Z"})))
    return DenseHamiltonian(matrix, model_name="xy_chain", n_sites=n_sites, terms=tuple(terms))


def xxz_chain(
    n_sites: int,
    coupling: float = 1.0,
    anisotropy: float = 1.0,
    field: float = 0.0,
    periodic: bool = False,
) -> DenseHamiltonian:
    """Return the dense XXZ spin-chain Hamiltonian."""

    matrix = heisenberg_chain(
        n_sites=n_sites,
        jx=coupling,
        jy=coupling,
        jz=coupling * anisotropy,
        field=field,
        periodic=periodic,
    )
    matrix.model_name = "xxz_chain"
    return matrix


def j1_j2_heisenberg_chain(
    n_sites: int,
    j1: float = 1.0,
    j2: float = 0.5,
    field: float = 0.0,
    periodic: bool = False,
) -> DenseHamiltonian:
    """Return a dense frustrated J1-J2 Heisenberg spin-chain Hamiltonian."""

    validate_positive_int(n_sites, "n_sites")
    matrix = np.zeros((2**n_sites, 2**n_sites), dtype=complex)
    terms: list[PauliTerm] = []
    for coefficient, bonds in (
        (j1, nearest_neighbor_bonds(n_sites, periodic)),
        (j2, next_nearest_neighbor_bonds(n_sites, periodic)),
    ):
        for i, k in bonds:
            for operator, label in ((PAULI_X, "X"), (PAULI_Y, "Y"), (PAULI_Z, "Z")):
                matrix += coefficient * two_site_operator(operator, operator, i, k, n_sites)
                terms.append(PauliTerm(coefficient, pauli_labels(n_sites, {i: label, k: label})))
    for i in range(n_sites):
        matrix += field * local_operator(PAULI_Z, i, n_sites)
        terms.append(PauliTerm(field, pauli_labels(n_sites, {i: "Z"})))
    return DenseHamiltonian(
        matrix, model_name="j1_j2_heisenberg_chain", n_sites=n_sites, terms=tuple(terms)
    )


def heisenberg_ladder(
    n_rungs: int,
    leg_coupling: float = 1.0,
    rung_coupling: float = 1.0,
    field: float = 0.0,
    periodic: bool = False,
) -> DenseHamiltonian:
    """Return a dense two-leg Heisenberg ladder with ``2 * n_rungs`` spins."""

    validate_positive_int(n_rungs, "n_rungs")
    n_sites = 2 * n_rungs
    matrix = np.zeros((2**n_sites, 2**n_sites), dtype=complex)
    terms: list[PauliTerm] = []
    leg_bonds = nearest_neighbor_bonds(n_rungs, periodic)
    ladder_bonds = [(i, k, leg_coupling) for i, k in leg_bonds]
    ladder_bonds.extend((n_rungs + i, n_rungs + k, leg_coupling) for i, k in leg_bonds)
    ladder_bonds.extend((r, n_rungs + r, rung_coupling) for r in range(n_rungs))
    for i, k, coefficient in ladder_bonds:
        for operator, label in ((PAULI_X, "X"), (PAULI_Y, "Y"), (PAULI_Z, "Z")):
            matrix += coefficient * two_site_operator(operator, operator, i, k, n_sites)
            terms.append(PauliTerm(coefficient, pauli_labels(n_sites, {i: label, k: label})))
    for i in range(n_sites):
        matrix += field * local_operator(PAULI_Z, i, n_sites)
        terms.append(PauliTerm(field, pauli_labels(n_sites, {i: "Z"})))
    return DenseHamiltonian(
        matrix, model_name="heisenberg_ladder", n_sites=n_sites, terms=tuple(terms)
    )
