from __future__ import annotations

import numpy as np

from quantum_lattice_models.models import (
    heisenberg_chain,
    ssh_edge_state_localizations,
    ssh_model,
    tight_binding_chain,
    transverse_field_ising,
)
from quantum_lattice_models.spectra import eigensystem, eigenvalues


def assert_hermitian(matrix: np.ndarray) -> None:
    assert np.allclose(matrix, matrix.conj().T)


def test_ising_two_site_shape() -> None:
    H = transverse_field_ising(n_sites=2, j=1.0, h=0.5, periodic=False)
    assert H.shape == (4, 4)


def test_spin_hamiltonians_are_hermitian() -> None:
    assert_hermitian(transverse_field_ising(n_sites=4, j=1.0, h=0.5))
    assert_hermitian(heisenberg_chain(n_sites=4, jx=1.0, jy=0.7, jz=1.2, field=0.1))


def test_tight_binding_shapes_and_hermiticity() -> None:
    chain = tight_binding_chain(n_sites=8, hopping=1.0, onsite=0.25)
    ssh = ssh_model(n_cells=6, t1=0.5, t2=1.0)
    assert chain.shape == (8, 8)
    assert ssh.shape == (12, 12)
    assert_hermitian(chain)
    assert_hermitian(ssh)


def test_eigenvalues_are_real_for_hermitian_models() -> None:
    models = [
        transverse_field_ising(n_sites=3, j=1.0, h=0.5),
        heisenberg_chain(n_sites=3),
        ssh_model(n_cells=4),
        tight_binding_chain(n_sites=5),
    ]
    for H in models:
        values = eigenvalues(H)
        assert np.allclose(values.imag, 0.0)


def test_open_ssh_topological_regime_has_edge_localized_low_energy_states() -> None:
    H = ssh_model(n_cells=12, t1=0.25, t2=1.0, periodic=False)
    values, vectors = eigensystem(H)
    near_zero = np.argsort(np.abs(values))[:2]
    localizations = ssh_edge_state_localizations(vectors[:, near_zero], n_cells=12, edge_cells=2)
    assert np.max(np.abs(values[near_zero])) < 1e-2
    assert np.all(localizations > 0.65)
