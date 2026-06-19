from __future__ import annotations

import numpy as np

from quantum_lattice_models.models import (
    bose_hubbard_chain,
    fermi_hubbard_chain,
    ssh_model,
    tight_binding_chain,
    transverse_field_ising,
)


def test_zero_field_open_ising_analytic_limit() -> None:
    values = np.linalg.eigvalsh(transverse_field_ising(n_sites=3, j=1.0, h=0.0))
    assert np.allclose(values, [-2.0, -2.0, 0.0, 0.0, 0.0, 0.0, 2.0, 2.0])


def test_uniform_open_tight_binding_chain_analytic_spectrum() -> None:
    n_sites = 5
    hopping = 0.7
    values = np.linalg.eigvalsh(tight_binding_chain(n_sites=n_sites, hopping=hopping))
    modes = np.arange(1, n_sites + 1)
    expected = np.sort(-2 * hopping * np.cos(modes * np.pi / (n_sites + 1)))
    assert np.allclose(values, expected)


def test_decoupled_ssh_dimer_analytic_spectrum() -> None:
    values = np.linalg.eigvalsh(ssh_model(n_cells=3, t1=0.8, t2=0.0))
    assert np.allclose(values, [-0.8, -0.8, -0.8, 0.8, 0.8, 0.8])


def test_noninteracting_single_site_hubbard_limits() -> None:
    bose = bose_hubbard_chain(
        n_sites=1,
        hopping=0.0,
        interaction=0.0,
        chemical_potential=0.0,
        max_occupancy=3,
    )
    fermi = fermi_hubbard_chain(
        n_sites=1,
        hopping=0.0,
        interaction=0.0,
        chemical_potential=0.0,
    )
    assert np.allclose(bose, 0.0)
    assert np.allclose(fermi, 0.0)
