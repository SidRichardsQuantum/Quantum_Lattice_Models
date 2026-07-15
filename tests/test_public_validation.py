from __future__ import annotations

import numpy as np
import pytest

from quantum_lattice_models.lattice import Lattice, custom_tight_binding
from quantum_lattice_models.observables import expectation, inverse_participation_ratio
from quantum_lattice_models.operators import (
    PAULI_X,
    kron_all,
    local_operator,
    pauli_string_matrix,
    two_site_operator,
)
from quantum_lattice_models.spectra import (
    density_of_states,
    lowest_eigenvalues,
    spectral_gap,
)


def test_model_namespaces_preserve_flat_compatibility_exports() -> None:
    from quantum_lattice_models.models import heisenberg_chain, ssh_model
    from quantum_lattice_models.models.benchmarks import graphene_lattice
    from quantum_lattice_models.models.particles import ssh_model as namespaced_ssh
    from quantum_lattice_models.models.periodic import ssh_unit_cell
    from quantum_lattice_models.models.spin import heisenberg_chain as namespaced_heisenberg

    assert namespaced_ssh is ssh_model
    assert namespaced_heisenberg is heisenberg_chain
    assert graphene_lattice(1, 1).shape == (2, 2)
    assert ssh_unit_cell().n_orbitals == 2


def test_operator_construction_and_validation() -> None:
    assert np.allclose(pauli_string_matrix(("X", "I")), np.kron(PAULI_X, np.eye(2)))
    with pytest.raises(ValueError, match="At least one operator"):
        kron_all(())
    with pytest.raises(ValueError, match="Unknown Pauli"):
        pauli_string_matrix(("A",))
    with pytest.raises(IndexError, match="site must satisfy"):
        local_operator(PAULI_X, site=2, n_sites=2)
    with pytest.raises(ValueError, match="distinct sites"):
        two_site_operator(PAULI_X, PAULI_X, 0, 0, 2)


def test_observable_and_spectral_validation() -> None:
    with pytest.raises(ValueError, match="operator shape"):
        expectation(np.ones(2), np.eye(3))
    with pytest.raises(ValueError, match="nonzero norm"):
        inverse_participation_ratio(np.zeros(2))
    with pytest.raises(ValueError, match="k must be positive"):
        lowest_eigenvalues(np.eye(2), k=0)
    with pytest.raises(ValueError, match="bins must be positive"):
        density_of_states(np.eye(2), bins=0)
    assert spectral_gap(np.diag([0.0, 0.0, 2.0])) == 2.0
    assert spectral_gap(np.array([[1.0]])) == 0.0


def test_lattice_input_validation_and_nonhermitian_mode() -> None:
    with pytest.raises(ValueError, match="n_sites is required"):
        Lattice()
    with pytest.raises(ValueError, match="less than n_sites"):
        Lattice(n_sites=2, bonds=[(0, 2)])
    with pytest.raises(ValueError, match="2D or 3D"):
        Lattice(positions=[(0.0,), (1.0,)])
    with pytest.raises(ValueError, match="either lattice or"):
        custom_tight_binding(
            n_sites=2,
            lattice=Lattice(n_sites=2),
        )

    directed = custom_tight_binding(
        n_sites=2,
        bonds=[(0, 1, 2j)],
        hermitian=False,
    )
    assert directed[0, 1] == 2j
    assert directed[1, 0] == 0
