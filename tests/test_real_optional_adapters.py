from __future__ import annotations

import pytest

from quantum_lattice_models import (
    SpinInteraction,
    create_graph_spin_spec,
    create_model_spec,
    to_netket_graph,
    to_openfermion,
    to_qiskit_sparse_pauli,
    to_quspin_hamiltonian,
    to_qutip_qobj,
)
from quantum_lattice_models.export import to_pennylane_terms
from quantum_lattice_models.spin import transverse_field_ising


def test_real_openfermion_adapter() -> None:
    pytest.importorskip("openfermion")
    model = create_model_spec("fermi_hubbard_chain", parameters={"n_sites": 2})
    assert to_openfermion(model).terms


def test_real_qiskit_adapter() -> None:
    pytest.importorskip("qiskit")
    model = create_graph_spin_spec(
        2,
        interactions=(SpinInteraction(0, 1, "Z", "Z", 1.0),),
    )
    assert to_qiskit_sparse_pauli(model).num_qubits == 2


def test_real_qutip_adapter() -> None:
    pytest.importorskip("qutip")
    model = create_model_spec("tight_binding_chain", parameters={"n_sites": 3})
    assert to_qutip_qobj(model).shape == (3, 3)


def test_real_netket_adapter() -> None:
    pytest.importorskip("netket")
    model = create_model_spec("tight_binding_chain", parameters={"n_sites": 3})
    assert to_netket_graph(model.lattice).n_nodes == 3


def test_real_quspin_adapter() -> None:
    pytest.importorskip("quspin")
    model = create_graph_spin_spec(
        2,
        interactions=(SpinInteraction(0, 1, "Z", "Z", 1.0),),
    )
    assert to_quspin_hamiltonian(model).Ns == 4


def test_real_pennylane_adapter() -> None:
    pytest.importorskip("pennylane")
    assert to_pennylane_terms(transverse_field_ising(2)) is not None
