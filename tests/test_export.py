from __future__ import annotations

import pytest

from quantum_lattice_models.export import to_pennylane_terms
from quantum_lattice_models.models import transverse_field_ising


def test_pennylane_export_if_installed() -> None:
    qml = pytest.importorskip("pennylane")
    H = transverse_field_ising(n_sites=2, j=1.0, h=0.5)
    qml_H = to_pennylane_terms(H)
    assert isinstance(qml_H, qml.Hamiltonian)
    assert len(qml_H.coeffs) == len(H.terms)
