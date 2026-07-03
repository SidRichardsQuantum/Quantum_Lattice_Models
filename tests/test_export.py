from __future__ import annotations

import sys
import types

import pytest

from quantum_lattice_models.export import to_pennylane_terms
from quantum_lattice_models.models import transverse_field_ising


class _FakeObservable:
    def __init__(self, name: str, wire: int) -> None:
        self.name = name
        self.wire = wire
        self.factors = ((name, wire),)

    def __matmul__(self, other: _FakeObservable) -> _FakeObservable:
        product = _FakeObservable("Product", self.wire)
        product.factors = self.factors + other.factors
        return product


class _FakeHamiltonian:
    def __init__(self, coeffs, observables) -> None:
        self.coeffs = coeffs
        self.observables = observables


def test_pennylane_export_builds_terms_with_optional_backend(monkeypatch) -> None:
    fake_qml = types.ModuleType("pennylane")
    fake_qml.Hamiltonian = _FakeHamiltonian
    fake_qml.Identity = lambda wire: _FakeObservable("I", wire)
    fake_qml.PauliX = lambda wire: _FakeObservable("X", wire)
    fake_qml.PauliY = lambda wire: _FakeObservable("Y", wire)
    fake_qml.PauliZ = lambda wire: _FakeObservable("Z", wire)
    monkeypatch.setitem(sys.modules, "pennylane", fake_qml)

    H = transverse_field_ising(n_sites=2, j=1.0, h=0.5)
    qml_H = to_pennylane_terms(H)

    assert isinstance(qml_H, _FakeHamiltonian)
    assert len(qml_H.coeffs) == len(H.terms)
    assert len(qml_H.observables) == len(H.terms)


def test_pennylane_export_reports_missing_optional_backend(monkeypatch) -> None:
    monkeypatch.setitem(sys.modules, "pennylane", None)

    with pytest.raises(ImportError, match="PennyLane is not installed"):
        to_pennylane_terms(transverse_field_ising(n_sites=2))
