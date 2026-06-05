"""Optional exports to external quantum software."""

from __future__ import annotations

import numpy as np

from quantum_lattice_models.types import DenseHamiltonian


def to_pennylane_terms(lattice_hamiltonian: DenseHamiltonian):
    """Convert a spin-chain Hamiltonian with Pauli metadata to a PennyLane Hamiltonian.

    PennyLane is optional. Install with ``pip install quantum-lattice-models[pennylane]``
    or ``pip install pennylane`` before calling this function.
    """

    try:
        import pennylane as qml
    except ImportError as exc:
        raise ImportError(
            "PennyLane is not installed. Install it with "
            "`pip install quantum-lattice-models[pennylane]` or `pip install pennylane`."
        ) from exc

    if not isinstance(lattice_hamiltonian, DenseHamiltonian) or not lattice_hamiltonian.terms:
        raise TypeError(
            "to_pennylane_terms expects a spin-chain DenseHamiltonian returned by this package."
        )

    coeffs: list[float | complex] = []
    observables = []
    for term in lattice_hamiltonian.terms:
        active_ops = [(wire, label) for wire, label in enumerate(term.operators) if label != "I"]
        if not active_ops:
            coeffs.append(np.real_if_close(term.coefficient).item())
            observables.append(qml.Identity(0))
            continue

        observable = _qml_pauli(qml, active_ops[0][1], active_ops[0][0])
        for wire, label in active_ops[1:]:
            observable = observable @ _qml_pauli(qml, label, wire)
        coeffs.append(np.real_if_close(term.coefficient).item())
        observables.append(observable)

    return qml.Hamiltonian(coeffs, observables)


def _qml_pauli(qml, label: str, wire: int):
    if label == "X":
        return qml.PauliX(wire)
    if label == "Y":
        return qml.PauliY(wire)
    if label == "Z":
        return qml.PauliZ(wire)
    raise ValueError(f"Unsupported Pauli label {label!r}.")
