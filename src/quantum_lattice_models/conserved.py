"""Conserved-quantity records and commutator diagnostics."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models.analysis import AnalysisResult


@dataclass(frozen=True)
class ConservedQuantity:
    """Named operator with optional sector and convention metadata."""

    name: str
    operator: np.ndarray | sp.spmatrix
    sector_value: float | int | None = None
    convention: str = ""
    metadata: dict[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name:
            raise ValueError("Conserved quantity name must be nonempty.")
        if self.operator.ndim != 2 or self.operator.shape[0] != self.operator.shape[1]:
            raise ValueError("Conserved quantity operator must be square.")


def commutator_diagnostic(
    hamiltonian: np.ndarray | sp.spmatrix,
    quantity: ConservedQuantity | np.ndarray | sp.spmatrix,
    *,
    tolerance: float = 1e-10,
    name: str | None = None,
) -> AnalysisResult:
    """Measure ``[H, Q]`` and report whether ``Q`` is conserved."""

    operator = quantity.operator if isinstance(quantity, ConservedQuantity) else quantity
    quantity_name = quantity.name if isinstance(quantity, ConservedQuantity) else (name or "Q")
    if hamiltonian.shape != operator.shape:
        raise ValueError("Hamiltonian and conserved-quantity operator shapes must match.")
    commutator = hamiltonian @ operator - operator @ hamiltonian
    if sp.issparse(commutator):
        frobenius = float(np.sqrt(np.sum(np.abs(commutator.data) ** 2)))
        maximum = float(np.max(np.abs(commutator.data), initial=0.0))
    else:
        values = np.asarray(commutator)
        frobenius = float(np.linalg.norm(values))
        maximum = float(np.max(np.abs(values), initial=0.0))
    return AnalysisResult(
        analysis="commutator_diagnostic",
        values={
            "frobenius_norm": np.asarray([frobenius]),
            "maximum_absolute_entry": np.asarray([maximum]),
        },
        parameters={"quantity": quantity_name, "tolerance": tolerance},
        solver={"method": "explicit commutator", "exact": True},
        metadata={"conserved": maximum <= tolerance, "dimension": hamiltonian.shape[0]},
    )


def sector_compatibility(
    hamiltonian: np.ndarray | sp.spmatrix,
    quantity: ConservedQuantity | np.ndarray | sp.spmatrix,
    *,
    tolerance: float = 1e-10,
) -> bool:
    """Return whether a Hamiltonian preserves a proposed sector quantity."""

    return bool(
        commutator_diagnostic(hamiltonian, quantity, tolerance=tolerance).metadata["conserved"]
    )
