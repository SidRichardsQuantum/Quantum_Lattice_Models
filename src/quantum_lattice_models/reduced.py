"""Generic mappings between reduced rows and a full discrete basis."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp


@dataclass(frozen=True)
class ReducedBasisMapping:
    """Explicit row-to-full-state mapping for a symmetry-reduced basis."""

    kind: str
    full_dimension: int
    states: tuple[int, ...]
    quantum_numbers: dict[str, int | float] = field(default_factory=dict)
    labels: tuple[str, ...] = ()
    metadata: dict[str, object] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.kind:
            raise ValueError("Reduced-basis kind must be nonempty.")
        if not isinstance(self.full_dimension, int) or self.full_dimension < 1:
            raise ValueError("full_dimension must be a positive integer.")
        if not self.states:
            raise ValueError("Reduced-basis mappings require at least one state.")
        if len(set(self.states)) != len(self.states):
            raise ValueError("Reduced-basis full states must be unique.")
        if any(
            not isinstance(state, int) or not 0 <= state < self.full_dimension
            for state in self.states
        ):
            raise ValueError("Reduced-basis states must index the full basis.")
        if self.labels and len(self.labels) != len(self.states):
            raise ValueError("Reduced-basis labels must align with states.")

    @property
    def dimension(self) -> int:
        return len(self.states)

    @property
    def state_to_index(self) -> dict[int, int]:
        return {state: index for index, state in enumerate(self.states)}

    def embed(self, state: np.ndarray) -> np.ndarray:
        vector = np.asarray(state, dtype=complex).reshape(-1)
        if vector.size != self.dimension:
            raise ValueError("Reduced state length must match mapping dimension.")
        full = np.zeros(self.full_dimension, dtype=complex)
        full[np.asarray(self.states, dtype=int)] = vector
        return full

    def project(self, state: np.ndarray) -> np.ndarray:
        vector = np.asarray(state, dtype=complex).reshape(-1)
        if vector.size != self.full_dimension:
            raise ValueError("Full state length must match mapping full_dimension.")
        return vector[np.asarray(self.states, dtype=int)]

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "full_dimension": self.full_dimension,
            "dimension": self.dimension,
            "states": list(self.states),
            "quantum_numbers": dict(self.quantum_numbers),
            "labels": list(self.labels),
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> ReducedBasisMapping:
        return cls(
            kind=str(data["kind"]),
            full_dimension=int(data["full_dimension"]),
            states=tuple(int(state) for state in data["states"]),
            quantum_numbers={
                str(name): value for name, value in dict(data.get("quantum_numbers", {})).items()
            },
            labels=tuple(str(label) for label in data.get("labels", ())),
            metadata=dict(data.get("metadata", {})),
        )


def reduced_operator(
    operator: np.ndarray | sp.spmatrix,
    mapping: ReducedBasisMapping,
) -> np.ndarray | sp.csr_matrix:
    """Project a full-basis operator into an explicit reduced basis."""

    if operator.shape != (mapping.full_dimension, mapping.full_dimension):
        raise ValueError("Operator shape must match mapping full_dimension.")
    states = np.asarray(mapping.states, dtype=int)
    if sp.issparse(operator):
        return operator.tocsr()[states][:, states]
    return np.asarray(operator)[np.ix_(states, states)]


def reduced_expectation(
    state: np.ndarray,
    operator: np.ndarray | sp.spmatrix,
    mapping: ReducedBasisMapping,
) -> complex:
    """Evaluate a full- or reduced-basis operator on a reduced state."""

    vector = np.asarray(state, dtype=complex).reshape(-1)
    if vector.size != mapping.dimension:
        raise ValueError("State length must match mapping dimension.")
    if operator.shape == (mapping.full_dimension, mapping.full_dimension):
        operator = reduced_operator(operator, mapping)
    elif operator.shape != (mapping.dimension, mapping.dimension):
        raise ValueError("Operator must use either the full or reduced basis.")
    return complex(np.vdot(vector, operator @ vector))
