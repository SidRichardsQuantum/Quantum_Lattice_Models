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
    components: tuple[tuple[tuple[int, complex], ...], ...] = ()
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
        if self.components:
            if len(self.components) != len(self.states):
                raise ValueError("Reduced-basis components must align with states.")
            support_sets: list[set[int]] = []
            for row in self.components:
                if not row:
                    raise ValueError("Each reduced-basis row requires at least one component.")
                support = {state for state, _coefficient in row}
                if len(support) != len(row):
                    raise ValueError("Reduced-basis row components must have unique states.")
                if any(
                    not isinstance(state, int) or not 0 <= state < self.full_dimension
                    for state in support
                ):
                    raise ValueError("Reduced-basis components must index the full basis.")
                norm = sum(abs(complex(coefficient)) ** 2 for _, coefficient in row)
                if not np.isclose(norm, 1.0):
                    raise ValueError("Reduced-basis component rows must be normalized.")
                support_sets.append(support)
            if any(
                left & right
                for index, left in enumerate(support_sets)
                for right in support_sets[index + 1 :]
            ):
                raise ValueError("Reduced-basis component rows must have disjoint support.")

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
        if self.components:
            for amplitude, row in zip(vector, self.components, strict=True):
                for full_state, coefficient in row:
                    full[full_state] += coefficient * amplitude
        else:
            full[np.asarray(self.states, dtype=int)] = vector
        return full

    def project(self, state: np.ndarray) -> np.ndarray:
        vector = np.asarray(state, dtype=complex).reshape(-1)
        if vector.size != self.full_dimension:
            raise ValueError("Full state length must match mapping full_dimension.")
        if not self.components:
            return np.asarray(vector[np.asarray(self.states, dtype=int)], dtype=complex)
        return np.asarray(
            [
                sum(
                    np.conjugate(coefficient) * vector[full_state]
                    for full_state, coefficient in row
                )
                for row in self.components
            ],
            dtype=complex,
        )

    def transformation(self) -> sp.csr_matrix:
        """Return the full-by-reduced isometry represented by this mapping."""

        if not self.components:
            rows = np.asarray(self.states, dtype=int)
            columns = np.arange(self.dimension, dtype=int)
            values = np.ones(self.dimension, dtype=complex)
        else:
            entries = [
                (full_state, column, coefficient)
                for column, row in enumerate(self.components)
                for full_state, coefficient in row
            ]
            rows = np.asarray([entry[0] for entry in entries], dtype=int)
            columns = np.asarray([entry[1] for entry in entries], dtype=int)
            values = np.asarray([entry[2] for entry in entries], dtype=complex)
        return sp.coo_matrix(
            (values, (rows, columns)),
            shape=(self.full_dimension, self.dimension),
        ).tocsr()

    def to_dict(self) -> dict[str, object]:
        return {
            "kind": self.kind,
            "full_dimension": self.full_dimension,
            "dimension": self.dimension,
            "states": list(self.states),
            "quantum_numbers": dict(self.quantum_numbers),
            "labels": list(self.labels),
            "components": [
                [
                    {
                        "state": full_state,
                        "coefficient": {
                            "__complex__": [
                                complex(coefficient).real,
                                complex(coefficient).imag,
                            ]
                        },
                    }
                    for full_state, coefficient in row
                ]
                for row in self.components
            ],
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> ReducedBasisMapping:
        states = _object_sequence(data["states"], "Reduced-basis states")
        quantum_numbers = _quantum_number_mapping(data.get("quantum_numbers", {}))
        labels = _object_sequence(data.get("labels", ()), "Reduced-basis labels")
        raw_components = _object_sequence(data.get("components", ()), "Reduced-basis components")
        metadata = _object_mapping(data.get("metadata", {}), "Reduced-basis metadata")
        return cls(
            kind=str(data["kind"]),
            full_dimension=_integer_value(data["full_dimension"], "Reduced-basis full_dimension"),
            states=tuple(_integer_value(state, "Reduced-basis state") for state in states),
            quantum_numbers={str(name): value for name, value in quantum_numbers.items()},
            labels=tuple(str(label) for label in labels),
            components=tuple(_component_row(row) for row in raw_components),
            metadata=metadata,
        )


def reduced_operator(
    operator: np.ndarray | sp.spmatrix,
    mapping: ReducedBasisMapping,
) -> np.ndarray | sp.csr_matrix:
    """Project a full-basis operator into an explicit reduced basis."""

    if operator.shape != (mapping.full_dimension, mapping.full_dimension):
        raise ValueError("Operator shape must match mapping full_dimension.")
    if mapping.components:
        transformation = mapping.transformation()
        if sp.issparse(operator):
            return sp.csr_matrix(transformation.getH() @ operator @ transformation)
        dense_transformation = transformation.toarray()
        return np.asarray(
            dense_transformation.conj().T @ np.asarray(operator) @ dense_transformation
        )
    states = np.asarray(mapping.states, dtype=int)
    if sp.issparse(operator):
        sparse_operator = sp.csr_matrix(operator)
        return sparse_operator[states][:, states]
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


def _component_coefficient(value: object) -> complex:
    if isinstance(value, (int, float, complex)) and not isinstance(value, bool):
        return complex(value)
    if isinstance(value, dict) and set(value) == {"__complex__"}:
        parts = value["__complex__"]
        if isinstance(parts, list) and len(parts) == 2:
            return complex(float(parts[0]), float(parts[1]))
    raise ValueError("Reduced-basis component coefficient must be numeric.")


def _integer_value(value: object, name: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError(f"{name} must be an integer.")
    return value


def _object_sequence(value: object, name: str) -> list[object] | tuple[object, ...]:
    if not isinstance(value, (list, tuple)):
        raise ValueError(f"{name} must be a list.")
    return value


def _object_mapping(value: object, name: str) -> dict[str, object]:
    if not isinstance(value, dict) or not all(isinstance(key, str) for key in value):
        raise ValueError(f"{name} must be an object with string keys.")
    return dict(value)


def _quantum_number_mapping(value: object) -> dict[str, int | float]:
    mapping = _object_mapping(value, "Reduced-basis quantum_numbers")
    if not all(
        isinstance(item, (int, float)) and not isinstance(item, bool) for item in mapping.values()
    ):
        raise ValueError("Reduced-basis quantum numbers must be numeric.")
    return {name: item for name, item in mapping.items() if isinstance(item, (int, float))}


def _component_row(value: object) -> tuple[tuple[int, complex], ...]:
    row = _object_sequence(value, "Reduced-basis component row")
    components = []
    for value in row:
        component = _object_mapping(value, "Reduced-basis component")
        if "state" not in component or "coefficient" not in component:
            raise ValueError("Reduced-basis components require state and coefficient.")
        components.append(
            (
                _integer_value(component["state"], "Reduced-basis component state"),
                _component_coefficient(component["coefficient"]),
            )
        )
    return tuple(components)
