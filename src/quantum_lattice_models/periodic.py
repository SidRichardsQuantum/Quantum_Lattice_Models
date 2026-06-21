"""Portable periodic unit cells, Bloch Hamiltonians, and band structures."""

from __future__ import annotations

import json
from collections.abc import Callable, Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

PERIODIC_SCHEMA_VERSION = "1.0"


@dataclass(frozen=True)
class CellBond:
    """A directed hopping from one orbital to an orbital in a displaced cell."""

    source: int
    target: int
    displacement: tuple[int, ...]
    value: complex

    def to_dict(self) -> dict[str, object]:
        return {
            "source": self.source,
            "target": self.target,
            "displacement": list(self.displacement),
            "value": {"__complex__": [complex(self.value).real, complex(self.value).imag]},
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> CellBond:
        value = data.get("value")
        if not isinstance(value, dict) or set(value) != {"__complex__"}:
            raise ValueError("Periodic bond values require complex-number encoding.")
        parts = value["__complex__"]
        if not isinstance(parts, list) or len(parts) != 2:
            raise ValueError("Periodic bond values require [real, imaginary].")
        return cls(
            source=int(data["source"]),
            target=int(data["target"]),
            displacement=tuple(int(item) for item in data["displacement"]),
            value=complex(float(parts[0]), float(parts[1])),
        )


@dataclass(frozen=True)
class PeriodicLatticeSpec:
    """Portable translationally invariant single-particle lattice model."""

    primitive_vectors: tuple[tuple[float, ...], ...]
    orbital_positions: tuple[tuple[float, ...], ...]
    bonds: tuple[CellBond, ...] = ()
    onsite: tuple[complex, ...] = ()
    orbital_labels: tuple[str, ...] = ()
    sublattice_labels: tuple[str, ...] = ()
    units: dict[str, str] = field(default_factory=dict)
    conventions: dict[str, str] = field(default_factory=dict)
    references: tuple[str, ...] = ()
    provenance: tuple[dict[str, object], ...] = ()
    metadata: dict[str, object] = field(default_factory=dict)
    schema_version: str = PERIODIC_SCHEMA_VERSION

    @property
    def dimension(self) -> int:
        return len(self.primitive_vectors)

    @property
    def n_orbitals(self) -> int:
        return len(self.orbital_positions)

    def validate(self) -> None:
        if self.schema_version != PERIODIC_SCHEMA_VERSION:
            raise ValueError(f"Unsupported periodic schema_version {self.schema_version!r}.")
        if not self.primitive_vectors:
            raise ValueError("primitive_vectors must contain at least one vector.")
        dimension = len(self.primitive_vectors)
        if any(len(vector) != dimension for vector in self.primitive_vectors):
            raise ValueError("primitive_vectors must form a square coordinate matrix.")
        vectors = np.asarray(self.primitive_vectors, dtype=float)
        if not np.all(np.isfinite(vectors)) or abs(np.linalg.det(vectors)) < 1e-12:
            raise ValueError("primitive_vectors must be finite and linearly independent.")
        if not self.orbital_positions:
            raise ValueError("orbital_positions must contain at least one orbital.")
        if any(len(position) != dimension for position in self.orbital_positions):
            raise ValueError("orbital positions must match the lattice dimension.")
        if not np.all(np.isfinite(np.asarray(self.orbital_positions, dtype=float))):
            raise ValueError("orbital positions must be finite.")
        if self.onsite and len(self.onsite) != self.n_orbitals:
            raise ValueError("onsite must contain one value per orbital.")
        for name, labels in (
            ("orbital_labels", self.orbital_labels),
            ("sublattice_labels", self.sublattice_labels),
        ):
            if labels and len(labels) != self.n_orbitals:
                raise ValueError(f"{name} must contain one value per orbital.")
        for bond in self.bonds:
            if not 0 <= bond.source < self.n_orbitals or not 0 <= bond.target < self.n_orbitals:
                raise ValueError("Periodic bond orbital indices are out of range.")
            if len(bond.displacement) != dimension:
                raise ValueError("Periodic bond displacements must match the lattice dimension.")
            if not np.isfinite(complex(bond.value).real) or not np.isfinite(
                complex(bond.value).imag
            ):
                raise ValueError("Periodic bond values must be finite.")

    def reciprocal_vectors(self) -> np.ndarray:
        """Return reciprocal primitive vectors as rows with ``a_i · b_j = 2πδ_ij``."""

        self.validate()
        return 2.0 * np.pi * np.linalg.inv(np.asarray(self.primitive_vectors)).T

    def bloch_hamiltonian(
        self,
        momentum: Sequence[float],
        *,
        gauge: str = "cell",
    ) -> np.ndarray:
        """Return ``H(k)`` using cell-periodic or orbital-position phases."""

        self.validate()
        k = np.asarray(momentum, dtype=float)
        if k.shape != (self.dimension,):
            raise ValueError("momentum must contain one component per lattice dimension.")
        if gauge not in {"cell", "orbital"}:
            raise ValueError("gauge must be 'cell' or 'orbital'.")
        matrix = np.diag(np.asarray(self.onsite or (0.0,) * self.n_orbitals, dtype=complex))
        vectors = np.asarray(self.primitive_vectors)
        positions = np.asarray(self.orbital_positions)
        for bond in self.bonds:
            translation = np.asarray(bond.displacement) @ vectors
            if gauge == "orbital":
                translation = translation + positions[bond.target] - positions[bond.source]
            amplitude = complex(bond.value) * np.exp(1j * np.dot(k, translation))
            matrix[bond.source, bond.target] += amplitude
            if bond.source != bond.target or any(bond.displacement):
                matrix[bond.target, bond.source] += np.conjugate(amplitude)
        return matrix

    def bands(
        self,
        path: MomentumPath | np.ndarray,
        *,
        eigenvectors: bool = True,
        gauge: str = "cell",
    ) -> BandStructure:
        """Evaluate sorted bands and optional eigenvectors along a momentum path."""

        momenta = path.momenta if isinstance(path, MomentumPath) else np.asarray(path, dtype=float)
        if momenta.ndim != 2 or momenta.shape[1] != self.dimension:
            raise ValueError("Band momenta must have shape (n_points, dimension).")
        energies = []
        vectors = []
        for momentum in momenta:
            values, states = np.linalg.eigh(self.bloch_hamiltonian(momentum, gauge=gauge))
            energies.append(values)
            if eigenvectors:
                vectors.append(states)
        distances = _path_distances(momenta)
        return BandStructure(
            momenta=momenta,
            distances=distances,
            energies=np.asarray(energies),
            eigenvectors=np.asarray(vectors) if eigenvectors else None,
            labels=path.labels if isinstance(path, MomentumPath) else (),
            label_indices=path.label_indices if isinstance(path, MomentumPath) else (),
            metadata={"gauge": gauge, "n_orbitals": self.n_orbitals},
        )

    def to_dict(self) -> dict[str, object]:
        self.validate()
        return {
            "schema_version": self.schema_version,
            "primitive_vectors": [list(vector) for vector in self.primitive_vectors],
            "orbital_positions": [list(position) for position in self.orbital_positions],
            "bonds": [bond.to_dict() for bond in self.bonds],
            "onsite": [_complex_dict(value) for value in self.onsite],
            "orbital_labels": list(self.orbital_labels),
            "sublattice_labels": list(self.sublattice_labels),
            "units": dict(self.units),
            "conventions": dict(self.conventions),
            "references": list(self.references),
            "provenance": _json_value(self.provenance),
            "metadata": _json_value(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> PeriodicLatticeSpec:
        allowed = {
            "schema_version",
            "primitive_vectors",
            "orbital_positions",
            "bonds",
            "onsite",
            "orbital_labels",
            "sublattice_labels",
            "units",
            "conventions",
            "references",
            "provenance",
            "metadata",
        }
        unknown = sorted(set(data) - allowed)
        if unknown:
            raise ValueError(
                f"Periodic specification contains unknown fields: {', '.join(unknown)}."
            )
        spec = cls(
            schema_version=str(data.get("schema_version", PERIODIC_SCHEMA_VERSION)),
            primitive_vectors=tuple(
                tuple(float(item) for item in vector) for vector in data["primitive_vectors"]
            ),
            orbital_positions=tuple(
                tuple(float(item) for item in position) for position in data["orbital_positions"]
            ),
            bonds=tuple(CellBond.from_dict(item) for item in data.get("bonds", [])),
            onsite=tuple(_decode_complex(item) for item in data.get("onsite", [])),
            orbital_labels=tuple(str(item) for item in data.get("orbital_labels", [])),
            sublattice_labels=tuple(str(item) for item in data.get("sublattice_labels", [])),
            units={str(key): str(value) for key, value in dict(data.get("units", {})).items()},
            conventions={
                str(key): str(value) for key, value in dict(data.get("conventions", {})).items()
            },
            references=tuple(str(item) for item in data.get("references", [])),
            provenance=tuple(dict(item) for item in data.get("provenance", [])),
            metadata=dict(data.get("metadata", {})),
        )
        spec.validate()
        return spec

    def save(self, path: str | Path) -> Path:
        output = Path(path)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(json.dumps(self.to_dict(), indent=2, sort_keys=True) + "\n")
        return output


@dataclass(frozen=True)
class MomentumPath:
    """Sampled momentum path with labeled high-symmetry vertices."""

    momenta: np.ndarray
    labels: tuple[str, ...] = ()
    label_indices: tuple[int, ...] = ()

    @classmethod
    def from_vertices(
        cls,
        vertices: Iterable[Sequence[float]],
        *,
        labels: Sequence[str] = (),
        points_per_segment: int = 64,
    ) -> MomentumPath:
        vertices_array = np.asarray(tuple(vertices), dtype=float)
        if vertices_array.ndim != 2 or len(vertices_array) < 2:
            raise ValueError("Momentum paths require at least two vertices.")
        if points_per_segment < 2:
            raise ValueError("points_per_segment must be at least 2.")
        pieces = []
        indices = [0]
        for index in range(len(vertices_array) - 1):
            piece = np.linspace(
                vertices_array[index],
                vertices_array[index + 1],
                points_per_segment,
            )
            if index:
                piece = piece[1:]
            pieces.append(piece)
            indices.append(indices[-1] + len(piece) - (0 if index == 0 else 0))
        momenta = np.concatenate(pieces)
        indices = [index * (points_per_segment - 1) for index in range(len(vertices_array))]
        if labels and len(labels) != len(vertices_array):
            raise ValueError("Momentum path labels must match the number of vertices.")
        return cls(momenta, tuple(labels), tuple(indices))


@dataclass(frozen=True)
class BandStructure:
    """Band energies and optional eigenvectors sampled over momentum."""

    momenta: np.ndarray
    distances: np.ndarray
    energies: np.ndarray
    eigenvectors: np.ndarray | None = None
    labels: tuple[str, ...] = ()
    label_indices: tuple[int, ...] = ()
    metadata: dict[str, object] = field(default_factory=dict)

    def to_dict(self, *, include_eigenvectors: bool = False) -> dict[str, object]:
        result = {
            "momenta": self.momenta.tolist(),
            "distances": self.distances.tolist(),
            "energies": self.energies.tolist(),
            "labels": list(self.labels),
            "label_indices": list(self.label_indices),
            "metadata": _json_value(self.metadata),
        }
        if include_eigenvectors and self.eigenvectors is not None:
            result["eigenvectors"] = _json_value(self.eigenvectors)
        return result

    def to_analysis_result(
        self,
        *,
        periodic: PeriodicLatticeSpec | None = None,
        parameters: dict[str, object] | None = None,
    ):
        """Return this band structure as a portable ``AnalysisResult``."""

        from quantum_lattice_models.analysis import band_structure_result

        return band_structure_result(self, periodic=periodic, parameters=parameters)


def load_periodic_lattice(path: str | Path) -> PeriodicLatticeSpec:
    data = json.loads(Path(path).read_text())
    if not isinstance(data, dict):
        raise ValueError("Periodic lattice JSON must contain an object.")
    return PeriodicLatticeSpec.from_dict(data)


def ssh_unit_cell(t1: complex = 1.0, t2: complex = 1.0) -> PeriodicLatticeSpec:
    return PeriodicLatticeSpec(
        primitive_vectors=((1.0,),),
        orbital_positions=((0.0,), (0.5,)),
        bonds=(CellBond(0, 1, (0,), -t1), CellBond(1, 0, (1,), -t2)),
        orbital_labels=("A", "B"),
        sublattice_labels=("A", "B"),
        conventions={"bloch_basis": "A, B", "bond_values": "Hamiltonian matrix elements"},
        metadata={"family": "ssh"},
    )


def rice_mele_unit_cell(
    t1: complex = 1.0,
    t2: complex = 1.0,
    staggered_potential: float = 0.0,
) -> PeriodicLatticeSpec:
    return PeriodicLatticeSpec(
        primitive_vectors=((1.0,),),
        orbital_positions=((0.0,), (0.5,)),
        bonds=(CellBond(0, 1, (0,), -t1), CellBond(1, 0, (1,), -t2)),
        onsite=(staggered_potential, -staggered_potential),
        orbital_labels=("A", "B"),
        sublattice_labels=("A", "B"),
        metadata={"family": "rice_mele"},
    )


def square_unit_cell(hopping: complex = 1.0, onsite: complex = 0.0) -> PeriodicLatticeSpec:
    return PeriodicLatticeSpec(
        primitive_vectors=((1.0, 0.0), (0.0, 1.0)),
        orbital_positions=((0.0, 0.0),),
        bonds=(CellBond(0, 0, (1, 0), -hopping), CellBond(0, 0, (0, 1), -hopping)),
        onsite=(onsite,),
        orbital_labels=("site",),
        metadata={"family": "square"},
    )


def honeycomb_unit_cell(hopping: complex = 1.0) -> PeriodicLatticeSpec:
    return PeriodicLatticeSpec(
        primitive_vectors=((1.5, np.sqrt(3.0) / 2.0), (1.5, -np.sqrt(3.0) / 2.0)),
        orbital_positions=((0.0, 0.0), (1.0, 0.0)),
        bonds=(
            CellBond(0, 1, (0, 0), -hopping),
            CellBond(0, 1, (-1, 0), -hopping),
            CellBond(0, 1, (0, -1), -hopping),
        ),
        orbital_labels=("A", "B"),
        sublattice_labels=("A", "B"),
        metadata={"family": "honeycomb"},
    )


def kagome_unit_cell(hopping: complex = 1.0) -> PeriodicLatticeSpec:
    return PeriodicLatticeSpec(
        primitive_vectors=((1.0, 0.0), (0.5, np.sqrt(3.0) / 2.0)),
        orbital_positions=((0.0, 0.0), (0.5, 0.0), (0.25, np.sqrt(3.0) / 4.0)),
        bonds=(
            CellBond(0, 1, (0, 0), -hopping),
            CellBond(0, 2, (0, 0), -hopping),
            CellBond(1, 2, (0, 0), -hopping),
            CellBond(0, 1, (-1, 0), -hopping),
            CellBond(0, 2, (0, -1), -hopping),
            CellBond(1, 2, (1, -1), -hopping),
        ),
        orbital_labels=("A", "B", "C"),
        sublattice_labels=("A", "B", "C"),
        metadata={"family": "kagome"},
    )


def haldane_unit_cell(
    t1: float = 1.0,
    t2: float = 0.1,
    phi: float = np.pi / 2,
    sublattice_potential: float = 0.0,
) -> PeriodicLatticeSpec:
    phase_a = -t2 * np.exp(1j * phi)
    phase_b = -t2 * np.exp(-1j * phi)
    nearest = (
        CellBond(0, 1, (0, 0), -t1),
        CellBond(0, 1, (-1, 0), -t1),
        CellBond(0, 1, (0, -1), -t1),
    )
    next_nearest = tuple(
        CellBond(sublattice, sublattice, displacement, phase)
        for sublattice, phase in ((0, phase_a), (1, phase_b))
        for displacement in ((1, 0), (0, 1), (-1, 1))
    )
    return PeriodicLatticeSpec(
        primitive_vectors=((1.5, np.sqrt(3.0) / 2.0), (1.5, -np.sqrt(3.0) / 2.0)),
        orbital_positions=((0.0, 0.0), (1.0, 0.0)),
        bonds=nearest + next_nearest,
        onsite=(sublattice_potential, -sublattice_potential),
        orbital_labels=("A", "B"),
        sublattice_labels=("A", "B"),
        metadata={"family": "haldane", "phi": phi},
    )


def standard_momentum_path(
    lattice: PeriodicLatticeSpec,
    *,
    points_per_segment: int = 64,
) -> MomentumPath:
    """Return a conventional path for supported one- and two-dimensional cells."""

    family = str(lattice.metadata.get("family", ""))
    reciprocal = lattice.reciprocal_vectors()
    if lattice.dimension == 1:
        return MomentumPath.from_vertices(
            ((-np.pi,), (0.0,), (np.pi,)),
            labels=("−π", "Γ", "π"),
            points_per_segment=points_per_segment,
        )
    gamma = np.zeros(2)
    if family in {"honeycomb", "kagome", "haldane"}:
        m_point = 0.5 * reciprocal[0]
        k_point = (2.0 * reciprocal[0] + reciprocal[1]) / 3.0
        return MomentumPath.from_vertices(
            (gamma, k_point, m_point, gamma),
            labels=("Γ", "K", "M", "Γ"),
            points_per_segment=points_per_segment,
        )
    x_point = 0.5 * reciprocal[0]
    m_point = 0.5 * (reciprocal[0] + reciprocal[1])
    return MomentumPath.from_vertices(
        (gamma, x_point, m_point, gamma),
        labels=("Γ", "X", "M", "Γ"),
        points_per_segment=points_per_segment,
    )


def bloch_function(
    source: PeriodicLatticeSpec | Callable[[np.ndarray], np.ndarray],
) -> Callable[[np.ndarray], np.ndarray]:
    if isinstance(source, PeriodicLatticeSpec):
        return source.bloch_hamiltonian
    return source


def _path_distances(momenta: np.ndarray) -> np.ndarray:
    distances = np.zeros(len(momenta), dtype=float)
    if len(momenta) > 1:
        distances[1:] = np.cumsum(np.linalg.norm(np.diff(momenta, axis=0), axis=1))
    return distances


def _complex_dict(value: complex) -> dict[str, list[float]]:
    number = complex(value)
    return {"__complex__": [number.real, number.imag]}


def _decode_complex(value: object) -> complex:
    if not isinstance(value, dict) or set(value) != {"__complex__"}:
        raise ValueError("Complex values require {'__complex__': [real, imaginary]}.")
    parts = value["__complex__"]
    if not isinstance(parts, list) or len(parts) != 2:
        raise ValueError("Complex values require [real, imaginary].")
    return complex(float(parts[0]), float(parts[1]))


def _json_value(value: object) -> object:
    if isinstance(value, np.ndarray):
        return _json_value(value.tolist())
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, complex):
        return _complex_dict(value)
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    return value
