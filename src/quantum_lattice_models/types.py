"""Shared types for Hamiltonian metadata."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np
import scipy.sparse as sp

if TYPE_CHECKING:
    from quantum_lattice_models.specs import ModelSpec


@dataclass(frozen=True)
class PauliTerm:
    """A product of single-site Pauli operators with a scalar coefficient."""

    coefficient: complex
    operators: tuple[str, ...]


@dataclass(frozen=True)
class HamiltonianResult:
    """Hamiltonian matrix together with portable construction metadata."""

    matrix: np.ndarray | sp.csr_matrix
    model: ModelSpec
    basis: str
    representation: str
    metadata: dict[str, object]

    def __post_init__(self) -> None:
        if self.representation not in {"dense", "sparse"}:
            raise ValueError("representation must be 'dense' or 'sparse'.")
        if self.representation == "sparse" and not sp.issparse(self.matrix):
            raise ValueError("Sparse Hamiltonian results require a SciPy sparse matrix.")
        if self.representation == "dense" and sp.issparse(self.matrix):
            raise ValueError("Dense Hamiltonian results require a NumPy array.")
        shape = self.matrix.shape
        if len(shape) != 2 or shape[0] != shape[1]:
            raise ValueError("Hamiltonian result matrices must be square.")

    @property
    def shape(self) -> tuple[int, int]:
        """Return the matrix shape."""

        return self.matrix.shape

    def to_metadata(self) -> dict[str, object]:
        """Return JSON-compatible metadata for persistence and inspection."""

        return {
            "model": self.model.to_dict(),
            "basis": self.basis,
            "representation": self.representation,
            "matrix_shape": list(self.matrix.shape),
            "matrix_dtype": str(self.matrix.dtype),
            "metadata": dict(self.metadata),
        }


class DenseHamiltonian(np.ndarray):
    """Dense NumPy array with optional spin-chain Pauli-term metadata."""

    model_name: str | None
    n_sites: int | None
    terms: tuple[PauliTerm, ...]

    def __new__(
        cls,
        matrix: np.ndarray,
        *,
        model_name: str | None = None,
        n_sites: int | None = None,
        terms: tuple[PauliTerm, ...] = (),
    ) -> DenseHamiltonian:
        obj = np.asarray(matrix, dtype=complex).view(cls)
        obj.model_name = model_name
        obj.n_sites = n_sites
        obj.terms = tuple(terms)
        return obj

    def __array_finalize__(self, obj: object) -> None:
        if obj is None:
            return
        self.model_name = getattr(obj, "model_name", None)
        self.n_sites = getattr(obj, "n_sites", None)
        self.terms = getattr(obj, "terms", ())


class LatticeHamiltonian(np.ndarray):
    """Dense NumPy array with optional lattice-model metadata."""

    model_name: str | None
    basis: str | None
    lattice_shape: tuple[int, ...] | None
    metadata: dict[str, object]

    def __new__(
        cls,
        matrix: np.ndarray,
        *,
        model_name: str | None = None,
        basis: str | None = None,
        lattice_shape: tuple[int, ...] | None = None,
        metadata: dict[str, object] | None = None,
    ) -> LatticeHamiltonian:
        obj = np.asarray(matrix, dtype=complex).view(cls)
        obj.model_name = model_name
        obj.basis = basis
        obj.lattice_shape = lattice_shape
        obj.metadata = dict(metadata or {})
        return obj

    def __array_finalize__(self, obj: object) -> None:
        if obj is None:
            return
        self.model_name = getattr(obj, "model_name", None)
        self.basis = getattr(obj, "basis", None)
        self.lattice_shape = getattr(obj, "lattice_shape", None)
        self.metadata = getattr(obj, "metadata", {})
