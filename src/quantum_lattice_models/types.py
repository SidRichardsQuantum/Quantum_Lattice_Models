"""Shared types for Hamiltonian metadata."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class PauliTerm:
    """A product of single-site Pauli operators with a scalar coefficient."""

    coefficient: complex
    operators: tuple[str, ...]


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
