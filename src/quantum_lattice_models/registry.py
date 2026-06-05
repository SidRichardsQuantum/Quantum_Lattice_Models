"""Structured registry for available lattice model builders."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass, field

from quantum_lattice_models import hubbard, spin, tight_binding, topological


@dataclass(frozen=True)
class ModelInfo:
    """Metadata describing a model builder."""

    name: str
    category: str
    basis: str
    dimension: str
    return_type: str
    description: str
    builder: Callable[..., object] | None = None
    defaults: dict[str, object] = field(default_factory=dict)


def _info(
    name: str,
    category: str,
    basis: str,
    dimension: str,
    return_type: str,
    description: str,
    builder: Callable[..., object],
    defaults: dict[str, object],
) -> ModelInfo:
    return ModelInfo(name, category, basis, dimension, return_type, description, builder, defaults)


MODEL_REGISTRY: dict[str, ModelInfo] = {
    "transverse_field_ising": _info(
        "transverse_field_ising",
        "spin",
        "qubit",
        "2**n_sites",
        "DenseHamiltonian",
        "TFIM spin chain",
        spin.transverse_field_ising,
        {"n_sites": 4, "j": 1.0, "h": 0.5},
    ),
    "longitudinal_field_ising": _info(
        "longitudinal_field_ising",
        "spin",
        "qubit",
        "2**n_sites",
        "DenseHamiltonian",
        "Ising chain with transverse and longitudinal fields",
        spin.longitudinal_field_ising,
        {"n_sites": 4, "j": 1.0, "h_x": 0.5, "h_z": 0.1},
    ),
    "next_nearest_neighbor_ising": _info(
        "next_nearest_neighbor_ising",
        "spin",
        "qubit",
        "2**n_sites",
        "DenseHamiltonian",
        "Frustrated Ising chain with next-nearest-neighbor ZZ couplings",
        spin.next_nearest_neighbor_ising,
        {"n_sites": 5, "j1": 1.0, "j2": 0.25, "h": 0.5},
    ),
    "heisenberg_chain": _info(
        "heisenberg_chain",
        "spin",
        "qubit",
        "2**n_sites",
        "DenseHamiltonian",
        "Anisotropic Heisenberg chain",
        spin.heisenberg_chain,
        {"n_sites": 4, "jx": 1.0, "jy": 1.0, "jz": 1.0},
    ),
    "xy_chain": _info(
        "xy_chain",
        "spin",
        "qubit",
        "2**n_sites",
        "DenseHamiltonian",
        "Anisotropic XY chain",
        spin.xy_chain,
        {"n_sites": 4, "coupling": 1.0, "anisotropy": 0.3},
    ),
    "xxz_chain": _info(
        "xxz_chain",
        "spin",
        "qubit",
        "2**n_sites",
        "DenseHamiltonian",
        "XXZ spin chain",
        spin.xxz_chain,
        {"n_sites": 4, "coupling": 1.0, "anisotropy": 0.7},
    ),
    "j1_j2_heisenberg_chain": _info(
        "j1_j2_heisenberg_chain",
        "spin",
        "qubit",
        "2**n_sites",
        "DenseHamiltonian",
        "Frustrated J1-J2 Heisenberg chain",
        spin.j1_j2_heisenberg_chain,
        {"n_sites": 5, "j1": 1.0, "j2": 0.4},
    ),
    "heisenberg_ladder": _info(
        "heisenberg_ladder",
        "spin",
        "qubit",
        "2**(2*n_rungs)",
        "DenseHamiltonian",
        "Two-leg Heisenberg ladder",
        spin.heisenberg_ladder,
        {"n_rungs": 2, "leg_coupling": 1.0, "rung_coupling": 0.7},
    ),
    "bose_hubbard_chain": _info(
        "bose_hubbard_chain",
        "hubbard",
        "truncated boson occupation",
        "(max_occupancy+1)**n_sites",
        "LatticeHamiltonian",
        "Truncated Bose-Hubbard chain",
        hubbard.bose_hubbard_chain,
        {"n_sites": 3, "hopping": 0.6, "interaction": 1.5, "max_occupancy": 2},
    ),
    "bose_hubbard_chain_sparse": _info(
        "bose_hubbard_chain_sparse",
        "hubbard",
        "truncated boson occupation",
        "(max_occupancy+1)**n_sites",
        "scipy.sparse.csr_matrix",
        "Sparse truncated Bose-Hubbard chain",
        hubbard.bose_hubbard_chain_sparse,
        {"n_sites": 3, "hopping": 0.6, "interaction": 1.5, "max_occupancy": 2},
    ),
    "fermi_hubbard_chain": _info(
        "fermi_hubbard_chain",
        "hubbard",
        "spinful fermion occupation",
        "2**(2*n_sites)",
        "LatticeHamiltonian",
        "Spinful Fermi-Hubbard chain",
        hubbard.fermi_hubbard_chain,
        {"n_sites": 3, "hopping": 0.5, "interaction": 3.0},
    ),
    "fermi_hubbard_chain_sparse": _info(
        "fermi_hubbard_chain_sparse",
        "hubbard",
        "spinful fermion occupation",
        "2**(2*n_sites)",
        "scipy.sparse.csr_matrix",
        "Sparse spinful Fermi-Hubbard chain",
        hubbard.fermi_hubbard_chain_sparse,
        {"n_sites": 3, "hopping": 0.5, "interaction": 3.0},
    ),
    "ssh_model": _info(
        "ssh_model",
        "tight_binding",
        "single particle",
        "2*n_cells",
        "LatticeHamiltonian",
        "SSH chain",
        tight_binding.ssh_model,
        {"n_cells": 8, "t1": 0.5, "t2": 1.0},
    ),
    "rice_mele_model": _info(
        "rice_mele_model",
        "tight_binding",
        "single particle",
        "2*n_cells",
        "LatticeHamiltonian",
        "Rice-Mele chain",
        tight_binding.rice_mele_model,
        {"n_cells": 8, "hopping": 1.0, "dimerization": 0.25, "staggering": 0.5},
    ),
    "tight_binding_chain": _info(
        "tight_binding_chain",
        "tight_binding",
        "single particle",
        "n_sites",
        "LatticeHamiltonian",
        "Generic 1D tight-binding chain",
        tight_binding.tight_binding_chain,
        {"n_sites": 8, "hopping": 1.0},
    ),
    "tight_binding_chain_sparse": _info(
        "tight_binding_chain_sparse",
        "tight_binding",
        "single particle",
        "n_sites",
        "scipy.sparse.csr_matrix",
        "Sparse generic 1D tight-binding chain",
        tight_binding.tight_binding_chain_sparse,
        {"n_sites": 32, "hopping": 1.0},
    ),
    "square_lattice_tight_binding": _info(
        "square_lattice_tight_binding",
        "tight_binding",
        "single particle",
        "n_rows*n_cols",
        "LatticeHamiltonian",
        "Rectangular square-lattice tight-binding model",
        tight_binding.square_lattice_tight_binding,
        {"n_rows": 3, "n_cols": 4},
    ),
    "square_lattice_tight_binding_sparse": _info(
        "square_lattice_tight_binding_sparse",
        "tight_binding",
        "single particle",
        "n_rows*n_cols",
        "scipy.sparse.csr_matrix",
        "Sparse rectangular square-lattice tight-binding model",
        tight_binding.square_lattice_tight_binding_sparse,
        {"n_rows": 8, "n_cols": 8},
    ),
    "aubry_andre_harper_chain": _info(
        "aubry_andre_harper_chain",
        "tight_binding",
        "single particle",
        "n_sites",
        "LatticeHamiltonian",
        "Aubry-Andre-Harper quasiperiodic chain",
        tight_binding.aubry_andre_harper_chain,
        {"n_sites": 16, "hopping": 1.0, "potential": 1.5},
    ),
    "triangular_lattice_tight_binding": _info(
        "triangular_lattice_tight_binding",
        "tight_binding",
        "single particle",
        "n_rows*n_cols",
        "LatticeHamiltonian",
        "Triangular-lattice tight-binding model",
        tight_binding.triangular_lattice_tight_binding,
        {"n_rows": 3, "n_cols": 3},
    ),
    "triangular_lattice_tight_binding_sparse": _info(
        "triangular_lattice_tight_binding_sparse",
        "tight_binding",
        "single particle",
        "n_rows*n_cols",
        "scipy.sparse.csr_matrix",
        "Sparse triangular-lattice tight-binding model",
        tight_binding.triangular_lattice_tight_binding_sparse,
        {"n_rows": 8, "n_cols": 8},
    ),
    "kagome_lattice_tight_binding": _info(
        "kagome_lattice_tight_binding",
        "tight_binding",
        "single particle",
        "3*n_rows*n_cols",
        "LatticeHamiltonian",
        "Kagome tight-binding lattice",
        tight_binding.kagome_lattice_tight_binding,
        {"n_rows": 2, "n_cols": 3},
    ),
    "kagome_lattice_tight_binding_sparse": _info(
        "kagome_lattice_tight_binding_sparse",
        "tight_binding",
        "single particle",
        "3*n_rows*n_cols",
        "scipy.sparse.csr_matrix",
        "Sparse kagome tight-binding lattice",
        tight_binding.kagome_lattice_tight_binding_sparse,
        {"n_rows": 4, "n_cols": 4},
    ),
    "harper_hofstadter_square_lattice": _info(
        "harper_hofstadter_square_lattice",
        "topological",
        "single particle",
        "n_rows*n_cols",
        "LatticeHamiltonian",
        "Harper-Hofstadter square lattice",
        topological.harper_hofstadter_square_lattice,
        {"n_rows": 4, "n_cols": 4, "flux": 0.25},
    ),
    "harper_hofstadter_square_lattice_sparse": _info(
        "harper_hofstadter_square_lattice_sparse",
        "topological",
        "single particle",
        "n_rows*n_cols",
        "scipy.sparse.csr_matrix",
        "Sparse Harper-Hofstadter square lattice",
        topological.harper_hofstadter_square_lattice_sparse,
        {"n_rows": 8, "n_cols": 8, "flux": 0.25},
    ),
    "kitaev_chain_bdg": _info(
        "kitaev_chain_bdg",
        "topological",
        "Nambu single particle",
        "2*n_sites",
        "LatticeHamiltonian",
        "Kitaev-chain Bogoliubov-de Gennes matrix",
        topological.kitaev_chain_bdg,
        {"n_sites": 8, "hopping": 1.0, "pairing": 0.5},
    ),
    "haldane_honeycomb_lattice": _info(
        "haldane_honeycomb_lattice",
        "topological",
        "single particle",
        "2*n_rows*n_cols",
        "LatticeHamiltonian",
        "Finite Haldane honeycomb lattice",
        topological.haldane_honeycomb_lattice,
        {"n_rows": 3, "n_cols": 3, "t1": 1.0, "t2": 0.18},
    ),
    "haldane_honeycomb_lattice_sparse": _info(
        "haldane_honeycomb_lattice_sparse",
        "topological",
        "single particle",
        "2*n_rows*n_cols",
        "scipy.sparse.csr_matrix",
        "Sparse finite Haldane honeycomb lattice",
        topological.haldane_honeycomb_lattice_sparse,
        {"n_rows": 4, "n_cols": 4, "t1": 1.0, "t2": 0.18},
    ),
}


def list_models(category: str | None = None) -> tuple[str, ...]:
    """Return registered model names, optionally filtered by category."""

    if category is None:
        return tuple(sorted(MODEL_REGISTRY))
    return tuple(sorted(name for name, info in MODEL_REGISTRY.items() if info.category == category))


def get_model_info(name: str) -> ModelInfo:
    """Return metadata for a registered model."""

    try:
        return MODEL_REGISTRY[name]
    except KeyError as exc:
        raise KeyError(f"Unknown model {name!r}.") from exc


def model_table() -> list[dict[str, str]]:
    """Return registry metadata as dictionaries for docs, CLIs, or notebooks."""

    return [
        {
            "name": info.name,
            "category": info.category,
            "basis": info.basis,
            "dimension": info.dimension,
            "return_type": info.return_type,
            "description": info.description,
            "defaults": repr(info.defaults),
        }
        for info in MODEL_REGISTRY.values()
    ]
