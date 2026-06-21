"""Structured registry for available lattice model builders."""

from __future__ import annotations

import inspect
from collections.abc import Callable
from dataclasses import dataclass, field, replace
from typing import get_args, get_type_hints

import numpy as np

from quantum_lattice_models import benchmarks, hubbard, lattice, spin, tight_binding, topological


@dataclass(frozen=True)
class ParameterInfo:
    """Typed metadata for one registered model parameter."""

    name: str
    type: type
    default: object = None
    required: bool = False
    description: str = ""
    cli_name: str | None = None
    minimum: float | None = None
    choices: tuple[object, ...] = ()
    multiple: bool = False


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
    parameters: tuple[ParameterInfo, ...] = ()
    validation_status: str = "tested"


@dataclass(frozen=True)
class ModelPreset:
    """Named parameter set for a canonical phase or reference limit."""

    name: str
    model: str
    description: str
    parameters: dict[str, object]


_PARAMETER_DESCRIPTIONS = {
    "n_sites": "Number of lattice or spin sites.",
    "n_cells": "Number of unit cells.",
    "n_rows": "Number of lattice rows.",
    "n_cols": "Number of lattice columns.",
    "n_rungs": "Number of ladder rungs.",
    "hopping": "Nearest-neighbor hopping amplitude.",
    "onsite": "Uniform onsite potential.",
    "j": "Nearest-neighbor Ising coupling.",
    "jx": "XX coupling.",
    "jy": "YY coupling.",
    "jz": "ZZ coupling.",
    "j1": "Nearest-neighbor coupling.",
    "j2": "Next-nearest-neighbor coupling.",
    "h": "Transverse field strength.",
    "h_x": "Transverse x-field strength.",
    "h_z": "Longitudinal z-field strength.",
    "field": "Uniform z-field strength.",
    "magnetization": "Total Pauli-Z eigenvalue for a conserved spin sector.",
    "coupling": "Overall interaction strength.",
    "anisotropy": "Interaction anisotropy.",
    "leg_coupling": "Spin-ladder leg coupling.",
    "rung_coupling": "Spin-ladder rung coupling.",
    "t1": "Intracell or nearest-neighbor hopping.",
    "t2": "Intercell or next-nearest-neighbor hopping.",
    "dimerization": "Alternating hopping offset.",
    "staggering": "Staggered onsite potential.",
    "flux": "Magnetic flux per plaquette.",
    "potential": "Quasiperiodic potential strength.",
    "beta": "Quasiperiodic inverse wavelength.",
    "phase": "Quasiperiodic phase offset.",
    "phi": "Complex next-nearest-neighbor hopping phase.",
    "pairing": "Superconducting pairing amplitude.",
    "chemical_potential": "Chemical potential.",
    "interaction": "Onsite interaction strength.",
    "max_occupancy": "Maximum boson occupancy per site.",
    "sublattice_potential": "Staggered honeycomb sublattice potential.",
    "periodic": "Enable periodic chain boundaries.",
    "periodic_x": "Enable periodic boundaries along x.",
    "periodic_y": "Enable periodic boundaries along y.",
    "hermitian": "Add Hermitian-conjugate reverse bonds.",
    "bonds": "Custom bond as source,target or source,target,value.",
    "disorder": "Disorder amplitude.",
    "seed": "Reproducible random seed.",
    "power": "Power-law decay exponent.",
    "diagonal_hopping": "Diagonal hopping amplitude.",
    "mass": "Sublattice or leg mass term.",
    "base_hopping": "Sawtooth-chain base hopping.",
    "tooth_hopping": "Sawtooth-chain tooth hopping.",
}

_EXCLUDED_CLI_PARAMETERS = {"lattice", "positions", "model_name"}
_POSITIVE_PARAMETERS = {"n_sites", "n_cells", "n_rows", "n_cols", "n_rungs", "max_occupancy"}


def _info(
    category: str,
    basis: str,
    dimension: str,
    return_type: str,
    description: str,
    builder: Callable[..., object],
    defaults: dict[str, object],
    validation_status: str = "tested",
) -> ModelInfo:
    parameters = _parameter_schema(builder, defaults)
    return ModelInfo(
        builder.__name__,
        category,
        basis,
        dimension,
        return_type,
        description,
        builder,
        defaults,
        parameters,
        validation_status,
    )


def _parameter_schema(
    builder: Callable[..., object],
    defaults: dict[str, object],
) -> tuple[ParameterInfo, ...]:
    signature = inspect.signature(builder)
    try:
        hints = get_type_hints(builder)
    except (NameError, TypeError):
        hints = {}
    parameters = []
    for name, parameter in signature.parameters.items():
        if name in _EXCLUDED_CLI_PARAMETERS:
            continue
        default = defaults.get(
            name,
            None if parameter.default is inspect.Parameter.empty else parameter.default,
        )
        parameter_type = _parameter_type(name, hints.get(name), default)
        parameters.append(
            ParameterInfo(
                name=name,
                type=parameter_type,
                default=default,
                required=parameter.default is inspect.Parameter.empty and name not in defaults,
                description=_PARAMETER_DESCRIPTIONS.get(name, name.replace("_", " ").capitalize()),
                cli_name="--bond" if name == "bonds" else f"--{name.replace('_', '-')}",
                minimum=1 if name in _POSITIVE_PARAMETERS else None,
                multiple=name == "bonds",
            )
        )
    return tuple(parameters)


def _parameter_type(name: str, annotation: object, default: object) -> type:
    if name == "bonds":
        return str
    if annotation is not None:
        candidates = get_args(annotation)
        for candidate in candidates or (annotation,):
            if candidate in (bool, int, float, complex, str):
                return candidate
    if isinstance(default, bool):
        return bool
    if isinstance(default, int):
        return int
    if isinstance(default, complex):
        return complex
    if isinstance(default, float):
        return float
    return float


def validate_parameter(parameter: ParameterInfo, value: object) -> None:
    """Validate a value against registry parameter constraints."""

    if value is None:
        if parameter.required:
            raise ValueError(f"{parameter.name} is required.")
        return
    if parameter.name not in {"bonds", "onsite"} and not _matches_parameter_type(
        parameter.type, value
    ):
        raise ValueError(
            f"{parameter.name} must have type {parameter.type.__name__}; "
            f"received {type(value).__name__}."
        )
    if parameter.name == "onsite" and not _matches_onsite_value(value):
        raise ValueError("onsite must be a numeric scalar or numeric sequence.")
    if parameter.minimum is not None and value < parameter.minimum:
        raise ValueError(f"{parameter.name} must be at least {parameter.minimum}.")
    if parameter.choices and value not in parameter.choices:
        raise ValueError(f"{parameter.name} must be one of {parameter.choices!r}.")


def _matches_parameter_type(expected: type, value: object) -> bool:
    if expected is bool:
        return isinstance(value, bool)
    if expected is int:
        return isinstance(value, int) and not isinstance(value, bool)
    if expected is float:
        return isinstance(value, (int, float)) and not isinstance(value, bool)
    if expected is complex:
        return isinstance(value, (int, float, complex)) and not isinstance(value, bool)
    return isinstance(value, expected)


def _matches_onsite_value(value: object) -> bool:
    if _matches_parameter_type(float, value):
        return True
    if not isinstance(value, (tuple, list)):
        return False
    return all(_matches_parameter_type(float, item) for item in value)


MODEL_REGISTRY: dict[str, ModelInfo] = {
    "transverse_field_ising": _info(
        "spin",
        "qubit",
        "2**n_sites",
        "DenseHamiltonian",
        "TFIM spin chain",
        spin.transverse_field_ising,
        {"n_sites": 4, "j": 1.0, "h": 0.5},
    ),
    "longitudinal_field_ising": _info(
        "spin",
        "qubit",
        "2**n_sites",
        "DenseHamiltonian",
        "Ising chain with transverse and longitudinal fields",
        spin.longitudinal_field_ising,
        {"n_sites": 4, "j": 1.0, "h_x": 0.5, "h_z": 0.1},
    ),
    "next_nearest_neighbor_ising": _info(
        "spin",
        "qubit",
        "2**n_sites",
        "DenseHamiltonian",
        "Frustrated Ising chain with next-nearest-neighbor ZZ couplings",
        spin.next_nearest_neighbor_ising,
        {"n_sites": 5, "j1": 1.0, "j2": 0.25, "h": 0.5},
    ),
    "heisenberg_chain": _info(
        "spin",
        "qubit",
        "2**n_sites",
        "DenseHamiltonian",
        "Anisotropic Heisenberg chain",
        spin.heisenberg_chain,
        {"n_sites": 4, "jx": 1.0, "jy": 1.0, "jz": 1.0},
    ),
    "xy_chain": _info(
        "spin",
        "qubit",
        "2**n_sites",
        "DenseHamiltonian",
        "Anisotropic XY chain",
        spin.xy_chain,
        {"n_sites": 4, "coupling": 1.0, "anisotropy": 0.3},
    ),
    "xxz_chain": _info(
        "spin",
        "qubit",
        "2**n_sites",
        "DenseHamiltonian",
        "XXZ spin chain",
        spin.xxz_chain,
        {"n_sites": 4, "coupling": 1.0, "anisotropy": 0.7},
    ),
    "j1_j2_heisenberg_chain": _info(
        "spin",
        "qubit",
        "2**n_sites",
        "DenseHamiltonian",
        "Frustrated J1-J2 Heisenberg chain",
        spin.j1_j2_heisenberg_chain,
        {"n_sites": 5, "j1": 1.0, "j2": 0.4},
    ),
    "heisenberg_ladder": _info(
        "spin",
        "qubit",
        "2**(2*n_rungs)",
        "DenseHamiltonian",
        "Two-leg Heisenberg ladder",
        spin.heisenberg_ladder,
        {"n_rungs": 2, "leg_coupling": 1.0, "rung_coupling": 0.7},
    ),
    "bose_hubbard_chain": _info(
        "hubbard",
        "truncated boson occupation",
        "(max_occupancy+1)**n_sites",
        "LatticeHamiltonian",
        "Truncated Bose-Hubbard chain",
        hubbard.bose_hubbard_chain,
        {"n_sites": 3, "hopping": 0.6, "interaction": 1.5, "max_occupancy": 2},
    ),
    "bose_hubbard_chain_sparse": _info(
        "hubbard",
        "truncated boson occupation",
        "(max_occupancy+1)**n_sites",
        "scipy.sparse.csr_matrix",
        "Sparse truncated Bose-Hubbard chain",
        hubbard.bose_hubbard_chain_sparse,
        {"n_sites": 3, "hopping": 0.6, "interaction": 1.5, "max_occupancy": 2},
    ),
    "fermi_hubbard_chain": _info(
        "hubbard",
        "spinful fermion occupation",
        "2**(2*n_sites)",
        "LatticeHamiltonian",
        "Spinful Fermi-Hubbard chain",
        hubbard.fermi_hubbard_chain,
        {"n_sites": 3, "hopping": 0.5, "interaction": 3.0},
    ),
    "fermi_hubbard_chain_sparse": _info(
        "hubbard",
        "spinful fermion occupation",
        "2**(2*n_sites)",
        "scipy.sparse.csr_matrix",
        "Sparse spinful Fermi-Hubbard chain",
        hubbard.fermi_hubbard_chain_sparse,
        {"n_sites": 3, "hopping": 0.5, "interaction": 3.0},
    ),
    "custom_tight_binding": _info(
        "user",
        "single particle",
        "n_sites",
        "LatticeHamiltonian",
        "User-defined graph tight-binding model",
        lattice.custom_tight_binding,
        {"n_sites": 3, "bonds": ((0, 1), (1, 2)), "hopping": 1.0},
    ),
    "custom_tight_binding_sparse": _info(
        "user",
        "single particle",
        "n_sites",
        "scipy.sparse.csr_matrix",
        "Sparse user-defined graph tight-binding model",
        lattice.custom_tight_binding_sparse,
        {"n_sites": 16, "bonds": ((0, 1), (1, 2)), "hopping": 1.0},
    ),
    "ssh_model": _info(
        "tight_binding",
        "single particle",
        "2*n_cells",
        "LatticeHamiltonian",
        "SSH chain",
        tight_binding.ssh_model,
        {"n_cells": 8, "t1": 0.5, "t2": 1.0},
    ),
    "rice_mele_model": _info(
        "tight_binding",
        "single particle",
        "2*n_cells",
        "LatticeHamiltonian",
        "Rice-Mele chain",
        tight_binding.rice_mele_model,
        {"n_cells": 8, "hopping": 1.0, "dimerization": 0.25, "staggering": 0.5},
    ),
    "tight_binding_chain": _info(
        "tight_binding",
        "single particle",
        "n_sites",
        "LatticeHamiltonian",
        "Generic 1D tight-binding chain",
        tight_binding.tight_binding_chain,
        {"n_sites": 8, "hopping": 1.0},
    ),
    "tight_binding_chain_sparse": _info(
        "tight_binding",
        "single particle",
        "n_sites",
        "scipy.sparse.csr_matrix",
        "Sparse generic 1D tight-binding chain",
        tight_binding.tight_binding_chain_sparse,
        {"n_sites": 32, "hopping": 1.0},
    ),
    "square_lattice_tight_binding": _info(
        "tight_binding",
        "single particle",
        "n_rows*n_cols",
        "LatticeHamiltonian",
        "Rectangular square-lattice tight-binding model",
        tight_binding.square_lattice_tight_binding,
        {"n_rows": 3, "n_cols": 4},
    ),
    "square_lattice_tight_binding_sparse": _info(
        "tight_binding",
        "single particle",
        "n_rows*n_cols",
        "scipy.sparse.csr_matrix",
        "Sparse rectangular square-lattice tight-binding model",
        tight_binding.square_lattice_tight_binding_sparse,
        {"n_rows": 8, "n_cols": 8},
    ),
    "aubry_andre_harper_chain": _info(
        "tight_binding",
        "single particle",
        "n_sites",
        "LatticeHamiltonian",
        "Aubry-Andre-Harper quasiperiodic chain",
        tight_binding.aubry_andre_harper_chain,
        {"n_sites": 16, "hopping": 1.0, "potential": 1.5},
    ),
    "triangular_lattice_tight_binding": _info(
        "tight_binding",
        "single particle",
        "n_rows*n_cols",
        "LatticeHamiltonian",
        "Triangular-lattice tight-binding model",
        tight_binding.triangular_lattice_tight_binding,
        {"n_rows": 3, "n_cols": 3},
    ),
    "triangular_lattice_tight_binding_sparse": _info(
        "tight_binding",
        "single particle",
        "n_rows*n_cols",
        "scipy.sparse.csr_matrix",
        "Sparse triangular-lattice tight-binding model",
        tight_binding.triangular_lattice_tight_binding_sparse,
        {"n_rows": 8, "n_cols": 8},
    ),
    "kagome_lattice_tight_binding": _info(
        "tight_binding",
        "single particle",
        "3*n_rows*n_cols",
        "LatticeHamiltonian",
        "Kagome tight-binding lattice",
        tight_binding.kagome_lattice_tight_binding,
        {"n_rows": 2, "n_cols": 3},
    ),
    "kagome_lattice_tight_binding_sparse": _info(
        "tight_binding",
        "single particle",
        "3*n_rows*n_cols",
        "scipy.sparse.csr_matrix",
        "Sparse kagome tight-binding lattice",
        tight_binding.kagome_lattice_tight_binding_sparse,
        {"n_rows": 4, "n_cols": 4},
    ),
    "harper_hofstadter_square_lattice": _info(
        "topological",
        "single particle",
        "n_rows*n_cols",
        "LatticeHamiltonian",
        "Harper-Hofstadter square lattice",
        topological.harper_hofstadter_square_lattice,
        {"n_rows": 4, "n_cols": 4, "flux": 0.25},
    ),
    "harper_hofstadter_square_lattice_sparse": _info(
        "topological",
        "single particle",
        "n_rows*n_cols",
        "scipy.sparse.csr_matrix",
        "Sparse Harper-Hofstadter square lattice",
        topological.harper_hofstadter_square_lattice_sparse,
        {"n_rows": 8, "n_cols": 8, "flux": 0.25},
    ),
    "kitaev_chain_bdg": _info(
        "topological",
        "Nambu single particle",
        "2*n_sites",
        "LatticeHamiltonian",
        "Kitaev-chain Bogoliubov-de Gennes matrix",
        topological.kitaev_chain_bdg,
        {"n_sites": 8, "hopping": 1.0, "pairing": 0.5},
    ),
    "haldane_honeycomb_lattice": _info(
        "topological",
        "single particle",
        "2*n_rows*n_cols",
        "LatticeHamiltonian",
        "Finite Haldane honeycomb lattice",
        topological.haldane_honeycomb_lattice,
        {"n_rows": 3, "n_cols": 3, "t1": 1.0, "t2": 0.18},
    ),
    "haldane_honeycomb_lattice_sparse": _info(
        "topological",
        "single particle",
        "2*n_rows*n_cols",
        "scipy.sparse.csr_matrix",
        "Sparse finite Haldane honeycomb lattice",
        topological.haldane_honeycomb_lattice_sparse,
        {"n_rows": 4, "n_cols": 4, "t1": 1.0, "t2": 0.18},
    ),
}

for _builder, _category, _basis, _dimension, _description, _defaults in (
    (
        benchmarks.anderson_chain,
        "tight_binding",
        "single particle",
        "n_sites",
        "Anderson-disordered tight-binding chain",
        {"n_sites": 16, "hopping": 1.0, "disorder": 2.0, "seed": 0},
    ),
    (
        benchmarks.long_range_tight_binding_chain,
        "tight_binding",
        "single particle",
        "n_sites",
        "Power-law long-range tight-binding chain",
        {"n_sites": 12, "hopping": 1.0, "power": 3.0},
    ),
    (
        benchmarks.creutz_ladder,
        "topological",
        "single particle",
        "2*n_cells",
        "Creutz ladder",
        {"n_cells": 8, "hopping": 1.0, "diagonal_hopping": 1.0, "flux": np.pi},
    ),
    (
        benchmarks.sawtooth_chain,
        "tight_binding",
        "single particle",
        "2*n_cells",
        "Sawtooth flat-band chain",
        {"n_cells": 8, "base_hopping": 1.0, "tooth_hopping": np.sqrt(2.0)},
    ),
    (
        benchmarks.lieb_lattice,
        "tight_binding",
        "single particle",
        "3*n_rows*n_cols",
        "Finite Lieb flat-band lattice",
        {"n_rows": 3, "n_cols": 3, "hopping": 1.0},
    ),
    (
        benchmarks.xyz_chain,
        "spin",
        "qubit",
        "2**n_sites",
        "XYZ spin chain",
        {"n_sites": 4, "jx": 1.0, "jy": 0.8, "jz": 0.6},
    ),
    (
        benchmarks.random_field_heisenberg_chain,
        "spin",
        "qubit",
        "2**n_sites",
        "Random-field Heisenberg chain",
        {"n_sites": 4, "coupling": 1.0, "disorder": 1.0, "seed": 0},
    ),
):
    MODEL_REGISTRY[_builder.__name__] = _info(
        _category,
        _basis,
        _dimension,
        "LatticeHamiltonian" if _basis == "single particle" else "DenseHamiltonian",
        _description,
        _builder,
        _defaults,
    )

MODEL_PRESETS: dict[str, ModelPreset] = {
    "ssh_trivial": ModelPreset(
        "ssh_trivial",
        "ssh_model",
        "Trivial SSH phase with stronger intracell hopping.",
        {"n_cells": 8, "t1": 1.0, "t2": 0.5, "periodic": False},
    ),
    "ssh_topological": ModelPreset(
        "ssh_topological",
        "ssh_model",
        "Topological SSH phase with open-boundary edge states.",
        {"n_cells": 8, "t1": 0.5, "t2": 1.0, "periodic": False},
    ),
    "ising_zero_field": ModelPreset(
        "ising_zero_field",
        "transverse_field_ising",
        "Classical zero-field Ising reference limit.",
        {"n_sites": 4, "j": 1.0, "h": 0.0, "periodic": False},
    ),
    "ising_critical": ModelPreset(
        "ising_critical",
        "transverse_field_ising",
        "Finite-chain transverse-field Ising critical reference J=h.",
        {"n_sites": 8, "j": 1.0, "h": 1.0, "periodic": False},
    ),
    "xxz_heisenberg": ModelPreset(
        "xxz_heisenberg",
        "xxz_chain",
        "Isotropic Heisenberg point of the XXZ chain.",
        {"n_sites": 6, "coupling": 1.0, "anisotropy": 1.0, "field": 0.0},
    ),
    "aubry_andre_localized": ModelPreset(
        "aubry_andre_localized",
        "aubry_andre_harper_chain",
        "Localized Aubry-Andre regime with potential greater than 2t.",
        {"n_sites": 34, "hopping": 1.0, "potential": 3.0},
    ),
    "kitaev_particle_hole": ModelPreset(
        "kitaev_particle_hole",
        "kitaev_chain_bdg",
        "Kitaev BdG particle-hole symmetry reference parameters.",
        {"n_sites": 8, "hopping": 1.0, "pairing": 0.5, "chemical_potential": 0.3},
    ),
}

for _sparse_spin_builder, _description, _defaults, _dimension in (
    (
        spin.transverse_field_ising_sparse,
        "Sparse TFIM spin chain",
        {"n_sites": 4, "j": 1.0, "h": 0.5},
        "2**n_sites",
    ),
    (
        spin.longitudinal_field_ising_sparse,
        "Sparse Ising chain with transverse and longitudinal fields",
        {"n_sites": 4, "j": 1.0, "h_x": 0.5, "h_z": 0.1},
        "2**n_sites",
    ),
    (
        spin.next_nearest_neighbor_ising_sparse,
        "Sparse frustrated Ising chain",
        {"n_sites": 5, "j1": 1.0, "j2": 0.25, "h": 0.5},
        "2**n_sites",
    ),
    (
        spin.heisenberg_chain_sparse,
        "Sparse anisotropic Heisenberg chain",
        {"n_sites": 4, "jx": 1.0, "jy": 1.0, "jz": 1.0},
        "2**n_sites",
    ),
    (
        spin.xy_chain_sparse,
        "Sparse anisotropic XY chain",
        {"n_sites": 4, "coupling": 1.0, "anisotropy": 0.3},
        "2**n_sites",
    ),
    (
        spin.xxz_chain_sparse,
        "Sparse XXZ spin chain",
        {"n_sites": 4, "coupling": 1.0, "anisotropy": 0.7},
        "2**n_sites",
    ),
    (
        spin.j1_j2_heisenberg_chain_sparse,
        "Sparse frustrated J1-J2 Heisenberg chain",
        {"n_sites": 5, "j1": 1.0, "j2": 0.4},
        "2**n_sites",
    ),
    (
        spin.heisenberg_ladder_sparse,
        "Sparse two-leg Heisenberg ladder",
        {"n_rungs": 2, "leg_coupling": 1.0, "rung_coupling": 0.7},
        "2**(2*n_rungs)",
    ),
):
    MODEL_REGISTRY[_sparse_spin_builder.__name__] = _info(
        "spin",
        "qubit",
        _dimension,
        "scipy.sparse.csr_matrix",
        _description,
        _sparse_spin_builder,
        _defaults,
    )

for _sector_builder, _description, _defaults in (
    (
        spin.heisenberg_chain_sector_sparse,
        "Fixed-magnetization anisotropic Heisenberg chain",
        {
            "n_sites": 6,
            "magnetization": 0,
            "jx": 1.0,
            "jy": 1.0,
            "jz": 1.0,
        },
    ),
    (
        spin.xxz_chain_sector_sparse,
        "Fixed-magnetization XXZ spin chain",
        {"n_sites": 6, "magnetization": 0, "coupling": 1.0, "anisotropy": 1.0},
    ),
):
    MODEL_REGISTRY[_sector_builder.__name__] = _info(
        "spin",
        "fixed magnetization qubit",
        "comb(n_sites,(n_sites-magnetization)//2)",
        "SpinSectorHamiltonian",
        _description,
        _sector_builder,
        _defaults,
    )

for _validated_model in (
    "transverse_field_ising",
    "tight_binding_chain",
    "tight_binding_chain_sparse",
    "ssh_model",
    "bose_hubbard_chain",
    "bose_hubbard_chain_sparse",
    "fermi_hubbard_chain",
    "fermi_hubbard_chain_sparse",
    "kitaev_chain_bdg",
    "haldane_honeycomb_lattice",
    "haldane_honeycomb_lattice_sparse",
    "heisenberg_chain_sector_sparse",
    "xxz_chain_sector_sparse",
):
    MODEL_REGISTRY[_validated_model] = replace(
        MODEL_REGISTRY[_validated_model],
        validation_status="validated",
    )


def list_models(
    category: str | None = None,
    *,
    basis: str | None = None,
    sparse: bool | None = None,
    validation_status: str | None = None,
) -> tuple[str, ...]:
    """Return registered model names filtered by discovery metadata."""

    return tuple(
        sorted(
            name
            for name, info in MODEL_REGISTRY.items()
            if (category is None or info.category == category)
            and (basis is None or info.basis == basis)
            and (sparse is None or supports_sparse(name) is sparse)
            and (validation_status is None or info.validation_status == validation_status)
        )
    )


def supports_sparse(name: str) -> bool:
    """Return whether a registered model can construct a sparse Hamiltonian."""

    base = name.removesuffix("_sparse")
    return f"{base}_sparse" in MODEL_REGISTRY


def list_presets(model: str | None = None) -> tuple[str, ...]:
    """Return named presets, optionally restricted to one model."""

    return tuple(
        sorted(
            name for name, preset in MODEL_PRESETS.items() if model is None or preset.model == model
        )
    )


def get_preset(name: str) -> ModelPreset:
    """Return one named model preset."""

    try:
        return MODEL_PRESETS[name]
    except KeyError as exc:
        raise KeyError(f"Unknown preset {name!r}.") from exc


def get_model_info(name: str) -> ModelInfo:
    """Return metadata for a registered model."""

    try:
        return MODEL_REGISTRY[name]
    except KeyError as exc:
        raise KeyError(f"Unknown model {name!r}.") from exc


def register_model(
    name: str,
    *,
    category: str,
    basis: str,
    dimension: str,
    return_type: str,
    description: str,
    builder: Callable[..., object],
    defaults: dict[str, object] | None = None,
    parameters: tuple[ParameterInfo, ...] | None = None,
    overwrite: bool = False,
    validation_status: str = "unvalidated",
) -> ModelInfo:
    """Register a model builder for discovery by notebooks, docs, and the CLI."""

    if not name:
        raise ValueError("name must be nonempty.")
    if name in MODEL_REGISTRY and not overwrite:
        raise ValueError(f"Model {name!r} is already registered.")
    registered_defaults = dict(defaults or {})
    info = ModelInfo(
        name=name,
        category=category,
        basis=basis,
        dimension=dimension,
        return_type=return_type,
        description=description,
        builder=builder,
        defaults=registered_defaults,
        parameters=parameters or _parameter_schema(builder, registered_defaults),
        validation_status=validation_status,
    )
    MODEL_REGISTRY[name] = info
    return info


def unregister_model(name: str) -> ModelInfo:
    """Remove and return a registered model."""

    try:
        return MODEL_REGISTRY.pop(name)
    except KeyError as exc:
        raise KeyError(f"Unknown model {name!r}.") from exc


def model_table() -> list[dict[str, object]]:
    """Return registry metadata as dictionaries for docs, CLIs, or notebooks."""

    return [
        {
            "name": info.name,
            "category": info.category,
            "basis": info.basis,
            "dimension": info.dimension,
            "return_type": info.return_type,
            "description": info.description,
            "defaults": dict(info.defaults),
            "parameters": tuple(parameter.name for parameter in info.parameters),
            "supports_sparse": supports_sparse(info.name),
            "validation_status": info.validation_status,
        }
        for info in sorted(MODEL_REGISTRY.values(), key=lambda item: item.name)
    ]
