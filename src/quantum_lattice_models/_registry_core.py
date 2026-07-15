"""Registry metadata and parameter-validation mechanics.

The built-in model catalog lives in :mod:`quantum_lattice_models.registry`;
this module contains the reusable, catalog-independent registry machinery.
"""

from __future__ import annotations

import inspect
from collections.abc import Callable
from dataclasses import dataclass, field
from typing import cast, get_args, get_type_hints


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
    "n_up": "Fixed number of spin-up fermions.",
    "n_down": "Fixed number of spin-down fermions.",
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


def model_info(
    category: str,
    basis: str,
    dimension: str,
    return_type: str,
    description: str,
    builder: Callable[..., object],
    defaults: dict[str, object],
    validation_status: str = "tested",
) -> ModelInfo:
    """Build catalog metadata for one model builder."""

    parameters = parameter_schema(builder, defaults)
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


def parameter_schema(
    builder: Callable[..., object],
    defaults: dict[str, object],
) -> tuple[ParameterInfo, ...]:
    """Infer CLI-facing parameter metadata from a builder signature."""

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
                return cast(type, candidate)
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
    if parameter.minimum is not None and cast(float, value) < parameter.minimum:
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
