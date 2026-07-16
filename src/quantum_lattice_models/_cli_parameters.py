"""Registry-driven CLI parameter parsing and validation."""

from __future__ import annotations

import argparse

from quantum_lattice_models.lattice import Bond
from quantum_lattice_models.registry import (
    MODEL_REGISTRY,
    ParameterInfo,
    validate_parameter,
)
from quantum_lattice_models.specs import LatticeSpec


def add_model_arguments(
    parser: argparse.ArgumentParser,
    *,
    required: bool = True,
) -> None:
    parser.add_argument(
        "--model",
        required=required,
        choices=tuple(sorted(MODEL_REGISTRY)),
    )
    add_parameter_arguments(parser)


def add_parameter_arguments(parser: argparse.ArgumentParser) -> None:
    for parameter in all_cli_parameters():
        kwargs: dict[str, object] = {
            "dest": parameter.name,
            "default": None,
            "help": parameter.description,
        }
        if parameter.type is bool:
            kwargs["action"] = argparse.BooleanOptionalAction
        else:
            kwargs["type"] = parameter.type
        if parameter.multiple:
            kwargs["action"] = "append"
        parser.add_argument(parameter.cli_name, **kwargs)


def all_cli_parameters() -> tuple[ParameterInfo, ...]:
    merged: dict[str, ParameterInfo] = {}
    for info in MODEL_REGISTRY.values():
        for parameter in info.parameters:
            existing = merged.get(parameter.name)
            if existing is None or parameter.type is complex:
                merged[parameter.name] = parameter
    return tuple(
        sorted(
            merged.values(),
            key=lambda parameter: parameter.cli_name or parameter.name,
        )
    )


def parameter_values(
    model_name: str,
    cli_values: dict[str, object],
) -> dict[str, object]:
    info = MODEL_REGISTRY[model_name]
    kwargs = dict(info.defaults)
    accepted = {parameter.name for parameter in info.parameters}
    known_parameters = {parameter.name for parameter in all_cli_parameters()}
    unsupported = sorted(
        name for name in known_parameters - accepted if cli_values.get(name) is not None
    )
    if unsupported:
        options = ", ".join(f"--{name.replace('_', '-')}" for name in unsupported)
        raise ValueError(f"Model {model_name!r} does not accept: {options}.")
    for parameter in info.parameters:
        value = cli_values.get(parameter.name)
        if value is None:
            continue
        if parameter.name == "bonds":
            value = tuple(parse_bond(item) for item in value)
        validate_parameter(parameter, value)
        kwargs[parameter.name] = value
    return kwargs


def explicit_parameter_values(
    model_name: str,
    cli_values: dict[str, object],
) -> dict[str, object]:
    info = MODEL_REGISTRY[model_name]
    accepted = {parameter.name for parameter in info.parameters}
    known = {parameter.name for parameter in all_cli_parameters()}
    unsupported = sorted(name for name in known - accepted if cli_values.get(name) is not None)
    if unsupported:
        options = ", ".join(f"--{name.replace('_', '-')}" for name in unsupported)
        raise ValueError(f"Model {model_name!r} does not accept: {options}.")
    values: dict[str, object] = {}
    for parameter in info.parameters:
        value = cli_values.get(parameter.name)
        if value is None:
            continue
        if parameter.name == "bonds":
            value = tuple(parse_bond(item) for item in value)
        validate_parameter(parameter, value)
        values[parameter.name] = value
    return values


def lattice_spec_from_parameters(
    model_name: str,
    parameters: dict[str, object],
) -> LatticeSpec | None:
    if model_name.removesuffix("_sparse") != "custom_tight_binding":
        return None
    bonds = tuple(parameters.get("bonds", ()))
    n_sites = parameters.get("n_sites")
    if n_sites is None:
        if not bonds:
            raise ValueError("Custom tight-binding specifications require n_sites or bonds.")
        n_sites = max(max(int(bond[0]), int(bond[1])) for bond in bonds) + 1
    lattice_bonds = tuple(
        Bond(
            int(bond[0]),
            int(bond[1]),
            None if len(bond) == 2 else complex(bond[2]),
        )
        for bond in bonds
    )
    return LatticeSpec(n_sites=int(n_sites), bonds=lattice_bonds)


def parse_bond(value: str) -> tuple[object, ...]:
    parts = [part.strip() for part in value.split(",")]
    if len(parts) == 2:
        return int(parts[0]), int(parts[1])
    if len(parts) == 3:
        return int(parts[0]), int(parts[1]), complex(parts[2])
    raise argparse.ArgumentTypeError("Bond must be source,target or source,target,value.")
