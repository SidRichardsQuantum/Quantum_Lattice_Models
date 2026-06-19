"""Command-line interface for small lattice model workflows."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from quantum_lattice_models.lattice import Bond
from quantum_lattice_models.plotting import plot_spectrum
from quantum_lattice_models.registry import MODEL_REGISTRY, ParameterInfo, validate_parameter
from quantum_lattice_models.specs import (
    LatticeSpec,
    create_model_spec,
    load_model,
)
from quantum_lattice_models.spectra import eigenvalues


def main(argv: list[str] | None = None) -> int:
    """Run the CLI."""

    parser = argparse.ArgumentParser(prog="quantum-lattice")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("models", help="List available registered models")

    create_parser = subparsers.add_parser("create", help="Create a portable model JSON file")
    create_parser.add_argument("model", choices=tuple(sorted(MODEL_REGISTRY)))
    _add_parameter_arguments(create_parser)
    create_parser.add_argument(
        "--representation",
        choices=("dense", "sparse"),
        default=None,
        help="Requested Hamiltonian representation.",
    )
    create_parser.add_argument("--output", required=True, help="Output model JSON path")

    inspect_parser = subparsers.add_parser("inspect", help="Inspect a model JSON file")
    inspect_parser.add_argument("path", help="Model JSON path")

    validate_parser = subparsers.add_parser("validate", help="Validate a model JSON file")
    validate_parser.add_argument("path", help="Model JSON path")

    spectrum_parser = subparsers.add_parser("spectrum", help="Print eigenvalues for a model")
    _add_model_arguments(spectrum_parser)

    plot_parser = subparsers.add_parser("plot", help="Save a spectrum plot for a model")
    _add_model_arguments(plot_parser)
    plot_parser.add_argument("--output", default="spectrum.png", help="Output PNG path")

    args = parser.parse_args(argv)
    if args.command == "models":
        for name, info in sorted(MODEL_REGISTRY.items()):
            print(f"{name}\t{info.category}\t{info.dimension}\t{info.description}")
        return 0

    if args.command == "create":
        parameters = _parameter_values(args.model, vars(args))
        lattice = _lattice_spec_from_parameters(args.model, parameters)
        if lattice is not None:
            parameters.pop("n_sites", None)
            parameters.pop("bonds", None)
        spec = create_model_spec(
            args.model,
            parameters=parameters,
            lattice=lattice,
            representation=args.representation,
        )
        output = spec.save(args.output)
        print(output)
        return 0

    if args.command == "inspect":
        print(json.dumps(load_model(args.path).summary(), indent=2, sort_keys=True, default=str))
        return 0

    if args.command == "validate":
        spec = load_model(args.path)
        spec.validate()
        print(f"valid\t{spec.family}\t{spec.schema_version}")
        return 0

    H = _build_model(args)
    if args.command == "spectrum":
        values = np.real_if_close(eigenvalues(H)).real
        for value in values:
            print(f"{value:.12g}")
        return 0

    if args.command == "plot":
        import matplotlib.pyplot as plt

        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        ax = plot_spectrum(H)
        ax.figure.tight_layout()
        ax.figure.savefig(output, dpi=160)
        plt.close(ax.figure)
        print(output)
        return 0

    parser.error(f"Unhandled command {args.command!r}")
    return 2


def _add_model_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--model",
        required=True,
        choices=tuple(sorted(MODEL_REGISTRY)),
    )
    _add_parameter_arguments(parser)


def _add_parameter_arguments(parser: argparse.ArgumentParser) -> None:
    for parameter in _all_cli_parameters():
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


def _all_cli_parameters() -> tuple[ParameterInfo, ...]:
    merged: dict[str, ParameterInfo] = {}
    for info in MODEL_REGISTRY.values():
        for parameter in info.parameters:
            existing = merged.get(parameter.name)
            if existing is None or parameter.type is complex:
                merged[parameter.name] = parameter
    return tuple(
        sorted(merged.values(), key=lambda parameter: parameter.cli_name or parameter.name)
    )


def _build_model(args: argparse.Namespace):
    info = MODEL_REGISTRY[args.model]
    if info.builder is None:
        raise ValueError(f"Registered model {args.model!r} does not define a builder.")

    kwargs = _parameter_values(args.model, vars(args))
    return info.builder(**kwargs)


def _parameter_values(model_name: str, cli_values: dict[str, object]) -> dict[str, object]:
    info = MODEL_REGISTRY[model_name]
    kwargs = dict(info.defaults)
    accepted = {parameter.name for parameter in info.parameters}
    known_parameters = {parameter.name for parameter in _all_cli_parameters()}
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
            value = tuple(_parse_bond(item) for item in value)
        validate_parameter(parameter, value)
        kwargs[parameter.name] = value
    return kwargs


def _lattice_spec_from_parameters(
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


def _parse_bond(value: str) -> tuple[object, ...]:
    parts = [part.strip() for part in value.split(",")]
    if len(parts) == 2:
        return int(parts[0]), int(parts[1])
    if len(parts) == 3:
        return int(parts[0]), int(parts[1]), complex(parts[2])
    raise argparse.ArgumentTypeError("Bond must be source,target or source,target,value.")


if __name__ == "__main__":
    raise SystemExit(main())
