"""Command-line interface for small lattice model workflows."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from quantum_lattice_models.comparison import compare_models
from quantum_lattice_models.diagnostics import inspect_model
from quantum_lattice_models.interchange import (
    export_graphml,
    export_lattice_csv,
    import_graphml,
    import_lattice_csv,
)
from quantum_lattice_models.lattice import Bond
from quantum_lattice_models.plotting import plot_spectrum
from quantum_lattice_models.registry import (
    MODEL_REGISTRY,
    ParameterInfo,
    get_preset,
    list_models,
    list_presets,
    model_table,
    validate_parameter,
)
from quantum_lattice_models.specs import (
    LatticeSpec,
    create_model_from_preset,
    create_model_spec,
    load_model,
)
from quantum_lattice_models.spectra import eigenvalues
from quantum_lattice_models.storage import save_hamiltonian


def main(argv: list[str] | None = None) -> int:
    """Run the CLI."""

    parser = argparse.ArgumentParser(prog="quantum-lattice")
    subparsers = parser.add_subparsers(dest="command", required=True)

    models_parser = subparsers.add_parser("models", help="List available registered models")
    models_parser.add_argument("--category")
    models_parser.add_argument("--basis")
    models_parser.add_argument("--sparse", action=argparse.BooleanOptionalAction, default=None)
    models_parser.add_argument("--validation-status")
    models_parser.add_argument("--json", action="store_true")

    create_parser = subparsers.add_parser("create", help="Create a portable model JSON file")
    create_parser.add_argument("model", nargs="?", choices=tuple(sorted(MODEL_REGISTRY)))
    create_parser.add_argument("--preset", choices=list_presets())
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
    inspect_parser.add_argument("--json", action="store_true")

    validate_parser = subparsers.add_parser("validate", help="Validate a model JSON file")
    validate_parser.add_argument("path", help="Model JSON path")
    validate_parser.add_argument("--json", action="store_true")

    presets_parser = subparsers.add_parser("presets", help="List named model presets")
    presets_parser.add_argument("--model", choices=tuple(sorted(MODEL_REGISTRY)))
    presets_parser.add_argument("--json", action="store_true")

    dry_run_parser = subparsers.add_parser(
        "dry-run", help="Inspect model resources without constructing a matrix"
    )
    dry_run_parser.add_argument("--model", choices=tuple(sorted(MODEL_REGISTRY)))
    dry_run_parser.add_argument("--preset", choices=list_presets())
    dry_run_parser.add_argument("--representation", choices=("dense", "sparse"), default=None)
    dry_run_parser.add_argument("--json", action="store_true")
    _add_parameter_arguments(dry_run_parser)

    compare_parser = subparsers.add_parser("compare", help="Compare two portable model JSON files")
    compare_parser.add_argument("left")
    compare_parser.add_argument("right")
    compare_parser.add_argument("--max-dimension", type=int, default=2048)
    compare_parser.add_argument("--json", action="store_true")

    import_parser = subparsers.add_parser("import", help="Import a custom lattice model")
    import_parser.add_argument("path", help="Site CSV or GraphML path")
    import_parser.add_argument("--format", choices=("csv", "graphml"), required=True)
    import_parser.add_argument("--bonds", help="Bond CSV path for CSV imports")
    import_parser.add_argument("--metadata", help="Metadata JSON sidecar for CSV imports")
    import_parser.add_argument("--output", required=True, help="Output model JSON path")

    lattice_export_parser = subparsers.add_parser(
        "export-lattice", help="Export lattice data from a model JSON file"
    )
    lattice_export_parser.add_argument("path", help="Portable model JSON path")
    lattice_export_parser.add_argument("--format", choices=("csv", "graphml"), required=True)
    lattice_export_parser.add_argument("--sites", help="Output site CSV path")
    lattice_export_parser.add_argument("--bonds", help="Output bond CSV path")
    lattice_export_parser.add_argument("--metadata", help="Output metadata JSON path")
    lattice_export_parser.add_argument("--output", help="Output GraphML path")

    spectrum_parser = subparsers.add_parser(
        "spectrum", help="Print eigenvalues for a model file or direct model invocation"
    )
    spectrum_parser.add_argument("path", nargs="?", help="Portable model JSON path")
    spectrum_parser.add_argument("--json", action="store_true")
    _add_model_arguments(spectrum_parser, required=False)

    export_parser = subparsers.add_parser(
        "export", help="Build and export a Hamiltonian from a model JSON file"
    )
    export_parser.add_argument("path", help="Portable model JSON path")
    export_parser.add_argument("--format", choices=("npy", "npz"), default="npz")
    export_parser.add_argument("--output", help="Output matrix path")

    plot_parser = subparsers.add_parser("plot", help="Save a spectrum plot for a model")
    _add_model_arguments(plot_parser)
    plot_parser.add_argument("--output", default="spectrum.png", help="Output PNG path")

    args = parser.parse_args(argv)
    if args.command == "models":
        names = list_models(
            args.category,
            basis=args.basis,
            sparse=args.sparse,
            validation_status=args.validation_status,
        )
        rows = [row for row in model_table() if row["name"] in names]
        if args.json:
            _print_json(rows)
        else:
            for row in rows:
                print(
                    f"{row['name']}\t{row['category']}\t{row['basis']}\t"
                    f"{row['supports_sparse']}\t{row['validation_status']}\t"
                    f"{row['description']}"
                )
        return 0

    if args.command == "create":
        if (args.model is None) == (args.preset is None):
            raise ValueError("Create requires exactly one of model or --preset.")
        model_name = args.model if args.model is not None else get_preset(args.preset).model
        if args.preset is not None:
            overrides = _explicit_parameter_values(model_name, vars(args))
            spec = create_model_from_preset(
                args.preset,
                parameters=overrides,
                representation=args.representation,
            )
            output = spec.save(args.output)
            print(output)
            return 0
        parameters = _parameter_values(model_name, vars(args))
        lattice = _lattice_spec_from_parameters(model_name, parameters)
        if lattice is not None:
            parameters.pop("n_sites", None)
            parameters.pop("bonds", None)
        spec = create_model_spec(
            model_name,
            parameters=parameters,
            lattice=lattice,
            representation=args.representation,
        )
        output = spec.save(args.output)
        print(output)
        return 0

    if args.command == "inspect":
        summary = load_model(args.path).summary()
        _print_json(summary)
        return 0

    if args.command == "validate":
        spec = load_model(args.path)
        spec.validate()
        if args.json:
            _print_json(
                {"valid": True, "family": spec.family, "schema_version": spec.schema_version}
            )
        else:
            print(f"valid\t{spec.family}\t{spec.schema_version}")
        return 0

    if args.command == "presets":
        rows = [
            {
                "name": name,
                "model": get_preset(name).model,
                "description": get_preset(name).description,
                "parameters": get_preset(name).parameters,
            }
            for name in list_presets(args.model)
        ]
        if args.json:
            _print_json(rows)
        else:
            for row in rows:
                print(f"{row['name']}\t{row['model']}\t{row['description']}")
        return 0

    if args.command == "dry-run":
        if (args.model is None) == (args.preset is None):
            raise ValueError("Dry-run requires exactly one of --model or --preset.")
        model_name = args.model if args.model is not None else get_preset(args.preset).model
        parameters = dict(get_preset(args.preset).parameters) if args.preset is not None else {}
        parameters.update(_explicit_parameter_values(model_name, vars(args)))
        report = inspect_model(
            model_name,
            representation=args.representation,
            **parameters,
        )
        if args.json:
            _print_json(report.to_dict())
        else:
            _print_key_values(report.to_dict())
        return 0

    if args.command == "compare":
        comparison = compare_models(
            load_model(args.left),
            load_model(args.right),
            max_dimension=args.max_dimension,
        )
        if args.json:
            _print_json(comparison.to_dict())
        else:
            _print_key_values(comparison.to_dict())
        return 0

    if args.command == "import":
        if args.format == "csv":
            if args.bonds is None:
                raise ValueError("CSV imports require --bonds.")
            lattice = import_lattice_csv(
                args.path,
                args.bonds,
                metadata_path=args.metadata,
            )
        else:
            lattice = import_graphml(args.path)
        spec = create_model_spec(
            "custom_tight_binding",
            lattice=lattice,
            parameters={"hopping": 1.0, "onsite": 0.0, "hermitian": True},
        )
        print(spec.save(args.output))
        return 0

    if args.command == "export-lattice":
        spec = load_model(args.path)
        if spec.lattice is None:
            raise ValueError("Model specification does not contain portable lattice data.")
        if args.format == "csv":
            if args.sites is None or args.bonds is None:
                raise ValueError("CSV exports require --sites and --bonds.")
            outputs = export_lattice_csv(
                spec.lattice,
                args.sites,
                args.bonds,
                metadata_path=args.metadata,
            )
            for output in outputs:
                print(output)
        else:
            if args.output is None:
                raise ValueError("GraphML exports require --output.")
            print(export_graphml(spec.lattice, args.output))
        return 0

    if args.command == "spectrum":
        H = _build_source(args)
        values = np.real_if_close(eigenvalues(H)).real
        if args.json:
            _print_json({"eigenvalues": values.tolist()})
        else:
            for value in values:
                print(f"{value:.12g}")
        return 0

    if args.command == "export":
        spec = load_model(args.path)
        output = args.output or str(Path(args.path).with_suffix(f".{args.format}"))
        saved = save_hamiltonian(spec.build_result(), output, format=args.format)
        print(saved)
        return 0

    if args.command == "plot":
        H = _build_model(args)
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


def _add_model_arguments(parser: argparse.ArgumentParser, *, required: bool = True) -> None:
    parser.add_argument(
        "--model",
        required=required,
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
    if args.model is None:
        raise ValueError("Direct model invocation requires --model.")
    info = MODEL_REGISTRY[args.model]
    if info.builder is None:
        raise ValueError(f"Registered model {args.model!r} does not define a builder.")

    kwargs = _parameter_values(args.model, vars(args))
    built = info.builder(**kwargs)
    return getattr(built, "matrix", built)


def _build_source(args: argparse.Namespace):
    if args.path is not None:
        if args.model is not None:
            raise ValueError("Use either a model JSON path or --model, not both.")
        supplied_parameters = [
            parameter.name
            for parameter in _all_cli_parameters()
            if getattr(args, parameter.name, None) is not None
        ]
        if supplied_parameters:
            raise ValueError("Model parameters cannot be combined with a model JSON path.")
        return load_model(args.path).hamiltonian()
    return _build_model(args)


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


def _explicit_parameter_values(model_name: str, cli_values: dict[str, object]) -> dict[str, object]:
    info = MODEL_REGISTRY[model_name]
    accepted = {parameter.name for parameter in info.parameters}
    known = {parameter.name for parameter in _all_cli_parameters()}
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
            value = tuple(_parse_bond(item) for item in value)
        validate_parameter(parameter, value)
        values[parameter.name] = value
    return values


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


def _print_json(value: object) -> None:
    print(json.dumps(value, indent=2, sort_keys=True, default=_json_default))


def _json_default(value: object) -> object:
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, tuple):
        return list(value)
    return str(value)


def _print_key_values(value: dict[str, object]) -> None:
    for key in sorted(value):
        item = value[key]
        if isinstance(item, (dict, list, tuple)):
            item = json.dumps(item, sort_keys=True, default=_json_default)
        print(f"{key}\t{item}")


if __name__ == "__main__":
    raise SystemExit(main())
