"""Command-line interface for small lattice model workflows."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from quantum_lattice_models.analysis import (
    load_analysis_result,
    save_analysis_result,
    spectrum_result,
    topology_result,
)
from quantum_lattice_models.comparison import compare_models
from quantum_lattice_models.diagnostics import diagnose_matrix, inspect_model
from quantum_lattice_models.intake import (
    adapter_capability_report,
    describe_model,
    lint_model,
)
from quantum_lattice_models.interchange import (
    export_graphml,
    export_lattice_csv,
    import_graphml,
    import_lattice_csv,
)
from quantum_lattice_models.lattice import Bond
from quantum_lattice_models.periodic import (
    haldane_unit_cell,
    honeycomb_unit_cell,
    kagome_unit_cell,
    load_periodic_lattice,
    rice_mele_unit_cell,
    square_unit_cell,
    ssh_unit_cell,
    standard_momentum_path,
)
from quantum_lattice_models.plotting import (
    plot_analysis_result,
    plot_band_structure,
    plot_spectrum,
)
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
from quantum_lattice_models.storage import (
    EXPORT_ARTIFACTS,
    export_hamiltonian_artifact,
    import_hamiltonian,
    load_hamiltonian,
    save_hamiltonian,
)
from quantum_lattice_models.topology import chern_number, winding_number, zak_phase
from quantum_lattice_models.visual_export import (
    export_band_data,
    export_periodic_svg,
)

_PERIODIC_BUILDERS = {
    "ssh": ssh_unit_cell,
    "rice-mele": rice_mele_unit_cell,
    "square": square_unit_cell,
    "honeycomb": honeycomb_unit_cell,
    "kagome": kagome_unit_cell,
    "haldane": haldane_unit_cell,
}


def main(argv: list[str] | None = None) -> int:
    """Run the CLI."""

    parser = argparse.ArgumentParser(prog="quantum-lattice")
    subparsers = parser.add_subparsers(dest="command", required=True)

    models_parser = subparsers.add_parser("models", help="List available registered models")
    models_parser.add_argument("--category")
    models_parser.add_argument("--basis")
    models_parser.add_argument("--sparse", action=argparse.BooleanOptionalAction, default=None)
    models_parser.add_argument("--validation-status")
    models_parser.add_argument(
        "--include-aliases",
        action="store_true",
        help="Include compatibility builder aliases such as names ending in _sparse.",
    )
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

    inspect_parser = subparsers.add_parser(
        "inspect", help="Inspect a model JSON or portable Hamiltonian file"
    )
    inspect_parser.add_argument("path", help="Model JSON or portable NPY/NPZ path")
    inspect_parser.add_argument("--json", action="store_true")

    validate_parser = subparsers.add_parser(
        "validate", help="Validate a model JSON or portable Hamiltonian file"
    )
    validate_parser.add_argument("path", help="Model JSON or portable NPY/NPZ path")
    validate_parser.add_argument("--json", action="store_true")

    describe_parser = subparsers.add_parser(
        "describe", help="Describe the physical content of a portable model or Hamiltonian"
    )
    describe_parser.add_argument("path", help="Model JSON or portable NPY/NPZ path")
    describe_parser.add_argument("--json", action="store_true")

    lint_parser = subparsers.add_parser(
        "lint", help="Report model-intake errors, warnings, and suggested checks"
    )
    lint_parser.add_argument("path", help="Model JSON or portable NPY/NPZ path")
    lint_parser.add_argument("--json", action="store_true")

    capability_parser = subparsers.add_parser(
        "adapter-capabilities",
        help="Report semantics preserved or lost by an export adapter",
    )
    capability_parser.add_argument("path", help="Model JSON or portable NPY/NPZ path")
    capability_parser.add_argument(
        "target",
        choices=(
            "ase",
            "csv",
            "dot",
            "graphml",
            "json",
            "netket",
            "networkx",
            "openfermion",
            "pennylane",
            "plot-data-json",
            "qiskit",
            "quspin",
            "qutip",
            "svg",
            "xyz",
            "yaml",
        ),
    )
    capability_parser.add_argument("--json", action="store_true")

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

    matrix_import_parser = subparsers.add_parser(
        "import-matrix", help="Import an external Hamiltonian matrix with metadata"
    )
    matrix_import_parser.add_argument("path", help="External NPY or NPZ matrix path")
    matrix_import_parser.add_argument(
        "--metadata", required=True, help="External matrix metadata JSON path"
    )
    matrix_import_parser.add_argument(
        "--allow-non-hermitian",
        action="store_true",
        help="Allow importing a non-Hermitian square matrix.",
    )
    matrix_import_parser.add_argument(
        "--format", choices=("npy", "npz"), default="npz", help="Portable output format"
    )
    matrix_import_parser.add_argument("--output", required=True, help="Portable matrix output path")

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
    spectrum_parser.add_argument("--result-output")
    _add_model_arguments(spectrum_parser, required=False)

    export_parser = subparsers.add_parser(
        "export", help="Export artifacts from a model JSON or portable matrix file"
    )
    export_parser.add_argument("path", help="Portable model JSON or NPY/NPZ path")
    export_parser.add_argument(
        "--artifact",
        choices=EXPORT_ARTIFACTS,
        default="matrix",
        help="Artifact to export; the default preserves matrix export behavior.",
    )
    export_parser.add_argument("--format", choices=("npy", "npz"), default="npz")
    export_parser.add_argument("--output", help="Output file or bundle directory")
    export_parser.add_argument(
        "--analysis",
        action="append",
        default=[],
        help="Analysis-result JSON/NPZ path to include in a bundle.",
    )

    plot_parser = subparsers.add_parser("plot", help="Save a spectrum plot for a model")
    _add_model_arguments(plot_parser)
    plot_parser.add_argument("--output", default="spectrum.png", help="Output PNG path")

    periodic_parser = subparsers.add_parser(
        "create-periodic", help="Create a portable periodic unit-cell JSON file"
    )
    periodic_parser.add_argument("family", choices=tuple(_PERIODIC_BUILDERS))
    periodic_parser.add_argument("--t1", type=float, default=None)
    periodic_parser.add_argument("--t2", type=float, default=None)
    periodic_parser.add_argument("--hopping", type=float, default=None)
    periodic_parser.add_argument("--phi", type=float, default=None)
    periodic_parser.add_argument("--sublattice-potential", type=float, default=None)
    periodic_parser.add_argument("--output", required=True)

    bands_parser = subparsers.add_parser(
        "bands", help="Calculate bands from a periodic unit-cell JSON file"
    )
    bands_parser.add_argument("path")
    bands_parser.add_argument("--points-per-segment", type=int, default=64)
    bands_parser.add_argument("--data-output")
    bands_parser.add_argument("--plot-output")
    bands_parser.add_argument("--result-output")
    bands_parser.add_argument("--json", action="store_true")

    topology_parser = subparsers.add_parser(
        "topology", help="Calculate a topological invariant for a periodic model"
    )
    topology_parser.add_argument("path")
    topology_parser.add_argument("invariant", choices=("zak", "winding", "chern"))
    topology_parser.add_argument("--band", type=int, default=0)
    topology_parser.add_argument("--n-points", type=int, default=401)
    topology_parser.add_argument("--mesh", type=int, default=31)
    topology_parser.add_argument("--occupied-bands", type=int, default=1)
    topology_parser.add_argument("--result-output")
    topology_parser.add_argument("--json", action="store_true")

    periodic_svg_parser = subparsers.add_parser(
        "export-periodic-svg", help="Export a repeated periodic unit-cell diagram"
    )
    periodic_svg_parser.add_argument("path")
    periodic_svg_parser.add_argument("--repeats-x", type=int, default=3)
    periodic_svg_parser.add_argument("--repeats-y", type=int, default=3)
    periodic_svg_parser.add_argument("--output", required=True)

    result_inspect_parser = subparsers.add_parser(
        "inspect-result", help="Inspect a portable analysis-result file"
    )
    result_inspect_parser.add_argument("path")

    result_export_parser = subparsers.add_parser(
        "export-result", help="Convert a portable analysis result between JSON and NPZ"
    )
    result_export_parser.add_argument("path")
    result_export_parser.add_argument("--format", choices=("json", "npz"), required=True)
    result_export_parser.add_argument("--output", required=True)

    result_plot_parser = subparsers.add_parser(
        "plot-result", help="Regenerate a plot from a portable analysis result"
    )
    result_plot_parser.add_argument("path")
    result_plot_parser.add_argument("--output", required=True)

    args = parser.parse_args(argv)
    if args.command == "models":
        names = list_models(
            args.category,
            basis=args.basis,
            sparse=args.sparse,
            validation_status=args.validation_status,
            include_aliases=args.include_aliases,
        )
        rows = [
            row for row in model_table(include_aliases=args.include_aliases) if row["name"] in names
        ]
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

    if args.command in {"describe", "lint", "adapter-capabilities"}:
        source = _load_portable_source(args.path)
        if args.command == "describe":
            report = describe_model(source).to_dict()
        elif args.command == "lint":
            report = lint_model(source).to_dict()
        else:
            report = adapter_capability_report(source, args.target).to_dict()
        if args.json:
            _print_json(report)
        else:
            _print_key_values(report)
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
        if _is_hamiltonian_path(args.path):
            result = load_hamiltonian(args.path)
            diagnostics = diagnose_matrix(result.matrix)
            summary = {
                **result.model.summary(),
                "matrix": {
                    "shape": list(diagnostics.shape),
                    "sparse": diagnostics.sparse,
                    "nonzero_entries": diagnostics.nonzero_entries,
                    "density": diagnostics.density,
                    "storage_bytes": diagnostics.storage_bytes,
                    "hermitian": diagnostics.hermitian,
                },
            }
        else:
            summary = load_model(args.path).summary()
        _print_json(summary)
        return 0

    if args.command == "validate":
        if _is_hamiltonian_path(args.path):
            result = load_hamiltonian(args.path)
            result.model.validate()
            family = result.model.family
            schema_version = result.model.schema_version
        else:
            spec = load_model(args.path)
            spec.validate()
            family = spec.family
            schema_version = spec.schema_version
        if args.json:
            _print_json({"valid": True, "family": family, "schema_version": schema_version})
        else:
            print(f"valid\t{family}\t{schema_version}")
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

    if args.command == "import-matrix":
        result = import_hamiltonian(
            args.path,
            metadata_path=args.metadata,
            require_hermitian=not args.allow_non_hermitian,
        )
        print(save_hamiltonian(result, args.output, format=args.format))
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
        if args.result_output:
            model = None
            if args.path is not None and not _is_hamiltonian_path(args.path):
                model = load_model(args.path)
            elif args.path is not None:
                model = load_hamiltonian(args.path).model
            result = spectrum_result(H, model=model)
            print(save_analysis_result(result, args.result_output))
        if args.json:
            _print_json({"eigenvalues": values.tolist()})
        else:
            for value in values:
                print(f"{value:.12g}")
        return 0

    if args.command == "export":
        if _is_hamiltonian_path(args.path):
            source = load_hamiltonian(args.path)
        else:
            spec = load_model(args.path)
            source = spec if args.artifact in {"model", "lattice"} else spec.build_result()
        output = args.output or str(_default_export_path(args.path, args.artifact, args.format))
        outputs = export_hamiltonian_artifact(
            source,
            output,
            artifact=args.artifact,
            format=args.format,
            analyses=tuple(load_analysis_result(path) for path in args.analysis),
        )
        for exported in outputs:
            print(exported)
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

    if args.command == "create-periodic":
        kwargs = _periodic_parameters(args)
        periodic = _PERIODIC_BUILDERS[args.family](**kwargs)
        print(periodic.save(args.output))
        return 0

    if args.command == "bands":
        periodic = load_periodic_lattice(args.path)
        path = standard_momentum_path(
            periodic,
            points_per_segment=args.points_per_segment,
        )
        bands = periodic.bands(path, eigenvectors=False)
        analysis = bands.to_analysis_result(
            periodic=periodic,
            parameters={"points_per_segment": args.points_per_segment},
        )
        if args.result_output:
            print(save_analysis_result(analysis, args.result_output))
        if args.data_output:
            print(export_band_data(bands, args.data_output))
        if args.plot_output:
            import matplotlib.pyplot as plt

            output = Path(args.plot_output)
            output.parent.mkdir(parents=True, exist_ok=True)
            ax = plot_band_structure(bands)
            ax.figure.tight_layout()
            ax.figure.savefig(output, dpi=160)
            plt.close(ax.figure)
            print(output)
        if args.json:
            _print_json(bands.to_dict())
        elif not args.data_output and not args.plot_output:
            for distance, energies in zip(bands.distances, bands.energies, strict=True):
                print("\t".join([f"{distance:.12g}", *(f"{value:.12g}" for value in energies)]))
        return 0

    if args.command == "topology":
        periodic = load_periodic_lattice(args.path)
        if args.invariant == "zak":
            value = zak_phase(periodic, band=args.band, n_points=args.n_points)
            parameters = {"band": args.band, "n_points": args.n_points}
            solver = {"method": "discrete Wilson loop", "gauge_tolerant": True}
        elif args.invariant == "winding":
            value = winding_number(periodic, n_points=args.n_points)
            parameters = {"n_points": args.n_points}
            solver = {"method": "unwrapped off-diagonal phase", "gauge_tolerant": True}
        else:
            value = chern_number(
                periodic,
                occupied_bands=args.occupied_bands,
                mesh=(args.mesh, args.mesh),
            )
            parameters = {
                "occupied_bands": args.occupied_bands,
                "mesh": [args.mesh, args.mesh],
            }
            solver = {"method": "Fukui-Hatsugai-Suzuki", "gauge_tolerant": True}
        if args.result_output:
            result = topology_result(
                args.invariant,
                value,
                periodic=periodic,
                parameters=parameters,
                solver=solver,
            )
            print(save_analysis_result(result, args.result_output))
        if args.json:
            _print_json({"invariant": args.invariant, "value": value})
        else:
            print(f"{value:.12g}")
        return 0

    if args.command == "export-periodic-svg":
        periodic = load_periodic_lattice(args.path)
        print(
            export_periodic_svg(
                periodic,
                args.output,
                repeats=(args.repeats_x, args.repeats_y),
            )
        )
        return 0

    if args.command == "inspect-result":
        _print_json(load_analysis_result(args.path).summary())
        return 0

    if args.command == "export-result":
        result = load_analysis_result(args.path)
        print(save_analysis_result(result, args.output, format=args.format))
        return 0

    if args.command == "plot-result":
        import matplotlib.pyplot as plt

        result = load_analysis_result(args.path)
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        ax = plot_analysis_result(result)
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
        if _is_hamiltonian_path(args.path):
            return load_hamiltonian(args.path).matrix
        return load_model(args.path).hamiltonian()
    return _build_model(args)


def _is_hamiltonian_path(path: str | Path) -> bool:
    return Path(path).suffix.lower() in {".npy", ".npz"}


def _load_portable_source(path: str | Path):
    return load_hamiltonian(path) if _is_hamiltonian_path(path) else load_model(path)


def _default_export_path(path: str | Path, artifact: str, output_format: str) -> Path:
    source = Path(path)
    if artifact == "matrix":
        return source.with_suffix(f".{output_format}")
    if artifact == "bundle":
        return source.with_suffix(".bundle")
    return source.with_suffix(f".{artifact}.json")


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


def _periodic_parameters(args: argparse.Namespace) -> dict[str, object]:
    supplied = {
        name: getattr(args, name)
        for name in ("t1", "t2", "hopping", "phi", "sublattice_potential")
        if getattr(args, name) is not None
    }
    accepted = {
        "ssh": {"t1", "t2"},
        "rice-mele": {"t1", "t2", "sublattice_potential"},
        "square": {"hopping"},
        "honeycomb": {"hopping"},
        "kagome": {"hopping"},
        "haldane": {"t1", "t2", "phi", "sublattice_potential"},
    }[args.family]
    unsupported = sorted(set(supplied) - accepted)
    if unsupported:
        raise ValueError(
            f"Periodic family {args.family!r} does not accept: {', '.join(unsupported)}."
        )
    return {name: value for name, value in supplied.items() if name in accepted}


if __name__ == "__main__":
    raise SystemExit(main())
