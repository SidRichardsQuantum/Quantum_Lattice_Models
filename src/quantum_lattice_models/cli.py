"""Command-line interface for small lattice model workflows."""

from __future__ import annotations

import argparse
import inspect
from pathlib import Path

import numpy as np

from quantum_lattice_models.plotting import plot_lattice_spectrum
from quantum_lattice_models.registry import MODEL_REGISTRY
from quantum_lattice_models.spectra import eigenvalues


def main(argv: list[str] | None = None) -> int:
    """Run the CLI."""

    parser = argparse.ArgumentParser(prog="quantum-lattice")
    subparsers = parser.add_subparsers(dest="command", required=True)

    subparsers.add_parser("models", help="List available registered models")

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
        ax = plot_lattice_spectrum(H)
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
    parser.add_argument("--n-sites", type=int, default=None)
    parser.add_argument("--n-cells", type=int, default=None)
    parser.add_argument("--n-rows", type=int, default=None)
    parser.add_argument("--n-cols", type=int, default=None)
    parser.add_argument("--n-rungs", type=int, default=None)
    parser.add_argument("--hopping", type=complex, default=None)
    parser.add_argument("--onsite", type=float, default=None)
    parser.add_argument("--j", type=float, default=None)
    parser.add_argument("--jx", type=float, default=None)
    parser.add_argument("--jy", type=float, default=None)
    parser.add_argument("--jz", type=float, default=None)
    parser.add_argument("--j1", type=float, default=None)
    parser.add_argument("--j2", type=float, default=None)
    parser.add_argument("--h", type=float, default=None)
    parser.add_argument("--h-x", type=float, default=None)
    parser.add_argument("--h-z", type=float, default=None)
    parser.add_argument("--field", type=float, default=None)
    parser.add_argument("--coupling", type=float, default=None)
    parser.add_argument("--anisotropy", type=float, default=None)
    parser.add_argument("--leg-coupling", type=float, default=None)
    parser.add_argument("--rung-coupling", type=float, default=None)
    parser.add_argument("--t1", type=float, default=None)
    parser.add_argument("--t2", type=float, default=None)
    parser.add_argument("--dimerization", type=float, default=None)
    parser.add_argument("--staggering", type=float, default=None)
    parser.add_argument("--flux", type=float, default=None)
    parser.add_argument("--potential", type=float, default=None)
    parser.add_argument("--beta", type=float, default=None)
    parser.add_argument("--phase", type=float, default=None)
    parser.add_argument("--chemical-potential", type=float, default=None)
    parser.add_argument("--interaction", type=float, default=None)
    parser.add_argument("--max-occupancy", type=int, default=None)
    parser.add_argument("--periodic", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--periodic-x", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--periodic-y", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument("--hermitian", action=argparse.BooleanOptionalAction, default=None)
    parser.add_argument(
        "--bond",
        action="append",
        default=None,
        help="Custom model bond as source,target or source,target,value",
    )


def _build_model(args: argparse.Namespace):
    info = MODEL_REGISTRY[args.model]
    if info.builder is None:
        raise ValueError(f"Registered model {args.model!r} does not define a builder.")

    kwargs = dict(info.defaults)
    cli_values = vars(args)
    signature = inspect.signature(info.builder)
    accepted = set(signature.parameters)
    for name in accepted:
        arg_name = name.replace("-", "_")
        if arg_name in cli_values and cli_values[arg_name] is not None:
            kwargs[name] = cli_values[arg_name]
    if "bonds" in accepted and args.bond is not None:
        kwargs["bonds"] = tuple(_parse_bond(value) for value in args.bond)
    return info.builder(**kwargs)


def _parse_bond(value: str) -> tuple[object, ...]:
    parts = [part.strip() for part in value.split(",")]
    if len(parts) == 2:
        return int(parts[0]), int(parts[1])
    if len(parts) == 3:
        return int(parts[0]), int(parts[1]), complex(parts[2])
    raise argparse.ArgumentTypeError("Bond must be source,target or source,target,value.")


if __name__ == "__main__":
    raise SystemExit(main())
