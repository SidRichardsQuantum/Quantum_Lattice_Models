"""Command-line interface for small lattice model workflows."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

from quantum_lattice_models.models import (
    bose_hubbard_chain,
    haldane_honeycomb_lattice,
    harper_hofstadter_square_lattice,
    heisenberg_chain,
    rice_mele_model,
    ssh_model,
    tight_binding_chain,
    transverse_field_ising,
)
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
        choices=(
            "bose_hubbard_chain",
            "haldane_honeycomb_lattice",
            "harper_hofstadter_square_lattice",
            "heisenberg_chain",
            "rice_mele_model",
            "ssh_model",
            "tight_binding_chain",
            "transverse_field_ising",
        ),
    )
    parser.add_argument("--n-sites", type=int, default=4)
    parser.add_argument("--n-cells", type=int, default=6)
    parser.add_argument("--n-rows", type=int, default=3)
    parser.add_argument("--n-cols", type=int, default=3)
    parser.add_argument("--hopping", type=float, default=1.0)
    parser.add_argument("--onsite", type=float, default=0.0)
    parser.add_argument("--j", type=float, default=1.0)
    parser.add_argument("--h", type=float, default=0.5)
    parser.add_argument("--field", type=float, default=0.0)
    parser.add_argument("--t1", type=float, default=0.5)
    parser.add_argument("--t2", type=float, default=1.0)
    parser.add_argument("--dimerization", type=float, default=0.25)
    parser.add_argument("--staggering", type=float, default=0.5)
    parser.add_argument("--flux", type=float, default=0.25)
    parser.add_argument("--interaction", type=float, default=1.0)
    parser.add_argument("--max-occupancy", type=int, default=2)
    parser.add_argument("--periodic", action="store_true")


def _build_model(args: argparse.Namespace):
    if args.model == "transverse_field_ising":
        return transverse_field_ising(args.n_sites, j=args.j, h=args.h, periodic=args.periodic)
    if args.model == "heisenberg_chain":
        return heisenberg_chain(args.n_sites, field=args.field, periodic=args.periodic)
    if args.model == "ssh_model":
        return ssh_model(args.n_cells, t1=args.t1, t2=args.t2, periodic=args.periodic)
    if args.model == "rice_mele_model":
        return rice_mele_model(
            args.n_cells,
            hopping=args.hopping,
            dimerization=args.dimerization,
            staggering=args.staggering,
            periodic=args.periodic,
        )
    if args.model == "tight_binding_chain":
        return tight_binding_chain(
            args.n_sites, hopping=args.hopping, onsite=args.onsite, periodic=args.periodic
        )
    if args.model == "harper_hofstadter_square_lattice":
        return harper_hofstadter_square_lattice(
            args.n_rows, args.n_cols, hopping=args.hopping, flux=args.flux
        )
    if args.model == "haldane_honeycomb_lattice":
        return haldane_honeycomb_lattice(args.n_rows, args.n_cols)
    if args.model == "bose_hubbard_chain":
        return bose_hubbard_chain(
            args.n_sites,
            hopping=args.hopping,
            interaction=args.interaction,
            max_occupancy=args.max_occupancy,
            periodic=args.periodic,
        )
    raise ValueError(f"Unsupported model {args.model!r}.")


if __name__ == "__main__":
    raise SystemExit(main())
