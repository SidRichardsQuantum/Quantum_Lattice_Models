"""Command-line entry point for quantum lattice model workflows."""

from quantum_lattice_models._cli_app import main

__all__ = ["main"]


if __name__ == "__main__":
    raise SystemExit(main())
