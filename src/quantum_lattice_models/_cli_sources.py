"""CLI source construction and portable-path helpers."""

from __future__ import annotations

import argparse
from pathlib import Path

from quantum_lattice_models._cli_parameters import (
    all_cli_parameters,
    parameter_values,
)
from quantum_lattice_models.registry import MODEL_REGISTRY
from quantum_lattice_models.specs import load_model
from quantum_lattice_models.storage import load_hamiltonian


def build_model(args: argparse.Namespace):
    if args.model is None:
        raise ValueError("Direct model invocation requires --model.")
    info = MODEL_REGISTRY[args.model]
    if info.builder is None:
        raise ValueError(f"Registered model {args.model!r} does not define a builder.")
    built = info.builder(**parameter_values(args.model, vars(args)))
    return getattr(built, "matrix", built)


def build_source(args: argparse.Namespace):
    if args.path is not None:
        if args.model is not None:
            raise ValueError("Use either a model JSON path or --model, not both.")
        supplied_parameters = [
            parameter.name
            for parameter in all_cli_parameters()
            if getattr(args, parameter.name, None) is not None
        ]
        if supplied_parameters:
            raise ValueError("Model parameters cannot be combined with a model JSON path.")
        if is_hamiltonian_path(args.path):
            return load_hamiltonian(args.path).matrix
        return load_model(args.path).hamiltonian()
    return build_model(args)


def is_hamiltonian_path(path: str | Path) -> bool:
    return Path(path).suffix.lower() in {".npy", ".npz"}


def load_portable_source(path: str | Path):
    return load_hamiltonian(path) if is_hamiltonian_path(path) else load_model(path)


def default_export_path(
    path: str | Path,
    artifact: str,
    output_format: str,
) -> Path:
    source = Path(path)
    if artifact == "matrix":
        return source.with_suffix(f".{output_format}")
    if artifact == "bundle":
        return source.with_suffix(".bundle")
    return source.with_suffix(f".{artifact}.json")
