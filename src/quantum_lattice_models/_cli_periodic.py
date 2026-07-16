"""Periodic-model CLI argument normalization."""

from __future__ import annotations

import argparse


def periodic_parameters(args: argparse.Namespace) -> dict[str, object]:
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
            f"Periodic family {args.family!r} does not accept: " f"{', '.join(unsupported)}."
        )
    return {name: value for name, value in supplied.items() if name in accepted}
