"""Collect and check small deterministic performance regression metrics."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from time import perf_counter

from quantum_lattice_models import (
    matrix_storage_bytes,
    solve_eigenpairs,
    transverse_field_ising_parity_sector,
    transverse_field_ising_sparse,
)

BASELINE = Path("benchmarks/performance_baseline.json")


def collect_metrics() -> dict[str, float]:
    """Return representative construction, reduction, and solver metrics."""

    start = perf_counter()
    full = transverse_field_ising_sparse(12, j=1.0, h=0.7, periodic=True)
    tfim_seconds = perf_counter() - start

    start = perf_counter()
    parity = transverse_field_ising_parity_sector(
        10,
        parity=1,
        j=1.0,
        h=0.7,
        periodic=True,
    )
    parity_seconds = perf_counter() - start

    solver_matrix = transverse_field_ising_sparse(8, j=1.0, h=0.7)
    start = perf_counter()
    solver = solve_eigenpairs(solver_matrix, k=4)
    solver_seconds = perf_counter() - start

    return {
        "tfim_sparse_12_seconds": tfim_seconds,
        "tfim_sparse_12_storage_bytes": float(matrix_storage_bytes(full)),
        "parity_sector_10_seconds": parity_seconds,
        "parity_sector_10_storage_bytes": float(matrix_storage_bytes(parity.matrix)),
        "parity_dimension_reduction": (2**10) / parity.basis.dimension,
        "sparse_eigenpairs_256_seconds": solver_seconds,
        "sparse_eigenpairs_max_residual": float(solver.values["residuals"].max(initial=0.0)),
    }


def check_metrics(
    metrics: dict[str, float],
    baseline: dict[str, dict[str, float]],
) -> list[str]:
    """Return human-readable threshold failures."""

    failures = []
    for name, maximum in baseline.get("maximum", {}).items():
        value = metrics[name]
        if value > maximum:
            failures.append(f"{name}={value:.6g} exceeds maximum {maximum:.6g}")
    for name, minimum in baseline.get("minimum", {}).items():
        value = metrics[name]
        if value < minimum:
            failures.append(f"{name}={value:.6g} is below minimum {minimum:.6g}")
    if metrics["sparse_eigenpairs_max_residual"] > 1e-8:
        failures.append("sparse eigensolver residual exceeds 1e-8")
    return failures


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--output")
    parser.add_argument("--baseline", default=str(BASELINE))
    args = parser.parse_args(argv)

    metrics = collect_metrics()
    rendered = json.dumps(metrics, indent=2, sort_keys=True) + "\n"
    print(rendered, end="")
    if args.output:
        output = Path(args.output)
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered)
    if args.check:
        baseline = json.loads(Path(args.baseline).read_text())
        failures = check_metrics(metrics, baseline)
        if failures:
            raise SystemExit("\n".join(failures))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
