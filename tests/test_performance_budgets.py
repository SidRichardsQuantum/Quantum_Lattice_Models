from __future__ import annotations

import json
from pathlib import Path

from quantum_lattice_models import (
    matrix_storage_bytes,
    transverse_field_ising_parity_sector,
    transverse_field_ising_sparse,
)


def test_representative_sparse_storage_and_reduction_budgets() -> None:
    baseline = json.loads(Path("benchmarks/performance_baseline.json").read_text())
    maximum = baseline["maximum"]
    minimum = baseline["minimum"]

    full = transverse_field_ising_sparse(12, periodic=True)
    parity = transverse_field_ising_parity_sector(10, parity=1, periodic=True)

    assert matrix_storage_bytes(full) <= maximum["tfim_sparse_12_storage_bytes"]
    assert matrix_storage_bytes(parity.matrix) <= maximum["parity_sector_10_storage_bytes"]
    assert (2**10) / parity.basis.dimension >= minimum["parity_dimension_reduction"]
