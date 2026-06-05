from __future__ import annotations

import subprocess
import sys

import numpy as np

from quantum_lattice_models.geometry import (
    honeycomb_lattice_positions,
    kagome_lattice_positions,
    square_lattice_positions,
    triangular_lattice_positions,
)
from quantum_lattice_models.registry import get_model_info, list_models, model_table


def test_lattice_position_shapes() -> None:
    assert square_lattice_positions(2, 3).shape == (6, 2)
    assert triangular_lattice_positions(2, 3).shape == (6, 2)
    assert honeycomb_lattice_positions(2, 3).shape == (12, 2)
    assert kagome_lattice_positions(2, 3).shape == (18, 2)
    assert np.allclose(square_lattice_positions(1, 2), [[0, 0], [1, 0]])


def test_model_registry_helpers() -> None:
    names = list_models()
    assert "ssh_model" in names
    assert "xxz_chain" in names
    assert "kitaev_chain_bdg" in names
    assert "haldane_honeycomb_lattice_sparse" in names
    assert get_model_info("ssh_model").basis == "single particle"
    assert get_model_info("ssh_model").builder is not None
    assert get_model_info("ssh_model").defaults["n_cells"] == 8
    assert any(row["name"] == "bose_hubbard_chain" for row in model_table())


def test_cli_models_and_spectrum() -> None:
    models = subprocess.run(
        [sys.executable, "-m", "quantum_lattice_models.cli", "models"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "ssh_model" in models.stdout

    spectrum = subprocess.run(
        [
            sys.executable,
            "-m",
            "quantum_lattice_models.cli",
            "spectrum",
            "--model",
            "tight_binding_chain",
            "--n-sites",
            "3",
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    assert len(spectrum.stdout.strip().splitlines()) == 3
