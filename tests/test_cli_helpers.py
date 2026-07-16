from __future__ import annotations

import argparse
import json

import numpy as np
import pytest

from quantum_lattice_models import _cli_output, _cli_sources


def test_cli_json_rendering_is_deterministic(capsys: pytest.CaptureFixture[str]) -> None:
    _cli_output.print_json(
        {
            "tuple": (2, 1),
            "scalar": np.int64(3),
            "complex": 1.0 - 2.0j,
        }
    )

    rendered = json.loads(capsys.readouterr().out)
    assert rendered == {
        "complex": {"imag": -2.0, "real": 1.0},
        "scalar": 3,
        "tuple": [2, 1],
    }
    assert _cli_output.json_default(object()).startswith("<object object at ")


def test_cli_key_value_rendering_sorts_and_serializes_nested_values(
    capsys: pytest.CaptureFixture[str],
) -> None:
    _cli_output.print_key_values({"z": (1, 2), "a": {"value": 3}})

    assert capsys.readouterr().out.splitlines() == [
        'a\t{"value": 3}',
        "z\t[1, 2]",
    ]


def test_cli_source_selection_rejects_ambiguous_inputs() -> None:
    with pytest.raises(ValueError, match="requires --model"):
        _cli_sources.build_model(argparse.Namespace(model=None))

    with pytest.raises(ValueError, match="not both"):
        _cli_sources.build_source(argparse.Namespace(path="model.json", model="ssh_model"))

    with pytest.raises(ValueError, match="cannot be combined"):
        _cli_sources.build_source(argparse.Namespace(path="model.json", model=None, n_sites=4))


def test_cli_path_helpers_distinguish_models_matrices_and_artifacts() -> None:
    assert _cli_sources.is_hamiltonian_path("matrix.NPZ")
    assert _cli_sources.is_hamiltonian_path("matrix.npy")
    assert not _cli_sources.is_hamiltonian_path("model.json")
    assert _cli_sources.default_export_path("model.json", "matrix", "npz").name == "model.npz"
    assert _cli_sources.default_export_path("model.json", "bundle", "npz").name == "model.bundle"
    assert (
        _cli_sources.default_export_path("model.json", "metadata", "npz").name
        == "model.metadata.json"
    )
