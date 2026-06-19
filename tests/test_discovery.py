from __future__ import annotations

import json

import numpy as np
import pytest

from quantum_lattice_models import (
    compare_models,
    create_model_from_preset,
    create_model_spec,
    get_model_info,
    get_preset,
    inspect_model,
    list_models,
    list_presets,
    supports_sparse,
)
from quantum_lattice_models.cli import main


def test_model_filters_and_sparse_capability() -> None:
    spin_models = list_models("spin")
    sparse_spin = list_models("spin", sparse=True)
    validated = list_models(validation_status="validated")

    assert "xxz_chain" in spin_models
    assert "xxz_chain" in sparse_spin
    assert "ssh_model" in validated
    assert "xy_chain" not in validated
    assert supports_sparse("xxz_chain")
    assert not supports_sparse("ssh_model")
    assert all(
        get_model_info(name).basis == "single particle"
        for name in list_models(basis="single particle")
    )


def test_named_presets_create_reproducible_specs() -> None:
    assert "ssh_topological" in list_presets("ssh_model")
    preset = get_preset("ssh_topological")
    spec = create_model_from_preset("ssh_topological", parameters={"n_cells": 5})

    assert preset.model == "ssh_model"
    assert spec.parameters["n_cells"] == 5
    assert spec.parameters["t1"] < spec.parameters["t2"]
    assert spec.provenance["preset"] == "ssh_topological"


def test_dry_run_reports_resources_without_building() -> None:
    report = inspect_model("transverse_field_ising", n_sites=12)
    sparse = inspect_model(
        "transverse_field_ising",
        n_sites=12,
        representation="sparse",
    )

    assert report.dimension == 4096
    assert report.dense_memory_bytes == 4096**2 * np.dtype(np.complex128).itemsize
    assert report.supports_sparse
    assert any("exponentially" in warning for warning in report.warnings)
    assert sparse.representation == "sparse"
    with pytest.raises(ValueError, match="does not support sparse"):
        inspect_model("ssh_model", representation="sparse")


def test_model_comparison_parameters_matrices_spectra_and_gaps() -> None:
    left = create_model_spec(
        "ssh_model",
        parameters={"n_cells": 4, "t1": 0.5, "t2": 1.0},
    )
    right = create_model_spec(
        "ssh_model",
        parameters={"n_cells": 4, "t1": 1.0, "t2": 0.5},
    )
    comparison = compare_models(left, right)

    assert comparison.same_basis
    assert comparison.same_shape
    assert set(comparison.parameter_differences) == {"t1", "t2"}
    assert comparison.matrix_frobenius_difference > 0
    assert comparison.spectrum_maximum_difference > 0
    assert comparison.left_gap is not None
    assert comparison.right_gap is not None


def test_model_comparison_skips_large_matrices() -> None:
    left = create_model_spec("transverse_field_ising", parameters={"n_sites": 12})
    right = create_model_spec("transverse_field_ising", parameters={"n_sites": 12, "h": 1.0})
    comparison = compare_models(left, right, max_dimension=100)

    assert comparison.matrix_frobenius_difference is None
    assert any("skipped" in warning for warning in comparison.warnings)


def test_discovery_cli_json_and_terminal_are_deterministic(tmp_path, capsys) -> None:
    assert main(["models", "--category", "topological", "--json"]) == 0
    models = json.loads(capsys.readouterr().out)
    assert models == sorted(models, key=lambda row: row["name"])
    assert all(row["category"] == "topological" for row in models)

    assert main(["presets", "--model", "ssh_model", "--json"]) == 0
    presets = json.loads(capsys.readouterr().out)
    assert [row["name"] for row in presets] == ["ssh_topological", "ssh_trivial"]

    assert (
        main(
            [
                "dry-run",
                "--model",
                "tight_binding_chain",
                "--n-sites",
                "20",
                "--json",
            ]
        )
        == 0
    )
    report = json.loads(capsys.readouterr().out)
    assert report["dimension"] == 20

    left = create_model_from_preset("ssh_topological").save(tmp_path / "left.json")
    right = create_model_from_preset("ssh_trivial").save(tmp_path / "right.json")
    assert main(["compare", str(left), str(right), "--json"]) == 0
    comparison = json.loads(capsys.readouterr().out)
    assert comparison["parameter_differences"]["t1"]["left"] == 0.5

    assert main(["spectrum", str(left), "--json"]) == 0
    spectrum = json.loads(capsys.readouterr().out)
    assert len(spectrum["eigenvalues"]) == 16


def test_create_from_preset_cli(tmp_path, capsys) -> None:
    output = tmp_path / "preset.json"
    assert (
        main(
            [
                "create",
                "--preset",
                "ssh_topological",
                "--n-cells",
                "3",
                "--output",
                str(output),
            ]
        )
        == 0
    )
    capsys.readouterr()
    data = json.loads(output.read_text())
    assert data["family"] == "ssh_model"
    assert data["parameters"]["n_cells"] == 3
