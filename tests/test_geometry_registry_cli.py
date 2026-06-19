from __future__ import annotations

import subprocess
import sys

import numpy as np
import pytest

from quantum_lattice_models.cli import main
from quantum_lattice_models.geometry import (
    honeycomb_lattice_positions,
    kagome_lattice_positions,
    square_lattice_positions,
    triangular_lattice_positions,
)
from quantum_lattice_models.registry import (
    ParameterInfo,
    get_model_info,
    list_models,
    model_table,
    register_model,
    supports_sparse,
    unregister_model,
)


def test_lattice_position_shapes() -> None:
    assert square_lattice_positions(2, 3).shape == (6, 2)
    assert triangular_lattice_positions(2, 3).shape == (6, 2)
    assert honeycomb_lattice_positions(2, 3).shape == (12, 2)
    assert kagome_lattice_positions(2, 3).shape == (18, 2)
    assert np.allclose(square_lattice_positions(1, 2), [[0, 0], [1, 0]])


def test_model_registry_helpers() -> None:
    names = list_models()
    assert "ssh_model" in names
    assert "custom_tight_binding" in names
    assert "xxz_chain" in names
    assert "kitaev_chain_bdg" in names
    assert "haldane_honeycomb_lattice_sparse" in names
    assert get_model_info("ssh_model").basis == "single particle"
    assert get_model_info("ssh_model").builder is not None
    assert get_model_info("ssh_model").defaults["n_cells"] == 8
    n_cells = next(
        parameter
        for parameter in get_model_info("ssh_model").parameters
        if parameter.name == "n_cells"
    )
    assert isinstance(n_cells, ParameterInfo)
    assert n_cells.type is int
    assert n_cells.minimum == 1
    assert n_cells.cli_name == "--n-cells"
    assert all(get_model_info(name).name == name for name in names)
    assert any(row["name"] == "bose_hubbard_chain" for row in model_table())
    assert supports_sparse("tight_binding_chain")
    assert get_model_info("ssh_model").validation_status == "validated"


def test_register_and_unregister_model() -> None:
    def builder(n_sites: int = 2) -> np.ndarray:
        return np.zeros((n_sites, n_sites))

    info = register_model(
        "test_custom_registry_model",
        category="test",
        basis="single particle",
        dimension="n_sites",
        return_type="ndarray",
        description="Temporary test model",
        builder=builder,
        defaults={"n_sites": 2},
    )

    try:
        assert info.name in list_models("test")
        assert get_model_info(info.name).builder is builder
        try:
            register_model(
                info.name,
                category="test",
                basis="single particle",
                dimension="n_sites",
                return_type="ndarray",
                description="duplicate",
                builder=builder,
            )
        except ValueError:
            pass
        else:
            raise AssertionError("duplicate model registration should fail")
    finally:
        removed = unregister_model(info.name)

    assert removed.name == info.name
    assert info.name not in list_models()
    assert info.validation_status == "unvalidated"


def test_cli_models_and_spectrum(capsys: pytest.CaptureFixture[str]) -> None:
    assert main(["models"]) == 0
    assert "ssh_model" in capsys.readouterr().out

    assert main(["spectrum", "--model", "tight_binding_chain", "--n-sites", "3"]) == 0
    assert len(capsys.readouterr().out.strip().splitlines()) == 3

    assert main(["spectrum", "--model", "xxz_chain", "--n-sites", "3"]) == 0
    assert len(capsys.readouterr().out.strip().splitlines()) == 8

    assert (
        main(
            [
                "spectrum",
                "--model",
                "kitaev_chain_bdg",
                "--n-sites",
                "3",
                "--pairing",
                "0.25j",
            ]
        )
        == 0
    )
    assert len(capsys.readouterr().out.strip().splitlines()) == 6


def test_cli_custom_tight_binding_bonds(capsys: pytest.CaptureFixture[str]) -> None:
    assert (
        main(
            [
                "spectrum",
                "--model",
                "custom_tight_binding",
                "--n-sites",
                "3",
                "--bond",
                "0,1",
                "--bond",
                "1,2",
            ]
        )
        == 0
    )
    assert len(capsys.readouterr().out.strip().splitlines()) == 3


def test_cli_rejects_parameters_not_in_model_schema() -> None:
    with pytest.raises(ValueError, match="does not accept: --flux"):
        main(
            [
                "spectrum",
                "--model",
                "ssh_model",
                "--n-cells",
                "3",
                "--flux",
                "0.2",
            ]
        )


def test_cli_module_entrypoint_smoke() -> None:
    result = subprocess.run(
        [sys.executable, "-m", "quantum_lattice_models.cli", "models"],
        check=True,
        capture_output=True,
        text=True,
    )
    assert "ssh_model" in result.stdout
