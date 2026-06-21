from __future__ import annotations

import json

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from quantum_lattice_models import (
    PeriodicLatticeSpec,
    chern_number,
    create_model_spec,
    export_band_data,
    export_interaction_plot_data,
    export_interaction_svg,
    export_lattice_plot_data,
    export_lattice_svg,
    export_periodic_svg,
    haldane_unit_cell,
    honeycomb_unit_cell,
    load_periodic_lattice,
    repeat_unit_cell,
    rice_mele_unit_cell,
    square_unit_cell,
    ssh_unit_cell,
    standard_momentum_path,
    winding_number,
    zak_phase,
)
from quantum_lattice_models.cli import main
from quantum_lattice_models.plotting import plot_band_structure


def test_periodic_spec_round_trip_and_reciprocal_vectors(tmp_path) -> None:
    original = honeycomb_unit_cell()
    path = original.save(tmp_path / "honeycomb.json")
    restored = load_periodic_lattice(path)

    assert restored == original
    assert np.allclose(
        np.asarray(restored.primitive_vectors) @ restored.reciprocal_vectors().T,
        2 * np.pi * np.eye(2),
    )
    assert PeriodicLatticeSpec.from_dict(json.loads(path.read_text())) == original


def test_ssh_and_square_bloch_analytic_dispersions() -> None:
    ssh = ssh_unit_cell(t1=0.4, t2=1.2)
    for momentum in np.linspace(-np.pi, np.pi, 11):
        values = np.linalg.eigvalsh(ssh.bloch_hamiltonian((momentum,)))
        magnitude = abs(0.4 + 1.2 * np.exp(1j * momentum))
        assert np.allclose(values, (-magnitude, magnitude))

    square = square_unit_cell(hopping=0.7)
    for kx, ky in ((0.0, 0.0), (np.pi, 0.0), (0.3, -1.1)):
        expected = -2 * 0.7 * (np.cos(kx) + np.cos(ky))
        assert np.allclose(square.bloch_hamiltonian((kx, ky)), [[expected]])


def test_band_paths_plotting_and_exports(tmp_path) -> None:
    lattice = honeycomb_unit_cell()
    path = standard_momentum_path(lattice, points_per_segment=8)
    bands = lattice.bands(path)

    assert bands.energies.shape == (22, 2)
    assert bands.eigenvectors is not None
    ax = plot_band_structure(bands, zero_line=True)
    assert ax.get_title() == "Band structure"
    plt.close(ax.figure)

    json_path = export_band_data(bands, tmp_path / "bands.json")
    csv_path = export_band_data(bands, tmp_path / "bands.csv")
    assert json.loads(json_path.read_text())["labels"] == ["Γ", "K", "M", "Γ"]
    assert csv_path.read_text().startswith("point,distance,k0,k1,band_0,band_1")


def test_repeat_unit_cell_and_visual_exports(tmp_path) -> None:
    periodic = honeycomb_unit_cell()
    finite = repeat_unit_cell(periodic, (2, 2))

    assert finite.n_sites == 8
    assert len(finite.unit_cells) == 8
    assert finite.provenance[-1]["operation"] == "repeat_unit_cell"
    assert export_lattice_svg(finite, tmp_path / "finite.svg").read_text().startswith("<?xml")
    assert export_periodic_svg(periodic, tmp_path / "periodic.svg").exists()
    assert export_lattice_plot_data(finite, tmp_path / "plot.json").exists()
    assert export_lattice_plot_data(finite, tmp_path / "plot.csv").exists()

    model = create_model_spec("xxz_chain", parameters={"n_sites": 3})
    assert export_interaction_svg(model, tmp_path / "interactions.svg").exists()
    interaction_data = export_interaction_plot_data(model, tmp_path / "interactions.json")
    assert len(json.loads(interaction_data.read_text())["interactions"]) == 9


def test_zak_winding_and_chern_reference_phases() -> None:
    trivial = ssh_unit_cell(t1=1.2, t2=0.4)
    topological = ssh_unit_cell(t1=0.4, t2=1.2)

    assert winding_number(trivial) == 0
    assert abs(winding_number(topological)) == 1
    assert abs(zak_phase(trivial)) < 1e-6
    assert abs(abs(zak_phase(topological)) - np.pi) < 2e-3

    positive = chern_number(haldane_unit_cell(t2=0.2, phi=np.pi / 2), mesh=(21, 21))
    negative = chern_number(haldane_unit_cell(t2=0.2, phi=-np.pi / 2), mesh=(21, 21))
    assert abs(positive) == pytest.approx(1.0, abs=1e-8)
    assert negative == pytest.approx(-positive, abs=1e-8)


def test_rice_mele_breaks_chiral_winding_assumption() -> None:
    with pytest.raises(ValueError, match="equal diagonal"):
        winding_number(rice_mele_unit_cell(staggered_potential=0.2))


def test_periodic_cli_workflow(tmp_path, capsys) -> None:
    model = tmp_path / "ssh-periodic.json"
    bands = tmp_path / "bands.csv"
    plot = tmp_path / "bands.png"

    assert (
        main(
            [
                "create-periodic",
                "ssh",
                "--t1",
                "0.4",
                "--t2",
                "1.2",
                "--output",
                str(model),
            ]
        )
        == 0
    )
    capsys.readouterr()
    assert (
        main(
            [
                "bands",
                str(model),
                "--points-per-segment",
                "8",
                "--data-output",
                str(bands),
                "--plot-output",
                str(plot),
            ]
        )
        == 0
    )
    capsys.readouterr()
    assert bands.exists() and plot.exists()
    assert main(["topology", str(model), "winding", "--json"]) == 0
    result = json.loads(capsys.readouterr().out)
    assert abs(result["value"]) == 1
