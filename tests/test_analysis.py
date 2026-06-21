from __future__ import annotations

import json

import matplotlib
import numpy as np
import pytest

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from quantum_lattice_models import (
    AnalysisResult,
    chern_number_result,
    create_model_spec,
    load_analysis_result,
    save_analysis_result,
    spectrum_result,
    spin_observables_result,
    ssh_unit_cell,
    standard_momentum_path,
    winding_number_result,
    zak_phase_result,
)
from quantum_lattice_models.cli import main
from quantum_lattice_models.plotting import plot_analysis_result
from quantum_lattice_models.storage import export_hamiltonian_artifact


def test_analysis_result_json_and_npz_round_trip(tmp_path) -> None:
    model = create_model_spec("tight_binding_chain", parameters={"n_sites": 4})
    original = spectrum_result(model.hamiltonian(), model=model)

    json_path = save_analysis_result(original, tmp_path / "spectrum.json")
    npz_path = save_analysis_result(original, tmp_path / "spectrum.npz")
    restored_json = load_analysis_result(json_path)
    restored_npz = load_analysis_result(npz_path)

    assert restored_json.analysis == "spectrum"
    assert restored_json.source["model"]["family"] == "tight_binding_chain"
    assert np.array_equal(restored_json.values["eigenvalues"], original.values["eigenvalues"])
    assert np.array_equal(restored_npz.coordinates["index"], np.arange(4))
    assert restored_npz.plot["kind"] == "spectrum"


def test_band_topology_and_observable_result_producers() -> None:
    periodic = ssh_unit_cell(t1=0.4, t2=1.2)
    bands = periodic.bands(standard_momentum_path(periodic, points_per_segment=6))
    band_result = bands.to_analysis_result(periodic=periodic)
    zak = zak_phase_result(periodic, n_points=101)
    winding = winding_number_result(periodic, n_points=201)

    assert band_result.values["energies"].shape == (11, 2)
    assert band_result.values["eigenvectors"].shape == (11, 2, 2)
    assert band_result.source["kind"] == "periodic_lattice"
    assert abs(abs(float(zak.values["value"])) - np.pi) < 0.01
    assert abs(int(winding.values["value"])) == 1

    state = np.zeros(4, dtype=complex)
    state[0] = 1.0
    observables = spin_observables_result(state, 2)
    assert np.allclose(observables.values["site_magnetization_z"], [1.0, 1.0])
    assert observables.values["correlations"].shape == (2, 2)


def test_chern_result_and_plot_regeneration() -> None:
    from quantum_lattice_models import haldane_unit_cell

    result = chern_number_result(
        haldane_unit_cell(t2=0.2, phi=np.pi / 2),
        mesh=(15, 15),
    )
    assert abs(float(result.values["value"])) == pytest.approx(1.0)

    ax = plot_analysis_result(result)
    assert ax.get_title() == "Chern Number"
    plt.close(ax.figure)

    spectrum = AnalysisResult(
        analysis="spectrum",
        coordinates={"index": np.arange(3)},
        values={"eigenvalues": np.array([-1.0, 0.0, 1.0])},
        plot={"kind": "spectrum", "x": "index", "y": "eigenvalues"},
    )
    ax = plot_analysis_result(spectrum)
    assert ax.get_title() == "Spectrum"
    plt.close(ax.figure)


def test_analysis_results_integrate_with_hamiltonian_bundle(tmp_path) -> None:
    spec = create_model_spec("tight_binding_chain", parameters={"n_sites": 3})
    hamiltonian = spec.build_result()
    analysis = spectrum_result(hamiltonian.matrix, model=spec)

    outputs = export_hamiltonian_artifact(
        hamiltonian,
        tmp_path / "bundle",
        artifact="bundle",
        analyses=(analysis,),
    )

    analysis_path = tmp_path / "bundle" / "analyses" / "000-spectrum.npz"
    assert analysis_path in outputs
    assert load_analysis_result(analysis_path).analysis == "spectrum"
    manifest = json.loads((tmp_path / "bundle" / "manifest.json").read_text())
    assert manifest["analyses"] == ["analyses/000-spectrum.npz"]
    assert "analyses/000-spectrum.npz" in manifest["files"]


def test_analysis_result_cli_workflow(tmp_path, capsys) -> None:
    model = create_model_spec(
        "tight_binding_chain",
        parameters={"n_sites": 3},
    ).save(tmp_path / "model.json")
    result_path = tmp_path / "spectrum.json"

    assert main(["spectrum", str(model), "--result-output", str(result_path)]) == 0
    capsys.readouterr()
    assert result_path.exists()

    assert main(["inspect-result", str(result_path)]) == 0
    summary = json.loads(capsys.readouterr().out)
    assert summary["analysis"] == "spectrum"
    assert summary["value_shapes"] == {"eigenvalues": [3]}

    converted = tmp_path / "spectrum.npz"
    assert (
        main(
            [
                "export-result",
                str(result_path),
                "--format",
                "npz",
                "--output",
                str(converted),
            ]
        )
        == 0
    )
    capsys.readouterr()
    assert load_analysis_result(converted).analysis == "spectrum"

    plot = tmp_path / "spectrum.png"
    assert main(["plot-result", str(converted), "--output", str(plot)]) == 0
    capsys.readouterr()
    assert plot.exists()


def test_cli_bundle_accepts_analysis_records(tmp_path, capsys) -> None:
    model = create_model_spec(
        "tight_binding_chain",
        parameters={"n_sites": 2},
    ).save(tmp_path / "model.json")
    result_model = create_model_spec(
        "tight_binding_chain",
        parameters={"n_sites": 2},
    )
    result = save_analysis_result(
        spectrum_result(result_model.hamiltonian()),
        tmp_path / "result.json",
    )
    bundle = tmp_path / "bundle"

    assert (
        main(
            [
                "export",
                str(model),
                "--artifact",
                "bundle",
                "--analysis",
                str(result),
                "--output",
                str(bundle),
            ]
        )
        == 0
    )
    capsys.readouterr()
    assert (bundle / "analyses" / "000-spectrum.npz").exists()
