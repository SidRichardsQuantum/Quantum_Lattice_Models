from __future__ import annotations

import matplotlib
import numpy as np
import pytest
import scipy.sparse as sp

matplotlib.use("Agg")

import matplotlib.pyplot as plt

from quantum_lattice_models import (
    anderson_chain,
    berry_curvature,
    bond_current,
    bose_hubbard_chain,
    bose_hubbard_chain_sector,
    boson_site_occupations,
    creutz_ladder,
    evolve_state,
    export_analysis_pdf,
    export_analysis_svg,
    export_matrix_plot_data,
    finite_size_sweep,
    fixed_particle_number_basis,
    haldane_unit_cell,
    lieb_lattice,
    local_density_of_states,
    long_range_tight_binding_chain,
    mixed_spin_correlation_matrix,
    parameter_sweep,
    quench_dynamics,
    random_field_heisenberg_chain,
    reciprocal_space_data,
    sawtooth_chain,
    single_particle_occupations,
    solve_eigenpairs,
    thermal_observables,
    tight_binding_chain,
    total_spin_squared,
    two_parameter_sweep,
    wilson_loop,
    xyz_chain,
)
from quantum_lattice_models.plotting import plot_analysis_result
from quantum_lattice_models.spectra import spectral_gap


def test_structured_solver_residuals_and_dense_guard() -> None:
    matrix = sp.diags(np.arange(10, dtype=float), format="csr")
    result = solve_eigenpairs(matrix, k=3)

    assert np.allclose(result.values["eigenvalues"], [0, 1, 2])
    assert result.solver["converged"] is True
    assert result.solver["max_residual"] < 1e-10
    with pytest.raises(ValueError, match="densify"):
        solve_eigenpairs(sp.eye(600), dense_threshold=512)


def test_dynamics_quench_and_loschmidt_records() -> None:
    hamiltonian = np.diag([0.0, 1.0])
    state = np.array([1.0, 1.0]) / np.sqrt(2)
    times = np.linspace(0, np.pi, 7)
    result = evolve_state(
        hamiltonian,
        state,
        times,
        observables={"z": np.diag([1.0, -1.0])},
    )

    assert result.values["states"].shape == (7, 2)
    assert np.allclose(np.linalg.norm(result.values["states"], axis=1), 1)
    assert result.values["loschmidt_echo"][0] == pytest.approx(1.0)
    quench = quench_dynamics(np.diag([-1.0, 1.0]), hamiltonian, times)
    assert quench.analysis == "quench_dynamics"


def test_sweeps_and_thermal_analysis() -> None:
    sweep = parameter_sweep(
        tight_binding_chain,
        "hopping",
        [0.5, 1.0, 1.5],
        spectral_gap,
        builder_parameters={"n_sites": 4},
    )
    grid = two_parameter_sweep(
        tight_binding_chain,
        "hopping",
        [0.5, 1.0],
        "onsite",
        [0.0, 0.2],
        spectral_gap,
        builder_parameters={"n_sites": 3},
    )
    sizes = finite_size_sweep(
        tight_binding_chain,
        [2, 3, 4],
        spectral_gap,
        builder_parameters={"hopping": 1.0},
    )
    thermal = thermal_observables(np.diag([0.0, 1.0]), [0.5, 1.0, 2.0])

    assert sweep.values["measurement"].shape == (3,)
    assert grid.values["measurement"].shape == (2, 2)
    assert sizes.analysis == "finite_size_sweep"
    assert np.all(thermal.values["heat_capacity"] >= 0)


def test_berry_curvature_wilson_loop_and_reciprocal_plots(tmp_path) -> None:
    model = haldane_unit_cell(t2=0.2, phi=np.pi / 2)
    curvature = berry_curvature(model, mesh=(15, 15))
    loops = wilson_loop(model, np.linspace(0, 1, 9), loop_points=31)
    reciprocal = reciprocal_space_data(model)

    assert abs(curvature.metadata["chern_number"]) == pytest.approx(1.0)
    assert loops.values["phases"].shape == (9, 1)
    assert reciprocal.values["reciprocal_vectors"].shape == (2, 2)
    for result in (curvature, loops, reciprocal):
        ax = plot_analysis_result(result)
        plt.close(ax.figure)
    assert export_analysis_svg(curvature, tmp_path / "curvature.svg").exists()
    assert export_analysis_pdf(loops, tmp_path / "wilson.pdf").exists()


def test_fixed_particle_number_bose_hubbard_sector_matches_full_block() -> None:
    sector = bose_hubbard_chain_sector(
        n_sites=3,
        n_particles=2,
        hopping=0.7,
        interaction=1.2,
        max_occupancy=2,
    )
    full = bose_hubbard_chain(
        n_sites=3,
        hopping=0.7,
        interaction=1.2,
        max_occupancy=2,
    )
    embedded_indices = []
    for occupations in sector.basis.states:
        index = 0
        for occupation in occupations:
            index = index * 3 + occupation
        embedded_indices.append(index)
    block = full[np.ix_(embedded_indices, embedded_indices)]

    assert np.allclose(sector.matrix.toarray(), block)
    state = np.ones(sector.basis.dimension) / np.sqrt(sector.basis.dimension)
    assert boson_site_occupations(state, sector.basis).sum() == pytest.approx(2.0)
    assert fixed_particle_number_basis(2, 1).dimension == 2


def test_additional_observables_and_matrix_export(tmp_path) -> None:
    state = np.array([1.0, 1j]) / np.sqrt(2)
    assert np.allclose(single_particle_occupations(state), [0.5, 0.5])
    assert bond_current(state, 0, 1, 1.0) == pytest.approx(1.0)
    ldos = local_density_of_states(np.array([[0.0, -1.0], [-1.0, 0.0]]), [-1, 0, 1])
    assert ldos.shape == (3, 2)

    spin_state = np.zeros(4)
    spin_state[0] = 1
    assert mixed_spin_correlation_matrix(spin_state, 2, "X", "Y").shape == (2, 2)
    assert total_spin_squared(spin_state, 2) == pytest.approx(2.0)
    assert export_matrix_plot_data(np.eye(2), tmp_path / "matrix.json").exists()


def test_benchmark_model_shapes_reproducibility_and_flat_bands() -> None:
    assert np.allclose(anderson_chain(6, seed=2), anderson_chain(6, seed=2))
    assert long_range_tight_binding_chain(5).shape == (5, 5)
    assert creutz_ladder(4).shape == (8, 8)
    assert sawtooth_chain(4).shape == (8, 8)
    lieb_values = np.linalg.eigvalsh(lieb_lattice(2, 2))
    assert np.count_nonzero(np.isclose(lieb_values, 0.0)) >= 4
    assert xyz_chain(3).shape == (8, 8)
    assert np.allclose(
        random_field_heisenberg_chain(3, seed=4),
        random_field_heisenberg_chain(3, seed=4),
    )
