from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models.models import (
    Lattice,
    aubry_andre_harper_chain,
    bose_hubbard_chain,
    bose_hubbard_chain_sparse,
    custom_tight_binding,
    custom_tight_binding_sparse,
    fermi_hubbard_chain,
    fermi_hubbard_chain_sparse,
    haldane_honeycomb_lattice,
    haldane_honeycomb_lattice_sparse,
    harper_hofstadter_square_lattice,
    harper_hofstadter_square_lattice_sparse,
    heisenberg_chain,
    heisenberg_ladder,
    honeycomb_lattice_index,
    j1_j2_heisenberg_chain,
    kagome_lattice_index,
    kagome_lattice_tight_binding,
    kagome_lattice_tight_binding_sparse,
    kitaev_chain_bdg,
    longitudinal_field_ising,
    next_nearest_neighbor_ising,
    rice_mele_model,
    square_lattice_index,
    square_lattice_tight_binding,
    square_lattice_tight_binding_sparse,
    ssh_edge_state_localizations,
    ssh_model,
    tight_binding_chain,
    tight_binding_chain_sparse,
    transverse_field_ising,
    triangular_lattice_tight_binding,
    triangular_lattice_tight_binding_sparse,
    xxz_chain,
    xy_chain,
)
from quantum_lattice_models.spectra import (
    eigensystem,
    eigenvalues,
    ground_state,
    lowest_eigenvalues,
)


def assert_hermitian(matrix: np.ndarray) -> None:
    assert np.allclose(matrix, matrix.conj().T)


def test_ising_two_site_shape() -> None:
    H = transverse_field_ising(n_sites=2, j=1.0, h=0.5, periodic=False)
    assert H.shape == (4, 4)


def test_spin_hamiltonians_are_hermitian() -> None:
    assert_hermitian(transverse_field_ising(n_sites=4, j=1.0, h=0.5))
    assert_hermitian(longitudinal_field_ising(n_sites=4, j=1.0, h_x=0.5, h_z=0.2))
    assert_hermitian(next_nearest_neighbor_ising(n_sites=4, j1=1.0, j2=0.3, h=0.5))
    assert_hermitian(heisenberg_chain(n_sites=4, jx=1.0, jy=0.7, jz=1.2, field=0.1))
    assert_hermitian(xy_chain(n_sites=4, coupling=1.0, anisotropy=0.3, field=0.2))
    assert_hermitian(xxz_chain(n_sites=4, coupling=1.0, anisotropy=0.5, field=0.1))
    assert_hermitian(j1_j2_heisenberg_chain(n_sites=4, j1=1.0, j2=0.4, field=0.1))
    assert_hermitian(heisenberg_ladder(n_rungs=2, leg_coupling=1.0, rung_coupling=0.7))


def test_new_spin_chain_metadata() -> None:
    xy = xy_chain(n_sites=3, coupling=1.0, anisotropy=0.2, field=0.3)
    xxz = xxz_chain(n_sites=3, coupling=1.0, anisotropy=0.7)
    j1_j2 = j1_j2_heisenberg_chain(n_sites=4, j1=1.0, j2=0.5)
    ladder = heisenberg_ladder(n_rungs=2)
    assert xy.model_name == "xy_chain"
    assert xxz.model_name == "xxz_chain"
    assert j1_j2.model_name == "j1_j2_heisenberg_chain"
    assert ladder.model_name == "heisenberg_ladder"
    assert len(xy.terms) == 2 * (xy.n_sites - 1) + xy.n_sites
    assert len(j1_j2.terms) == 3 * (3 + 2) + 4
    assert ladder.n_sites == 4


def test_tight_binding_shapes_and_hermiticity() -> None:
    chain = tight_binding_chain(n_sites=8, hopping=1.0, onsite=0.25)
    ssh = ssh_model(n_cells=6, t1=0.5, t2=1.0)
    rice_mele = rice_mele_model(n_cells=5)
    square = square_lattice_tight_binding(n_rows=3, n_cols=4, hopping=1.0)
    hofstadter = harper_hofstadter_square_lattice(n_rows=3, n_cols=4, flux=0.2)
    aah = aubry_andre_harper_chain(n_sites=9, hopping=1.0, potential=0.7)
    kitaev = kitaev_chain_bdg(n_sites=5, pairing=0.4)
    haldane = haldane_honeycomb_lattice(n_rows=2, n_cols=3)
    triangular = triangular_lattice_tight_binding(n_rows=3, n_cols=3)
    kagome = kagome_lattice_tight_binding(n_rows=2, n_cols=2)
    assert chain.shape == (8, 8)
    assert ssh.shape == (12, 12)
    assert rice_mele.shape == (10, 10)
    assert square.shape == (12, 12)
    assert hofstadter.shape == (12, 12)
    assert aah.shape == (9, 9)
    assert kitaev.shape == (10, 10)
    assert haldane.shape == (12, 12)
    assert triangular.shape == (9, 9)
    assert kagome.shape == (12, 12)
    assert_hermitian(chain)
    assert_hermitian(ssh)
    assert_hermitian(rice_mele)
    assert_hermitian(square)
    assert_hermitian(hofstadter)
    assert_hermitian(aah)
    assert_hermitian(kitaev)
    assert_hermitian(haldane)
    assert_hermitian(triangular)
    assert_hermitian(kagome)
    assert ssh.model_name == "ssh_model"
    assert square.lattice_shape == (3, 4)
    assert hofstadter.metadata["flux"] == 0.2


def test_public_lattice_indices() -> None:
    assert square_lattice_index(row=2, col=1, n_cols=4) == 9
    assert honeycomb_lattice_index(row=1, col=2, sublattice=1, n_cols=3) == 11
    assert kagome_lattice_index(row=1, col=1, sublattice=2, n_cols=2) == 11


def test_custom_lattice_tight_binding_builder() -> None:
    lattice = Lattice(
        positions=[(0.0, 0.0), (1.0, 0.0), (0.5, 0.8)],
        bonds=[(0, 1), (1, 2, 0.5j)],
        metadata={"label": "triangle fragment"},
    )

    H = custom_tight_binding(lattice=lattice, hopping=2.0, onsite=[0.0, 0.1, 0.2])

    assert H.shape == (3, 3)
    assert H.model_name == "custom_tight_binding"
    assert H.basis == "single_particle"
    assert H.metadata["label"] == "triangle fragment"
    assert np.allclose(np.diag(H), [0.0, 0.1, 0.2])
    assert H[0, 1] == -2.0
    assert H[1, 0] == -2.0
    assert H[1, 2] == 0.5j
    assert H[2, 1] == -0.5j
    assert_hermitian(H)


def test_custom_lattice_infers_site_count_and_sparse_matches_dense() -> None:
    dense = custom_tight_binding(bonds=[(0, 1), (1, 3)], hopping=0.7, onsite=0.2)
    sparse = custom_tight_binding_sparse(bonds=[(0, 1), (1, 3)], hopping=0.7, onsite=0.2)

    assert dense.shape == (4, 4)
    assert sp.issparse(sparse)
    assert np.allclose(sparse.toarray(), dense)


def test_sparse_builders_match_dense_models() -> None:
    chain_dense = tight_binding_chain(n_sites=5, hopping=0.7, onsite=0.2, periodic=True)
    chain_sparse = tight_binding_chain_sparse(n_sites=5, hopping=0.7, onsite=0.2, periodic=True)
    square_dense = square_lattice_tight_binding(
        n_rows=2,
        n_cols=3,
        onsite=np.arange(6) / 10,
        periodic_x=True,
        periodic_y=True,
    )
    square_sparse = square_lattice_tight_binding_sparse(
        n_rows=2,
        n_cols=3,
        onsite=np.arange(6) / 10,
        periodic_x=True,
        periodic_y=True,
    )
    bose_dense = bose_hubbard_chain(n_sites=2, hopping=0.4, max_occupancy=2)
    bose_sparse = bose_hubbard_chain_sparse(n_sites=2, hopping=0.4, max_occupancy=2)
    fermi_dense = fermi_hubbard_chain(n_sites=2, hopping=0.4)
    fermi_sparse = fermi_hubbard_chain_sparse(n_sites=2, hopping=0.4)
    hofstadter_dense = harper_hofstadter_square_lattice(
        n_rows=2, n_cols=3, flux=0.2, periodic_x=True, periodic_y=True
    )
    hofstadter_sparse = harper_hofstadter_square_lattice_sparse(
        n_rows=2, n_cols=3, flux=0.2, periodic_x=True, periodic_y=True
    )
    haldane_dense = haldane_honeycomb_lattice(
        n_rows=2, n_cols=2, phi=0.4, periodic_x=True, periodic_y=True
    )
    haldane_sparse = haldane_honeycomb_lattice_sparse(
        n_rows=2, n_cols=2, phi=0.4, periodic_x=True, periodic_y=True
    )
    triangular_dense = triangular_lattice_tight_binding(
        n_rows=2, n_cols=3, periodic_x=True, periodic_y=True
    )
    triangular_sparse = triangular_lattice_tight_binding_sparse(
        n_rows=2, n_cols=3, periodic_x=True, periodic_y=True
    )
    kagome_dense = kagome_lattice_tight_binding(
        n_rows=2, n_cols=2, periodic_x=True, periodic_y=True
    )
    kagome_sparse = kagome_lattice_tight_binding_sparse(
        n_rows=2, n_cols=2, periodic_x=True, periodic_y=True
    )
    assert sp.issparse(chain_sparse)
    assert sp.issparse(square_sparse)
    assert np.allclose(chain_sparse.toarray(), chain_dense)
    assert np.allclose(square_sparse.toarray(), square_dense)
    assert np.allclose(bose_sparse.toarray(), bose_dense)
    assert np.allclose(fermi_sparse.toarray(), fermi_dense)
    assert np.allclose(hofstadter_sparse.toarray(), hofstadter_dense)
    assert np.allclose(haldane_sparse.toarray(), haldane_dense)
    assert np.allclose(triangular_sparse.toarray(), triangular_dense)
    assert np.allclose(kagome_sparse.toarray(), kagome_dense)
    assert np.allclose(eigenvalues(chain_sparse), eigenvalues(chain_dense))
    assert np.allclose(lowest_eigenvalues(chain_sparse, k=2), eigenvalues(chain_dense)[:2])
    energy, vector = ground_state(chain_sparse)
    assert np.isclose(energy, eigenvalues(chain_dense)[0])
    assert vector.shape == (chain_dense.shape[0],)


def test_hubbard_models_are_hermitian() -> None:
    bose = bose_hubbard_chain(n_sites=3, hopping=0.5, interaction=1.2, max_occupancy=2)
    fermi = fermi_hubbard_chain(n_sites=2, hopping=0.5, interaction=1.0)
    assert bose.shape == (27, 27)
    assert fermi.shape == (16, 16)
    assert_hermitian(bose)
    assert_hermitian(fermi)


def test_square_lattice_periodic_connectivity() -> None:
    H = square_lattice_tight_binding(
        n_rows=2,
        n_cols=3,
        hopping=2.0,
        periodic_x=True,
        periodic_y=False,
    )
    assert H[0, 1] == -2.0
    assert H[0, 2] == -2.0
    assert H[0, 3] == -2.0
    assert H[0, 5] == 0.0


def test_aubry_andre_harper_onsite_potential() -> None:
    H = aubry_andre_harper_chain(
        n_sites=4,
        hopping=0.0,
        potential=2.0,
        beta=0.25,
        phase=0.0,
    )
    assert np.allclose(np.diag(H), [2.0, 0.0, -2.0, 0.0])


def test_rice_mele_staggering_and_hoppings() -> None:
    H = rice_mele_model(n_cells=2, hopping=1.0, dimerization=0.2, staggering=0.7)
    assert np.allclose(np.diag(H), [0.7, -0.7, 0.7, -0.7])
    assert H[0, 1] == -1.2
    assert H[1, 2] == -0.8


def test_harper_hofstadter_vertical_phase() -> None:
    H = harper_hofstadter_square_lattice(n_rows=2, n_cols=2, hopping=1.0, flux=0.25)
    assert H[0, 2] == -1.0
    assert np.allclose(H[1, 3], -1j)


def test_bose_hubbard_single_site_diagonal() -> None:
    H = bose_hubbard_chain(
        n_sites=1,
        hopping=0.0,
        interaction=2.0,
        chemical_potential=0.5,
        max_occupancy=2,
    )
    assert np.allclose(np.diag(H), [0.0, -0.5, 1.0])


def test_bose_hubbard_conserves_total_number() -> None:
    H = bose_hubbard_chain(n_sites=2, hopping=0.5, interaction=1.0, max_occupancy=2)
    local_dim = 3
    for bra in range(H.shape[0]):
        bra_total = bra // local_dim + bra % local_dim
        for ket in range(H.shape[1]):
            if abs(H[bra, ket]) > 1e-12:
                ket_total = ket // local_dim + ket % local_dim
                assert bra_total == ket_total


def test_fermi_hubbard_single_site_diagonal() -> None:
    H = fermi_hubbard_chain(
        n_sites=1,
        hopping=0.0,
        interaction=3.0,
        chemical_potential=0.5,
    )
    assert np.allclose(np.diag(H), [0.0, -0.5, -0.5, 2.0])


def test_fermi_hubbard_hopping_sign_convention() -> None:
    H = fermi_hubbard_chain(n_sites=2, hopping=1.0, interaction=0.0)
    # Moving an up fermion from site 0 orbital 0 to site 1 orbital 2 crosses
    # one occupied down orbital in state |site0_up site0_down>, giving +t.
    source_state = 0b0011
    target_state = 0b0110
    assert H[target_state, source_state] == 1.0


def test_kitaev_bdg_particle_hole_symmetry() -> None:
    H = kitaev_chain_bdg(n_sites=4, hopping=1.0, chemical_potential=0.3, pairing=0.4)
    values = np.sort(np.linalg.eigvalsh(H))
    assert np.allclose(values, -values[::-1])


def test_haldane_complex_phase_changes_next_nearest_hopping() -> None:
    H = haldane_honeycomb_lattice(n_rows=2, n_cols=2, t1=0.0, t2=0.2, phi=np.pi / 2)
    a00 = honeycomb_lattice_index(0, 0, 0, n_cols=2)
    a01 = honeycomb_lattice_index(0, 1, 0, n_cols=2)
    assert np.isclose(H[a00, a01], -0.2j)
    assert np.isclose(H[a01, a00], 0.2j)


def test_eigenvalues_are_real_for_hermitian_models() -> None:
    models = [
        transverse_field_ising(n_sites=3, j=1.0, h=0.5),
        longitudinal_field_ising(n_sites=3),
        next_nearest_neighbor_ising(n_sites=4),
        heisenberg_chain(n_sites=3),
        xy_chain(n_sites=3),
        xxz_chain(n_sites=3),
        j1_j2_heisenberg_chain(n_sites=4),
        heisenberg_ladder(n_rungs=2),
        ssh_model(n_cells=4),
        rice_mele_model(n_cells=4),
        tight_binding_chain(n_sites=5),
        square_lattice_tight_binding(n_rows=2, n_cols=3),
        harper_hofstadter_square_lattice(n_rows=2, n_cols=3, flux=0.2),
        aubry_andre_harper_chain(n_sites=5),
        bose_hubbard_chain(n_sites=2, max_occupancy=2),
        fermi_hubbard_chain(n_sites=2),
        kitaev_chain_bdg(n_sites=3),
        haldane_honeycomb_lattice(n_rows=2, n_cols=2),
        triangular_lattice_tight_binding(n_rows=2, n_cols=3),
        kagome_lattice_tight_binding(n_rows=1, n_cols=2),
    ]
    for H in models:
        values = eigenvalues(H)
        assert np.allclose(values.imag, 0.0)


def test_open_ssh_topological_regime_has_edge_localized_low_energy_states() -> None:
    H = ssh_model(n_cells=12, t1=0.25, t2=1.0, periodic=False)
    values, vectors = eigensystem(H)
    near_zero = np.argsort(np.abs(values))[:2]
    localizations = ssh_edge_state_localizations(vectors[:, near_zero], n_cells=12, edge_cells=2)
    assert np.max(np.abs(values[near_zero])) < 1e-2
    assert np.all(localizations > 0.65)
