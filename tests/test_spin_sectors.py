from __future__ import annotations

from math import comb

import numpy as np
import pytest
import scipy.sparse as sp

from quantum_lattice_models import (
    create_model_spec,
    estimate_model_dimension,
    fixed_magnetization_basis,
    heisenberg_chain_sector,
    heisenberg_chain_sparse,
    load_hamiltonian,
    save_hamiltonian,
    xxz_chain_sector,
    xxz_chain_sparse,
)
from quantum_lattice_models.cli import main


def test_fixed_magnetization_basis_mapping_and_state_conversion() -> None:
    basis = fixed_magnetization_basis(n_sites=4, magnetization=0)

    assert basis.dimension == comb(4, 2)
    assert all(state.bit_count() == 2 for state in basis.states)
    assert basis.state_to_index[basis.states[2]] == 2

    reduced = np.arange(basis.dimension, dtype=complex)
    embedded = basis.embed(reduced)
    assert embedded.shape == (16,)
    assert np.allclose(basis.project(embedded), reduced)
    assert np.count_nonzero(embedded) == basis.dimension - 1


@pytest.mark.parametrize(
    ("n_sites", "magnetization", "periodic"),
    [(4, 0, False), (5, 1, False), (4, -2, True)],
)
def test_xxz_sector_matches_full_space_block(n_sites, magnetization, periodic) -> None:
    sector = xxz_chain_sector(
        n_sites=n_sites,
        magnetization=magnetization,
        coupling=0.7,
        anisotropy=1.3,
        field=0.2,
        periodic=periodic,
    )
    full = xxz_chain_sparse(
        n_sites=n_sites,
        coupling=0.7,
        anisotropy=1.3,
        field=0.2,
        periodic=periodic,
    )
    states = np.asarray(sector.basis.states)
    block = full[states][:, states]

    assert sp.issparse(sector.matrix)
    assert np.allclose(sector.matrix.toarray(), block.toarray())
    assert np.allclose(
        np.linalg.eigvalsh(sector.matrix.toarray()),
        np.linalg.eigvalsh(block.toarray()),
    )


def test_heisenberg_sector_matches_full_space_block() -> None:
    sector = heisenberg_chain_sector(
        n_sites=5,
        magnetization=1,
        jx=0.8,
        jy=0.8,
        jz=1.2,
        field=-0.15,
        periodic=True,
    )
    full = heisenberg_chain_sparse(
        n_sites=5,
        jx=0.8,
        jy=0.8,
        jz=1.2,
        field=-0.15,
        periodic=True,
    )
    states = np.asarray(sector.basis.states)

    assert np.allclose(sector.matrix.toarray(), full[states][:, states].toarray())
    assert sector.to_metadata()["sector"]["magnetization"] == 1


def test_sector_model_specs_dimension_metadata_and_storage(tmp_path) -> None:
    spec = create_model_spec(
        "xxz_chain_sector_sparse",
        parameters={"n_sites": 6, "magnetization": 0, "anisotropy": 0.5},
    )
    result = spec.build_result()

    assert result.shape == (20, 20)
    assert result.representation == "sparse"
    assert result.metadata["sector"]["basis_states"]
    assert estimate_model_dimension("xxz_chain_sector_sparse", n_sites=8, magnetization=0) == comb(
        8, 4
    )

    path = save_hamiltonian(result, tmp_path / "sector.npz")
    restored = load_hamiltonian(path)
    assert restored.metadata["sector"]["magnetization"] == 0
    assert np.allclose(restored.matrix.toarray(), result.matrix.toarray())


def test_sector_file_cli_workflow(tmp_path, capsys) -> None:
    model_path = tmp_path / "sector.json"
    assert (
        main(
            [
                "create",
                "xxz_chain_sector_sparse",
                "--n-sites",
                "6",
                "--magnetization",
                "0",
                "--output",
                str(model_path),
            ]
        )
        == 0
    )
    capsys.readouterr()

    assert main(["spectrum", str(model_path)]) == 0
    assert len(capsys.readouterr().out.strip().splitlines()) == comb(6, 3)


def test_sector_validation_failures() -> None:
    with pytest.raises(ValueError, match="same parity"):
        fixed_magnetization_basis(4, 1)
    with pytest.raises(ValueError, match=r"abs\(magnetization\)"):
        fixed_magnetization_basis(4, 6)
    with pytest.raises(ValueError, match="jx == jy"):
        heisenberg_chain_sector(4, 0, jx=1.0, jy=0.5)
    with pytest.raises(ValueError, match="same parity"):
        estimate_model_dimension("xxz_chain_sector_sparse", n_sites=5, magnetization=0)
