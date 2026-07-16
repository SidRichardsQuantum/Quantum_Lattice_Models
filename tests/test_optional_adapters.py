from __future__ import annotations

import runpy

import numpy as np
import pytest

from quantum_lattice_models import from_ase


def test_ase_import_preserves_coordinates_cell_and_periodicity() -> None:
    ase = pytest.importorskip("ase")
    atoms = ase.Atoms(
        symbols=("H", "He"),
        positions=((0.0, 0.0, 0.0), (1.25, 0.5, 0.0)),
        cell=np.diag((3.0, 4.0, 5.0)),
        pbc=(True, False, True),
    )

    lattice = from_ase(atoms)

    assert lattice.n_sites == 2
    assert lattice.site_labels == ("H", "He")
    assert lattice.positions == ((0.0, 0.0, 0.0), (1.25, 0.5, 0.0))
    assert lattice.units == {"position": "angstrom"}
    assert lattice.metadata["cell"] == [[3.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 5.0]]
    assert lattice.metadata["periodic_axes"] == [True, False, True]


def test_ase_structure_intake_case_study() -> None:
    pytest.importorskip("ase")
    namespace = runpy.run_path("case_studies/ase_structure_intake.py")
    assert namespace["model"].lattice.n_sites == 3
