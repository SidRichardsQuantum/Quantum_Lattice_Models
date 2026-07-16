"""ASE geometry-intake case study with explicit post-import bonds."""

from __future__ import annotations

from dataclasses import replace

from ase import Atoms

from quantum_lattice_models import Bond, create_model_spec, from_ase

atoms = Atoms(
    symbols=["H", "H", "H"],
    positions=[(0.0, 0.0, 0.0), (0.8, 0.0, 0.0), (1.6, 0.0, 0.0)],
    cell=(3.0, 3.0, 3.0),
    pbc=(False, False, False),
)
geometry = from_ase(atoms)

# ASE supplies structure data, not a hopping model. Add bonds and conventions
# explicitly before constructing a portable Hamiltonian.
lattice = replace(
    geometry,
    bonds=(Bond(0, 1), Bond(1, 2)),
    conventions={
        **geometry.conventions,
        "bond_source": "explicit case-study nearest-neighbor assignment",
    },
)
model = create_model_spec(
    "custom_tight_binding",
    lattice=lattice,
    parameters={"hopping": 1.0, "onsite": 0.0, "hermitian": True},
    provenance={"source": "ASE Atoms case study"},
)

print(model.summary())
