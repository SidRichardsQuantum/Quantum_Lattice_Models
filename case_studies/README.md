# Model-intake case studies

These examples exercise real optional third-party objects while documenting
where model semantics must be supplied explicitly.

## NetworkX weighted chain

`networkx_weighted_chain.py` imports a labeled `MultiDiGraph` with real and
complex edge values, creates a portable tight-binding model, and reports
GraphML adapter capabilities. Geometry, labels, directed edges, and values are
available to the graph adapter; the model basis and interaction interpretation
remain portable model metadata.

Run with:

```bash
python -m pip install -e ".[networkx]"
python case_studies/networkx_weighted_chain.py
```

## ASE structure intake

`ase_structure_intake.py` starts from an `ase.Atoms` object. ASE provides
species, coordinates, cell, and periodicity, but it does not define a
tight-binding Hamiltonian. The case study therefore adds bonds, hopping
conventions, and provenance explicitly before model construction.

Run with:

```bash
python -m pip install -e ".[ase]"
python case_studies/ase_structure_intake.py
```

This separation is intentional: atomic coordinates alone are not treated as a
physically complete lattice Hamiltonian.
