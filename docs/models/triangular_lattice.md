<!-- builders: triangular_lattice_tight_binding, triangular_lattice_tight_binding_sparse -->
# Triangular-Lattice Tight Binding

## Purpose and structure

The triangular lattice adds diagonal neighbors to the square-grid indexing,
giving six bulk nearest neighbors and geometric frustration.

![Finite lattice geometries](../diagrams/lattice_geometry_comparison.svg)

## Basis and scaling

The single-particle dimension is $N_rN_c$. Dense and CSR builders are
available, with independent periodic flags.

```python
from quantum_lattice_models import triangular_lattice_tight_binding

H = triangular_lattice_tight_binding(n_rows=4, n_cols=5, hopping=1.0)
```

## Parameters

{{PARAMETERS}}

## User notes

Use `triangular_lattice_positions` for geometry-aware plotting. Small periodic
dimensions may cause multiple geometric bonds to connect the same pair; matrix
entries accumulate those contributions.
