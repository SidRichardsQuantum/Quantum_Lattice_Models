<!-- builders: kagome_lattice_tight_binding, kagome_lattice_tight_binding_sparse -->
# Kagome-Lattice Tight Binding

## Purpose and structure

The kagome lattice has three sublattice sites per unit cell and triangular
connectivity. It is a standard flat-band and frustrated-lattice testbed.

![Finite lattice geometries](../diagrams/lattice_geometry_comparison.svg)

## Basis and scaling

The single-particle dimension is $3N_rN_c$. Dense and CSR builders are
available.

```python
from quantum_lattice_models import kagome_lattice_tight_binding_sparse

H = kagome_lattice_tight_binding_sparse(n_rows=4, n_cols=4)
```

## Parameters

{{PARAMETERS}}

## User notes

Use `kagome_lattice_positions` and `kagome_lattice_index` for plotting and
indexing. Sublattice labels are not automatically attached to sparse matrices.
