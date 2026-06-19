<!-- builders: square_lattice_tight_binding, square_lattice_tight_binding_sparse -->
# Square-Lattice Tight Binding

## Purpose and structure

This model places one orbital at each site of an $N_r\times N_c$ rectangular
grid with row-major indexing and nearest-neighbor hopping.

![Finite lattice geometries](../diagrams/lattice_geometry_comparison.svg)

## Basis and scaling

The single-particle dimension is $N_rN_c$. Dense and CSR builders are
available. Opposite edges can be reconnected independently along $x$ and $y$.

```python
from quantum_lattice_models import square_lattice_tight_binding_sparse

H = square_lattice_tight_binding_sparse(
    n_rows=20, n_cols=20, periodic_x=True
)
```

## Parameters

{{PARAMETERS}}

## User notes

Use `square_lattice_positions` for plotting. This model has real hopping and
no magnetic flux; use [Harper-Hofstadter](harper_hofstadter.md) for Peierls
phases.
