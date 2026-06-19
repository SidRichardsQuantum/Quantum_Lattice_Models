<!-- builders: haldane_honeycomb_lattice, haldane_honeycomb_lattice_sparse -->
# Haldane Honeycomb Lattice

## Purpose and structure

The finite Haldane model has two honeycomb sublattices, real nearest-neighbor
hopping $t_1$, oriented complex next-nearest-neighbor hopping
$t_2e^{\pm i\phi}$, and staggered onsite potential $\pm M$.

![Haldane honeycomb model](../diagrams/haldane_honeycomb.svg)

## Basis and scaling

The single-particle dimension is $2N_rN_c$. Dense and CSR builders are
available.

```python
from quantum_lattice_models import haldane_honeycomb_lattice

H = haldane_honeycomb_lattice(
    n_rows=4, n_cols=4, t1=1.0, t2=0.18, phi=1.5708
)
```

## Parameters

{{PARAMETERS}}

## Validation and cautions

Complex next-nearest-neighbor phases and Hermitian conjugates are tested.
Finite real-space spectra do not directly provide a Chern number; Bloch and
topological-analysis support is planned.
