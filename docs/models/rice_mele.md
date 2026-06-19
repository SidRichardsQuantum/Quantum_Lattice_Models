<!-- builders: rice_mele_model -->
# Rice-Mele Model

## Purpose and structure

The Rice-Mele chain extends SSH with alternating hoppings and staggered onsite
energies:

$$
t_{\rm intra}=t+\delta,\qquad t_{\rm inter}=t-\delta,
$$

with onsite values $+\Delta$ on $A$ and $-\Delta$ on $B$.

![SSH and Rice-Mele unit cells](../diagrams/ssh_rice_mele_unit_cell.svg)

## Basis and use

The single-particle dimension is $2N_c$.

```python
from quantum_lattice_models import rice_mele_model

H = rice_mele_model(
    n_cells=10, hopping=1.0, dimerization=0.25, staggering=0.5
)
```

## Parameters

{{PARAMETERS}}

## User notes

Nonzero staggering breaks the SSH sublattice symmetry and generally shifts
edge-state energies away from zero.
