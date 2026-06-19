<!-- builders: ssh_model -->
# SSH Model

## Purpose and structure

The Su-Schrieffer-Heeger chain has $A$ and $B$ sites in each unit cell, with
intracell hopping $t_1$ and intercell hopping $t_2$:

$$
H=-\sum_m(t_1|m,A\rangle\langle m,B|
+t_2|m+1,A\rangle\langle m,B|+\mathrm{h.c.}).
$$

![SSH unit cell](../diagrams/ssh_rice_mele_unit_cell.svg)

For open boundaries and $|t_1|<|t_2|$, finite chains support near-zero
edge-localized states.

## Basis and use

The single-particle dimension is $2N_c$.

```python
from quantum_lattice_models import ssh_model

H = ssh_model(n_cells=12, t1=0.4, t2=1.0)
```

## Parameters

{{PARAMETERS}}

## Validation and cautions

The decoupled-dimer spectrum and topological edge localization are tested.
This is a finite real-space model; winding numbers require future Bloch
support.

Related: [Rice-Mele model](rice_mele.md), [Kitaev BdG chain](kitaev_bdg.md).
