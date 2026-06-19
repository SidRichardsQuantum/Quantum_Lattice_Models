<!-- builders: j1_j2_heisenberg_chain -->
# J1-J2 Heisenberg Chain

## Purpose and structure

This frustrated spin chain combines nearest- and next-nearest-neighbor
Heisenberg exchange:

$$
H=J_1\sum_i\mathbf P_i\cdot\mathbf P_{i+1}
+J_2\sum_i\mathbf P_i\cdot\mathbf P_{i+2}+g\sum_iZ_i,
$$

with $\mathbf P=(X,Y,Z)$. Competition between $J_1$ and $J_2$ creates a
compact frustration benchmark.

![Spin-chain coupling comparison](../diagrams/spin_chain_couplings.svg)

## Basis and use

```python
from quantum_lattice_models import j1_j2_heisenberg_chain

H = j1_j2_heisenberg_chain(n_sites=6, j1=1.0, j2=0.4)
```

The dense computational-basis dimension is $2^N$.

## Parameters

{{PARAMETERS}}

## User notes

Periodic next-nearest-neighbor bonds are defined only when enough distinct
sites exist. The package uses Pauli products directly.
