<!-- builders: xxz_chain -->
# XXZ Chain

## Purpose and Hamiltonian

The XXZ chain is the $J_x=J_y$ specialization of the anisotropic Heisenberg
chain:

$$
H=J\sum_i(X_iX_{i+1}+Y_iY_{i+1}+\Delta Z_iZ_{i+1})
+g\sum_iZ_i.
$$

It is useful for anisotropy, magnetization, gap, and conserved-$S^z$
benchmarks.

## Basis and use

The dense computational-basis matrix has dimension $2^N$.

```python
from quantum_lattice_models import xxz_chain

H = xxz_chain(n_sites=6, coupling=1.0, anisotropy=0.7)
```

## Parameters

{{PARAMETERS}}

## User notes

`xxz_chain` delegates to `heisenberg_chain`; its field term therefore has a
positive sign. Fixed-magnetization sectors are planned but not yet available.

Related: [Heisenberg chain](heisenberg_chain.md), [XY chain](xy_chain.md).
