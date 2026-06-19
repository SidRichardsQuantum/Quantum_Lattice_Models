<!-- builders: xxz_chain, xxz_chain_sparse, xxz_chain_sector_sparse -->
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

The full computational-basis matrix has dimension $2^N$. Dense and CSR sparse
builders are available. A fixed total Pauli-$Z$ magnetization sector
$M=\sum_i Z_i$ has dimension $\binom{N}{(N-M)/2}$.

```python
from quantum_lattice_models import xxz_chain, xxz_chain_sector

H = xxz_chain(n_sites=6, coupling=1.0, anisotropy=0.7)
sector = xxz_chain_sector(n_sites=10, magnetization=0, anisotropy=0.7)
H_sector = sector.matrix
```

## Parameters

{{PARAMETERS}}

## User notes

`xxz_chain` delegates to `heisenberg_chain`; its field term therefore has a
positive sign. `magnetization` is the total Pauli-$Z$ eigenvalue, so it must
have the same parity as `n_sites` and satisfy $|M|\le N$. Sector basis mappings
support projection from and embedding into the full computational basis.

Related: [Heisenberg chain](heisenberg_chain.md), [XY chain](xy_chain.md).
