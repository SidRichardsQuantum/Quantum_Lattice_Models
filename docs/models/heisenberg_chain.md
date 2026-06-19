<!-- builders: heisenberg_chain, heisenberg_chain_sparse, heisenberg_chain_sector_sparse -->
# Heisenberg Chain

## Purpose and structure

The anisotropic nearest-neighbor Heisenberg chain supports independent $XX$,
$YY$, and $ZZ$ couplings and a uniform longitudinal field. It is useful for
spin correlations, symmetry studies, and exact-diagonalization benchmarks.

$$
H=\sum_i(J_xX_iX_{i+1}+J_yY_iY_{i+1}+J_zZ_iZ_{i+1})
+g\sum_iZ_i.
$$

The package uses Pauli products directly and a positive sign for `field`.

## Basis and scaling

The full computational basis has dimension $2^N$. Dense and CSR sparse
builders are available. When $J_x=J_y$, total Pauli-$Z$ magnetization
$M=\sum_i Z_i$ is conserved and the fixed-sector dimension is

$$
\binom{N}{(N-M)/2}.
$$

Sector basis states retain their full computational-basis integer indices.

## Package use

```python
from quantum_lattice_models import heisenberg_chain, heisenberg_chain_sector

H = heisenberg_chain(n_sites=5, jx=1.0, jy=0.8, jz=1.2, field=0.1)
sector = heisenberg_chain_sector(
    n_sites=6,
    magnetization=0,
    jx=1.0,
    jy=1.0,
    jz=1.2,
)
H_sector = sector.matrix
```

## Parameters

{{PARAMETERS}}

## Validation and cautions

Hermiticity and real spectra are tested. The isotropic limit is
$J_x=J_y=J_z$. Fixed-magnetization construction rejects $J_x\ne J_y$ because
that anisotropy does not conserve total $Z$ magnetization. Reduced matrices
are tested against the corresponding full-space blocks.

Related: [XXZ chain](xxz_chain.md), [Heisenberg ladder](heisenberg_ladder.md).
