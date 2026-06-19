<!-- builders: heisenberg_chain -->
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

The computational basis has dimension $2^N$. The builder returns a dense
matrix with Pauli-term metadata.

## Package use

```python
from quantum_lattice_models import heisenberg_chain

H = heisenberg_chain(n_sites=5, jx=1.0, jy=0.8, jz=1.2, field=0.1)
```

## Parameters

{{PARAMETERS}}

## Validation and cautions

Hermiticity and real spectra are tested. The isotropic limit is
$J_x=J_y=J_z$. Sparse and fixed-magnetization construction remain roadmap
work.

Related: [XXZ chain](xxz_chain.md), [Heisenberg ladder](heisenberg_ladder.md).
