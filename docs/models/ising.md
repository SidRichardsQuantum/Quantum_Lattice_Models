<!-- builders: transverse_field_ising, longitudinal_field_ising, next_nearest_neighbor_ising, transverse_field_ising_parity_sector_sparse -->
# Ising Chains

## Purpose and structure

These spin-$\tfrac12$ chains combine Ising $ZZ$ interactions with transverse
and optional longitudinal fields. The next-nearest-neighbor variant adds
frustration. They are standard small-system benchmarks for phase-transition
intuition, quantum simulation, and variational algorithms. The transverse-field
model also supports explicit global spin-flip parity sectors.

![Ising and related spin-chain couplings](../diagrams/spin_chain_couplings.svg)

## Hamiltonians

$$
H_{\rm TFIM}=-J\sum_i Z_iZ_{i+1}-h\sum_iX_i,
$$

$$
H_{\rm long}=-J\sum_iZ_iZ_{i+1}-h_x\sum_iX_i-h_z\sum_iZ_i,
$$

$$
H_{\rm NNN}=-J_1\sum_iZ_iZ_{i+1}-J_2\sum_iZ_iZ_{i+2}-h\sum_iX_i.
$$

The operators are Pauli matrices, not spin operators divided by two.

For the transverse-field model,

$$
P=\prod_i X_i,\qquad [H_{\rm TFIM},P]=0,
$$

so the Hamiltonian can be reduced to either parity eigenvalue $p=\pm1$.

## Basis and scaling

The full computational basis has dimension $2^N$. The ordinary builder returns
a dense `DenseHamiltonian`; the compatible `_sparse` builder returns a CSR
matrix.

Each spin-flip parity sector has dimension $2^{N-1}$. Its basis rows are the
normalized superpositions

$$
|s;p\rangle=\frac{|s\rangle+p|\bar{s}\rangle}{\sqrt{2}},
$$

where $|\bar{s}\rangle$ is the bitwise-complement state. The sector builder
returns a `SpinParitySectorHamiltonian` containing a CSR matrix and a portable
reduced-basis mapping with both components and coefficients.

## Package use

```python
from quantum_lattice_models import (
    transverse_field_ising,
    transverse_field_ising_parity_sector,
)

H = transverse_field_ising(n_sites=6, j=1.0, h=0.7, periodic=False)
even = transverse_field_ising_parity_sector(
    n_sites=6,
    parity=1,
    j=1.0,
    h=0.7,
    periodic=False,
)
H_even = even.matrix
```

```bash
quantum-lattice create transverse_field_ising --n-sites 6 --j 1 --h 0.7 --output ising.json
quantum-lattice create transverse_field_ising_parity_sector_sparse \
  --n-sites 6 --parity 1 --j 1 --h 0.7 --output ising-even.json
```

## Parameters

{{PARAMETERS}}

## Validation and cautions

The zero-field three-site spectrum is checked analytically. Parity mappings are
checked as orthonormal eigenbases of $P$, each reduced matrix is compared with
an explicit full-space isometric reduction, and the two sector spectra are
verified to recombine into the full spectrum.

Dense memory grows as $4^N$. A parity sector halves the matrix dimension but
does not remove exponential scaling, so inspect the estimated dimension and
storage before increasing $N$. A nonzero longitudinal $Z$ field breaks the
global spin-flip symmetry and is not supported by the parity-sector builder.

Related: [XY chain](xy_chain.md), [XXZ chain](xxz_chain.md).
