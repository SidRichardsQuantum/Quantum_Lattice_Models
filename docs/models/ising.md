<!-- builders: transverse_field_ising, longitudinal_field_ising, next_nearest_neighbor_ising -->
# Ising Chains

## Purpose and structure

These spin-$\tfrac12$ chains combine Ising $ZZ$ interactions with transverse
and optional longitudinal fields. The next-nearest-neighbor variant adds
frustration. They are standard small-system benchmarks for phase-transition
intuition, quantum simulation, and variational algorithms.

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

## Basis and scaling

The computational basis has dimension $2^N$. Builders currently return dense
`DenseHamiltonian` arrays with Pauli-term metadata.

## Package use

```python
from quantum_lattice_models import transverse_field_ising

H = transverse_field_ising(n_sites=6, j=1.0, h=0.7, periodic=False)
```

```bash
quantum-lattice create transverse_field_ising --n-sites 6 --j 1 --h 0.7 --output ising.json
```

## Parameters

{{PARAMETERS}}

## Validation and cautions

The zero-field three-site spectrum is checked analytically. Dense memory grows
as $4^N$; inspect the estimated dimension and memory before increasing $N$.

Related: [XY chain](xy_chain.md), [XXZ chain](xxz_chain.md).
