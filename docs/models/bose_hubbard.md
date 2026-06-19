<!-- builders: bose_hubbard_chain, bose_hubbard_chain_sparse -->
# Bose-Hubbard Chain

## Purpose and Hamiltonian

The truncated Bose-Hubbard chain describes hopping bosons with onsite
interaction and chemical potential:

$$
H=-t\sum_{\langle i,j\rangle}(a_i^\dagger a_j+a_j^\dagger a_i)
+\frac U2\sum_i n_i(n_i-1)-\mu\sum_i n_i.
$$

![Truncated boson basis](../diagrams/hubbard_basis.svg)

## Basis and scaling

Each site uses occupations $0,\ldots,n_{\max}$, giving dimension
$(n_{\max}+1)^N$. Dense and CSR builders are available.

## Package use

```python
from quantum_lattice_models import bose_hubbard_chain_sparse

H = bose_hubbard_chain_sparse(
    n_sites=4, hopping=0.6, interaction=1.5, max_occupancy=2
)
```

## Parameters

{{PARAMETERS}}

## Validation and cautions

Single-site energies and particle-number conservation are tested. Sparse
storage does not remove exponential basis growth.
