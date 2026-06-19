<!-- builders: tight_binding_chain, tight_binding_chain_sparse -->
# Tight-Binding Chain

## Purpose and Hamiltonian

The generic one-dimensional single-particle model is

$$
H=-t\sum_{\langle i,j\rangle}(|i\rangle\langle j|+\mathrm{h.c.})
+\sum_i\epsilon_i|i\rangle\langle i|.
$$

`onsite` may be a scalar or one value per site. Dense and CSR builders are
available.

## Package use

```python
from quantum_lattice_models import tight_binding_chain_sparse

H = tight_binding_chain_sparse(
    n_sites=128, hopping=1.0, onsite=0.2, periodic=False
)
```

## Parameters

{{PARAMETERS}}

## Validation and cautions

The open uniform-chain spectrum is checked against
$E_m=-2t\cos[m\pi/(N+1)]$. For `n_sites=2`, periodic construction includes
both directed wrap contributions, consistent with the package bond convention.
