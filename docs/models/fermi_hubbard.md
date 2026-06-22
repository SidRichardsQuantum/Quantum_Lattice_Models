<!-- builders: fermi_hubbard_chain, fermi_hubbard_chain_sparse -->
# Fermi-Hubbard Chain

## Purpose and Hamiltonian

The spinful Fermi-Hubbard chain is

$$
H=-t\sum_{\langle i,j\rangle,\sigma}
(c_{i\sigma}^\dagger c_{j\sigma}+\mathrm{h.c.})
+U\sum_i n_{i\uparrow}n_{i\downarrow}
-\mu\sum_{i,\sigma}n_{i\sigma}.
$$

The orbital order is
$(0\uparrow,0\downarrow,1\uparrow,1\downarrow,\ldots)$.

![Fermionic hopping sign](../diagrams/fermi_hubbard_sign.svg)

## Basis and scaling

There are $2N$ binary orbitals and dimension $2^{2N}$. Dense and CSR builders
are available. Fixed spin-resolved particle numbers reduce the dimension to
$\binom{N}{N_\uparrow}\binom{N}{N_\downarrow}$ while preserving explicit
full-basis occupation-state mappings.

## Package use

```python
from quantum_lattice_models import fermi_hubbard_chain_sparse

H = fermi_hubbard_chain_sparse(n_sites=4, hopping=1.0, interaction=4.0)

from quantum_lattice_models import fermi_hubbard_chain_sector

sector = fermi_hubbard_chain_sector(
    n_sites=8, n_up=4, n_down=4, interaction=4.0
)
```

## Parameters

{{PARAMETERS}}

## Validation and cautions

Single-site energies and explicit fermionic parity-sign cases are tested.
Sector matrices are validated against full-space blocks for open and periodic
chains. The full occupation basis grows faster than the spin-chain basis.
