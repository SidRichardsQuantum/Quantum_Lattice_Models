<!-- builders: kitaev_chain_bdg -->
# Kitaev BdG Chain

## Purpose and structure

The spinless p-wave superconducting chain combines hopping, chemical
potential, and nearest-neighbor pairing. The package returns

$$
\mathcal H_{\rm BdG}=
\begin{pmatrix}A&D\\-D^*&-A^*\end{pmatrix}
$$

in the Nambu basis
$(c_0,\ldots,c_{N-1},c_0^\dagger,\ldots,c_{N-1}^\dagger)^T$.

![Kitaev Nambu chain](../diagrams/kitaev_nambu_chain.svg)

## Basis and use

The BdG dimension is $2N$, not $2^N$.

```python
from quantum_lattice_models import kitaev_chain_bdg

H = kitaev_chain_bdg(n_sites=10, hopping=1.0, pairing=0.5)
```

## Parameters

{{PARAMETERS}}

## Validation and cautions

Hermiticity and spectral particle-hole pairing $E\leftrightarrow-E$ are
tested. The returned matrix is a single-particle BdG representation, not the
many-body Hamiltonian.
