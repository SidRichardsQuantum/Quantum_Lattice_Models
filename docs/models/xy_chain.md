<!-- builders: xy_chain -->
# XY Chain

## Purpose and Hamiltonian

The anisotropic XY chain isolates $XX$ and $YY$ exchange:

$$
H=-J\sum_i\left[\frac{1+\gamma}{2}X_iX_{i+1}
+\frac{1-\gamma}{2}Y_iY_{i+1}\right]-g\sum_iZ_i.
$$

`coupling` is $J$, `anisotropy` is $\gamma$, and `field` is $g$.
$\gamma=0$ gives equal $XX$ and $YY$ couplings.

## Basis and use

The dense computational-basis matrix has dimension $2^N$.

```python
from quantum_lattice_models import xy_chain

H = xy_chain(n_sites=5, coupling=1.0, anisotropy=0.3, field=0.2)
```

## Parameters

{{PARAMETERS}}

## User notes

The field sign differs from the general Heisenberg builder: this model uses
$-g\sum Z_i$. Check conventions when comparing parameterizations.

Related: [Ising chains](ising.md), [XXZ chain](xxz_chain.md).
