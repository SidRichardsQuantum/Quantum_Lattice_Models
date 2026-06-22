<!-- builders: heisenberg_ladder -->
# Heisenberg Ladder

## Purpose and structure

The two-leg ladder couples two Heisenberg chains through rungs:

$$
H=J_{\rm leg}\sum_{\ell,r}\mathbf P_{\ell,r}\cdot\mathbf P_{\ell,r+1}
+J_{\rm rung}\sum_r\mathbf P_{0,r}\cdot\mathbf P_{1,r}
+g\sum_{\ell,r}Z_{\ell,r}.
$$

![Two-leg Heisenberg ladder](../diagrams/heisenberg_ladder.svg)

## Basis and scaling

$R$ rungs contain $2R$ spins, so the dense matrix dimension is $2^{2R}$.
This grows particularly quickly. Fixed total Pauli-$Z$ magnetization sectors
use dimension $\binom{2R}{(2R-M)/2}$.

## Package use

```python
from quantum_lattice_models import heisenberg_ladder

H = heisenberg_ladder(n_rungs=3, leg_coupling=1.0, rung_coupling=0.7)

from quantum_lattice_models import heisenberg_ladder_sector

sector = heisenberg_ladder_sector(n_rungs=6, magnetization=0)
```

## Parameters

{{PARAMETERS}}

## User notes

`periodic=True` closes each leg but does not change rung connectivity. Use
memory diagnostics before increasing the rung count.
