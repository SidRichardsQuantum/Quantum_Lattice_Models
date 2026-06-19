<!-- builders: harper_hofstadter_square_lattice, harper_hofstadter_square_lattice_sparse -->
# Harper-Hofstadter Square Lattice

## Purpose and gauge

This square lattice adds magnetic Peierls phases in Landau gauge. Horizontal
hoppings are real; a vertical hopping in column $x$ carries

$$
-t\exp(2\pi ifx),
$$

where $f$ is flux per plaquette in flux-quantum units.

![Hofstadter plaquette phases](../diagrams/hofstadter_plaquette.svg)

## Basis and use

The single-particle dimension is $N_rN_c$. Dense and CSR builders are
available.

```python
from quantum_lattice_models import harper_hofstadter_square_lattice

H = harper_hofstadter_square_lattice(n_rows=8, n_cols=8, flux=0.25)
```

## Parameters

{{PARAMETERS}}

## Validation and cautions

Representative vertical phases and dense/sparse equivalence are tested. A
finite flux sweep is not the infinite-system Hofstadter butterfly.
