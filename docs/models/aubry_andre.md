<!-- builders: aubry_andre_harper_chain -->
# Aubry-André-Harper Chain

## Purpose and structure

This quasiperiodic chain adds

$$
V_i=\lambda\cos(2\pi\beta i+\phi)
$$

to the onsite terms of a tight-binding chain. Irrational $\beta$ produces
quasiperiodicity and supports localization studies.

![Aubry-André localization](../diagrams/aubry_andre_localization.svg)

## Package use

```python
from quantum_lattice_models import aubry_andre_harper_chain

H = aubry_andre_harper_chain(
    n_sites=34, hopping=1.0, potential=2.0, phase=0.0
)
```

## Parameters

{{PARAMETERS}}

## User notes

This builder currently returns a dense single-particle matrix. Use inverse
participation ratios to compare extended and localized eigenstates. Finite
rational approximants depend on system size and boundary conditions.
