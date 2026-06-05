# Usage

Install locally:

```bash
pip install -e ".[dev]"
```

## Transverse-Field Ising Chain

```python
from quantum_lattice_models.models import transverse_field_ising
from quantum_lattice_models.spectra import eigenvalues, ground_energy, spectral_gap

H = transverse_field_ising(n_sites=4, j=1.0, h=0.5, periodic=False)

print(H.shape)
print(eigenvalues(H))
print(ground_energy(H))
print(spectral_gap(H))
```

The convention is:

```text
H = -J sum_i Z_i Z_{i+1} - h sum_i X_i
```

## Heisenberg Chain

```python
from quantum_lattice_models.models import heisenberg_chain

H = heisenberg_chain(
    n_sites=4,
    jx=1.0,
    jy=1.0,
    jz=1.0,
    field=0.0,
    periodic=False,
)
```

The convention is:

```text
H = sum_i (Jx X_i X_{i+1} + Jy Y_i Y_{i+1} + Jz Z_i Z_{i+1})
    + field sum_i Z_i
```

## SSH Model

```python
import numpy as np

from quantum_lattice_models.models import ssh_edge_state_localizations, ssh_model
from quantum_lattice_models.spectra import eigensystem

H = ssh_model(n_cells=6, t1=0.5, t2=1.0, periodic=False)
values, vectors = eigensystem(H)

near_zero = np.argsort(abs(values))[:2]
edge_weights = ssh_edge_state_localizations(vectors[:, near_zero], n_cells=6, edge_cells=1)
print(values[near_zero])
print(edge_weights)
```

`ssh_model` returns a single-particle matrix with dimension `2 * n_cells`.

## Generic Tight-Binding Chain

```python
from quantum_lattice_models.models import tight_binding_chain

H = tight_binding_chain(n_sites=8, hopping=1.0, onsite=0.0, periodic=False)
```

The `onsite` argument can be a scalar or a length-`n_sites` iterable.

## Observables

```python
from quantum_lattice_models.models import transverse_field_ising
from quantum_lattice_models.observables import correlation_zz, magnetization_z
from quantum_lattice_models.spectra import eigensystem

H = transverse_field_ising(n_sites=4, j=1.0, h=0.5)
values, vectors = eigensystem(H)
ground_state = vectors[:, 0]

print(magnetization_z(ground_state, n_sites=4))
print(correlation_zz(ground_state, n_sites=4, i=0, j=1))
```

## Plotting

```python
import matplotlib.pyplot as plt

from quantum_lattice_models.models import transverse_field_ising
from quantum_lattice_models.plotting import plot_spectrum

H = transverse_field_ising(n_sites=5, j=1.0, h=0.7)
ax = plot_spectrum(H)
ax.figure.tight_layout()
plt.show()
```

## Optional PennyLane Export

```python
from quantum_lattice_models.export import to_pennylane_terms
from quantum_lattice_models.models import transverse_field_ising

H = transverse_field_ising(n_sites=3, j=1.0, h=0.5)
qml_H = to_pennylane_terms(H)
```

If PennyLane is not installed, `to_pennylane_terms` raises an `ImportError` with installation guidance.

## Command-Line Examples

Each example saves one PNG under `images/`:

```bash
python examples/ising_spectrum.py
python examples/heisenberg_density.py
python examples/ssh_edge_state.py
python examples/tight_binding_spectrum.py
```
