# Usage

Install locally:

```bash
pip install -e ".[dev]"
```

Optional notebook support:

```bash
python -m pip install -e ".[notebooks]"
python -m ipykernel install --user --name quantum-lattice-models --display-name "Quantum Lattice Models"
```

Notebook workflows, numbered in the recommended learning order:

```text
notebooks/01_ising_spin_chains.ipynb
notebooks/02_spin_observables_and_correlations.ipynb
notebooks/03_spin_chain_model_comparison.ipynb
notebooks/04_heisenberg_ladder_spectrum.ipynb
notebooks/05_hubbard_exact_diagonalization.ipynb
notebooks/06_ssh_rice_mele_comparison.ipynb
notebooks/07_boundary_conditions_and_finite_size.ipynb
notebooks/08_aubry_andre_localization.ipynb
notebooks/09_kitaev_bdg_symmetry.ipynb
notebooks/10_hofstadter_flux_sweep.ipynb
notebooks/11_haldane_kagome_lattices.ipynb
notebooks/12_custom_lattice_workflow.ipynb
notebooks/13_hamiltonian_structure_gallery.ipynb
notebooks/14_sparse_dense_scaling.ipynb
notebooks/15_pennylane_export.ipynb
notebooks/16_model_registry_and_cli.ipynb
notebooks/17_cli_plot_walkthrough.ipynb
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

$$
H=-J\sum_i Z_iZ_{i+1}-h\sum_iX_i.
$$

## Ising Variants

```python
from quantum_lattice_models.models import (
    longitudinal_field_ising,
    next_nearest_neighbor_ising,
)

H_lfi = longitudinal_field_ising(n_sites=4, j=1.0, h_x=0.5, h_z=0.1)
H_nnn = next_nearest_neighbor_ising(n_sites=5, j1=1.0, j2=0.25, h=0.5)
```

The longitudinal-field model adds $-h_z\sum_iZ_i$.
The next-nearest-neighbor model adds $-J_2\sum_iZ_iZ_{i+2}$.

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

$$
H=\sum_i\left(
J_xX_iX_{i+1}+J_yY_iY_{i+1}+J_zZ_iZ_{i+1}
\right)+g\sum_iZ_i.
$$

## XY Chain

```python
from quantum_lattice_models.models import xy_chain

H = xy_chain(
    n_sites=5,
    coupling=1.0,
    anisotropy=0.3,
    field=0.2,
    periodic=False,
)
```

The convention is:

$$
H=-J\sum_i\left[
\frac{1+\gamma}{2}X_iX_{i+1}
+\frac{1-\gamma}{2}Y_iY_{i+1}
\right]-g\sum_iZ_i.
$$

## XXZ Chain

```python
from quantum_lattice_models.models import xxz_chain

H = xxz_chain(n_sites=5, coupling=1.0, anisotropy=0.7, field=0.0)
```

The XXZ chain is a convenience wrapper around the anisotropic Heisenberg builder
with $J_x=J_y=J$ and $J_z=J\Delta$, where $J$ is `coupling` and
$\Delta$ is `anisotropy`.

## J1-J2 Heisenberg Chain

```python
from quantum_lattice_models.models import j1_j2_heisenberg_chain

H = j1_j2_heisenberg_chain(n_sites=6, j1=1.0, j2=0.4, field=0.0)
```

This model adds next-nearest-neighbor Heisenberg couplings to the nearest-neighbor
chain, which gives a compact frustrated spin-system testbed.

## Two-Leg Heisenberg Ladder

```python
from quantum_lattice_models.models import heisenberg_ladder

H = heisenberg_ladder(n_rungs=3, leg_coupling=1.0, rung_coupling=0.7)
```

Sites are ordered as the top leg first, then the bottom leg. Rung $r$ connects
sites $r$ and $N_r+r$, where $N_r$ is `n_rungs`.

## Hubbard Models

```python
from quantum_lattice_models.models import bose_hubbard_chain, fermi_hubbard_chain

H_bose = bose_hubbard_chain(
    n_sites=3,
    hopping=1.0,
    interaction=2.0,
    chemical_potential=0.0,
    max_occupancy=2,
)

H_fermi = fermi_hubbard_chain(
    n_sites=3,
    hopping=1.0,
    interaction=4.0,
    chemical_potential=0.0,
)
```

The Bose-Hubbard builder uses the local basis
$\{|0\rangle,\ldots,|n_{\max}\rangle\}$, where $n_{\max}$ is
`max_occupancy`.
The Fermi-Hubbard builder uses orbital order
$(0\uparrow,0\downarrow,1\uparrow,1\downarrow,\ldots)$.

Sparse variants are available for larger small-system workflows:

```python
from quantum_lattice_models.models import (
    bose_hubbard_chain_sparse,
    fermi_hubbard_chain_sparse,
    haldane_honeycomb_lattice_sparse,
    harper_hofstadter_square_lattice_sparse,
    kagome_lattice_tight_binding_sparse,
    square_lattice_tight_binding_sparse,
    tight_binding_chain_sparse,
    triangular_lattice_tight_binding_sparse,
)

H_sparse = tight_binding_chain_sparse(n_sites=128, hopping=1.0)
```

Sparse lattice builders assemble SciPy sparse matrices directly instead of
constructing an intermediate dense matrix. Dense and sparse variants share the
same construction path and are tested for equivalent results, including
periodic boundaries and complex hopping phases.

Sparse spectral helpers are available from `quantum_lattice_models.spectra`:

```python
from quantum_lattice_models.models import fermi_hubbard_chain_sparse
from quantum_lattice_models.spectra import ground_state, lowest_eigenvalues

H = fermi_hubbard_chain_sparse(n_sites=4, hopping=1.0, interaction=4.0)
print(lowest_eigenvalues(H, k=3))
energy, state = ground_state(H)
```

## Kitaev Chain

```python
from quantum_lattice_models.models import kitaev_chain_bdg

H = kitaev_chain_bdg(n_sites=8, hopping=1.0, chemical_potential=0.0, pairing=0.5)
```

This returns a $2N\times2N$ Bogoliubov-de Gennes matrix for $N$ sites.

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

`ssh_model` returns a single-particle matrix with dimension $2N_c$, where
$N_c$ is `n_cells`.

## Rice-Mele Model

```python
from quantum_lattice_models.models import rice_mele_model

H = rice_mele_model(
    n_cells=8,
    hopping=1.0,
    dimerization=0.25,
    staggering=0.5,
)
```

The Rice-Mele model extends SSH with staggered onsite potentials.

## Generic Tight-Binding Chain

```python
from quantum_lattice_models.models import tight_binding_chain

H = tight_binding_chain(n_sites=8, hopping=1.0, onsite=0.0, periodic=False)
```

The `onsite` argument can be a scalar or a length-$N$ iterable for $N$ sites.

## User-Defined Lattices and Models

Use `Lattice` when the model geometry is easier to express as sites and bonds
than as one of the built-in lattice families:

```python
from quantum_lattice_models import Lattice, TightBindingModel

lattice = Lattice(
    positions=[(0.0, 0.0), (1.0, 0.0), (0.5, 0.8)],
    bonds=[(0, 1), (1, 2, 0.25j), (2, 0)],
)

H = TightBindingModel(lattice).hamiltonian(
    hopping=1.0,
    onsite=[0.0, 0.1, 0.0],
)
```

Two-item bonds use `-hopping` as the matrix element.
Three-item bonds use the third value directly:

```python
from quantum_lattice_models import custom_tight_binding

H = custom_tight_binding(
    n_sites=4,
    bonds=[(0, 1), (1, 2, -0.5), (2, 3, 0.2j)],
    hopping=1.0,
)
```

Sparse custom models use the same inputs:

```python
from quantum_lattice_models import custom_tight_binding_sparse

H_sparse = custom_tight_binding_sparse(
    n_sites=128,
    bonds=[(i, i + 1) for i in range(127)],
    hopping=1.0,
)
```

Register a custom builder when you want it to appear in model discovery helpers
and CLI model choices:

```python
from quantum_lattice_models.registry import register_model


def three_site_chain(hopping: float = 1.0):
    return custom_tight_binding(n_sites=3, bonds=[(0, 1), (1, 2)], hopping=hopping)


register_model(
    "three_site_chain",
    category="user",
    basis="single particle",
    dimension="3",
    return_type="LatticeHamiltonian",
    description="Three-site user-defined chain",
    builder=three_site_chain,
    defaults={"hopping": 1.0},
)
```

## Visualization Helpers

The plotting helpers understand dense and sparse Hamiltonians. Custom
Hamiltonians built from `Lattice` can carry positions in their metadata, so
graph and state plots do not need duplicate coordinate arguments:

```python
import numpy as np

from quantum_lattice_models import Lattice, custom_tight_binding
from quantum_lattice_models.plotting import (
    plot_hamiltonian_matrix,
    plot_lattice_graph,
    plot_lattice_state,
    plot_parameter_sweep,
)
from quantum_lattice_models.spectra import eigensystem

lattice = Lattice(
    positions=[(0.0, 0.0), (1.0, 0.0), (0.5, 0.8)],
    bonds=[(0, 1), (1, 2, 0.25j), (2, 0)],
)
H = custom_tight_binding(lattice=lattice)

plot_lattice_graph(H, show_colorbar=True)
plot_hamiltonian_matrix(H, mode="phase")

values, vectors = eigensystem(H)
plot_lattice_state(H, vectors[:, 0])
```

`plot_spectrum(..., highlight_gap=True, zero_line=True)` can emphasize the
ground/first-excited pair and add an energy-zero reference. Lattice graph plots
also accept per-site `sublattices` and `unit_cells` arrays for color coding and
dashed unit-cell outlines. Phase colorbars use explicit `-π`, `0`, and `π`
ticks.

Parameter sweeps work with any builder that accepts one parameter value:

```python
from quantum_lattice_models.models import ssh_model

plot_parameter_sweep(
    lambda t1: ssh_model(n_cells=8, t1=t1, t2=1.0),
    np.linspace(0.2, 1.8, 32),
    parameter_name="Intracell hopping",
    title="SSH finite-chain spectrum",
)
```

## Square-Lattice Tight Binding

```python
from quantum_lattice_models.models import square_lattice_tight_binding

H = square_lattice_tight_binding(
    n_rows=3,
    n_cols=4,
    hopping=1.0,
    onsite=0.0,
    periodic_x=False,
    periodic_y=False,
)
```

Sites are ordered row-major, so matrix index $rN_c+c$ represents lattice
coordinate $(r,c)$, where $r$ is the row, $c$ is the column, and $N_c$ is
`n_cols`.

## Harper-Hofstadter Square Lattice

```python
from quantum_lattice_models.models import harper_hofstadter_square_lattice

H = harper_hofstadter_square_lattice(
    n_rows=4,
    n_cols=4,
    hopping=1.0,
    flux=0.25,
)
```

Landau gauge is used: horizontal hoppings are real and vertical hoppings carry
phase $\exp(2\pi i f x)$, where $f$ is `flux` and $x$ is the column index.

## Aubry-Andre-Harper Chain

```python
from quantum_lattice_models.models import aubry_andre_harper_chain

H = aubry_andre_harper_chain(
    n_sites=16,
    hopping=1.0,
    potential=1.5,
    beta=(5**0.5 - 1) / 2,
    phase=0.0,
)
```

The onsite potential is
$V_i=\lambda\cos(2\pi\beta i+\phi)$, where $\lambda$ is `potential`,
$\beta$ is `beta`, and $\phi$ is `phase`.

## Haldane, Triangular, and Kagome Lattices

```python
from quantum_lattice_models.models import (
    haldane_honeycomb_lattice,
    kagome_lattice_tight_binding,
    triangular_lattice_tight_binding,
)

H_haldane = haldane_honeycomb_lattice(n_rows=2, n_cols=3, t1=1.0, t2=0.1)
H_triangular = triangular_lattice_tight_binding(n_rows=3, n_cols=3)
H_kagome = kagome_lattice_tight_binding(n_rows=2, n_cols=2)
```

The Haldane model uses two sublattices per honeycomb unit cell.
The kagome builder uses three sublattices per unit cell.

## Geometry Helpers

```python
from quantum_lattice_models.geometry import (
    honeycomb_lattice_positions,
    kagome_lattice_positions,
    square_lattice_positions,
    triangular_lattice_positions,
)

positions = kagome_lattice_positions(n_rows=2, n_cols=3)
```

These helpers return coordinate arrays compatible with `plot_lattice_graph`.

## Model Registry

```python
from quantum_lattice_models.registry import get_model_info, list_models, model_table

print(list_models())
print(get_model_info("ssh_model"))
rows = model_table()
```

The registry records model category, basis, dimension scaling, return type,
short description, and typed `ParameterInfo` entries. Parameter metadata
includes defaults, CLI names, descriptions, repeatability, and basic
constraints. The CLI options are generated from this metadata.

## Computational Diagnostics

```python
from quantum_lattice_models import (
    diagnose_matrix,
    estimate_dense_memory,
    estimate_model_dimension,
    has_particle_hole_symmetric_spectrum,
)
from quantum_lattice_models.models import kitaev_chain_bdg

dimension = estimate_model_dimension("heisenberg_ladder", n_rungs=4)
print(dimension, estimate_dense_memory(dimension))

H = kitaev_chain_bdg(n_sites=8)
print(diagnose_matrix(H))
print(has_particle_hole_symmetric_spectrum(H))
```

`estimate_dense_memory` returns bytes for a square dense matrix. Sparse storage
reported by `diagnose_matrix` includes the CSR data, index, and pointer arrays.
See [VALIDATION.md](VALIDATION.md) for the automated scientific checks.
For equations, parameter definitions, basis conventions, and model-specific
cautions, browse the [model reference](docs/models/index.md). Shared numerical
conventions are documented in [THEORY.md](THEORY.md).

## Portable Model Specifications

```python
from quantum_lattice_models import create_model_spec, load_model

spec = create_model_spec(
    "ssh_model",
    parameters={"n_cells": 12, "t1": 0.4, "t2": 1.0},
)
spec.save("ssh.json")

loaded = load_model("ssh.json")
loaded.validate()
H_dense = loaded.hamiltonian()
```

When a registered sparse counterpart exists, the same specification can build
it without changing the saved model parameters:

```python
H_sparse = loaded.hamiltonian(sparse=True)
```

`LatticeSpec` stores custom site coordinates, complex bonds, labels, unit-cell
membership, boundary conditions, and metadata. Complex numbers use an explicit
JSON representation so they survive round trips without relying on
implementation-specific encoding.

## Command-Line Interface

After installation, use the `quantum-lattice` entry point:

```bash
quantum-lattice models
quantum-lattice models --category spin --sparse --json
quantum-lattice presets --model ssh_model
quantum-lattice dry-run --preset ssh_topological --n-cells 20 --json
quantum-lattice create ssh_model --n-cells 12 --t1 0.4 --output ssh.json
quantum-lattice create --preset ssh_topological --n-cells 12 --output ssh.json
quantum-lattice inspect ssh.json
quantum-lattice validate ssh.json
quantum-lattice compare ssh.json other-ssh.json --json
quantum-lattice spectrum --model ssh_model --n-cells 8
quantum-lattice plot --model harper_hofstadter_square_lattice --n-rows 4 --n-cols 4 --flux 0.25 --output images/hofstadter_cli.png
quantum-lattice spectrum --model custom_tight_binding --n-sites 3 --bond 0,1 --bond 1,2,0.25j
```

Model choices and defaults come from the model registry. Custom bonds use
`source,target` or `source,target,value`; the optional value accepts complex
numbers. `create` writes a versioned JSON specification, `inspect` reports its
model and resource summary, and `validate` checks the schema, registered
parameters, geometry, and requested dense or sparse representation.

`models` supports category, basis, sparse-capability, and validation-status
filters. `dry-run` reports dimensions, dense-memory requirements,
representation, basis, sparse availability, validation status, and warnings
without constructing a matrix. Add `--json` for deterministic
machine-readable output. Presets are transparent parameter dictionaries and
can be overridden with normal model options.

## Observables

```python
from quantum_lattice_models import (
    bipartite_entanglement_entropy,
    site_magnetization_z,
    spin_correlation_matrix,
    static_spin_structure_factor,
)
from quantum_lattice_models.models import transverse_field_ising
from quantum_lattice_models.spectra import eigensystem

H = transverse_field_ising(n_sites=4, j=1.0, h=0.5)
values, vectors = eigensystem(H)
ground_state = vectors[:, 0]

print(site_magnetization_z(ground_state, n_sites=4))
print(spin_correlation_matrix(ground_state, n_sites=4, connected=True))
print(static_spin_structure_factor(ground_state, n_sites=4, momenta=[0.0, 3.14159]))
print(bipartite_entanglement_entropy(ground_state, n_sites=4, subsystem=[0, 1]))
```

Pass `basis=sector.basis` to the same functions when the state comes from
`xxz_chain_sector` or `heisenberg_chain_sector`. Correlations and reduced
density matrices operate directly on the reduced state vector.

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

Additional plotting helpers include:

- `plot_spectrum`
- `plot_density`
- `plot_site_probabilities`
- `plot_hofstadter_butterfly`
- `plot_parameter_sweep`
- `plot_lattice_graph`
- `plot_lattice_state`
- `plot_hamiltonian_matrix`
- `apply_plot_style`

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
python examples/rice_mele_spectrum.py
python examples/hofstadter_butterfly.py
python examples/bose_hubbard_spectrum.py
python examples/haldane_spectrum.py
python examples/kagome_graph.py
python examples/hamiltonian_matrix.py
```
