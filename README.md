# Quantum Lattice Models

<p align="center">

<a href="https://pypi.org/project/quantum-lattice-models/">
<img src="https://img.shields.io/pypi/v/quantum-lattice-models?style=flat-square" alt="PyPI Version">
</a>

<a href="https://pypi.org/project/quantum-lattice-models/">
<img src="https://img.shields.io/pypi/pyversions/quantum-lattice-models?style=flat-square" alt="Python Versions">
</a>

<a href="https://github.com/SidRichardsQuantum/Quantum_Lattice_Models/actions/workflows/tests.yml">
<img src="https://img.shields.io/github/actions/workflow/status/SidRichardsQuantum/Quantum_Lattice_Models/tests.yml?label=tests&style=flat-square" alt="Tests">
</a>

<a href="https://sidrichardsquantum.github.io/Quantum_Lattice_Models/">
<img src="https://img.shields.io/github/actions/workflow/status/SidRichardsQuantum/Quantum_Lattice_Models/pages.yml?label=docs&style=flat-square" alt="Docs">
</a>

<a href="LICENSE">
<img src="https://img.shields.io/github/license/SidRichardsQuantum/Quantum_Lattice_Models?style=flat-square" alt="License">
</a>

<a href="https://github.com/sponsors/SidRichardsQuantum">
<img src="https://img.shields.io/badge/sponsor-GitHub-ea4aaa?style=flat-square&logo=githubsponsors" alt="Sponsor">
</a>

</p>

Quantum Lattice Models is a lightweight, package-first Python library for constructing, analyzing, plotting, and exporting small lattice Hamiltonians used in physics workflows and quantum algorithm research prototypes.

PyPI: [https://pypi.org/project/quantum-lattice-models/](https://pypi.org/project/quantum-lattice-models/)  

This repository is organized as an installable package first.
The real logic lives in `src/quantum_lattice_models/`; notebooks, scripts, and examples should stay thin and import the public package API.
`quantum_lattice_models.models` remains the compatibility import surface, while implementations are split across focused modules such as `spin`, `tight_binding`, `hubbard`, and `topological`.

## Implemented Models

- Transverse-field Ising spin chain
- Longitudinal-field Ising spin chain
- Next-nearest-neighbor Ising spin chain
- Anisotropic Heisenberg spin chain
- XY spin chain
- XXZ spin chain
- Frustrated J1-J2 Heisenberg spin chain
- Two-leg Heisenberg spin ladder
- Truncated Bose-Hubbard chain
- Spinful Fermi-Hubbard chain
- Kitaev-chain Bogoliubov-de Gennes matrix
- Su-Schrieffer-Heeger single-particle tight-binding model
- Rice-Mele single-particle chain
- Generic one-dimensional single-particle tight-binding chain
- Square-lattice single-particle tight-binding model
- Harper-Hofstadter square-lattice model
- Aubry-Andre-Harper quasiperiodic tight-binding chain
- Haldane honeycomb-lattice model
- Triangular-lattice single-particle tight-binding model
- Kagome-lattice single-particle tight-binding model

Spin-chain Hamiltonians are dense qubit-space matrices
Tight-binding Hamiltonians are single-particle matrices.
This distinction is intentional and explicit.

## Why Lattice Models Matter

Small lattice Hamiltonians are useful because they are concrete, inspectable testbeds.
They connect physics intuition to numerical linear algebra, and they give quantum algorithm researchers controlled problems for VQE, QPE, QSVT, spectral transforms, quantum walks, and simulation workflows.

This package does not claim quantum advantage.
It provides honest small-system tools for exact diagonalization, prototyping, teaching, and notebook-first experiments.

## Installation

From a local checkout:

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[dev]"
```

Minimal runtime dependencies are `numpy`, `scipy`, and `matplotlib`.

PennyLane export is optional:

```bash
pip install -e ".[pennylane]"
```

Notebook support is optional:

```bash
python -m pip install -e ".[notebooks]"
python -m ipykernel install --user --name quantum-lattice-models --display-name "Quantum Lattice Models"
```

## Quickstart

```python
from quantum_lattice_models.models import transverse_field_ising
from quantum_lattice_models.spectra import ground_energy, spectral_gap

H = transverse_field_ising(n_sites=4, j=1.0, h=0.5, periodic=False)

print(H.shape)
print(ground_energy(H))
print(spectral_gap(H))
```

```python
from quantum_lattice_models.models import ssh_model, ssh_edge_state_localizations
from quantum_lattice_models.spectra import eigensystem

H = ssh_model(n_cells=8, t1=0.5, t2=1.0, periodic=False)
values, vectors = eigensystem(H)
weights = ssh_edge_state_localizations(vectors, n_cells=8, edge_cells=2)
```

## Repository Structure

```text
src/quantum_lattice_models/  Package source
tests/                       Pytest test suite
examples/                    Command-line examples that save plots
notebooks/                   Thin-client exploratory notebooks
README.md                    Project overview
CHANGELOG.md                 Release notes
USAGE.md                     API examples
THEORY.md                    Model and method notes
RESULTS.md                   Generated results
```

Key package modules:

```text
spin.py                      Dense spin-chain and ladder builders
tight_binding.py             Single-particle tight-binding builders
hubbard.py                   Bose-Hubbard and Fermi-Hubbard builders
topological.py               Haldane, Hofstadter, and Kitaev builders
geometry.py                  Coordinate helpers for plotting
registry.py                  Structured model metadata
cli.py                       quantum-lattice command-line entry point
models.py                    Backwards-compatible re-export layer
```

## Notebooks as Thin Clients

Notebooks should import from `quantum_lattice_models` rather than defining their own model logic.
A notebook can choose parameters, run spectra, plot results, and tell a story.
The package should remain the source of truth.

Current notebooks:

- `notebooks/ising_spin_chains.ipynb`
- `notebooks/ssh_rice_mele_comparison.ipynb`
- `notebooks/hofstadter_flux_sweep.ipynb`
- `notebooks/hubbard_exact_diagonalization.ipynb`
- `notebooks/haldane_kagome_lattices.ipynb`
- `notebooks/model_registry_and_cli.ipynb`
- `notebooks/kitaev_bdg_symmetry.ipynb`
- `notebooks/heisenberg_ladder_spectrum.ipynb`
- `notebooks/sparse_dense_scaling.ipynb`
- `notebooks/cli_plot_walkthrough.ipynb`

## Development

Use the virtual environment for examples, notebooks, tests, and packaging commands.
The standard local checks are:

```bash
make format
make lint
make test
```

The `Makefile` runs Black one file at a time to avoid multi-file formatter stalls observed in some Codespace environments.

## Limitations / Truth Contract

- Dense spin-chain matrices scale as `2**n_sites` by `2**n_sites`.
- These tools are for small systems, education, exact diagonalization, and research prototypes.
- SSH, Rice-Mele, square, Harper-Hofstadter, Haldane, triangular, kagome, and generic tight-binding builders return single-particle matrices, not many-body Fock-space Hamiltonians.
- The Bose-Hubbard builder uses a truncated local occupation basis.
- The Fermi-Hubbard builder uses a dense occupation-number basis with explicit fermionic signs.
- The Kitaev-chain builder returns a Bogoliubov-de Gennes matrix, not a many-body Hamiltonian.
- Sparse builders are available for selected tight-binding and Hubbard chains, but exact diagonalization remains a small-system workflow.
- PennyLane is optional and only used when explicitly installed.
- The project is a backend for experiments, not a benchmark suite proving speedup or quantum advantage.

---

## Support development

If this repository is useful for research, learning, or experimentation, you can support continued development via GitHub Sponsors:

https://github.com/sponsors/SidRichardsQuantum

Sponsorship helps support ongoing work on open-source implementations of quantum algorithms, including improvements to documentation, reproducible workflows, and example notebooks.

Support helps maintain and expand practical tooling for variational quantum methods, quantum simulation workflows, and related experimentation.

## Citation

Sid Richards (2026)

Unified Variational and Phase-Estimation Quantum Simulation Suite

## Author

Sid Richards

- LinkedIn: [sid-richards-21374b30b](https://www.linkedin.com/in/sid-richards-21374b30b/)
- GitHub: [SidRichardsQuantum](https://github.com/SidRichardsQuantum)

## License

MIT. See [LICENSE](LICENSE).
