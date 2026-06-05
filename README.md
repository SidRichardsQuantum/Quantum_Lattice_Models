# Quantum Lattice Models

<p align="center">

<a href="https://pypi.org/project/quantum-lattice-models/">
<img src="https://img.shields.io/pypi/v/quantum-lattice-models?style=flat-square" alt="PyPI Version">
</a>

<a href="https://pypi.org/project/quantum-lattice-models/">
<img src="https://img.shields.io/pypi/pyversions/quantum-lattice-models?style=flat-square" alt="Python Versions">
</a>

<a href="https://github.com/SidRichardsQuantum/Variational_Quantum_Eigensolver/actions/workflows/tests.yml">
<img src="https://img.shields.io/github/actions/workflow/status/SidRichardsQuantum/Variational_Quantum_Eigensolver/tests.yml?label=tests&style=flat-square" alt="Tests">
</a>

<a href="https://sidrichardsquantum.github.io/Variational_Quantum_Eigensolver/">
<img src="https://img.shields.io/github/actions/workflow/status/SidRichardsQuantum/Variational_Quantum_Eigensolver/pages.yml?label=docs&style=flat-square" alt="Docs">
</a>

<a href="LICENSE">
<img src="https://img.shields.io/github/license/SidRichardsQuantum/Variational_Quantum_Eigensolver?style=flat-square" alt="License">
</a>

<a href="https://github.com/sponsors/SidRichardsQuantum">
<img src="https://img.shields.io/badge/sponsor-GitHub-ea4aaa?style=flat-square&logo=githubsponsors" alt="Sponsor">
</a>

</p>

Quantum Lattice Models is a lightweight, package-first Python library for constructing, analyzing, plotting, and exporting small lattice Hamiltonians used in physics workflows and quantum algorithm research prototypes.

PyPI: placeholder  
Website: placeholder

This repository is organized as an installable package first.
The real logic lives in `src/quantum_lattice_models/`; notebooks, scripts, and examples should stay thin and import the public package API.

## Implemented Models

- Transverse-field Ising spin chain
- Anisotropic Heisenberg spin chain
- Su-Schrieffer-Heeger single-particle tight-binding model
- Generic one-dimensional single-particle tight-binding chain

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
README.md                    Project overview
USAGE.md                     API examples
THEORY.md                    Model and method notes
RESULTS.md                   Placeholder for generated results
```

## Notebooks as Thin Clients

Notebooks should import from `quantum_lattice_models` rather than defining their own model logic.
A notebook can choose parameters, run spectra, plot results, and tell a story.
The package should remain the source of truth.

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
- The SSH and generic tight-binding builders return single-particle matrices, not many-body Fock-space Hamiltonians.
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
