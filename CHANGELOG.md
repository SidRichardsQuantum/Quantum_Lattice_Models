# Changelog

## v0.1.1

This release expands Quantum Lattice Models from a small set of dense toy
Hamiltonians into a broader package-first toolkit for spin chains, Hubbard
models, topological lattices, sparse workflows, examples, notebooks, and Pages
documentation.

### Added

- Spin models:
  - Longitudinal-field Ising chain
  - Next-nearest-neighbor Ising chain
  - XY chain
  - XXZ chain
  - J1-J2 Heisenberg chain
  - Two-leg Heisenberg ladder
- Hubbard and superconducting models:
  - Truncated Bose-Hubbard chain
  - Spinful Fermi-Hubbard chain with explicit fermionic signs
  - Kitaev-chain Bogoliubov-de Gennes matrix
- Single-particle and topological lattice models:
  - Rice-Mele chain
  - Square-lattice tight-binding model
  - Harper-Hofstadter square lattice
  - Aubry-Andre-Harper chain
  - Haldane honeycomb lattice
  - Triangular-lattice tight-binding model
  - Kagome-lattice tight-binding model
- Sparse builders for selected tight-binding, Hubbard, and topological models.
- Sparse spectral helpers:
  - `lowest_eigenvalues`
  - `ground_state`
- `LatticeHamiltonian` metadata array type for non-Pauli lattice matrices.
- Geometry helpers for square, triangular, honeycomb, and kagome plotting.
- Model registry with category, basis, dimension scaling, return type, defaults,
  and builder references.
- `quantum-lattice` command-line interface for listing models, printing spectra,
  and saving spectrum plots.
- Plotting helpers for spectra, densities, site probabilities, Hofstadter sweeps,
  and lattice graph connectivity.
- Command-line examples and generated plot gallery artifacts.
- Executed notebooks with saved outputs and rendered Pages HTML:
  - Ising spin chains
  - SSH and Rice-Mele comparison
  - Hofstadter flux sweep
  - Hubbard exact diagonalization
  - Haldane, triangular, and kagome lattices
  - Kitaev BdG symmetry
  - Heisenberg ladder spectra
  - Sparse vs dense scaling
  - Registry and CLI workflows

### Changed

- Split model implementations into focused modules:
  - `spin.py`
  - `tight_binding.py`
  - `hubbard.py`
  - `topological.py`
  - `_model_utils.py`
- Kept `quantum_lattice_models.models` as a backwards-compatible re-export
  layer.
- Made spectral utilities sparse-aware for sparse inputs.
- Expanded `README.md`, `USAGE.md`, `THEORY.md`, `RESULTS.md`, and GitHub Pages
  documentation.
- Updated notebook numerical outputs to use readable tables, rounded arrays, and
  labeled scalar summaries.

### Fixed

- Corrected repository links and badges to point at
  `SidRichardsQuantum/Quantum_Lattice_Models`.
- Added CI workflow definitions for tests and GitHub Pages deployment.

### Notes

- Dense spin-chain matrices still scale as `2**n_sites`.
- Hubbard models remain small-system exact diagonalization tools even when sparse
  builders are used.
- SSH, Rice-Mele, tight-binding, Haldane, Hofstadter, triangular, kagome, and AAH
  builders return single-particle matrices.
- The Kitaev-chain builder returns a BdG matrix, not a many-body Hamiltonian.
