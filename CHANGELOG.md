# Changelog

## v0.1.3

This release reorganizes the project documentation into a numbered learning
curriculum, expands the executed notebook coverage, and improves plot
readability for spectra, lattice structure, and complex Hamiltonians.

### Added

- Seven notebook workflows covering:
  - Spin observables and correlations
  - XY, XXZ, and J1-J2 spin-chain comparisons
  - Boundary conditions and finite-size effects
  - Aubry-Andre localization
  - Custom lattice construction and registration
  - Hamiltonian magnitude and phase structure
  - PennyLane export validation
- `plot_spectrum` options to highlight the ground-state gap and draw an
  energy-zero reference line.
- Sublattice coloring, legends, and dashed unit-cell outlines for
  `plot_lattice_graph`.
- Explicit $-\pi$, $0$, and $\pi$ ticks for phase colorbars.
- A colorblind-friendly plotting palette.
- A Hamiltonian matrix example showing the magnitude and phase of a finite
  Haldane model.
- Theory diagrams for spin-chain couplings, Heisenberg ladders, Hubbard basis
  and fermionic-sign conventions, the Kitaev Nambu representation, SSH and
  Rice-Mele unit cells and boundary conditions, Aubry-Andre localization,
  two-dimensional lattice geometries, Hofstadter flux phases, and Haldane
  hoppings.
- A notebook-summary script and Makefile targets for rebuilding example and
  documentation assets.

### Changed

- Renamed and ordered the notebook collection as a 17-part learning curriculum.
- Rebuilt the executed notebooks, rendered HTML, output figures, results
  summary, usage guide, and Pages notebook index around the new curriculum.
- Standardized documentation equations on dollar-delimited LaTeX and defined
  theory variables at first use.
- Expanded the example gallery with clearer parameter captions and refreshed
  Ising, SSH, Hofstadter, kagome, and other model figures.
- Increased the Hofstadter example to an $8\times8$ lattice and 121 flux
  values.
- Moved reproducible example images to the top-level `images/` directory and
  updated the Pages workflow to stage them with the documentation site.
- Updated lattice and matrix plots to distinguish sublattices, unit cells,
  hopping magnitudes, and hopping phases more clearly.

### Fixed

- Masked zero matrix elements in phase plots so absent couplings are not
  displayed as zero-phase hoppings.
- Standardized phase color scales to the full $[-\pi,\pi]$ interval.

### Notes

- The notebook numbering defines the recommended reading order.
- Generated notebook HTML and figures remain published under `docs/`.

---

## v0.1.2

This release adds user-defined lattice models, makes model registration usable
at runtime, and expands the plotting and command-line interfaces for generic
lattice workflows.

### Added

- `Bond` and `Lattice` containers for finite user-defined lattice geometries.
- `TightBindingModel` plus dense and sparse `custom_tight_binding` builders.
- Support for per-bond complex matrix elements, Hermitian conjugate completion,
  onsite terms, inferred site counts, and lattice position metadata.
- Runtime `register_model` and `unregister_model` helpers.
- Registry entries for dense and sparse custom tight-binding models.
- Plotting helpers:
  - `plot_lattice_state` for probability- and phase-resolved lattice states
  - `plot_hamiltonian_matrix` for real, imaginary, magnitude, and phase views
  - `plot_parameter_sweep` for generic one-parameter spectra
  - `apply_plot_style` for consistent axes styling
- CLI `--bond` arguments for custom graph models.
- `.gitattributes` rules for LF text normalization, binary assets, and generated
  documentation classification.
- Tag-triggered PyPI publishing through GitHub Actions and PyPI Trusted
  Publishing, with version, formatting, lint, test, build, and distribution
  checks before upload.

### Changed

- The CLI now builds every supported model through `MODEL_REGISTRY`, uses each
  builder's registered defaults, and exposes all registered models as choices.
- `plot_lattice_graph` can read positions from Hamiltonian metadata and encode
  hopping phase or magnitude with edge color and width.
- Existing plotting helpers now apply a consistent visual style.
- `plot_hofstadter_butterfly` now delegates to the generic parameter-sweep
  implementation.
- The Hofstadter, kagome, and SSH examples use the generalized plotting API.
- Regenerated the example plot gallery with the updated plotting style.
- Exported the custom lattice types, builders, and registry mutation helpers
  from the package's public API.
- Expanded `README.md` and `USAGE.md` with custom model, registration,
  visualization, and sparse-builder examples.
- Expanded `.gitignore` coverage for Python tooling, coverage reports, virtual
  environments, editors, operating-system files, and local generated outputs.

### Fixed

- Wrapped long notebook source lines so full-repository Black and Ruff checks
  pass without changing the saved notebook output.

### Notes

- Two-item bonds use `-hopping` as their matrix element; three-item bonds use
  the supplied value directly.
- Custom tight-binding builders remain single-particle Hamiltonians rather than
  many-body Fock-space models.
- The optional PennyLane export test is exercised when the `pennylane` extra is
  installed.

---

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

- Dense spin-chain matrices still scale as $2^N$ for $N$ sites.
- Hubbard models remain small-system exact diagonalization tools even when sparse
  builders are used.
- SSH, Rice-Mele, tight-binding, Haldane, Hofstadter, triangular, kagome, and AAH
  builders return single-particle matrices.
- The Kitaev-chain builder returns a BdG matrix, not a many-body Hamiltonian.
