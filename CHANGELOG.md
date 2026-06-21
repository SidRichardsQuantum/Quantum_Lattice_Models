# Changelog

## v0.1.7

This release adds a portable periodic-lattice layer, generic Bloch and band
analysis, deterministic visual-data exports, spatial lattice transformations,
reference topological invariants, and reproducible analysis-result records.

### Added

- Versioned periodic unit cells with primitive and reciprocal vectors, orbital
  positions, cell-displacement bonds, onsite terms, labels, and JSON
  persistence.
- Generic cell- and orbital-gauge Bloch Hamiltonians, momentum paths, band
  energies, eigenvectors, plotting, and CSV or JSON export.
- Periodic SSH, Rice-Mele, square, honeycomb, kagome, and Haldane constructors.
- Finite supercell expansion with independently selectable periodic axes.
- Deterministic SVG and plot-data exports for lattices, periodic diagrams,
  physical interaction graphs, and band structures.
- Impurity and domain-wall potentials, twisted boundaries, and distance-based
  power-law bond generation.
- Discretized Zak phases, two-band chiral winding numbers, and
  Fukui-Hatsugai-Suzuki Chern numbers.
- CLI workflows for periodic model creation, bands, topology, and SVG export.
- Analytic dispersion and topological-phase validation for SSH, square, and
  Haldane models.
- Versioned `AnalysisResult` records carrying source identity, parameters,
  coordinates, numerical arrays, solver details, warnings, declarative plot
  metadata, package version, provenance, and optional creation time.
- Deterministic readable JSON and compressed self-contained NPZ
  analysis-result persistence.
- Portable result producers for spectra, band structures, Zak phases, winding
  and Chern numbers, and spin magnetization or correlation observables.
- Hamiltonian bundles with optional `analyses/` records and manifest entries.
- The `inspect-result`, `export-result`, and `plot-result` commands, plus
  `--result-output` support for spectrum, band, and topology workflows.
- Plot regeneration from stored spectrum, band, scalar-topology, and observable
  results without reconstructing the source Hamiltonian.
- Structured dense and sparse eigensolver selection with sparse-to-dense
  safeguards, residuals, convergence metadata, degeneracy groups, and memory
  estimates.
- Gauge-invariant Berry-curvature meshes, occupied-subspace Wilson loops,
  reciprocal-space records, and corresponding heatmap or reciprocal diagrams.
- Dense and sparse time evolution, quench workflows, expectation-value time
  series, and Loschmidt amplitudes and echoes.
- One- and two-parameter sweeps, finite-size studies, extrema metadata, and
  reusable line or heatmap result plots.
- Fixed-particle-number Bose-Hubbard bases, sparse sector Hamiltonians,
  embedding, and sector-aware site occupations.
- Stable canonical thermal reference calculations for partition functions,
  internal and free energy, entropy, heat capacity, and thermal expectations.
- Single-particle occupations, bond currents, mixed-axis spin correlations,
  total-spin diagnostics, and broadened local density of states.
- Analysis SVG/PDF rendering and matrix real, imaginary, magnitude, phase, and
  sparsity data export.
- Anderson, long-range tight-binding, Creutz-ladder, sawtooth, Lieb-lattice,
  XYZ-chain, and random-field Heisenberg benchmark models.

### Changed

- Updated project version metadata to `0.1.7`.

## v0.1.6

This release expands portable interchange, adds fixed-magnetization spin
sectors and reduced-basis observables, and introduces registry-driven discovery,
presets, resource inspection, and model comparison workflows.

### Added

- Paired site/bond CSV lattice import and export with metadata sidecars.
- Optional NetworkX and GraphML round-trip adapters preserving directed bonds,
  duplicate edges, complex matrix elements, labels, and portable metadata.
- Immutable lattice transformations for relabeling, subgraphs, vacancies,
  bond addition and removal, boundary conditions, and reproducible onsite or
  bond disorder.
- Standardized lattice units, conventions, references, and ordered provenance,
  plus model convention metadata.
- External dense NPY and dense or CSR NPZ Hamiltonian import with explicit
  basis, optional geometry, units, conventions, references, and provenance
  metadata.
- Strict imported-matrix validation for shape, numeric finite values,
  Hermiticity, and optional basis dimension.
- Portable `external_matrix` results that support persistence, inspection,
  spectra, validation, and re-export without claiming builder reconstruction.
- The `quantum-lattice import-matrix` command.
- Fixed total Pauli-$Z$ magnetization bases with explicit computational-basis
  state lists, state-to-index mappings, projection, and embedding helpers.
- Reduced sparse XXZ and magnetization-conserving Heisenberg-chain builders.
- Registered sector builders with binomial dimension estimates and portable
  sector metadata in `HamiltonianResult` and NPZ persistence.
- Full-space block and spectrum validation for open, periodic, positive, zero,
  and negative magnetization sectors.
- Site-resolved and total Pauli-$Z$ magnetization.
- Full and connected $X$, $Y$, and $Z$ correlation matrices and static spin
  structure factors.
- Reduced density matrices and bipartite von Neumann entropy for full and
  fixed-magnetization state vectors.
- Product-state, Bell-state, and reduced-basis embedding validation for spin
  observables and entanglement.
- Model discovery filters for category, basis, sparse capability, and
  validation status.
- Named presets for canonical phases and analytic reference limits.
- Construction-free dry-run reports with dimensions, dense-memory estimates,
  representations, basis information, and warnings.
- Deterministic JSON output for discovery, presets, validation, spectra, dry
  runs, and model comparisons.
- Model comparison summaries covering parameters, matrices, spectra, and
  spectral gaps with configurable construction safeguards.
- Selective model, lattice, matrix, metadata, and deterministic directory-bundle
  export through the public API and `quantum-lattice export --artifact`.
- Portable `LocalDegreeOfFreedom`, `BasisIndexMapping`, and `InteractionTerm`
  records for explicit physical sites, local state spaces, basis ordering, and
  onsite or two-body operators.
- Automatically populated physical-system data for Ising, Heisenberg, XXZ,
  SSH, and custom tight-binding model specifications.
- Interaction-graph plotting directly from portable model specifications,
  including spin axes, hopping links, coupling strengths, site labels, and
  onsite terms.
- Physical-system records for Bose-Hubbard, spinful Fermi-Hubbard, and
  Kitaev-chain BdG specifications, including truncated boson factors,
  spin-resolved fermionic modes, particle/hole ordering, interactions, chemical
  potentials, hopping, and pairing.
- Portable `graph_spin` specifications created with
  `create_graph_spin_spec`, preserving arbitrary two-site Pauli interactions,
  local fields, positions, labels, and dense or sparse reconstruction.
- Multi-mode interaction diagrams that offset spin, fermionic, or Nambu
  components sharing one physical site.

### Changed

- Updated the project version metadata to `0.1.6`.
- Merged the redundant CLI plot walkthrough into the registry notebook,
  preserving contiguous curriculum numbering while combining discovery,
  portable files, spectra, plotting, and artifact bundles.
- Replaced five one-call spectrum examples with workflow examples for portable
  physical interaction graphs, deterministic model bundles, and external
  matrix import; retained the distinct density, state, sweep, geometry, and
  matrix-inspection examples.
- Expanded the executed notebook curriculum from 16 to 21 workflows with
  dedicated notebooks for physical-system records, model and matrix
  interchange, graph-spin models, fixed-magnetization sectors, and lattice
  import and transformations.
- Fixed notebook 16 to avoid Python 3.12-only nested f-string syntax, retaining
  compatibility with the package's Python 3.10 minimum.
- Reorganized the roadmap around three explicit scope layers: core model data
  and interchange, lightweight reference analysis, and optional ecosystem
  adapters.
- Refocused the near-term roadmap on periodic unit-cell geometry, Bloch
  Hamiltonians, topological analysis, dynamics, additional symmetry sectors,
  sweeps, and thermal reference tools.
- Expanded `inspect`, `validate`, `spectrum`, and `export` file workflows to
  accept portable Hamiltonian files where applicable.
- Added schema compatibility documentation and an end-to-end importing and
  transformation guide.

### Fixed

- Portable Hamiltonian loading now rejects matrix shape, representation, basis,
  and external-dimension metadata mismatches.
- Hamiltonian metadata JSON now normalizes NumPy arrays, NumPy scalars, tuples,
  and complex values before persistence.

## v0.1.5

This release completes the first portable Hamiltonian workflow and introduces
shared sparse infrastructure for the package's spin-$1/2$ models.

### Added

- Explicit model and lattice schema migration behavior, including support for
  loading legacy unversioned specifications and clear rejection of unsupported
  future versions.
- Stricter specification validation for fields, types, coordinates, bonds,
  boundary conditions, and portable metadata values.
- `HamiltonianResult` for carrying dense or sparse matrix data with its model
  specification, basis, representation, and construction metadata.
- Self-contained NPZ persistence for dense and CSR sparse Hamiltonians.
- Dense NPY persistence with deterministic JSON metadata sidecars.
- Public `save_hamiltonian`, `load_hamiltonian`, and `metadata_path` helpers.
- File-oriented `quantum-lattice spectrum MODEL.json` and
  `quantum-lattice export MODEL.json` workflows.
- `SpinInteraction`, `SpinField`, and graph-spin builders supporting arbitrary
  two-site Pauli-axis interactions and site-resolved fields.
- Sparse variants of the Ising, XY, XXZ, Heisenberg, J1-J2, and Heisenberg
  ladder builders.
- Dense/sparse equivalence, persistence round-trip, schema migration, strict
  validation, and file-oriented CLI tests.

### Changed

- Updated the project version metadata to `0.1.5`.
- Existing dense spin builders now use the same graph interaction
  specification and sparse assembly backend as their sparse counterparts.
- Registered sparse spin builders for model discovery, portable
  reconstruction, and CLI use.
- Updated the README and roadmap to describe the portable result, persistence,
  CLI, and sparse spin workflows.

## v0.1.4

This release adds a capability-based development roadmap and removes several
sources of redundant package maintenance and dense-memory use.

### Added

- `ROADMAP.md` describing planned work for:
  - Portable lattice and model specifications
  - Validation, serialization, and import/export
  - Metadata-preserving dense and sparse Hamiltonians
  - Registry-driven parameter schemas and file-oriented CLI workflows
  - Periodic lattices, Bloch Hamiltonians, and band structures
  - General observables, topological analysis, and dynamics
  - Scientific verification and documented physical conventions
  - Reproducible experiment and result records
  - Parameter sweeps and finite-size studies
  - Symmetry-reduced Hilbert spaces
  - Solver interfaces and computational safeguards
  - Optional ecosystem adapters and plugin discovery
  - API stability, testing, typing, documentation, and benchmarking
- Completion criteria, implementation priorities, engineering work, and
  explicit non-goals for future development.
- Expanded roadmap coverage for sparse and graph-based spin construction,
  symmetry sectors, spin observables, entanglement, dynamics, and prioritized
  spin-$1/2$ and spin-1 model families.
- Added prioritized lattice candidates including plain honeycomb, Lieb,
  sawtooth, Creutz-ladder, Anderson-disordered, long-range, Chern-insulator,
  dice, pyrochlore, quasicrystal, and multi-orbital models.
- Expanded roadmap coverage for thermal and response analysis, reproducible
  disorder ensembles, defects and interfaces, model presets and comparison,
  dry-run resource inspection, batch and machine-readable CLI workflows, and
  task-oriented user recipes.
- Typed `ParameterInfo` registry metadata with parameter defaults,
  descriptions, CLI names, repeatability, and basic constraints.
- Registry-generated CLI arguments, including previously unavailable options
  such as `--pairing`, `--phi`, and `--sublattice-potential`.
- Public computational diagnostics for registered-model dimensions, dense
  memory estimates, matrix density and storage, Hermiticity, and particle-hole
  spectral symmetry.
- `VALIDATION.md`, documenting the initial scientific validation matrix and its
  corresponding automated tests.
- Analytic-limit tests for the zero-field Ising chain, uniform open
  tight-binding chain, decoupled SSH dimers, and noninteracting single-site
  Hubbard models.
- Versioned `ModelSpec` and `LatticeSpec` dataclasses for portable registered
  model parameters and finite-lattice geometry.
- Canonical JSON save/load support with explicit complex-number encoding.
- Dense or sparse Hamiltonian reconstruction from loaded model
  specifications.
- CLI `create`, `inspect`, and `validate` commands for reusable model JSON
  files.
- Round-trip, schema validation, custom-lattice, dense/sparse reconstruction,
  and file-oriented CLI tests.
- A concise Markdown reference page for every implemented model family,
  including equations, structure, basis, scaling, package use, generated
  registry parameter tables, validation notes, and computational cautions.
- Generated model-reference HTML pages, a general theory webpage, MathJax
  rendering, a website model index, and freshness tests for generated
  documentation.

### Changed

- Updated the project version metadata to `0.1.4`.
- Added the roadmap to the documented repository structure.
- Reworked square, triangular, kagome, Harper-Hofstadter, Haldane, and custom
  tight-binding sparse builders to assemble SciPy sparse matrices directly.
  Their dense counterparts now reuse the same sparse construction path, keeping
  dense and sparse model semantics synchronized without allocating a dense
  matrix inside sparse APIs.
- Changed the top-level `quantum_lattice_models` exports to import from focused
  implementation modules directly. `quantum_lattice_models.models` remains the
  backwards-compatible re-export module.
- Derived built-in `ModelInfo.name` values from their registered builder names,
  removing duplicated model-name declarations from the registry.
- Consolidated internal dense conversion for NumPy and SciPy matrices into one
  shared helper used by spectral and plotting utilities.
- Migrated the CLI, examples, notebooks, rendered notebook pages, tests, and
  usage documentation to use `plot_spectrum` consistently.
- Standardized notebook Markdown cells on the equation and variable conventions
  used by `THEORY.md`: dollar-delimited LaTeX, mathematical notation for
  physical variables and formulas, and code formatting only for API names,
  Python identifiers, and CLI options.
- Improved saved notebook outputs with indexed eigenvalue tables, explicit
  dense/sparse error comparisons, compact registry summaries, grouped model
  listings, and labeled relative artifact paths.
- Refocused `THEORY.md` on shared basis, operator, boundary, sparse-matrix,
  exact-diagonalization, observable, symmetry, and reproducibility conventions;
  model-specific material now lives under `docs/models/`.

### Removed

- Removed the redundant `plot_lattice_spectrum` alias. Use `plot_spectrum` for
  spin, lattice, dense, and sparse Hamiltonians.

### Fixed

- Sparse lattice builders no longer construct a full dense matrix before
  converting it to CSR format.
- Added stronger dense/sparse equivalence coverage for onsite terms, periodic
  boundaries, and complex hopping phases.
- Added compact validation tests for public operator, observable, spectral,
  lattice, diagnostics, registry, and CLI failure paths.
- Reduced test-suite process startup overhead while retaining an end-to-end CLI
  module smoke test.

### Notes

- The roadmap is organized by capability rather than assigning planned work to
  specific package versions.
- Roadmap items describe future direction and are not implemented by this
  release unless documented elsewhere.
- Removing `plot_lattice_spectrum` is an API cleanup; downstream imports should
  be changed to `quantum_lattice_models.plotting.plot_spectrum`.

---

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
