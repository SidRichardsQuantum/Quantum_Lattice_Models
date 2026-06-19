# Roadmap

Quantum Lattice Models helps users create or import lattice models and preserve
the information required to understand, reconstruct, validate, analyze, and
exchange them.

```text
create or import -> validate -> inspect -> build -> analyze -> export or save
```

The package includes spin, tight-binding, Hubbard, and BdG models because these
are useful lattice-model representations with distinct bases and conventions.
It also includes lightweight analysis tools that make models immediately
inspectable and scientifically verifiable. It is not intended to become a
general-purpose many-body solver or workflow-execution framework.

This roadmap describes intended direction, not a compatibility guarantee.
Priorities may change as the model schema and import workflows are exercised.

## Product Scope

Development is organized into three layers.

### Core model and data layer

This is the package's primary responsibility:

- Lattice geometry, sites, orbitals, labels, unit cells, and bonds
- Model parameters, basis conventions, boundary conditions, and units
- Dense, sparse, and symmetry-reduced Hamiltonian construction
- Portable specifications, validation, schema migration, and provenance
- Import and export of model, lattice, matrix, and metadata formats
- Reusable transformations such as disorder, defects, boundaries, and
  relabeling
- Registry and plugin contracts for built-in and third-party model families

Core features should preserve enough information for a model to be reconstructed
and interpreted in a different process or external package.

### Lightweight reference analysis

Small, dependency-light tools are in scope when they help users inspect,
validate, compare, or demonstrate a model:

- Exact and low-energy spectra, gaps, eigenstates, and density of states
- Matrix diagnostics, conserved-sector checks, and resource estimates
- Occupations, magnetization, correlations, structure factors, localization,
  reduced density matrices, and entanglement
- Band structures and topological invariants
- Time-independent evolution, compact quench workflows, and expectation-value
  time series
- Small parameter sweeps, finite-size studies, and thermal reference
  calculations
- Plotting and deterministic tabular summaries

These capabilities should remain transparent and modest in scope. Where
specialized packages provide stronger solvers or workflows, this package should
export a well-described model to them rather than duplicate their functionality.

### Optional ecosystem adapters

Integrations should be optional and translate through the stable model
representation:

- Graph and structure formats such as NetworkX, GraphML, CSV, ASE, or pymatgen
- Operator packages such as OpenFermion, PennyLane, and Qiskit
- Solver ecosystems such as QuSpin, NetKet, or QuTiP where conventions map
  explicitly
- Third-party model, importer, exporter, and analysis plugins

## Guiding Principles

- Treat geometry, parameters, basis conventions, symmetries, units, references,
  and provenance as first-class data.
- Use one validated representation across creation, import, construction,
  analysis, and export.
- Keep dense, sparse, and reduced-basis behavior consistent.
- Prefer explicit physical conventions over implicit defaults.
- Keep core dependencies lightweight and place integrations behind optional
  extras.
- Add named models when they exercise reusable infrastructure or provide a
  useful analytic, geometric, or topological benchmark.
- Keep notebooks, examples, and command-line workflows as clients of the
  public package API.
- Report computational scaling honestly and avoid surprising dense
  conversions.

## Current Foundation

### Implemented through v0.1.4

- Typed model registry and generated command-line parameters
- Dimension, memory, storage, Hermiticity, and particle-hole diagnostics
- Scientific validation matrix and analytic reference checks
- Versioned `ModelSpec` and `LatticeSpec` JSON files
- Dense and sparse model reconstruction
- CLI `create`, `inspect`, and `validate` workflows
- Generated model-reference documentation

### Implemented in v0.1.5

- Explicit schema migration behavior and stricter portable-file validation
- Metadata-preserving `HamiltonianResult`
- Dense NPY and dense/sparse NPZ persistence
- File-oriented CLI `spectrum` and `export` workflows
- Shared graph-based dense and sparse spin construction
- Sparse Ising, XY, XXZ, Heisenberg, J1-J2, and ladder builders

### Implemented after v0.1.5

- Fixed-magnetization bases and reduced XXZ and Heisenberg Hamiltonians
- Explicit reduced-basis mappings, dimension estimates, metadata, and
  persistence
- Site and total magnetization
- Same-axis full and connected spin correlation matrices
- Static spin structure factors
- Reduced density matrices and bipartite entanglement entropy
- Full-space and reduced-basis cross-validation
- Model filtering by category, basis, sparse support, and validation status
- Named canonical-phase and reference-limit presets
- Pre-build dimension, memory, representation, and warning reports
- Deterministic JSON CLI output
- Model parameter, matrix, spectrum, and gap comparisons

## Priority 1: Complete Model Import and Interchange

The package is currently stronger at creating models than importing them. This
is the highest-priority gap.

Implemented:

- Published schema, compatibility, metadata, and migration policy
- End-to-end CSV and GraphML import/export documentation
- Paired site/bond CSV with metadata sidecars
- Optional NetworkX and GraphML adapters
- Immutable transformations for relabeling, subgraphs, vacancies, bonds,
  boundaries, and reproducible disorder
- Standardized lattice/model units, conventions, references, and provenance

### Schema and metadata

- Document the JSON schema and complex-number encoding.
- Define compatibility and deprecation policies.
- Add explicit migrations when the schema changes.
- Complete lattice metadata for site, orbital, sublattice, and unit-cell
  labels.
- Record physical units, references, gauge choices, basis ordering, and
  provenance consistently.
- Add validation status and convention metadata to registry entries.

### Tabular and graph import/export

- Import and export site coordinates and labels as CSV.
- Import and export bond or edge tables, including complex amplitudes.
- Support paired coordinate and edge-list files.
- Add NetworkX and GraphML adapters behind an optional extra.
- Preserve labels, edge attributes, directedness, Hermiticity conventions, and
  user metadata during round trips.

### Matrix and model persistence

- Extend round-trip coverage for NPY and NPZ.
- Add optional YAML model specifications.
- Define how imported matrices are associated with basis and geometry metadata.
- Add deterministic summaries and machine-readable CLI output.
- Add file-oriented `import` and richer `export` commands.

### Completion criteria

A user should be able to:

1. Create a built-in model or import a custom lattice.
2. Attach parameters, geometry, conventions, labels, units, and provenance.
3. Validate and inspect the model before construction.
4. Save and reload it in another process.
5. Build an equivalent dense, sparse, or supported reduced Hamiltonian.
6. Run reference analysis or export it to another tool without reconstructing
   missing context manually.

## Priority 2: Lattice Transformations

Finite lattices are particularly useful when translation symmetry is broken.
Transformations should operate on portable lattice/model data rather than
creating unrelated builder functions.

- Reproducible onsite and bond disorder with explicit distributions and seeds
- Vacancies, impurities, modified bonds, and boundary potentials
- Relabeling, subgraph extraction, and repeated unit-cell construction
- Domain walls, interfaces, and spatially varying parameters
- Open, periodic, and twisted boundary transformations where meaningful
- Long-range bonds generated by distance or power-law rules
- Transformation provenance and deterministic serialization

Lightweight ensemble summaries may report means, variances, and realization
metadata, but distributed ensemble execution is outside the core package.

## Priority 3: Periodic Lattice Representation

Periodic geometry should be represented as data rather than encoded separately
inside every model builder.

- Primitive lattice vectors and reciprocal vectors
- Basis sites and orbitals within a unit cell
- Bonds with cell-displacement vectors
- Clear distinction between finite periodic real-space matrices and
  translationally invariant models
- Bloch Hamiltonians \(H(\mathbf{k})\)
- Standard and user-defined Brillouin-zone paths
- Band structures and k-resolved eigenvectors

The representation should support SSH, Rice-Mele, square, triangular,
honeycomb, kagome, Haldane, and related user-defined systems.

## Priority 4: Reference Analysis

Reference analysis should stay directly connected to model inspection,
validation, and interchange.

### Spectral and solver safeguards

- Prevent accidental dense conversion of large sparse matrices.
- Select sensible dense or sparse eigensolvers while allowing explicit
  control.
- Return convergence status, tolerances, residuals, and whether a result is
  exact or iterative.
- Add dry-run estimates before expensive construction or diagonalization.
- Add degeneracy and commutator-based conserved-quantity diagnostics.

This should remain a compact interface, not a replacement for specialized
many-body solver packages.

### Observables and entanglement

- Add mixed-axis spin correlation matrices.
- Add total-spin \(S^2\) diagnostics where conventions permit.
- Add site and orbital occupations for supported occupation bases.
- Add current operators for supported tight-binding models.
- Extend reduced-basis-aware observables as new sectors are introduced.

### Topological analysis

- Berry and Zak phases
- One-dimensional winding numbers
- Wilson loops where applicable
- Chern numbers for supported Bloch models
- Gauge-tolerant numerical routines with convergence controls
- Reference tests against known SSH and Haldane phases

### Lightweight dynamics

- Time-independent dense and sparse state evolution
- Parameter quenches and expectation-value time series
- Loschmidt amplitudes and echoes
- Structured results with times, tolerances, model metadata, and solver details

### Lightweight sweeps and thermal tools

- One- and two-parameter sweeps returning labeled coordinates
- Finite-size summaries and standard plots
- Gap-closing and extrema detection
- Canonical partition functions and thermal expectation values for small
  systems
- Free energy, entropy, heat capacity, and susceptibility
- Broadened density-of-states and simple single-particle spectral functions

Checkpointing, distributed execution, and general workflow orchestration should
be delegated to external tools.

## Priority 5: Symmetry-Reduced Bases

Reduced bases should be added where the conserved quantity is explicit and the
mapping to the full basis is testable.

- Extend fixed-magnetization support to compatible graph and ladder models.
- Add fixed-particle-number sectors for Bose-Hubbard models.
- Add particle-number and spin sectors for Fermi-Hubbard models.
- Consider spin-inversion or parity sectors where they offer a clear benefit.
- Consider translation sectors only after periodic geometry and momentum
  conventions stabilize.
- Preserve sector labels, basis-state mappings, and resource estimates.
- Validate sector spectra and observables against full-space blocks.

## Model Families

New named builders should provide a clear benchmark or exercise reusable
infrastructure.

Near-term candidates:

- Plain honeycomb or graphene lattice
- Lieb lattice and its flat band
- Sawtooth or delta chain
- Creutz ladder
- Anderson-disordered chains and lattices
- Configurable long-range tight-binding chains
- XYZ and random-field spin chains

Candidates after periodic and topological support:

- Checkerboard or related two-band Chern-insulator model
- Dice or \(T_3\) lattice
- Kitaev honeycomb spin model

Research-oriented geometries such as pyrochlore lattices, quasicrystal
approximants, or general multi-orbital structures should be added only with
documented conventions, reference checks, and a concrete interchange use case.

## Ecosystem Interoperability and Plugins

### Planned adapters

- NetworkX and GraphML
- ASE or pymatgen structure import, including selected CIF workflows
- OpenFermion operator export
- Expanded PennyLane and Qiskit export
- QuSpin, NetKet, or QuTiP adapters where basis mappings are explicit

### Plugin discovery

- Python entry points for third-party model packages
- Registration of builders, parameter schemas, importers, exporters, and
  lightweight analysis routines
- Compatibility checks for plugin and schema versions
- Explicit, diagnosable plugin loading

## Scientific Verification and Project Maturity

Every built-in model and format should document and test its conventions.

- Basis ordering and local-state conventions
- Operator normalization and signs
- Gauge choices and boundary conditions
- Parameter units and references
- Analytic limits and known spectra or eigenstates
- Dense/sparse/reduced equivalence where applicable
- Gauge-equivalent constructions where applicable
- Import/export round trips and malformed-input failures
- Supported Python-version policy
- Type-complete public API and `py.typed`
- Coverage expectations and representative benchmarks
- Installation and CLI tests from built wheels
- Reproducible release and generated-documentation checks

## Near-Term Backlog

Recommended implementation order:

1. Design periodic unit-cell and displacement-vector geometry.
2. Add Bloch Hamiltonians and band-structure workflows.
3. Add Berry/Zak phases and winding numbers, followed by Haldane Chern-number
   validation.
4. Add compact dense/sparse dynamics and quench helpers.
5. Add fixed-particle-number Hubbard sectors.
6. Add lightweight sweeps and finite-size summaries.
7. Add small-system thermal observables and broadened spectral helpers.
8. Add high-value benchmark models only when their conventions and validation
    are ready.

## Explicit Non-Goals

The project is not intended to:

- Compete with large-scale tensor-network, quantum Monte Carlo, or specialized
  many-body simulation frameworks.
- Provide distributed workflow orchestration, cluster scheduling, or a general
  resumable execution engine.
- Hide exponential Hilbert-space scaling from users.
- Infer a physically complete model from atomic coordinates alone.
- Promise lossless conversion to every external package or basis convention.
- Treat single-particle, BdG, spin, bosonic, and fermionic matrices as
  interchangeable.
- Claim quantum advantage from small exact-diagonalization examples.

New features should preserve these distinctions and state their basis,
conventions, dimensions, metadata behavior, and computational limits
explicitly.
