# Roadmap

This document lists future work only. Implemented capabilities belong in the
README, changelog, validation matrix, and model reference.

The package will continue to focus on portable lattice-model data, lightweight
reference analysis, and optional ecosystem adapters rather than becoming a
general-purpose many-body solver or workflow-execution framework.

This roadmap describes intended direction, not a compatibility guarantee.
Priorities may change as the package is exercised.

## Product Scope

Future development is organized into three layers.

### Core model and data layer

Future core work includes:

- A unified physical-system representation connecting geometry, local degrees
  of freedom, basis ordering, and interaction terms
- Periodic unit-cell geometry, reciprocal-space data, and Bloch Hamiltonians
- Additional symmetry-reduced Hamiltonian construction
- Richer matrix imports and optional model-specification formats
- Spatial transformations such as interfaces, twisted boundaries, and
  long-range bond generation
- Portable result and visual-artifact descriptions
- Stable plugin contracts for third-party model families and data adapters

Core features should preserve enough information for a model to be reconstructed
and interpreted in a different process or external package.

### Lightweight reference analysis

Small, dependency-light tools are in scope when they help users inspect,
validate, compare, or demonstrate a model:

- Band structures and topological invariants
- Solver safeguards, convergence reporting, and conserved-quantity diagnostics
- Additional occupations, correlations, currents, and entanglement support
- Time-independent evolution, compact quench workflows, and expectation-value
  time series
- Small parameter sweeps, finite-size studies, and thermal reference
  calculations
- Reproducible lattice diagrams, interaction graphs, matrix views, and
  analysis plots

These capabilities should remain transparent and modest in scope. Where
specialized packages provide stronger solvers or workflows, this package should
export a well-described model to them rather than duplicate their functionality.

### Optional ecosystem adapters

Integrations should be optional and translate through the stable model
representation:

- Structure formats such as ASE, pymatgen, and selected CIF workflows
- Visual and tabular formats such as SVG, Graphviz DOT, CSV, and plot-data JSON
- Operator packages such as OpenFermion and Qiskit
- Solver ecosystems such as QuSpin, NetKet, or QuTiP where conventions map
  explicitly
- Third-party model, importer, exporter, and analysis plugins

## Guiding Principles

- Treat geometry, local degrees of freedom, parameters, basis conventions,
  symmetries, units, references, and provenance as first-class data.
- Use one validated representation across creation, import, construction,
  analysis, and export.
- Keep physical data and declarative visual metadata portable; treat rendered
  images as derived artifacts.
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

## Priority 1: Extend the Unified Physical-System Representation

Portable models should explicitly connect physical structure to Hamiltonian
indices and local state spaces. This representation should work for spin,
tight-binding, Hubbard, bosonic, and BdG systems without treating their bases
as interchangeable.

- Add carefully scoped multi-site interaction records beyond current onsite and
  two-body terms.
- Add richer multi-orbital site semantics, orbital geometry, and local
  constraints.
- Add explicit conserved quantities and symmetry actions on local degrees.
- Define stable mappings for reduced bases without confusing sector rows with
  full-basis or local-factor indices.
- Preserve physical-system records through additional importers and ecosystem
  adapters.

### Completion criteria

A user should be able to:

1. Import or create a spin, tight-binding, Hubbard, bosonic, or BdG model and
   identify what each basis component physically represents.
2. Map geometry sites, orbitals, local degrees of freedom, and interaction
   terms to Hamiltonian indices without model-specific assumptions.
3. Preserve labels and conventions through construction, persistence,
   visualization, and export.
4. Use the same portable representation as the input to generic diagram,
   graph, and analysis tools.

## Priority 2: Extend Model Import and Interchange

- Add optional YAML model specifications.
- Add schema migrations only when future schema revisions require them.

## Priority 3: Additional Lattice Transformations

Finite lattices are particularly useful when translation symmetry is broken.
Transformations should operate on portable lattice/model data rather than
creating unrelated builder functions.

- General interfaces and continuously spatially varying parameters
- Geometry-aware defect templates beyond site removal
- Model-level transformations that update interactions alongside geometry

Lightweight ensemble summaries may report means, variances, and realization
metadata, but distributed ensemble execution is outside the core package.

## Priority 4: Visualization and Diagram Workflows

Visualization should be a lightweight client of portable model and result data,
not a parallel source of physical truth.

### Geometry and interaction diagrams

- Lattice geometry and interaction graphs for both particle and spin models
- Site, orbital, sublattice, unit-cell, and local-degree-of-freedom labels
- Bond and interaction styling by type, direction, phase, strength, range, or
  spin axis
- Unit-cell outlines, primitive vectors, boundary links, defects, impurities,
  interfaces, and disorder annotations
- Spin textures and local expectation-value arrows where state data is
  available
- Reusable layout controls that do not alter physical geometry

### Reciprocal-space and matrix views

- Reciprocal lattices, Brillouin zones, high-symmetry points, and momentum paths
- Hamiltonian magnitude, phase, sparsity, block structure, and basis labels
- Clear distinction between physical adjacency diagrams and matrix-connectivity
  diagrams

### Analysis plots and output formats

- Standard plots for spectra, bands, density of states, occupations,
  correlations, entanglement, dynamics, and parameter sweeps
- PDF figure output
- Graphviz DOT and styled NetworkX export
- Deterministic visual defaults suitable for documentation and regression
  testing

Physical data and declarative visual metadata may be portable. Rendered image
bytes should remain derived artifacts rather than being embedded in
`ModelSpec`.

## Priority 6: Extend Result and Visual-Artifact Records

Portable analysis records now cover spectra, bands, topology, and selected spin
observables with JSON or NPZ persistence, bundle integration, and plot
regeneration. Future extensions include:

- Solver residuals, convergence histories, and iterative-status conventions
- Generated table, figure, diagram, and plot-data filename references
- Multi-panel declarative plot specifications and styling validation
- Result-schema migrations when future revisions require them

## Priority 7: Reference Analysis

Reference analysis should stay directly connected to model inspection,
validation, and interchange.

### Spectral and solver safeguards

- Add commutator-based conserved-quantity diagnostics.
- Extend iterative solver reporting with backend-specific iteration counts.

This should remain a compact interface, not a replacement for specialized
many-body solver packages.

### Observables and entanglement

- Extend reduced-basis-aware observables as new sectors are introduced.

### Extended topological analysis

- Local Berry curvature and convergence reports
- Additional multi-band and symmetry-aware invariants

### Lightweight dynamics

- Time-dependent Hamiltonians and higher-order propagation controls
- Additional state-distance and dynamical-order diagnostics

### Lightweight sweeps and thermal tools

- Gap-closing and extrema detection
- Susceptibility helpers with explicit conjugate fields
- Broadened density-of-states and simple single-particle spectral functions

Checkpointing, distributed execution, and general workflow orchestration should
be delegated to external tools.

## Priority 8: Symmetry-Reduced Bases

Reduced bases should be added where the conserved quantity is explicit and the
mapping to the full basis is testable.

- Extend fixed-magnetization support to compatible graph and ladder models.
- Add particle-number and spin sectors for Fermi-Hubbard models.
- Consider spin-inversion or parity sectors where they offer a clear benefit.
- Consider translation sectors only after periodic geometry and momentum
  conventions stabilize.

New sectors must preserve labels, basis-state mappings, and resource estimates,
and their spectra and observables must be validated against full-space blocks.

## Model Families

New named builders should provide a clear benchmark or exercise reusable
infrastructure.

Near-term candidates:

- Plain honeycomb or graphene lattice
- Anderson-disordered two-dimensional lattices

Candidates after periodic and topological support:

- Checkerboard or related two-band Chern-insulator model
- Dice or \(T_3\) lattice
- Kitaev honeycomb spin model

Research-oriented geometries such as pyrochlore lattices, quasicrystal
approximants, or general multi-orbital structures should be added only with
documented conventions, reference checks, and a concrete interchange use case.

## Ecosystem Interoperability and Plugins

### Planned adapters

- ASE or pymatgen structure import, including selected CIF workflows
- XYZ and selected structure/coordinate formats where model semantics remain
  explicit
- SVG and Graphviz DOT diagram export
- Styled NetworkX graph export
- CSV and JSON export for observables, coordinates, and plot data
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

Future maturity work:

- Gauge-equivalent constructions where applicable
- Type-complete public API and `py.typed`
- Coverage expectations and representative benchmarks
- Installation and CLI tests from built wheels

## Near-Term Backlog

Recommended implementation order:

1. Add reciprocal-space and Berry-curvature diagrams.
2. Add compact dense/sparse dynamics and quench helpers with expectation-value
   plots.
3. Add fixed-particle-number Hubbard sectors.
4. Add lightweight sweeps, finite-size summaries, and machine-readable plot
    data.
5. Add small-system thermal observables and broadened spectral helpers.
6. Add high-value benchmark models only when their conventions, diagrams, and
    validation are ready.

## Explicit Non-Goals

The project is not intended to:

- Compete with large-scale tensor-network, quantum Monte Carlo, or specialized
  many-body simulation frameworks.
- Provide distributed workflow orchestration, cluster scheduling, or a general
  resumable execution engine.
- Hide exponential Hilbert-space scaling from users.
- Infer a physically complete model from atomic coordinates alone.
- Treat a visually plausible diagram as sufficient physical model data.
- Promise lossless conversion to every external package or basis convention.
- Treat single-particle, BdG, spin, bosonic, and fermionic matrices as
  interchangeable.
- Claim quantum advantage from small exact-diagonalization examples.

New features should preserve these distinctions and state their basis,
conventions, dimensions, metadata behavior, and computational limits
explicitly.
