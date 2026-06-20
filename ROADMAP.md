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

## Priority 3: Lattice Transformations

Finite lattices are particularly useful when translation symmetry is broken.
Transformations should operate on portable lattice/model data rather than
creating unrelated builder functions.

- Impurity and boundary-potential transformations
- Repeated unit-cell construction
- Domain walls, interfaces, and spatially varying parameters
- Twisted boundary transformations where meaningful
- Long-range bonds generated by distance or power-law rules

Lightweight ensemble summaries may report means, variances, and realization
metadata, but distributed ensemble execution is outside the core package.

## Priority 4: Periodic Lattice Representation

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

## Priority 5: Visualization and Diagram Workflows

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
- SVG, PNG, and PDF figure output
- Graphviz DOT and styled NetworkX export
- CSV and JSON export for underlying plot coordinates and tabular results
- Deterministic visual defaults suitable for documentation and regression
  testing

Physical data and declarative visual metadata may be portable. Rendered image
bytes should remain derived artifacts rather than being embedded in
`ModelSpec`.

## Priority 6: Result and Visual-Artifact Records

Add a compact portable record that connects numerical and visual outputs to the
model and analysis that produced them.

- Source model identity or embedded portable specification
- Analysis name, parameters, coordinate labels, and numerical results
- Solver method, tolerances, residuals, convergence status, and warnings
- Declarative plot or diagram specifications, including selections, labels,
  styling, and annotations
- Generated table, figure, diagram, and plot-data filenames
- Package version, creation time, references, and provenance
- Deterministic JSON persistence and bundle integration

These records should support reproducibility and exchange without becoming a
general workflow engine or embedding large rendered images in model data.

## Priority 7: Reference Analysis

Reference analysis should stay directly connected to model inspection,
validation, and interchange.

### Spectral and solver safeguards

- Prevent accidental dense conversion of large sparse matrices.
- Select sensible dense or sparse eigensolvers while allowing explicit
  control.
- Return convergence status, tolerances, residuals, and whether a result is
  exact or iterative.
- Add pre-diagonalization cost and solver-selection estimates.
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

## Priority 8: Symmetry-Reduced Bases

Reduced bases should be added where the conserved quantity is explicit and the
mapping to the full basis is testable.

- Extend fixed-magnetization support to compatible graph and ladder models.
- Add fixed-particle-number sectors for Bose-Hubbard models.
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

1. Design periodic unit-cell and displacement-vector geometry.
2. Add reusable lattice and interaction-graph SVG and
   plot-data export.
3. Add Bloch Hamiltonians, band structures, reciprocal-space diagrams, and
   Brillouin-zone paths.
4. Add portable result and declarative visual-artifact records integrated with
   export bundles.
5. Add Berry/Zak phases and winding numbers, followed by Haldane Chern-number
   validation and corresponding plots.
6. Add compact dense/sparse dynamics and quench helpers with expectation-value
   plots.
7. Add fixed-particle-number Hubbard sectors.
8. Add lightweight sweeps, finite-size summaries, and machine-readable plot
    data.
9. Add small-system thermal observables and broadened spectral helpers.
10. Add high-value benchmark models only when their conventions, diagrams, and
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
