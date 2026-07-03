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
- Additional symmetry-reduced Hamiltonian construction
- Richer matrix imports and optional model-specification formats
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
- Add symmetry actions on local degrees beyond the current named conserved
  quantities and commutator diagnostics.
- Extend shared reduced-basis mappings to future parity and translation
  sectors.
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

- Add schema migrations only when future schema revisions require them.
- Preserve richer orbital and interaction semantics through structure imports.
- Add portable adapter capability reports before translating a model.

## Model Intake and Authoring Experience

Users should be able to import, create, inspect, validate, and translate lattice
models through a small number of predictable workflows. Future work in this area
should make the package more useful as a model-intake layer rather than only a
collection of named Hamiltonian builders.

### Canonical physical-model summaries

- Provide a user-facing summary for imported or constructed models, including
  sites, orbitals, local degrees of freedom, basis ordering, terms, geometry,
  boundary conditions, symmetries, sectors, provenance, and reconstruction
  limits.
- Make the normalized physical-system representation easy to access from every
  supported model family, including spin, tight-binding, Hubbard, bosonic, and
  BdG systems.
- Clearly distinguish reconstructable model specifications, external imported
  matrices, lossy ecosystem translations, analysis-only result records, and
  rendered visual artifacts.

### Import diagnostics and model linting

- Report what each importer inferred, what metadata was missing, which
  conventions were assumed, and which downstream operations may be unsafe or
  lossy.
- Add model-linting checks for Hermiticity, duplicate or missing labels,
  inconsistent bond directions, invalid local Hilbert spaces, basis-dimension
  mismatches, advertised symmetry violations, and dense-construction resource
  risks.
- Suggest relevant validation checks after import, such as connected-component
  inspection, conserved-quantity commutators, sector block validation, or
  reference spectra where applicable.

### Creation recipes and graph-spin workflows

- Add higher-level recipe constructors for common lattice-model authoring
  tasks, while preserving explicit physical conventions and portable metadata.
- Keep graph-spin models first-class in import, export, visualization, sector
  construction, and NetworkX-style workflows.
- Support arbitrary graph-spin interactions with site labels, coordinates,
  edge-type metadata, Pauli-term labels, symmetry-sector metadata, and clear
  conversion limitations.

### Adapter capability reports

- Provide capability reports for importers and exporters before translation,
  including whether geometry, labels, complex hoppings, spin terms, sectors,
  units, references, provenance, and reconstruction metadata will be preserved.
- Expand adapter convention coverage so partial translations fail clearly or
  report loss explicitly.

## Priority 3: Additional Lattice Transformations

Finite lattices are particularly useful when translation symmetry is broken.
Transformations should operate on portable lattice/model data rather than
creating unrelated builder functions.

- Geometry-aware defect templates beyond current vacancies, subgraphs, and
  bond substitutions
- Extend model-level geometry/interaction transformations to interacting
  particle models.

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
- Reusable layout controls beyond the current declarative multi-panel records

### Reciprocal-space and matrix views

- Reciprocal lattices, Brillouin zones, high-symmetry points, and momentum paths
- Richer automatic symmetry-block discovery for matrix views
- Clear distinction between physical adjacency diagrams and matrix-connectivity
  diagrams

### Analysis plots and output formats

- Standard plots for spectra, bands, density of states, occupations,
  correlations, entanglement, dynamics, and parameter sweeps
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
- Broader styling validation and cross-panel axis/link conventions
- Result-schema migrations when future revisions require them

## Priority 7: Reference Analysis

Reference analysis should stay directly connected to model inspection,
validation, and interchange.

### Spectral and solver safeguards

- Extend iterative solver reporting with backend-specific iteration counts.

This should remain a compact interface, not a replacement for specialized
many-body solver packages.

### Observables and entanglement

- Extend reduced-basis-aware observables as new sectors are introduced.

### Extended topological analysis

- Local Berry curvature and convergence reports
- Additional multi-band and symmetry-aware invariants

### Lightweight dynamics

- Additional state-distance and dynamical-order diagnostics

### Lightweight sweeps and thermal tools

- Simple single-particle spectral functions beyond the current broadened local
  density of states

Checkpointing, distributed execution, and general workflow orchestration should
be delegated to external tools.

## Priority 8: Symmetry-Reduced Bases

Reduced bases should be added where the conserved quantity is explicit and the
mapping to the full basis is testable.

- Consider spin-inversion or parity sectors where they offer a clear benefit.
- Consider translation sectors only after periodic geometry and momentum
  conventions stabilize.

New sectors must preserve labels, basis-state mappings, and resource estimates,
and their spectra and observables must be validated against full-space blocks.

## Model Families

New named builders should provide a clear benchmark or exercise reusable
infrastructure.

Implemented benchmark families include graphene, two-dimensional Anderson,
checkerboard Chern-insulator, and dice/\(T_3\) lattices.

Remaining candidate:

- Kitaev honeycomb spin model

Research-oriented geometries such as pyrochlore lattices, quasicrystal
approximants, or general multi-orbital structures should be added only with
documented conventions, reference checks, and a concrete interchange use case.

## Ecosystem Interoperability and Plugins

### Planned adapters

YAML, ASE, XYZ, SVG, Graphviz DOT, CSV/JSON plot data, OpenFermion, Qiskit,
QuSpin, NetKet, and QuTiP adapters are available. Remaining adapter work:

- pymatgen and selected CIF workflows
- Styled NetworkX graph export
- Expanded PennyLane export
- Capability reports and broader convention coverage for existing adapters

### Plugin discovery

- Registration of importers, exporters, and lightweight analysis routines in
  addition to the current model-builder entry points
- Compatibility checks for plugin and schema versions
- Explicit, diagnosable plugin loading

## Scientific Verification and Project Maturity

Future maturity work:

- Gauge-equivalent constructions where applicable
- Increase annotation completeness now that the package ships `py.typed`.
- Raise coverage expectations as optional-adapter test environments are added.
- Add performance regression thresholds to representative benchmarks.

## Near-Term Backlog

Recommended implementation order:

1. Add parity sectors and extend transformations to interacting particle
   models.
2. Add adapter capability reports and selected pymatgen/CIF workflows.
3. Add solver iteration reporting and simple spectral-function records.
4. Improve annotation completeness, optional-adapter CI, and benchmark
   regression thresholds.
5. Add the Kitaev honeycomb model only with documented gauge and sector
   conventions.

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
