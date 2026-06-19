# Roadmap

Quantum Lattice Models aims to support a complete, reproducible workflow:

```text
create or import -> validate -> inspect -> build -> analyze -> export or save
```

The package already provides a useful collection of finite lattice-model
builders, exact-diagonalization helpers, plotting tools, a model registry, and
a command-line interface. The next priority is to make models portable,
reproducible, and interoperable rather than adding many more isolated builder
functions.

This roadmap describes intended direction, not a compatibility guarantee.
Priorities and implementation order may change as the data model and user
workflows are tested.

Implemented in v0.1.4: typed registry parameter metadata, registry-generated CLI
arguments, dimension and memory estimates, matrix storage and Hermiticity
diagnostics, particle-hole spectral checks, and an initial scientific validation
matrix. Versioned `ModelSpec` and `LatticeSpec` objects, canonical JSON model
files, dense/sparse reconstruction, and CLI `create`, `inspect`, and `validate`
commands are also available. Rich result containers, additional file formats,
and metadata-preserving sparse result containers remain future work.

## Guiding Principles

- Keep the package useful from both Python and the command line.
- Treat geometry, model parameters, basis conventions, and provenance as
  first-class data.
- Use one validated representation across creation, import, analysis, and
  export workflows.
- Keep dense and sparse Hamiltonians behaviorally consistent.
- Prefer explicit physical conventions over implicit defaults.
- Preserve lightweight core dependencies and place integrations behind
  optional extras.
- Keep notebooks and examples as thin clients of the public package API.

## Portable Model Specifications

A primary near-term goal is to make lattice models serializable and
reproducible.

### Model and lattice specifications

- Extend the versioned `LatticeSpec` representation with:
  - Site count and coordinates
  - Bonds and complex hopping amplitudes
  - Site, orbital, and sublattice labels
  - Unit-cell membership
  - Boundary conditions
  - User metadata
- Extend the versioned `ModelSpec` representation with:
  - Model family and parameters
  - Lattice specification
  - Basis and Hilbert-space conventions
  - Dense or sparse construction preferences
  - Units, references, and provenance
- Add a structured Hamiltonian result/container that associates matrix data
  with its model, basis, geometry, and construction metadata.
- Preserve equivalent metadata for dense and sparse matrices.

### Validation

- Validate site indices, coordinate dimensions, bond records, parameter
  types, boundary conditions, and model-specific requirements.
- Detect malformed, ambiguous, and non-Hermitian specifications when
  Hermiticity is required.
- Produce error messages that identify the invalid field and expected form.
- Provide schema checks and a migration path for future format changes.

### Serialization and import/export

- Continue evolving JSON as the canonical portable format with documented
  migrations.
- Support YAML through an optional dependency.
- Import and export lattice coordinates and edge lists as CSV.
- Save and load dense Hamiltonians with NumPy formats.
- Save and load sparse Hamiltonians with SciPy NPZ.
- Encode complex values explicitly and portably rather than relying on
  implementation-specific JSON behavior.
- Guarantee round-trip preservation of supported model information.

Example target API:

```python
from quantum_lattice_models import load_model

model = load_model("ssh.json")
model.validate()
H = model.hamiltonian(sparse=True)
model.save("experiment.json")
```

### Registry and parameter schemas

- Extend registry entries with typed parameter definitions.
- Record required parameters, defaults, descriptions, constraints, and CLI
  names in one place.
- Generate model-specific CLI options and documentation from registry
  metadata.
- Retain runtime registration while preparing for third-party plugins.

### File-oriented CLI

Extend commands centered on reusable model files:

```bash
quantum-lattice create ssh_model --n-cells 20 --output ssh.json
quantum-lattice import lattice.csv --output lattice.json
quantum-lattice inspect ssh.json
quantum-lattice validate ssh.json
quantum-lattice spectrum ssh.json
quantum-lattice export ssh.json --format npz
```

Existing direct model invocation should remain available during the
transition.

### Documentation and testing

- Add a create, save, load, analyze, and export tutorial.
- Document the model-file schema and complex-number encoding.
- Add round-trip tests for every supported format.
- Add compatibility tests for dense and sparse construction.
- Test malformed files and schema-version failures.
- Add public API reference documentation.

### Completion criteria

This capability is complete when a user can:

1. Create a built-in or custom lattice model.
2. Save its geometry, parameters, conventions, and metadata.
3. Load it in a new process without manually reconstructing Python objects.
4. Validate and inspect it.
5. Build an equivalent dense or sparse Hamiltonian.
6. Analyze it through the Python API or CLI.
7. Export the resulting matrix and relevant metadata.

## Spin-System Infrastructure

Spin support should become more scalable and general before the package adds a
large number of isolated named spin builders.

### Sparse spin construction

- Add direct sparse construction for:
  - Transverse- and longitudinal-field Ising chains
  - Next-nearest-neighbor and long-range Ising chains
  - XY, XYZ, XXZ, and Heisenberg chains
  - J1-J2 chains and Heisenberg ladders
  - Random-field variants
- Avoid constructing dense Kronecker products inside sparse APIs.
- Keep dense and sparse builders generated from the same interaction
  specification where practical.
- Cross-check sparse matrices and low-energy spectra against dense reference
  constructions.

### General spin-model representation

- Add a graph-based spin Hamiltonian builder that supports:
  - Arbitrary interaction bonds
  - XX, YY, ZZ, and mixed-axis couplings
  - Site-dependent x, y, and z fields
  - Complex couplings where physically meaningful
  - Optional Dzyaloshinskii-Moriya interactions
- Refactor named spin-chain and ladder builders to delegate to this common
  representation without changing their public convenience APIs.
- Preserve Pauli-term metadata for export and verification.
- Make interaction conventions and factors of $1/2$ explicit.

### Spin observables and entanglement

- Add site-resolved magnetization and total $S^z$.
- Add full one- and two-point correlation matrices.
- Add connected correlations and static spin structure factors.
- Add total-spin $S^2$ diagnostics where basis conventions permit.
- Add reduced density matrices and bipartite entanglement entropy for small
  systems.
- Make observables compatible with dense, sparse, and symmetry-reduced bases.

### Spin dynamics

- Add time-independent state evolution using dense exponentiation and sparse
  `expm_multiply` paths.
- Add parameter-quench workflows and expectation-value time series.
- Add Loschmidt amplitudes and echoes.
- Track magnetization, correlations, and entanglement during evolution.
- Return structured results containing times, tolerances, solver details, and
  model provenance.

### Additional spin models

Add these after the common sparse and graph-based infrastructure is stable:

- XYZ chain
- Random-field Heisenberg and XXZ chains
- Long-range power-law Ising models
- Chains and graphs with Dzyaloshinskii-Moriya interactions
- Spin-1 bilinear-biquadratic and AKLT chains
- Kitaev honeycomb spin model

Spin-1 and Kitaev honeycomb support require explicit local-basis, operator, and
physical-convention designs and are therefore lower priority than spin-$1/2$
chains and ladders.

## Periodic Lattices and Analysis

Represent translationally invariant systems and provide more complete physical
analysis.

### Periodic lattice geometry

- Add primitive lattice vectors.
- Add basis sites and orbitals within a unit cell.
- Represent bonds using cell-displacement vectors.
- Distinguish finite open systems from periodic real-space systems.
- Compute reciprocal lattice vectors.
- Provide standard and user-defined Brillouin-zone paths.
- Construct Bloch Hamiltonians \(H(\mathbf{k})\).
- Calculate and plot band structures.
- Support k-resolved eigenvectors and observables.

The periodic representation should work naturally with SSH, Rice-Mele,
Haldane, Hofstadter, square, triangular, honeycomb, and kagome systems.

### Additional finite lattice and model families

New lattice builders should exercise existing infrastructure or provide a
useful analytic or topological benchmark rather than only expanding the model
count.

Near-term candidates:

- Plain nearest-neighbor honeycomb/graphene lattice, independent of the Haldane
  model
- Lieb lattice with a three-site unit cell and flat band
- Sawtooth or delta chain as a compact frustrated flat-band model
- Creutz ladder with complex hopping, flux, topology, and edge states
- Anderson-disordered tight-binding chains and lattices with reproducible
  onsite disorder
- Configurable long-range power-law tight-binding chains

Candidates to add after Bloch-Hamiltonian and topological-analysis support:

- Checkerboard or related two-band Chern-insulator model
- Dice or $T_3$ lattice

Longer-term or optional research-oriented candidates:

- Pyrochlore lattice
- Penrose or other quasicrystal approximants
- General multi-orbital lattice models

The more specialized candidates should be introduced only with documented
geometry, conventions, reference limits, and a clear analysis workflow.

### General observables

- Site and orbital occupations
- Generic one- and two-point expectation values
- Connected correlations and correlation matrices
- Static structure factors
- Participation ratios and localization diagnostics
- Reduced density matrices and entanglement entropy for small systems
- Current operators for supported tight-binding models

### Topological analysis

- Berry and Zak phases
- One-dimensional winding numbers
- Wilson loops where applicable
- Chern numbers for supported two-dimensional Bloch models
- Gauge-tolerant numerical routines with documented convergence controls
- Reference tests against known SSH and Haldane phases

### Dynamics

- State evolution under time-independent Hamiltonians
- Expectation-value time series
- Dense and sparse evolution paths
- Parameter quench helpers
- Structured results containing times, parameters, tolerances, and
  provenance

### Symmetry and diagnostics

- Hermiticity and particle-hole checks
- Commutator-based conserved-quantity checks
- Degeneracy reporting with configurable tolerances
- Basis-dimension and memory estimates before matrix construction
- Clear warnings for physically or computationally risky configurations

## Scientific Verification and Conventions

Every built-in model should have an explicit, testable account of its physical
and numerical conventions.

- Record references and implemented conventions for every built-in model.
- Document basis ordering, operator normalization, hopping signs, gauge
  choices, boundary conditions, and parameter units.
- Add analytic-limit tests for small systems and special parameter values.
- Add regression fixtures for known spectra, eigenstates, observables, and
  topological invariants.
- Test gauge-equivalent constructions where applicable.
- Cross-check dense and sparse builders against the same reference results.
- Publish a model-validation matrix summarizing the checks applied to each
  builder.
- Include references and validation status in model registry metadata.

## Reproducible Experiments and Results

Model construction and analysis should produce reusable records rather than
disconnected arrays and figures.

- Add an `ExperimentSpec` or equivalent representation that combines a model
  specification with solver and analysis settings.
- Add structured result objects for spectra, eigenstates, observables,
  parameter sweeps, and dynamics.
- Record the model specification, package version, solver, tolerances, random
  seeds, and relevant environment metadata.
- Save and load results without losing labels, coordinates, units, or
  provenance.
- Export tabular results to CSV and richer numerical results to an appropriate
  portable scientific-data format.
- Support deterministic output naming and resumable long-running analyses.
- Make result summaries directly usable in notebooks, reports, and CLI output.

## Parameter Sweeps and Finite-Size Studies

Parameter studies should be first-class analysis workflows instead of custom
loops repeated across notebooks.

- Provide general one- and multi-parameter sweep APIs.
- Return labeled parameter coordinates with spectra and observables.
- Support finite-size sweeps and scaling summaries.
- Detect candidate gap closings, level crossings, and extrema.
- Track eigenstate continuity between nearby parameter values where feasible.
- Parallelize independent parameter points through an optional execution
  backend.
- Support checkpoints and resuming partially completed sweeps.
- Add standard plots for sweep curves, heat maps, phase diagrams, and
  finite-size scaling.

## Thermal, Spectral, and Response Analysis

Users should be able to move beyond zero-temperature eigenvalue inspection
without assembling common statistical-mechanics calculations by hand.

- Add canonical partition functions and thermal expectation values for small
  dense systems.
- Add free energy, entropy, heat capacity, and susceptibility helpers.
- Add local and total density of states with configurable broadening.
- Add local spectral functions and simple Green-function utilities for
  single-particle models.
- Add frequency-domain response functions for selected observables.
- Support exact full-spectrum paths and clearly labeled low-energy
  approximations.
- Record temperature grids, broadening, truncation, and numerical tolerances in
  structured results.
- Add analytic tests for noninteracting and few-level limits.

## Disorder, Defects, and Interfaces

Finite lattices are especially useful for studying broken translation symmetry.

- Add reproducible disorder specifications with explicit random seeds and
  distributions.
- Support onsite, bond, and correlated disorder.
- Add vacancies, impurities, modified bonds, and boundary potentials through
  reusable lattice transformations.
- Add domain walls, heterostructures, and interfaces between parameter regions.
- Add open-boundary edge and corner-state diagnostics.
- Provide disorder-ensemble workflows reporting means, variances, confidence
  intervals, and realization metadata.
- Integrate disorder ensembles with parameter sweeps, localization measures,
  and resumable result storage.

## User Workflows and Ergonomics

The package should make common workflows discoverable without requiring users
to inspect implementation signatures or write repetitive glue code.

- Add model search and filtering by category, basis, dimensionality, sparse
  support, and validation status.
- Add named, documented parameter presets for canonical phases and reference
  limits without hiding the underlying values.
- Add a dry-run or inspection API that reports basis, dimension, estimated
  memory, expected matrix representation, and warnings before construction.
- Add model-to-model comparison helpers for spectra, gaps, observables, and
  matrix differences.
- Add deterministic tabular summaries suitable for terminals, notebooks, CSV,
  and data-frame conversion without requiring a data-frame dependency.
- Add batch CLI execution from a manifest and machine-readable JSON output for
  scripting.
- Add progress reporting and cancellation points for sweeps and ensemble runs.
- Add task-oriented recipes for localization, edge states, finite-size scaling,
  quench dynamics, and dense-versus-sparse solver selection.
- Improve error messages with model-specific parameter context and suggested
  valid alternatives.

## Symmetry-Reduced Hilbert Spaces

Sparse storage alone does not remove exponential Hilbert-space growth.
Conserved-sector construction should reduce the represented basis whenever the
model permits it.

- Add fixed-magnetization sectors for compatible spin models.
- Add spin-inversion or parity sectors where they produce a clear practical
  benefit.
- Consider translation sectors only after periodic geometry and momentum
  conventions stabilize.
- Add fixed-particle-number sectors for Bose-Hubbard models.
- Add fixed particle-number and spin sectors for Fermi-Hubbard models.
- Provide explicit basis-state enumeration and state-to-index mappings.
- Preserve sector labels and basis information in Hamiltonian metadata.
- Make operators and observables aware of reduced bases.
- Validate reduced-sector spectra against the corresponding blocks of
  full-space Hamiltonians.
- Report sector dimensions and estimated memory use before construction.

## Solver Interface and Computational Safeguards

Analysis routines should expose numerical behavior and avoid surprising dense
conversions.

- Define a common solver interface for dense, sparse, and optional external
  eigensolvers.
- Select sensible solvers from matrix properties and requested results while
  allowing explicit user control.
- Return convergence status, residuals, tolerances, and iteration information.
- Prevent accidental conversion of large sparse matrices to dense arrays.
- Estimate memory and computational cost before expensive construction or
  diagonalization.
- Support reproducible initial vectors and random seeds where solvers use
  randomness.
- Distinguish exact results from iterative approximations in result metadata
  and user-facing output.

## Ecosystem Interoperability

Exchange models with established scientific and quantum-software ecosystems.

Integrations should be optional and should translate through the stable model
specification rather than introducing separate internal representations.

### Planned adapters

- NetworkX graph import and export
- GraphML import and export
- ASE and/or pymatgen structure import, including selected CIF workflows
- OpenFermion operator export
- Expanded PennyLane export
- Qiskit operator export
- QuSpin or NetKet adapters where basis and convention mappings are explicit

### Plugin discovery

- Support Python entry points for third-party model packages.
- Allow external packages to register builders, schemas, importers, exporters,
  and analysis routines.
- Provide compatibility checks for plugin and schema versions.
- Keep plugin loading explicit and diagnosable.

## API Stability and Project Maturity

The core model schema and public API should be exercised extensively before
they are treated as stable.

Expected requirements:

- Stable public model and lattice specifications
- Documented compatibility and deprecation policy
- Schema migration support
- Reliable dense and sparse metadata behavior
- Type-complete public API with `py.typed`
- API reference and task-oriented documentation
- Tested import/export round trips
- Supported Python-version policy
- Performance and memory guidance
- Reproducible release and provenance metadata

## Engineering Work

The following work should proceed alongside feature development:

- Reduce reliance on `numpy.ndarray` subclasses as the only metadata carrier.
- Ensure numerical operations do not silently discard essential metadata.
- Add static type checking and distribute a `py.typed` marker.
- Add coverage reporting and define coverage expectations for new modules.
- Test package installation and CLI behavior from built wheels.
- Add tests for supported minimum and maximum Python versions.
- Keep generated documentation and notebook artifacts reproducible.
- Add benchmark coverage for representative dense, sparse, and import/export
  workflows.

## Near-Term Backlog

Recommended implementation order:

1. Expand schema validation and define migration behavior for future versions.
2. Document and test model conventions and analytic reference cases.
3. Introduce metadata-preserving Hamiltonian and result containers.
4. Update existing builders to emit complete specifications and results.
5. Add NPY, NPZ, CSV, and optional YAML support.
6. Add file-oriented spectrum and export commands that consume model files.
7. Add a common sparse graph-based spin builder and migrate existing spin
   models onto it.
8. Add fixed-magnetization sectors for XXZ and Heisenberg models.
9. Add spin correlation matrices, structure factors, reduced density
    matrices, and entanglement entropy.
10. Add dense and sparse spin dynamics and quench workflows.
11. Add reusable parameter-sweep and finite-size-study APIs.
12. Add model inspection, dry-run resource estimates, presets, comparison
    helpers, and machine-readable CLI output.
13. Add reproducible disorder, defect, interface, and ensemble workflows.
14. Add thermal observables, broadened density-of-states, and basic response
    analysis.
15. Add plain honeycomb, Lieb, sawtooth, Creutz-ladder, and Anderson model
    builders with analytic reference tests.
16. Publish the end-to-end tutorial and schema reference.
17. Begin periodic-lattice and Bloch-Hamiltonian design only after the portable
    formats stabilize.
18. Add checkerboard/Chern-insulator and other advanced lattice models only
    after topological analysis is available.

## Explicit Non-Goals

The project is not currently intended to:

- Compete with large-scale tensor-network or many-body simulation frameworks.
- Hide exponential Hilbert-space scaling from users.
- Infer a physically complete model from atomic coordinates alone.
- Promise lossless conversion to every external package or basis convention.
- Treat single-particle, BdG, spin, bosonic, and fermionic matrices as
  interchangeable.
- Claim quantum advantage from small exact-diagonalization examples.

New features should preserve these distinctions and state their basis,
conventions, dimensions, and computational limits explicitly.
