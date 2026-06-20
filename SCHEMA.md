# Portable Model Schema

Quantum Lattice Models uses versioned JSON-compatible `ModelSpec` and
`LatticeSpec` records. JSON is the canonical interchange format; CSV and
GraphML adapters translate through the same lattice representation.

The current schema version is `1.0`.

## Compatibility policy

- Files with `schema_version: "1.0"` are validated strictly.
- Legacy unversioned files are treated as the pre-versioned format and migrated
  to `1.0`.
- Unknown fields are rejected to catch misspellings and unsupported data.
- Unsupported past or future versions are rejected until an explicit migration
  is implemented.
- Additive Python API changes do not automatically change the schema version.
  A schema version changes when serialized meaning or required migration
  behavior changes.
- Minor package releases may add optional fields while preserving existing
  `1.0` files.
- Removal or reinterpretation of fields requires a documented schema migration
  and package deprecation period where practical.

## Model record

| Field | Type | Meaning |
|---|---|---|
| `schema_version` | string | Portable schema version |
| `family` | string | Registered model builder family |
| `parameters` | object | Validated builder parameters |
| `lattice` | object or null | Optional `LatticeSpec` |
| `basis` | string or null | Registered Hilbert-space basis |
| `local_degrees` | object list | Indexed physical spins, orbitals, modes, or Nambu components |
| `basis_mappings` | object list | Mapping from local degrees to tensor factors, modes, or single-particle indices |
| `interactions` | object list | Portable onsite and two-body operator terms |
| `representation` | `dense` or `sparse` | Requested matrix representation |
| `units` | string map | Units keyed by parameter or quantity |
| `conventions` | string map | Signs, gauges, orderings, and interpretation rules |
| `references` | string list | DOI, URL, citation, or local reference identifiers |
| `provenance` | object | Source, generator, import, or environment information |
| `metadata` | object | Portable user-defined data |

The reserved `external_matrix` family represents imported Hamiltonians that do
not have a registered reconstruction builder. Its `basis` must be nonempty,
`parameters` must be empty, and the matrix itself must be retained in a
portable NPY/NPZ Hamiltonian file. The model record may carry optional lattice
geometry, units, conventions, references, provenance, and metadata like any
other model.

The reserved `graph_spin` family represents portable arbitrary spin-$1/2$
graphs. Its `parameters` contain only `n_sites`; its physical-system records
carry the Pauli interactions, fields, geometry, labels, and tensor-factor
mapping used for dense or sparse reconstruction.

## Physical-system records

`LocalDegreeOfFreedom` records contain:

- Contiguous `index` and physical `site`
- `kind`: `spin`, `fermion`, `boson`, `orbital`, or `nambu`
- Positive `local_dimension`
- Optional label, component, orbital, and portable metadata

`BasisIndexMapping` assigns every local degree exactly once to a nonnegative
`basis_index`. Its role distinguishes local tensor-factor order,
single-particle matrix-state order, and general mode order. A spin tensor
factor is not a Hamiltonian row index.

`InteractionTerm` currently represents one-site or two-site terms. It stores
the participating local-degree indices, one compatible operator per degree, a
finite complex coefficient, interaction kind, optional label and unit, and
portable metadata. Spin operators are `I`, `X`, `Y`, and `Z`; orbital,
fermionic, and bosonic records use identity, creation, annihilation, and number
operators where applicable. Truncated bosons additionally use `number_pair`
for \(n(n-1)\). Nambu records use explicit particle and hole components for
BdG matrix-state ordering and pairing links.

Multiple local degrees may reference the same physical site. For example,
Fermi-Hubbard has up/down fermionic modes per site and a Kitaev BdG
specification has particle/hole components per site. Their basis mappings and
indices remain distinct.

These fields are optional for compatibility with existing schema `1.0` files.
Newly created specifications populate them for supported model families.

## Artifact bundles

A directory bundle contains independently readable portable artifacts:

- `matrix.npz`, or `matrix.npy` with `matrix.npy.json`
- `model.json`
- `metadata.json`
- `lattice.json` when lattice data exists
- `manifest.json` listing the generated artifact files and matrix format

The matrix artifact remains the canonical round-trip input for
`load_hamiltonian`; the other files support selective interchange and
inspection.

Registry discovery reports whether a model has a sparse construction path and
its package validation status (`validated`, `tested`, or `unvalidated`).
Validation status describes package test coverage; it is not a claim that a
parameter choice is physically appropriate for a particular study.

## Lattice record

| Field | Type | Meaning |
|---|---|---|
| `schema_version` | string | Portable schema version |
| `n_sites` | positive integer | Number of indexed sites or orbitals |
| `positions` | list of 2D or 3D coordinates | One coordinate per site |
| `bonds` | bond-record list | Directed source, target, and optional value |
| `site_labels` | string list | One label per site |
| `orbital_labels` | string list | One orbital label per site |
| `sublattice_labels` | string list | One sublattice label per site |
| `unit_cells` | integer list | Unit-cell membership per site |
| `boundary_conditions` | string map | Axis to `open` or `periodic` |
| `units` | string map | Position, hopping, or other lattice-data units |
| `conventions` | string map | Edge direction, sign, gauge, and indexing rules |
| `references` | string list | Geometry or model references |
| `provenance` | object list | Ordered lattice creation/import/transformation records |
| `metadata` | object | Portable user-defined data |

Bond values are matrix elements when present. A missing value is distinct from
zero and allows a model builder to apply its default hopping convention.

## Complex-number encoding

JSON complex values use:

```json
{"__complex__": [0.25, -0.5]}
```

The first item is the real part and the second is the imaginary part.

## Metadata conventions

- Use `units` for physical units, not free-form explanatory text.
- Use `conventions` for basis ordering, hopping signs, gauge choices,
  directed-edge interpretation, and operator normalization.
- Use `references` for stable citations or URLs.
- Use `provenance` for how data was created, imported, or transformed.
- Use `metadata` only for information that does not fit a standardized field.
- All values must remain JSON-portable. NumPy scalar and complex values are
  encoded by the package.

## CSV representation

CSV interchange uses:

- A site table with coordinates and per-site labels
- A bond table with source, target, value-presence, real, and imaginary columns
- A JSON sidecar containing boundary conditions, units, conventions,
  references, provenance, metadata, and schema version

The sidecar defaults to `<sites-file>.json`.

## GraphML representation

GraphML uses integer site nodes, node attributes for coordinates and labels,
and edge attributes for complex matrix elements. A canonical lattice JSON
record is stored as a graph attribute so standardized metadata survives
round trips. NetworkX and GraphML support requires the optional `networkx`
extra.

## Validation

Validation checks:

- Contiguous site indexing and in-range bond endpoints
- Consistent 2D or 3D coordinates
- Complete label columns when labels are supplied
- Registered model parameters and basis
- Supported matrix representation
- Boundary-condition values
- String-valued units and conventions
- Portable references, provenance, and metadata

Malformed files report the invalid field and expected form.
