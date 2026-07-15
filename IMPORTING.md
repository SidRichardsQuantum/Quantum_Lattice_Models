# Importing and Transforming Lattices

This guide covers the end-to-end portable workflow:

```text
import -> validate -> transform -> save -> build -> analyze -> export
```

## External Hamiltonian matrices

Dense NPY and dense or CSR NPZ Hamiltonians can be imported with explicit JSON
metadata:

```json
{
  "basis": "single-particle site basis",
  "basis_dimension": 2,
  "lattice": {
    "n_sites": 2,
    "positions": [[0.0, 0.0], [1.0, 0.0]],
    "bonds": []
  },
  "local_degrees": [
    {
      "index": 0,
      "site": 0,
      "kind": "orbital",
      "local_dimension": 2,
      "label": "left"
    },
    {
      "index": 1,
      "site": 1,
      "kind": "orbital",
      "local_dimension": 2,
      "label": "right"
    }
  ],
  "basis_mappings": [
    {"local_degree": 0, "basis_index": 0, "role": "single_particle_state"},
    {"local_degree": 1, "basis_index": 1, "role": "single_particle_state"}
  ],
  "interactions": [
    {
      "degrees": [0, 1],
      "operators": ["create", "annihilate"],
      "coefficient": {"__complex__": [1.0, 0.0]},
      "kind": "hopping"
    }
  ],
  "units": {"energy": "eV"},
  "conventions": {"basis_ordering": "site index"},
  "references": [],
  "provenance": {"generator": "external-code"},
  "metadata": {}
}
```

```python
from quantum_lattice_models import import_hamiltonian, save_hamiltonian

result = import_hamiltonian("external.npy", metadata_path="external-metadata.json")
save_hamiltonian(result, "portable.npz")
```

```bash
quantum-lattice import-matrix external.npy \
  --metadata external-metadata.json \
  --output portable.npz
```

Imports reject non-square, nonnumeric, nonfinite, or non-Hermitian matrices.
Use `require_hermitian=False` in Python or `--allow-non-hermitian` in the CLI
only when the non-Hermitian interpretation is intentional. `basis_dimension`
is optional; when supplied, it must match the matrix dimension. Lattice site
count is not assumed to equal Hilbert-space dimension.
Optional `local_degrees`, `basis_mappings`, and `interactions` records preserve
the imported matrix's physical interpretation using the same validation as
created model specifications.

## CSV import

```python
from quantum_lattice_models import import_lattice_csv, create_model_spec

lattice = import_lattice_csv("sites.csv", "bonds.csv")
model = create_model_spec(
    "custom_tight_binding",
    lattice=lattice,
    parameters={"hopping": 1.0, "onsite": 0.0, "hermitian": True},
    units={"hopping": "eV"},
    conventions={"hopping_sign": "two-column bonds use -hopping"},
    provenance={"source": "experiment-42"},
)
model.save("model.json")
```

Equivalent CLI workflow:

```bash
quantum-lattice import sites.csv \
  --format csv \
  --bonds bonds.csv \
  --output model.json
```

Use `--metadata sites.csv.json` when the sidecar has a non-default name.

## CSV export

```python
from quantum_lattice_models import export_lattice_csv, load_model

model = load_model("model.json")
export_lattice_csv(model.lattice, "sites.csv", "bonds.csv")
```

```bash
quantum-lattice export-lattice model.json \
  --format csv \
  --sites sites.csv \
  --bonds bonds.csv
```

## Selective artifact export

The general export command accepts registered model JSON and portable
Hamiltonian NPY/NPZ inputs:

```bash
quantum-lattice export model.json --artifact matrix --format npz
quantum-lattice export model.json --artifact model --output model-only.json
quantum-lattice export model.json --artifact lattice --output lattice-only.json
quantum-lattice export model.json --artifact metadata --output metadata-only.json
quantum-lattice export model.json --artifact bundle --output model.bundle
```

Model and lattice exports do not construct a Hamiltonian. Metadata exports
describe a constructed or loaded `HamiltonianResult`. Bundles use deterministic
filenames and contain a self-contained matrix file, model JSON, result metadata,
optional lattice JSON, and a manifest. With `--format npy`, the required
`matrix.npy.json` sidecar is included.

## NetworkX and GraphML

Install the optional adapter:

```bash
pip install "quantum-lattice-models[networkx]"
```

```python
from quantum_lattice_models import (
    export_graphml,
    from_networkx,
    import_graphml,
    to_networkx,
)

graph = to_networkx(lattice)
restored = from_networkx(graph)

export_graphml(lattice, "lattice.graphml")
restored = import_graphml("lattice.graphml")
```

`MultiDiGraph` is used to preserve directed bonds, duplicate edges, and
explicit complex matrix elements.

## Transformations

Transformations return new immutable `LatticeSpec` values and append provenance
records:

```python
from quantum_lattice_models import (
    apply_onsite_disorder,
    lattice_subgraph,
    remove_lattice_sites,
    with_boundary_conditions,
)

fragment = lattice_subgraph(lattice, [0, 1, 4, 5])
vacancy = remove_lattice_sites(lattice, [3])
disordered = apply_onsite_disorder(lattice, strength=0.5, seed=123)
periodic = with_boundary_conditions(lattice, {"x": "periodic"})
```

Available transformations include:

- Site relabeling
- Induced subgraphs and vacancy removal
- Bond addition and removal
- Boundary-condition replacement
- Reproducible onsite and bond disorder

Onsite disorder values are stored under `lattice.metadata["onsite_potential"]`.
Bond disorder produces explicit bond values. Seeds and distribution parameters
are recorded in lattice provenance.

## Build and inspect

```python
result = model.build_result(sparse=True)
print(model.summary())
print(result.matrix.shape)
```

```bash
quantum-lattice validate model.json
quantum-lattice inspect model.json
quantum-lattice spectrum model.json
quantum-lattice export model.json --format npz
```

Use the intake commands before translating an imported model:

```bash
quantum-lattice describe model.json --json
quantum-lattice lint model.json --json
quantum-lattice adapter-capabilities model.json graphml --json
```

`describe` reports physical content and reconstruction limits. `lint` reports
missing or inconsistent metadata and suggests validation checks.
`adapter-capabilities` identifies which semantics a target can preserve before
an export is attempted.

See [SCHEMA.md](SCHEMA.md) for the compatibility policy and complete field
definitions.
