# Importing and Transforming Lattices

This guide covers the end-to-end portable workflow:

```text
import -> validate -> transform -> save -> build -> analyze -> export
```

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

See [SCHEMA.md](SCHEMA.md) for the compatibility policy and complete field
definitions.
