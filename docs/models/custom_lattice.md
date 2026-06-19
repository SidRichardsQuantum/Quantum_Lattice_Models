<!-- builders: custom_tight_binding, custom_tight_binding_sparse -->
# Custom Lattice and Tight Binding

## Purpose and structure

`Lattice`, `LatticeSpec`, and `TightBindingModel` support finite user-defined
graphs with coordinates, complex bonds, labels, unit cells, boundary metadata,
and onsite terms.

Two-item bonds `(i, j)` use `-hopping`. Three-item bonds `(i, j, value)` use
`value` directly. With `hermitian=True`, the conjugate reverse matrix element
is added automatically.

## Package use

```python
from quantum_lattice_models import Lattice, TightBindingModel

lattice = Lattice(
    positions=[(0, 0), (1, 0), (0.5, 0.8)],
    bonds=[(0, 1), (1, 2, 0.25j), (2, 0)],
)
H = TightBindingModel(lattice).hamiltonian(hopping=1.0)
```

```bash
quantum-lattice create custom_tight_binding \
  --n-sites 3 --bond 0,1 --bond 1,2,0.25j --output custom.json
```

## Parameters

{{PARAMETERS}}

## Validation and cautions

Bond indices, coordinate dimensions, and dense/sparse agreement are tested.
Set `hermitian=False` only when a directed or explicitly non-Hermitian model is
intended.
