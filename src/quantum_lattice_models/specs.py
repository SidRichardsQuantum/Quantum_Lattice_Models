"""Versioned portable specifications for lattices and registered models.

The stable import surface is intentionally small. Record definitions,
factories, physical inference, and migrations have separate internal module
boundaries so they can evolve without expanding this public module.
"""

from quantum_lattice_models._spec_factories import (
    create_graph_spin_spec,
    create_model_from_preset,
    create_model_spec,
    load_model,
)
from quantum_lattice_models._spec_migrations import migrate_spec_data
from quantum_lattice_models._spec_records import (
    EXTERNAL_MATRIX_FAMILY,
    GRAPH_SPIN_FAMILY,
    SPEC_SCHEMA_VERSION,
    LatticeSpec,
    ModelSpec,
)

for _public_value in (
    LatticeSpec,
    ModelSpec,
    create_model_spec,
    create_graph_spin_spec,
    create_model_from_preset,
    load_model,
    migrate_spec_data,
):
    _public_value.__module__ = __name__
del _public_value

__all__ = [
    "SPEC_SCHEMA_VERSION",
    "EXTERNAL_MATRIX_FAMILY",
    "GRAPH_SPIN_FAMILY",
    "LatticeSpec",
    "ModelSpec",
    "create_model_spec",
    "create_graph_spin_spec",
    "create_model_from_preset",
    "load_model",
    "migrate_spec_data",
]
