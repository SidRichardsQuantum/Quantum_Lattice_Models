"""Built-in model-to-physical-record inference boundary."""

from quantum_lattice_models._specification import (
    _infer_physical_system as infer_physical_system,
)
from quantum_lattice_models._specification import (
    _infer_symmetry_actions as infer_symmetry_actions,
)

__all__ = ["infer_physical_system", "infer_symmetry_actions"]
