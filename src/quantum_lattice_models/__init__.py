"""Small quantum lattice Hamiltonians for educational research prototypes.

The package uses dense matrices for spin-chain Hamiltonians. This is deliberate
and convenient for exact diagonalization on small systems, but it is not meant
for large-scale simulation.
"""

from quantum_lattice_models.models import (
    heisenberg_chain,
    ssh_model,
    tight_binding_chain,
    transverse_field_ising,
)

__all__ = [
    "heisenberg_chain",
    "ssh_model",
    "tight_binding_chain",
    "transverse_field_ising",
]
