"""Small quantum lattice Hamiltonians for educational research prototypes.

The package uses dense matrices for spin-chain Hamiltonians. This is deliberate
and convenient for exact diagonalization on small systems, but it is not meant
for large-scale simulation.
"""

from quantum_lattice_models.diagnostics import (
    MatrixDiagnostics,
    diagnose_matrix,
    estimate_dense_memory,
    estimate_model_dimension,
    has_particle_hole_symmetric_spectrum,
    is_hermitian,
    matrix_density,
    matrix_storage_bytes,
)
from quantum_lattice_models.geometry import (
    honeycomb_lattice_positions,
    kagome_lattice_positions,
    square_lattice_positions,
    triangular_lattice_positions,
)
from quantum_lattice_models.hubbard import (
    bose_hubbard_chain,
    bose_hubbard_chain_sparse,
    fermi_hubbard_chain,
    fermi_hubbard_chain_sparse,
)
from quantum_lattice_models.lattice import (
    Bond,
    Lattice,
    TightBindingModel,
    custom_tight_binding,
    custom_tight_binding_sparse,
)
from quantum_lattice_models.registry import (
    get_model_info,
    list_models,
    model_table,
    register_model,
    unregister_model,
)
from quantum_lattice_models.specs import (
    SPEC_SCHEMA_VERSION,
    LatticeSpec,
    ModelSpec,
    create_model_spec,
    load_model,
)
from quantum_lattice_models.spin import (
    heisenberg_chain,
    heisenberg_ladder,
    j1_j2_heisenberg_chain,
    longitudinal_field_ising,
    next_nearest_neighbor_ising,
    transverse_field_ising,
    xxz_chain,
    xy_chain,
)
from quantum_lattice_models.tight_binding import (
    aubry_andre_harper_chain,
    kagome_lattice_index,
    kagome_lattice_tight_binding,
    kagome_lattice_tight_binding_sparse,
    rice_mele_model,
    square_lattice_index,
    square_lattice_tight_binding,
    square_lattice_tight_binding_sparse,
    ssh_model,
    tight_binding_chain,
    tight_binding_chain_sparse,
    triangular_lattice_tight_binding,
    triangular_lattice_tight_binding_sparse,
)
from quantum_lattice_models.topological import (
    haldane_honeycomb_lattice,
    haldane_honeycomb_lattice_sparse,
    harper_hofstadter_square_lattice,
    harper_hofstadter_square_lattice_sparse,
    honeycomb_lattice_index,
    kitaev_chain_bdg,
)

__all__ = [
    "aubry_andre_harper_chain",
    "Bond",
    "bose_hubbard_chain",
    "bose_hubbard_chain_sparse",
    "custom_tight_binding",
    "custom_tight_binding_sparse",
    "fermi_hubbard_chain",
    "fermi_hubbard_chain_sparse",
    "haldane_honeycomb_lattice",
    "haldane_honeycomb_lattice_sparse",
    "harper_hofstadter_square_lattice",
    "harper_hofstadter_square_lattice_sparse",
    "heisenberg_chain",
    "heisenberg_ladder",
    "honeycomb_lattice_index",
    "j1_j2_heisenberg_chain",
    "kagome_lattice_index",
    "kagome_lattice_tight_binding",
    "kagome_lattice_tight_binding_sparse",
    "kitaev_chain_bdg",
    "Lattice",
    "longitudinal_field_ising",
    "next_nearest_neighbor_ising",
    "rice_mele_model",
    "square_lattice_index",
    "ssh_model",
    "square_lattice_tight_binding",
    "square_lattice_tight_binding_sparse",
    "tight_binding_chain",
    "tight_binding_chain_sparse",
    "TightBindingModel",
    "triangular_lattice_tight_binding",
    "triangular_lattice_tight_binding_sparse",
    "transverse_field_ising",
    "xxz_chain",
    "xy_chain",
    "honeycomb_lattice_positions",
    "kagome_lattice_positions",
    "square_lattice_positions",
    "triangular_lattice_positions",
    "get_model_info",
    "list_models",
    "model_table",
    "register_model",
    "unregister_model",
    "MatrixDiagnostics",
    "diagnose_matrix",
    "estimate_dense_memory",
    "estimate_model_dimension",
    "has_particle_hole_symmetric_spectrum",
    "is_hermitian",
    "matrix_density",
    "matrix_storage_bytes",
    "SPEC_SCHEMA_VERSION",
    "LatticeSpec",
    "ModelSpec",
    "create_model_spec",
    "load_model",
]
