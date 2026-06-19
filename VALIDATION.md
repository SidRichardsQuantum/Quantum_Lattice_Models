# Scientific Validation Matrix

Quantum Lattice Models uses small analytic limits, symmetry checks, and
dense/sparse cross-checks to verify its numerical conventions. These checks are
intended to catch implementation regressions; they do not replace convergence
studies or model-specific validation for research calculations.

## Automated reference checks

| Model or capability | Validation | Test coverage |
|---|---|---|
| Transverse-field Ising chain | At zero field, the open three-site spectrum matches the classical Ising energies and degeneracies. | `test_zero_field_open_ising_analytic_limit` |
| Uniform tight-binding chain | Open-chain eigenvalues match $E_m=-2t\cos[m\pi/(N+1)]$. | `test_uniform_open_tight_binding_chain_analytic_spectrum` |
| SSH chain | The decoupled-dimer limit $t_2=0$ gives repeated energies $\pm t_1$; the topological regime has near-zero edge-localized states. | `test_decoupled_ssh_dimer_analytic_spectrum`, `test_open_ssh_topological_regime_has_edge_localized_low_energy_states` |
| Bose-Hubbard model | Single-site interacting energies match the occupation formula; the noninteracting zero-chemical-potential limit vanishes; hopping conserves particle number. | `test_bose_hubbard_single_site_diagonal`, `test_noninteracting_single_site_hubbard_limits`, `test_bose_hubbard_conserves_total_number` |
| Fermi-Hubbard model | Single-site energies and a fermionic hopping-sign case are checked explicitly; the noninteracting single-site limit vanishes. | `test_fermi_hubbard_single_site_diagonal`, `test_fermi_hubbard_hopping_sign_convention`, `test_noninteracting_single_site_hubbard_limits` |
| Kitaev BdG chain | Eigenvalues occur in particle-hole pairs $E,-E$. | `test_kitaev_bdg_particle_hole_symmetry`, `test_hermiticity_and_particle_hole_diagnostics` |
| Haldane model | Complex next-nearest-neighbor hopping has the documented phase and Hermitian conjugate. | `test_haldane_complex_phase_changes_next_nearest_hopping` |
| Dense/sparse builders | Dense and sparse matrices agree for chain, custom, Hubbard, square, triangular, kagome, Hofstadter, and Haldane models, including periodic boundaries and complex phases. | `test_sparse_builders_match_dense_models`, `test_custom_lattice_infers_site_count_and_sparse_matches_dense` |
| Fixed-magnetization spin sectors | XXZ and magnetization-conserving Heisenberg sector matrices and spectra match the corresponding full computational-basis blocks for open and periodic chains. | `test_xxz_sector_matches_full_space_block`, `test_heisenberg_sector_matches_full_space_block` |
| Spin observables and entanglement | Product and Bell states verify site and total magnetization, full and connected correlations, structure factors, reduced density matrices, and von Neumann entropy. Reduced-basis results are cross-checked against full-space embeddings for $X$, $Y$, and $Z$ correlations. | `test_site_and_total_magnetization_for_product_state`, `test_bell_state_correlations_structure_factor_and_entropy`, `test_sector_observables_match_full_space_embedding`, `test_sector_reduced_density_matrix_and_entropy_match_embedding` |
| Discovery and comparison | Registry filters, sparse capability, validation status, named presets, resource-only dry runs, deterministic JSON output, and model parameter/matrix/spectrum/gap comparisons are checked without hidden construction in dry-run paths. | `test_model_filters_and_sparse_capability`, `test_named_presets_create_reproducible_specs`, `test_dry_run_reports_resources_without_building`, `test_model_comparison_parameters_matrices_spectra_and_gaps`, `test_discovery_cli_json_and_terminal_are_deterministic` |
| Hermiticity | Representative built-in Hamiltonians have real spectra and pass direct Hermiticity checks. | `test_spin_hamiltonians_are_hermitian`, `test_eigenvalues_are_real_for_hermitian_models`, `test_dense_and_sparse_matrix_diagnostics` |

## Computational diagnostics

The public diagnostics API provides:

- Registered-model dimension estimates
- Dense matrix memory estimates
- Dense and sparse storage measurements
- Matrix density summaries
- Hermiticity checks
- Particle-hole spectral-symmetry checks

Example:

```python
from quantum_lattice_models import (
    diagnose_matrix,
    estimate_dense_memory,
    estimate_model_dimension,
)
from quantum_lattice_models.models import tight_binding_chain_sparse

dimension = estimate_model_dimension("transverse_field_ising", n_sites=12)
memory = estimate_dense_memory(dimension)
summary = diagnose_matrix(tight_binding_chain_sparse(n_sites=128))
```

All tolerances used for scientific conclusions should be selected for the
problem scale and documented by downstream analyses.
