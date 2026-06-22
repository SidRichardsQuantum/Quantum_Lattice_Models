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
| Fermi-Hubbard and ladder sectors | Fixed `(N_up, N_down)` Fermi-Hubbard and fixed-magnetization Heisenberg-ladder matrices match full-space blocks, including periodic fermionic signs. Basis projection, embedding, registry dimensions, persistence metadata, and sector observables are checked. | `test_fermi_hubbard_sector_matches_full_space_block`, `test_fermi_sector_registry_spec_and_mapping`, `test_ladder_sector_and_commutator_diagnostics` |
| Conserved quantities and extended dynamics | Explicit commutator norms distinguish conserved and broken quantities. Time-dependent RK4 propagation preserves state norm, while gap-closing and Berry-curvature convergence records preserve machine-readable coordinates and metadata. | `test_ladder_sector_and_commutator_diagnostics`, `test_time_dependent_gap_and_topology_extensions` |
| Spatial and optional interchange adapters | Continuous onsite profiles and interfaces are deterministic. DOT and XYZ preserve their documented subsets, and YAML round-trips through the canonical model schema when PyYAML is installed. | `test_spatial_transformations_and_text_adapters`, `test_yaml_model_adapter` |
| Generic reduced mappings and graph sectors | Spin, boson, and fermion sectors round-trip through one explicit reduced-row mapping. Arbitrary conserving graph-spin sectors match full-space blocks, while leaking Pauli terms are rejected. | `test_generic_reduced_mapping_shared_by_all_sector_bases`, `test_graph_spin_sector_matches_full_block_and_rejects_leakage` |
| Particle-model transformations | Single-particle subgraphs and vacancies reproduce selected full-matrix blocks while remapping geometry and interactions; bond substitutions preserve portable provenance. | `test_particle_model_transformations_preserve_selected_matrix_blocks` |
| Extended visualization records | Spin textures, styled interaction graphs, matrix block metadata, multi-panel layouts, and lattice annotations are validated independently of rendered image bytes. | `test_visual_metadata_spin_texture_blocks_and_styled_graph` |
| Graphene, 2D Anderson, checkerboard, and dice benchmarks | Dense and sparse builders agree. Graphene has bipartite spectral symmetry, Anderson disorder is reproducible, the checkerboard model is Hermitian, and the dice lattice has the expected zero-energy flat-band subspace. | `test_new_benchmark_models_dense_sparse_and_reference_properties` |
| Spin observables and entanglement | Product and Bell states verify site and total magnetization, full and connected correlations, structure factors, reduced density matrices, and von Neumann entropy. Reduced-basis results are cross-checked against full-space embeddings for $X$, $Y$, and $Z$ correlations. | `test_site_and_total_magnetization_for_product_state`, `test_bell_state_correlations_structure_factor_and_entropy`, `test_sector_observables_match_full_space_embedding`, `test_sector_reduced_density_matrix_and_entropy_match_embedding` |
| Discovery and comparison | Registry filters, sparse capability, validation status, named presets, resource-only dry runs, deterministic JSON output, and model parameter/matrix/spectrum/gap comparisons are checked without hidden construction in dry-run paths. | `test_model_filters_and_sparse_capability`, `test_named_presets_create_reproducible_specs`, `test_dry_run_reports_resources_without_building`, `test_model_comparison_parameters_matrices_spectra_and_gaps`, `test_discovery_cli_json_and_terminal_are_deterministic` |
| Lattice interchange and transformations | CSV, NetworkX, and GraphML round trips preserve geometry, labels, complex bonds, and portable metadata. Immutable transformations verify reindexing, provenance, boundary replacement, and deterministic disorder. | `test_csv_lattice_round_trip_preserves_tables_and_metadata`, `test_networkx_and_graphml_round_trip`, `test_relabel_subgraph_and_vacancy_preserve_site_metadata`, `test_disorder_is_reproducible_and_records_provenance` |
| External Hamiltonian import and artifact export | Dense NPY and CSR NPZ imports are checked for metadata preservation, shape, finite numeric values, Hermiticity policy, basis dimension, portable persistence, CLI analysis, re-export, and stored-metadata consistency. Selective model, lattice, matrix, metadata, and bundle exports are checked for deterministic names and matrix round trips. | `test_import_external_dense_matrix_with_portable_metadata`, `test_import_external_sparse_npz_and_validation`, `test_import_matrix_cli_supports_analysis_and_reexport`, `test_load_rejects_matrix_metadata_mismatch`, `test_artifact_exports_are_deterministic_and_round_trippable`, `test_cli_selective_artifact_exports_preserve_default_behavior` |
| Portable physical-system data | Ising, Heisenberg, XXZ, SSH, and custom tight-binding specifications are checked for local-degree, basis-mapping, geometry, operator-axis, hopping, onsite-term, and JSON round-trip consistency. Invalid indices and operator/local-kind combinations are rejected. | `test_supported_models_include_portable_physical_system`, `test_spin_physical_terms_preserve_axes_coefficients_and_geometry`, `test_ssh_and_custom_tight_binding_physical_terms`, `test_physical_system_validation_rejects_invalid_indices_and_operators` |
| Interaction diagrams | Portable spin and hopping interaction records are rendered without reconstructing graph semantics from matrix entries; labels, axes, onsite fields, and coupling annotations are checked. | `test_interaction_graph_consumes_portable_spin_and_hopping_terms` |
| Hubbard, BdG, and graph-spin physical data | Truncated boson factors, spin-resolved fermionic modes, Hubbard interactions, chemical potentials, Nambu particle/hole ordering, pairing links, and arbitrary graph-spin reconstruction are checked through JSON round trips and dense/sparse construction. | `test_hubbard_and_bdg_specs_include_physical_records`, `test_hubbard_physical_terms_match_basis_conventions`, `test_kitaev_bdg_records_particle_hole_order_and_pairing`, `test_portable_graph_spin_spec_round_trip_and_construction` |
| Multi-mode diagrams | Fermi-Hubbard and BdG plots verify that multiple local degrees sharing one physical site receive distinct visual positions. | `test_interaction_graph_offsets_multiple_modes_on_one_site` |
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
