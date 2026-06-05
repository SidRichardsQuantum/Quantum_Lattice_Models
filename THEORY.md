# Theory

This document gives short, working notes for the models and methods implemented in the package. It is not a full textbook treatment.

## Transverse-Field Ising Model

The transverse-field Ising chain is a spin-1/2 model with nearest-neighbor `ZZ` interactions and a transverse `X` field:

```text
H = -J sum_i Z_i Z_{i+1} - h sum_i X_i
```

The interaction term favors aligned or anti-aligned computational-basis states depending on the sign of `J`. The transverse field mixes computational-basis states and drives quantum fluctuations. This makes the model a standard testbed for exact diagonalization, phase-transition intuition, VQE ansatz studies, and quantum simulation circuits.

The longitudinal-field variant adds `-h_z sum_i Z_i`, which breaks the simple
spin-flip symmetry of the transverse-field model. The next-nearest-neighbor
variant adds `-J2 sum_i Z_i Z_{i+2}`, giving a small frustrated Ising testbed.

## Heisenberg Chain

The anisotropic Heisenberg chain includes `XX`, `YY`, and `ZZ` couplings:

```text
H = sum_i (Jx X_i X_{i+1} + Jy Y_i Y_{i+1} + Jz Z_i Z_{i+1})
    + field sum_i Z_i
```

Special choices recover familiar limits. `Jx = Jy = Jz` gives the isotropic Heisenberg model. Unequal couplings give anisotropic variants. The model is useful for testing symmetry-aware methods, spin transport ideas, exact diagonalization, and variational workflows.

## XY and XXZ Spin Chains

The XY chain keeps `XX` and `YY` couplings with an anisotropy parameter and a transverse `Z` field:

```text
H = -J sum_i [((1 + gamma) / 2) X_i X_{i+1}
             + ((1 - gamma) / 2) Y_i Y_{i+1}]
    - field sum_i Z_i
```

The XXZ chain is the `Jx = Jy` specialization of the anisotropic Heisenberg chain:

```text
H = J sum_i (X_i X_{i+1} + Y_i Y_{i+1} + Delta Z_i Z_{i+1})
    + field sum_i Z_i
```

These models are useful for studying anisotropy, magnetization, spectral gaps, and small spin-system benchmarks with exact references.

## J1-J2 Heisenberg Chain

The J1-J2 chain adds next-nearest-neighbor Heisenberg interactions:

```text
H = J1 sum_i S_i . S_{i+1} + J2 sum_i S_i . S_{i+2}
    + field sum_i Z_i
```

In this package, the Pauli-matrix convention is used directly, so `S_i . S_j`
is represented as `X_i X_j + Y_i Y_j + Z_i Z_j`. The model is a compact
frustrated spin-chain testbed.

## Heisenberg Ladders

A two-leg Heisenberg ladder connects two spin chains with rung couplings. It is
a compact way to move beyond strictly one-dimensional chains while keeping the
dense Hilbert space small enough for exact diagonalization at low rung counts.

## Hubbard Models

The Bose-Hubbard model describes bosons hopping on a lattice with onsite
repulsion and chemical potential:

```text
H = -t sum_<ij> (a_i^\dagger a_j + a_j^\dagger a_i)
    + U/2 sum_i n_i (n_i - 1) - mu sum_i n_i
```

The implementation uses a truncated local occupation basis
`|0>, ..., |max_occupancy>`.

The spinful Fermi-Hubbard chain is implemented in occupation-number basis with
explicit fermionic signs:

```text
H = -t sum_<ij>,s (c_i,s^\dagger c_j,s + c_j,s^\dagger c_i,s)
    + U sum_i n_i,up n_i,down - mu sum_i,s n_i,s
```

Orbital order is `site0_up, site0_down, site1_up, site1_down, ...`.

## Kitaev Chain

The Kitaev chain is represented here as a Bogoliubov-de Gennes matrix for a
spinless p-wave superconducting chain. This is a single-particle doubled-Nambu
representation, not a many-body occupation-basis Hamiltonian.

## SSH Model

The Su-Schrieffer-Heeger model is a one-dimensional tight-binding model with two sites per unit cell and alternating hopping amplitudes `t1` and `t2`.

For open boundaries and `t1 < t2`, the finite chain supports near-zero-energy states localized near the two ends. In this package, the SSH Hamiltonian is a single-particle matrix of size `2 * n_cells`, not a many-body qubit Hamiltonian.

The SSH model is a compact way to study band structure, boundary modes, topology in simple lattice systems, and quantum walk dynamics.

The Rice-Mele model extends SSH with staggered onsite potentials, making it a
simple model for dimerized chains with inversion-symmetry breaking.

## Tight-Binding Hamiltonians

A one-dimensional tight-binding Hamiltonian describes hopping between neighboring sites plus onsite potentials. In single-particle form, it is a matrix whose diagonal entries are onsite energies and whose off-diagonal entries are hopping amplitudes.

The square-lattice tight-binding builder extends this to a rectangular two-dimensional lattice with row-major site ordering and optional periodic boundaries in each direction.

The Harper-Hofstadter square-lattice builder adds magnetic Peierls phases in
Landau gauge. Horizontal hoppings are real and vertical hoppings carry
`exp(2 pi i flux * col)`.

The Aubry-Andre-Harper chain adds a quasiperiodic onsite potential:

```text
V_i = lambda cos(2 pi beta i + phase)
```

It is useful for localization studies and for testing spectral algorithms on structured but nonuniform single-particle Hamiltonians.

The Haldane honeycomb-lattice builder provides a finite two-sublattice Chern
insulator model with complex next-nearest-neighbor hoppings. The triangular and
kagome builders provide additional non-square single-particle lattice geometries
useful for flat-band and frustrated-lattice examples.

These matrices are useful because they are easy to inspect, diagonalize, and plot while still representing real lattice physics.

## Exact Diagonalization

Exact diagonalization constructs the Hamiltonian matrix directly and computes its eigenvalues and eigenvectors using dense linear algebra. It gives direct access to spectra, ground states, gaps, observables, and eigenstate structure.

The cost grows quickly. For spin chains, the Hilbert-space dimension is `2**n_sites`, so dense exact diagonalization is only appropriate for small systems. The benefit is clarity: every matrix element is explicit.

Sparse builders and sparse eigensolvers are provided for selected tight-binding
and Hubbard workflows. They reduce memory pressure for matrix construction and
low-energy calculations, but they do not change the exponential basis scaling
of many-body Hubbard models.

## Quantum Algorithm Testbeds

These models are useful for quantum algorithm work because they are structured but still small enough to verify classically.

- VQE: prepare parameterized states and compare variational energies against exact ground energies.
- QPE: use known spectra to test phase-estimation workflows.
- QSVT: study polynomial transformations of Hamiltonian spectra on controlled examples.
- Quantum walks: use tight-binding matrices as graph or lattice walk generators.
- Quantum simulation: compare Trotterized or block-encoded dynamics against exact dense evolution.

The honest role of this package is to provide small, transparent reference problems. It does not claim speedups by itself.
