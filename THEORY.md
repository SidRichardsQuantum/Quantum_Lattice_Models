# Theory

This document gives short, working notes for the models and methods implemented in the package. It is not a full textbook treatment.

## Transverse-Field Ising Model

The transverse-field Ising chain is a spin-1/2 model with nearest-neighbor `ZZ` interactions and a transverse `X` field:

```text
H = -J sum_i Z_i Z_{i+1} - h sum_i X_i
```

The interaction term favors aligned or anti-aligned computational-basis states depending on the sign of `J`. The transverse field mixes computational-basis states and drives quantum fluctuations. This makes the model a standard testbed for exact diagonalization, phase-transition intuition, VQE ansatz studies, and quantum simulation circuits.

## Heisenberg Chain

The anisotropic Heisenberg chain includes `XX`, `YY`, and `ZZ` couplings:

```text
H = sum_i (Jx X_i X_{i+1} + Jy Y_i Y_{i+1} + Jz Z_i Z_{i+1})
    + field sum_i Z_i
```

Special choices recover familiar limits. `Jx = Jy = Jz` gives the isotropic Heisenberg model. Unequal couplings give anisotropic variants. The model is useful for testing symmetry-aware methods, spin transport ideas, exact diagonalization, and variational workflows.

## SSH Model

The Su-Schrieffer-Heeger model is a one-dimensional tight-binding model with two sites per unit cell and alternating hopping amplitudes `t1` and `t2`.

For open boundaries and `t1 < t2`, the finite chain supports near-zero-energy states localized near the two ends. In this package, the SSH Hamiltonian is a single-particle matrix of size `2 * n_cells`, not a many-body qubit Hamiltonian.

The SSH model is a compact way to study band structure, boundary modes, topology in simple lattice systems, and quantum walk dynamics.

## Tight-Binding Hamiltonians

A one-dimensional tight-binding Hamiltonian describes hopping between neighboring sites plus onsite potentials. In single-particle form, it is a matrix whose diagonal entries are onsite energies and whose off-diagonal entries are hopping amplitudes.

These matrices are useful because they are easy to inspect, diagonalize, and plot while still representing real lattice physics.

## Exact Diagonalization

Exact diagonalization constructs the Hamiltonian matrix directly and computes its eigenvalues and eigenvectors using dense linear algebra. It gives direct access to spectra, ground states, gaps, observables, and eigenstate structure.

The cost grows quickly. For spin chains, the Hilbert-space dimension is `2**n_sites`, so dense exact diagonalization is only appropriate for small systems. The benefit is clarity: every matrix element is explicit.

## Quantum Algorithm Testbeds

These models are useful for quantum algorithm work because they are structured but still small enough to verify classically.

- VQE: prepare parameterized states and compare variational energies against exact ground energies.
- QPE: use known spectra to test phase-estimation workflows.
- QSVT: study polynomial transformations of Hamiltonian spectra on controlled examples.
- Quantum walks: use tight-binding matrices as graph or lattice walk generators.
- Quantum simulation: compare Trotterized or block-encoded dynamics against exact dense evolution.

The honest role of this package is to provide small, transparent reference problems. It does not claim speedups by itself.
