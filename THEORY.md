# Theory and Shared Conventions

This document describes concepts shared across Quantum Lattice Models. Concise
Hamiltonians, variables, structures, package examples, and validation notes for
individual models live in the [model reference](docs/models/index.md).

## Represented systems

The package uses several distinct matrix representations:

- Spin-$\tfrac12$ models use the computational basis
  $|s_0s_1\ldots s_{N-1}\rangle$ and have dimension $2^N$.
- Single-particle tight-binding models use one basis state $|i\rangle$ per
  site or orbital; their dimension is the number of sites or orbitals.
- Bose-Hubbard models use truncated local occupations
  $|0\rangle,\ldots,|n_{\max}\rangle$ and dimension
  $(n_{\max}+1)^N$.
- Spinful Fermi-Hubbard models use binary occupation of $2N$ ordered
  spin-orbitals and dimension $2^{2N}$.
- Bogoliubov-de Gennes models use a doubled Nambu single-particle basis. They
  are not many-body occupation-basis Hamiltonians.

These representations are not interchangeable. Model reference pages state
the basis and dimension explicitly.

## Operators and normalization

Spin builders use Pauli matrices $X$, $Y$, and $Z$ directly. Consequently,
$X_iX_j+Y_iY_j+Z_iZ_j$ is four times the conventional
$\mathbf S_i\cdot\mathbf S_j$ expression for spin operators
$\mathbf S=\mathbf P/2$. Couplings in the package multiply Pauli products as
written in each model page.

Site index zero is the leftmost tensor factor in Kronecker-product spin
operators. Fermi-Hubbard orbitals are ordered

$$
(0\uparrow,0\downarrow,1\uparrow,1\downarrow,\ldots).
$$

## Tight-binding sign convention

A generic Hermitian hopping term is represented as

$$
H_{ij}=v,\qquad H_{ji}=v^*.
$$

Named tight-binding builders conventionally use $v=-t$ for a real hopping
parameter $t$. Custom three-item bond records use the supplied matrix element
directly; two-item records use `-hopping`.

Complex hoppings therefore encode phases without an additional implicit sign.
The relevant gauge convention is documented on each topological-model page.

## Boundary conditions

Open boundaries include only bonds that remain inside the finite geometry.
Periodic boundaries reconnect opposite ends or edges. Chain models use
`periodic`; two-dimensional finite lattices use independent `periodic_x` and
`periodic_y` flags.

Periodic finite real-space matrices are distinct from Bloch Hamiltonians
$H(\mathbf k)$. The periodic-lattice API constructs Bloch Hamiltonians in
documented cell or orbital gauges and supports reciprocal vectors, momentum
paths, bands, Berry curvature, Wilson loops, and reference invariants.

## Dense and sparse matrices

Dense matrices are NumPy arrays or metadata-bearing NumPy subclasses. Selected
lattice and Hubbard models also provide SciPy CSR builders. Matching dense and
sparse lattice builders share construction paths and are tested for numerical
equivalence.

Sparse storage reduces memory when the matrix contains few nonzero entries,
but it does not change exponential Hilbert-space growth in many-body models.
Use `estimate_model_dimension`, `estimate_dense_memory`, and `diagnose_matrix`
before constructing larger systems.

For XXZ and Heisenberg chains with $J_x=J_y$, fixed-magnetization builders
restrict the basis to a total Pauli-$Z$ eigenvalue

$$
M=\sum_i Z_i=N-2n_1.
$$

The reduced dimension is $\binom{N}{(N-M)/2}$. Reduced basis states are stored
as their integer indices in the full computational basis, allowing explicit
projection, embedding, and full-space block validation.

## Exact diagonalization

Exact diagonalization solves

$$
H|\psi_n\rangle=E_n|\psi_n\rangle.
$$

Hermitian matrices use Hermitian eigensolvers. Full dense diagonalization gives
all eigenpairs; sparse routines can request only low-energy eigenvalues or the
ground state. Iterative sparse results should not be interpreted as complete
spectra.

The package is intended for transparent small-system calculations and
reference workflows, not large-scale many-body simulation.

## Observables and states

State vectors are expected to use the same basis ordering as their
Hamiltonian. `expectation` computes $\langle\psi|O|\psi\rangle`;
spin helpers provide site-resolved and total $Z$ magnetization, same-axis
correlation matrices, connected correlations, and static structure factors.
These routines accept either the full computational basis or a
`FixedMagnetizationBasis`.

Reduced density matrices group amplitudes by subsystem and environment bit
patterns. Sector states are handled directly from their full-basis integer
labels rather than first allocating a vector of length $2^N$. Bipartite
entanglement uses the von Neumann entropy

$$
S_A=-\operatorname{Tr}(\rho_A\log \rho_A).
$$

The default logarithm base is two, so Bell-state entropy is one bit.
Localization helpers provide inverse participation ratios or SSH edge weights.

Probability plots normalize $|\psi_i|^2$ by default. Phase-resolved lattice
plots use the complex argument of each amplitude.

## Symmetries and diagnostics

Hermiticity is checked through $H=H^\dagger$. BdG particle-hole diagnostics
currently test spectral pairing $E\leftrightarrow-E$ rather than constructing
the full antiunitary symmetry operator. Numerical tolerances must be chosen for
the problem scale.

Further parity, translation-sector, commutator, degeneracy, and topological
diagnostics are listed in [ROADMAP.md](ROADMAP.md).

## Portable model specifications

`ModelSpec` records a registered model family, parameters, basis, requested
dense or sparse representation, and optional `LatticeSpec`. JSON encodes
complex values explicitly. Loading a specification validates it against the
registry before reconstruction.

`HamiltonianResult` preserves matrix construction metadata, including reduced
spin-sector basis information where applicable.

## Quantum-algorithm testbeds

Small lattice Hamiltonians provide classically verifiable inputs for VQE, QPE,
Hamiltonian simulation, QSVT, and quantum-walk prototypes. The package supplies
reference matrices and conventions; it does not claim quantum advantage.

## Model reference

Browse the [model index](docs/models/index.md) for per-model equations,
variables, structures, examples, limitations, and validation status.
