# Model Reference

- [Logical model capability matrix](capabilities.md)
- [Complete public API inventory](api.md)

These pages give concise, implementation-specific references for every model
family currently included in Quantum Lattice Models. Shared basis, operator,
boundary, sparse-matrix, and exact-diagonalization conventions are described in
[Theory and Shared Conventions](../../THEORY.md).

## Spin systems

- [Ising chains](ising.md)
- [Heisenberg chain](heisenberg_chain.md)
- [XY chain](xy_chain.md)
- [XXZ chain](xxz_chain.md)
- [J1-J2 Heisenberg chain](j1_j2_heisenberg.md)
- [Heisenberg ladder](heisenberg_ladder.md)

## Interacting occupation-basis systems

- [Bose-Hubbard chain](bose_hubbard.md)
- [Fermi-Hubbard chain](fermi_hubbard.md)

## One-dimensional single-particle and BdG systems

- [SSH model](ssh.md)
- [Rice-Mele model](rice_mele.md)
- [Tight-binding chain](tight_binding_chain.md)
- [Aubry-André-Harper chain](aubry_andre.md)
- [Kitaev BdG chain](kitaev_bdg.md)

## Two-dimensional finite lattices

- [Additional benchmark lattices](benchmark_lattices.md)
- [Square lattice](square_lattice.md)
- [Triangular lattice](triangular_lattice.md)
- [Kagome lattice](kagome_lattice.md)
- [Harper-Hofstadter lattice](harper_hofstadter.md)
- [Haldane honeycomb lattice](haldane.md)

## User-defined systems

- [Custom lattice and tight-binding models](custom_lattice.md)

Dense and sparse variants share one page because they represent the same
physical model. Parameter tables are generated from the package registry.
