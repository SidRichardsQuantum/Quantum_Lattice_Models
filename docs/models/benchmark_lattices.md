<!-- builders: graphene_lattice, graphene_lattice_sparse, anderson_square_lattice, anderson_square_lattice_sparse, checkerboard_chern_insulator, checkerboard_chern_insulator_sparse, dice_lattice, dice_lattice_sparse -->
# Additional Benchmark Lattices

## Graphene

`graphene_lattice` is the nearest-neighbor limit of the finite honeycomb
construction. With zero onsite terms it is bipartite and its finite spectrum
is symmetric around zero.

## Two-dimensional Anderson model

`anderson_square_lattice` adds reproducible uniform onsite disorder in
$[-W/2,W/2]$ to a square-lattice hopping Hamiltonian. The random seed and
boundary conditions are explicit parameters.

## Checkerboard Chern insulator

`checkerboard_chern_insulator` uses a two-orbital Qi-Wu-Zhang-type real-space
representation,

$$
H(\mathbf k)=t\sin k_x\,\sigma_x+t\sin k_y\,\sigma_y+
(m+t\cos k_x+t\cos k_y)\sigma_z.
$$

The orbital ordering is row-major unit cells with two orbitals per cell.

## Dice or $T_3$ lattice

`dice_lattice` has one hub and two rim sites per unit cell. Its bipartite
imbalance produces a zero-energy flat-band subspace in finite open systems.

## Package use

```python
from quantum_lattice_models import (
    anderson_square_lattice_sparse,
    checkerboard_chern_insulator,
    dice_lattice,
    graphene_lattice,
)

graphene = graphene_lattice(3, 4)
disordered = anderson_square_lattice_sparse(12, 12, disorder=4.0, seed=7)
chern = checkerboard_chern_insulator(4, 4, mass=1.0)
dice = dice_lattice(3, 4)
```

## Parameters

{{PARAMETERS}}

## Validation and cautions

Dense and sparse builders are cross-checked. The package validates graphene
spectral symmetry, Anderson reproducibility, checkerboard Hermiticity, and the
dice zero-energy flat-band subspace. Finite models do not by themselves imply
thermodynamic or bulk topological conclusions.
