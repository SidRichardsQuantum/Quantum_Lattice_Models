from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))

import matplotlib.pyplot as plt

from quantum_lattice_models.models import haldane_honeycomb_lattice
from quantum_lattice_models.plotting import plot_hamiltonian_matrix


def main() -> None:
    output_dir = Path("images")
    output_dir.mkdir(exist_ok=True)

    H = haldane_honeycomb_lattice(
        n_rows=2,
        n_cols=3,
        t1=1.0,
        t2=0.18,
        phi=1.5707963267948966,
        sublattice_potential=0.1,
    )
    figure, axes = plt.subplots(1, 2, figsize=(9.5, 4.0))
    plot_hamiltonian_matrix(H, ax=axes[0], mode="magnitude")
    axes[0].set_title("Haldane Hamiltonian magnitude")
    plot_hamiltonian_matrix(H, ax=axes[1], mode="phase")
    axes[1].set_title("Haldane hopping phase")
    figure.tight_layout()
    figure.savefig(output_dir / "haldane_matrix_structure.png", dpi=160)
    plt.close(figure)


if __name__ == "__main__":
    main()
