from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))

import matplotlib.pyplot as plt

from quantum_lattice_models.geometry import kagome_lattice_positions
from quantum_lattice_models.models import kagome_lattice_tight_binding
from quantum_lattice_models.plotting import plot_lattice_graph


def main() -> None:
    output_dir = Path("images")
    output_dir.mkdir(exist_ok=True)

    n_rows = 2
    n_cols = 3
    H = kagome_lattice_tight_binding(n_rows=n_rows, n_cols=n_cols)
    positions = kagome_lattice_positions(n_rows=n_rows, n_cols=n_cols)
    ax = plot_lattice_graph(H, positions, node_size=42, color="0.45")
    ax.set_title("Kagome lattice connectivity")
    ax.figure.tight_layout()
    ax.figure.savefig(output_dir / "kagome_graph.png", dpi=160)
    plt.close(ax.figure)


if __name__ == "__main__":
    main()
