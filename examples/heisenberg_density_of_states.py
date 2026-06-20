from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))

import matplotlib.pyplot as plt

from quantum_lattice_models.models import heisenberg_chain
from quantum_lattice_models.plotting import plot_density_of_states


def main() -> None:
    output_dir = Path("images")
    output_dir.mkdir(exist_ok=True)

    H = heisenberg_chain(n_sites=5, jx=1.0, jy=1.0, jz=1.0, field=0.2)
    ax = plot_density_of_states(H, bins=16)
    ax.figure.tight_layout()
    ax.figure.savefig(output_dir / "heisenberg_density_of_states.png", dpi=160)
    plt.close(ax.figure)


if __name__ == "__main__":
    main()
