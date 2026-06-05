from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))

import matplotlib.pyplot as plt

from quantum_lattice_models.models import bose_hubbard_chain
from quantum_lattice_models.plotting import plot_lattice_spectrum


def main() -> None:
    output_dir = Path("images")
    output_dir.mkdir(exist_ok=True)

    H = bose_hubbard_chain(
        n_sites=3,
        hopping=0.6,
        interaction=1.5,
        chemical_potential=0.2,
        max_occupancy=2,
    )
    ax = plot_lattice_spectrum(H, color="tab:green")
    ax.set_title("Truncated Bose-Hubbard spectrum")
    ax.figure.tight_layout()
    ax.figure.savefig(output_dir / "bose_hubbard_spectrum.png", dpi=160)
    plt.close(ax.figure)


if __name__ == "__main__":
    main()
