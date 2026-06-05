from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))

import matplotlib.pyplot as plt

from quantum_lattice_models.models import transverse_field_ising
from quantum_lattice_models.plotting import plot_spectrum


def main() -> None:
    output_dir = Path("images")
    output_dir.mkdir(exist_ok=True)

    H = transverse_field_ising(n_sites=5, j=1.0, h=0.7, periodic=False)
    ax = plot_spectrum(H)
    ax.figure.tight_layout()
    ax.figure.savefig(output_dir / "ising_spectrum.png", dpi=160)
    plt.close(ax.figure)


if __name__ == "__main__":
    main()
