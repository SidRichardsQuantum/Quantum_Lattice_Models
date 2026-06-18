from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))

import matplotlib.pyplot as plt

from quantum_lattice_models.models import rice_mele_model
from quantum_lattice_models.plotting import plot_lattice_spectrum


def main() -> None:
    output_dir = Path("images")
    output_dir.mkdir(exist_ok=True)

    H = rice_mele_model(n_cells=12, hopping=1.0, dimerization=0.35, staggering=0.4)
    ax = plot_lattice_spectrum(H, highlight_gap=True, zero_line=True)
    ax.set_title("Rice-Mele finite-chain spectrum")
    ax.figure.tight_layout()
    ax.figure.savefig(output_dir / "rice_mele_spectrum.png", dpi=160)
    plt.close(ax.figure)


if __name__ == "__main__":
    main()
