from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))

import matplotlib.pyplot as plt

from quantum_lattice_models.models import haldane_honeycomb_lattice
from quantum_lattice_models.plotting import plot_spectrum


def main() -> None:
    output_dir = Path("images")
    output_dir.mkdir(exist_ok=True)

    H = haldane_honeycomb_lattice(
        n_rows=3,
        n_cols=3,
        t1=1.0,
        t2=0.18,
        phi=1.5707963267948966,
        sublattice_potential=0.1,
    )
    ax = plot_spectrum(H, highlight_gap=True, zero_line=True, color="#E69F00")
    ax.set_title("Finite Haldane honeycomb spectrum")
    ax.figure.tight_layout()
    ax.figure.savefig(output_dir / "haldane_spectrum.png", dpi=160)
    plt.close(ax.figure)


if __name__ == "__main__":
    main()
