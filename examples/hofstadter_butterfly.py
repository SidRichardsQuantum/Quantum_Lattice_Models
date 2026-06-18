from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))

import matplotlib.pyplot as plt
import numpy as np

from quantum_lattice_models.models import harper_hofstadter_square_lattice
from quantum_lattice_models.plotting import plot_parameter_sweep


def main() -> None:
    output_dir = Path("images")
    output_dir.mkdir(exist_ok=True)

    fluxes = np.linspace(0.0, 1.0, 121)
    ax = plot_parameter_sweep(
        lambda flux: harper_hofstadter_square_lattice(n_rows=8, n_cols=8, flux=flux),
        fluxes,
        parameter_name="Flux per plaquette",
        title="Finite-lattice Hofstadter spectrum (8×8)",
        s=2.5,
        color="#CC79A7",
        alpha=0.5,
    )
    ax.figure.tight_layout()
    ax.figure.savefig(output_dir / "hofstadter_butterfly.png", dpi=160)
    plt.close(ax.figure)


if __name__ == "__main__":
    main()
