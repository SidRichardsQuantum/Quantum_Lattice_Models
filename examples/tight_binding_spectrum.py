from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))

import matplotlib.pyplot as plt

from quantum_lattice_models.models import tight_binding_chain
from quantum_lattice_models.plotting import plot_spectrum


def main() -> None:
    output_dir = Path("images")
    output_dir.mkdir(exist_ok=True)

    H = tight_binding_chain(n_sites=16, hopping=1.0, onsite=0.0, periodic=False)
    ax = plot_spectrum(H, highlight_gap=True, zero_line=True)
    ax.figure.tight_layout()
    ax.figure.savefig(output_dir / "tight_binding_spectrum.png", dpi=160)
    plt.close(ax.figure)


if __name__ == "__main__":
    main()
