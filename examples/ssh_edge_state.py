from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))

import matplotlib.pyplot as plt
import numpy as np

from quantum_lattice_models.models import ssh_model
from quantum_lattice_models.plotting import plot_lattice_state
from quantum_lattice_models.spectra import eigensystem


def main() -> None:
    output_dir = Path("images")
    output_dir.mkdir(exist_ok=True)

    H = ssh_model(n_cells=12, t1=0.4, t2=1.0, periodic=False)
    values, vectors = eigensystem(H)
    edge_state_index = int(np.argmin(np.abs(values)))
    positions = np.array(
        [(site / 2.0, 0.35 if site % 2 else -0.35) for site in range(24)], dtype=float
    )
    _, ax = plt.subplots(figsize=(9.0, 3.2))
    ax = plot_lattice_state(H, vectors[:, edge_state_index], positions=positions, ax=ax)
    ax.set_aspect("auto", adjustable="box")
    ax.set_ylim(-0.9, 0.9)
    ax.set_title("SSH edge-state probability and phase")
    ax.figure.tight_layout()
    ax.figure.savefig(output_dir / "ssh_edge_state.png", dpi=160)
    plt.close(ax.figure)


if __name__ == "__main__":
    main()
