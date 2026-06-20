from __future__ import annotations

import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))

import matplotlib.pyplot as plt

from quantum_lattice_models import (
    SpinField,
    SpinInteraction,
    create_graph_spin_spec,
    create_model_spec,
)
from quantum_lattice_models.plotting import plot_interaction_graph


def main() -> None:
    output_dir = Path("images")
    output_dir.mkdir(exist_ok=True)

    xxz = create_model_spec(
        "xxz_chain",
        parameters={"n_sites": 4, "coupling": 1.0, "anisotropy": 0.6, "field": 0.2},
    )
    fermi = create_model_spec(
        "fermi_hubbard_chain",
        parameters={"n_sites": 3, "hopping": 0.7, "interaction": 2.0},
    )
    graph_spin = create_graph_spin_spec(
        4,
        interactions=(
            SpinInteraction(0, 1, "X", "Y", 0.35),
            SpinInteraction(1, 2, "Z", "Z", -0.8),
            SpinInteraction(2, 3, "Y", "X", 0.25),
            SpinInteraction(3, 0, "Z", "X", 0.4),
        ),
        fields=(SpinField(0, "X", -0.3), SpinField(2, "Z", 0.2)),
        positions=((0.0, 0.0), (1.0, 0.0), (1.0, 1.0), (0.0, 1.0)),
        site_labels=("A", "B", "C", "D"),
    )

    figure, axes = plt.subplots(1, 3, figsize=(13.5, 4.0))
    plot_interaction_graph(xxz, ax=axes[0], show_coefficients=False)
    axes[0].set_title("XXZ spin chain")
    plot_interaction_graph(fermi, ax=axes[1], show_coefficients=False)
    axes[1].set_title("Fermi-Hubbard modes")
    plot_interaction_graph(graph_spin, ax=axes[2], show_coefficients=False)
    axes[2].set_title("User-defined spin graph")
    figure.tight_layout()
    figure.savefig(output_dir / "physical_interaction_graph.png", dpi=160)
    plt.close(figure)


if __name__ == "__main__":
    main()
