from __future__ import annotations

import json
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(".matplotlib-cache").resolve()))

import matplotlib.pyplot as plt
import numpy as np

from quantum_lattice_models import import_hamiltonian, save_hamiltonian
from quantum_lattice_models.plotting import plot_hamiltonian_matrix


def main() -> None:
    artifact_dir = Path("results/examples/external_matrix")
    image_dir = Path("images")
    artifact_dir.mkdir(parents=True, exist_ok=True)
    image_dir.mkdir(exist_ok=True)

    matrix = np.array(
        [
            [0.2, -1.0, 0.0],
            [-1.0, 0.0, 0.35j],
            [0.0, -0.35j, -0.2],
        ],
        dtype=complex,
    )
    source = artifact_dir / "external.npy"
    np.save(source, matrix, allow_pickle=False)

    metadata = {
        "basis": "single-particle site basis",
        "basis_dimension": 3,
        "lattice": {
            "n_sites": 3,
            "positions": [[0.0, 0.0], [1.0, 0.0], [0.5, 0.8]],
            "bonds": [],
            "site_labels": ["left", "right", "top"],
        },
        "local_degrees": [
            {
                "index": site,
                "site": site,
                "kind": "orbital",
                "local_dimension": 2,
                "label": label,
            }
            for site, label in enumerate(("left", "right", "top"))
        ],
        "basis_mappings": [
            {
                "local_degree": site,
                "basis_index": site,
                "role": "single_particle_state",
            }
            for site in range(3)
        ],
        "interactions": [
            {
                "degrees": [0],
                "operators": ["number"],
                "coefficient": {"__complex__": [0.2, 0.0]},
                "kind": "onsite",
            },
            {
                "degrees": [2],
                "operators": ["number"],
                "coefficient": {"__complex__": [-0.2, 0.0]},
                "kind": "onsite",
            },
            {
                "degrees": [0, 1],
                "operators": ["create", "annihilate"],
                "coefficient": {"__complex__": [-1.0, 0.0]},
                "kind": "hopping",
                "metadata": {"hermitian_conjugate": True},
            },
            {
                "degrees": [1, 2],
                "operators": ["create", "annihilate"],
                "coefficient": {"__complex__": [0.0, 0.35]},
                "kind": "hopping",
                "metadata": {"hermitian_conjugate": True},
            },
        ],
        "units": {"energy": "eV"},
        "conventions": {"basis_ordering": "left, right, top"},
        "provenance": {"generator": "external_matrix_import example"},
    }
    metadata_path = artifact_dir / "metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")

    imported = import_hamiltonian(source, metadata_path=metadata_path)
    portable = save_hamiltonian(imported, artifact_dir / "portable.npz")

    figure, axes = plt.subplots(1, 2, figsize=(8.5, 3.8))
    plot_hamiltonian_matrix(imported.matrix, ax=axes[0], mode="magnitude")
    plot_hamiltonian_matrix(imported.matrix, ax=axes[1], mode="phase")
    axes[0].set_title("Imported magnitude")
    axes[1].set_title("Imported phase")
    figure.tight_layout()
    figure.savefig(image_dir / "external_matrix_import.png", dpi=160)
    plt.close(figure)

    print("External matrix import")
    print(f"  source: {source}")
    print(f"  portable: {portable}")
    print(f"  basis: {imported.basis}")
    print(f"  local degrees: {len(imported.model.local_degrees)}")
    print(f"  interactions: {len(imported.model.interactions)}")


if __name__ == "__main__":
    main()
