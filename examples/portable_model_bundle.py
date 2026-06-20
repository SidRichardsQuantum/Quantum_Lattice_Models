from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from quantum_lattice_models import (
    create_model_spec,
    export_hamiltonian_artifact,
    load_hamiltonian,
)


def main() -> None:
    output_dir = Path("results/examples")
    output_dir.mkdir(parents=True, exist_ok=True)

    model = create_model_spec(
        "ssh_model",
        parameters={"n_cells": 6, "t1": 0.35, "t2": 1.0, "periodic": False},
        units={"t1": "eV", "t2": "eV"},
        provenance={"example": "portable_model_bundle"},
    )
    model_path = model.save(output_dir / "ssh_model.json")
    bundle_dir = output_dir / "ssh_model.bundle"
    outputs = export_hamiltonian_artifact(
        model.build_result(),
        bundle_dir,
        artifact="bundle",
        format="npz",
    )

    restored = load_hamiltonian(bundle_dir / "matrix.npz")
    manifest = json.loads((bundle_dir / "manifest.json").read_text())
    assert np.allclose(restored.matrix, model.hamiltonian())
    assert restored.model == model

    print("Portable model bundle")
    print(f"  model: {model_path}")
    print(f"  files: {', '.join(path.name for path in outputs)}")
    print(f"  manifest format: {manifest['format']}")
    print(f"  restored shape: {restored.shape}")


if __name__ == "__main__":
    main()
