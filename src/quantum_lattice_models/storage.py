"""Persistence helpers for Hamiltonian matrices and their metadata."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models.specs import ModelSpec
from quantum_lattice_models.types import HamiltonianResult


def save_hamiltonian(
    result: HamiltonianResult,
    path: str | Path,
    *,
    format: str | None = None,
) -> Path:
    """Save a Hamiltonian result as NPY plus JSON or as self-contained NPZ."""

    output_format = (format or Path(path).suffix.lstrip(".")).lower()
    if output_format not in {"npy", "npz"}:
        raise ValueError("Hamiltonian format must be 'npy' or 'npz'.")
    output = _with_suffix(Path(path), output_format)
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata = json.dumps(result.to_metadata(), sort_keys=True)

    if output_format == "npy":
        if sp.issparse(result.matrix):
            raise ValueError("NPY export requires a dense matrix; use NPZ for sparse matrices.")
        np.save(output, np.asarray(result.matrix), allow_pickle=False)
        serialized = json.dumps(result.to_metadata(), indent=2, sort_keys=True) + "\n"
        metadata_path(output).write_text(serialized)
        return output

    if sp.issparse(result.matrix):
        matrix = result.matrix.tocsr()
        np.savez_compressed(
            output,
            storage=np.array("csr"),
            data=matrix.data,
            indices=matrix.indices,
            indptr=matrix.indptr,
            shape=np.asarray(matrix.shape, dtype=np.int64),
            metadata=np.array(metadata),
        )
    else:
        np.savez_compressed(
            output,
            storage=np.array("dense"),
            matrix=np.asarray(result.matrix),
            metadata=np.array(metadata),
        )
    return output


def load_hamiltonian(path: str | Path) -> HamiltonianResult:
    """Load a Hamiltonian result saved by :func:`save_hamiltonian`."""

    source = Path(path)
    if source.suffix == ".npy":
        matrix = np.load(source, allow_pickle=False)
        metadata = json.loads(metadata_path(source).read_text())
    elif source.suffix == ".npz":
        with np.load(source, allow_pickle=False) as archive:
            metadata = json.loads(str(archive["metadata"].item()))
            storage = str(archive["storage"].item())
            if storage == "dense":
                matrix = archive["matrix"]
            elif storage == "csr":
                shape = tuple(int(value) for value in archive["shape"])
                matrix = sp.csr_matrix(
                    (archive["data"], archive["indices"], archive["indptr"]),
                    shape=shape,
                )
            else:
                raise ValueError(f"Unsupported matrix storage {storage!r}.")
    else:
        raise ValueError("Hamiltonian path must end in .npy or .npz.")

    model = ModelSpec.from_dict(metadata["model"])
    return HamiltonianResult(
        matrix=matrix,
        model=model,
        basis=str(metadata["basis"]),
        representation=str(metadata["representation"]),
        metadata=dict(metadata.get("metadata", {})),
    )


def metadata_path(matrix_path: str | Path) -> Path:
    """Return the deterministic JSON sidecar path for an NPY matrix."""

    path = Path(matrix_path)
    return path.with_suffix(path.suffix + ".json")


def _with_suffix(path: Path, output_format: str) -> Path:
    suffix = f".{output_format}"
    return path if path.suffix.lower() == suffix else path.with_suffix(suffix)
