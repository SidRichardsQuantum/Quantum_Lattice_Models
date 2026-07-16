"""Persistence helpers for Hamiltonian matrices and their metadata."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models.analysis import AnalysisResult, save_analysis_result
from quantum_lattice_models.diagnostics import is_hermitian
from quantum_lattice_models.physical import (
    BasisIndexMapping,
    InteractionTerm,
    LocalDegreeOfFreedom,
    SymmetryAction,
)
from quantum_lattice_models.specs import EXTERNAL_MATRIX_FAMILY, LatticeSpec, ModelSpec
from quantum_lattice_models.types import HamiltonianResult

_IMPORT_METADATA_FIELDS = {
    "basis",
    "basis_dimension",
    "lattice",
    "local_degrees",
    "basis_mappings",
    "interactions",
    "symmetry_actions",
    "units",
    "conventions",
    "references",
    "provenance",
    "metadata",
}
EXPORT_ARTIFACTS = ("matrix", "metadata", "model", "lattice", "bundle")


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
    representation = "sparse" if sp.issparse(matrix) else "dense"
    expected_shape = metadata.get("matrix_shape")
    if expected_shape is not None and list(matrix.shape) != expected_shape:
        raise ValueError("Stored matrix shape does not match its metadata.")
    if metadata.get("representation") != representation:
        raise ValueError("Stored matrix representation does not match its metadata.")
    if model.representation != representation:
        raise ValueError("Stored matrix representation does not match its model specification.")
    if metadata.get("basis") != model.basis:
        raise ValueError("Stored matrix basis does not match its model specification.")
    if model.family == EXTERNAL_MATRIX_FAMILY:
        expected_dimension = model.metadata.get("matrix_dimension")
        if expected_dimension != matrix.shape[0]:
            raise ValueError("Stored external matrix dimension does not match its metadata.")
    return HamiltonianResult(
        matrix=matrix,
        model=model,
        basis=str(metadata["basis"]),
        representation=representation,
        metadata=dict(metadata.get("metadata", {})),
    )


def import_hamiltonian(
    path: str | Path,
    *,
    metadata_path: str | Path,
    require_hermitian: bool = True,
) -> HamiltonianResult:
    """Import an external NPY or NPZ matrix with explicit portable metadata."""

    source = Path(path)
    metadata_source = Path(metadata_path)
    decoded = json.loads(metadata_source.read_text())
    if not isinstance(decoded, dict):
        raise ValueError("External matrix metadata must contain a JSON object.")
    unknown = sorted(set(decoded) - _IMPORT_METADATA_FIELDS)
    if unknown:
        raise ValueError(f"Unknown external matrix metadata fields: {', '.join(unknown)}.")

    basis = decoded.get("basis")
    if not isinstance(basis, str) or not basis.strip():
        raise ValueError("External matrix metadata requires a nonempty 'basis' string.")
    lattice_data = decoded.get("lattice")
    if lattice_data is not None and not isinstance(lattice_data, dict):
        raise ValueError("External matrix metadata 'lattice' must be an object or null.")
    lattice = None if lattice_data is None else LatticeSpec.from_dict(lattice_data)
    local_degrees = tuple(
        LocalDegreeOfFreedom.from_dict(record)
        for record in _metadata_records(decoded.get("local_degrees", []), "local_degrees")
    )
    basis_mappings = tuple(
        BasisIndexMapping.from_dict(record)
        for record in _metadata_records(decoded.get("basis_mappings", []), "basis_mappings")
    )
    interactions = tuple(
        InteractionTerm.from_dict(record)
        for record in _metadata_records(decoded.get("interactions", []), "interactions")
    )
    symmetry_actions = tuple(
        SymmetryAction.from_dict(record)
        for record in _metadata_records(decoded.get("symmetry_actions", []), "symmetry_actions")
    )

    matrix = _load_external_matrix(source)
    basis_dimension = decoded.get("basis_dimension")
    if basis_dimension is not None and (
        not isinstance(basis_dimension, int) or isinstance(basis_dimension, bool)
    ):
        raise ValueError("External matrix metadata 'basis_dimension' must be an integer.")
    _validate_external_matrix(
        matrix,
        basis_dimension=basis_dimension,
        require_hermitian=require_hermitian,
    )
    representation = "sparse" if sp.issparse(matrix) else "dense"
    dimension = int(matrix.shape[0])
    provenance = _object_mapping(decoded.get("provenance", {}), "provenance")
    provenance = {
        **provenance,
        "import": {
            "format": source.suffix.lower().lstrip("."),
            "source": str(source),
            "metadata_source": str(metadata_source),
        },
    }
    user_metadata = _object_mapping(decoded.get("metadata", {}), "metadata")
    model = ModelSpec(
        family=EXTERNAL_MATRIX_FAMILY,
        basis=basis,
        representation=representation,
        lattice=lattice,
        local_degrees=local_degrees,
        basis_mappings=basis_mappings,
        interactions=interactions,
        symmetry_actions=symmetry_actions,
        units=_string_mapping(decoded.get("units", {}), "units"),
        conventions=_string_mapping(decoded.get("conventions", {}), "conventions"),
        references=_string_list(decoded.get("references", []), "references"),
        provenance=provenance,
        metadata={**user_metadata, "matrix_dimension": dimension},
    )
    model.validate()
    return HamiltonianResult(
        matrix=matrix,
        model=model,
        basis=basis,
        representation=representation,
        metadata={
            "imported": True,
            "source_format": source.suffix.lower().lstrip("."),
            "hermitian": is_hermitian(matrix),
        },
    )


def export_hamiltonian_artifact(
    source: HamiltonianResult | ModelSpec,
    path: str | Path,
    *,
    artifact: str = "matrix",
    format: str = "npz",
    analyses: tuple[AnalysisResult, ...] = (),
) -> tuple[Path, ...]:
    """Export one portable artifact or a deterministic directory bundle."""

    if artifact not in EXPORT_ARTIFACTS:
        raise ValueError(f"Artifact must be one of: {', '.join(EXPORT_ARTIFACTS)}.")
    if artifact != "bundle" and analyses:
        raise ValueError("Analysis records can only be attached to bundle exports.")
    output = Path(path)
    if artifact == "model":
        model = source.model if isinstance(source, HamiltonianResult) else source
        return (model.save(_with_suffix(output, "json")),)
    if artifact == "lattice":
        model = source.model if isinstance(source, HamiltonianResult) else source
        if model.lattice is None:
            raise ValueError("Export source does not contain portable lattice data.")
        return (_write_json(_with_suffix(output, "json"), model.lattice.to_dict()),)
    if not isinstance(source, HamiltonianResult):
        raise ValueError(f"Artifact {artifact!r} requires a constructed Hamiltonian result.")
    result = source
    if artifact == "matrix":
        return (save_hamiltonian(result, output, format=format),)
    if artifact == "metadata":
        return (_write_json(_with_suffix(output, "json"), result.to_metadata()),)

    _prepare_bundle_directory(output)
    matrix = save_hamiltonian(result, output / f"matrix.{format}", format=format)
    model = result.model.save(output / "model.json")
    metadata = _write_json(output / "metadata.json", result.to_metadata())
    paths = [matrix]
    if format == "npy":
        paths.append(metadata_path(matrix))
    paths.extend((model, metadata))
    if result.model.lattice is not None:
        paths.append(_write_json(output / "lattice.json", result.model.lattice.to_dict()))
    if analyses:
        analysis_directory = output / "analyses"
        analysis_directory.mkdir(parents=True, exist_ok=True)
        for index, analysis in enumerate(analyses):
            name = _safe_analysis_name(analysis.analysis)
            paths.append(
                save_analysis_result(
                    analysis,
                    analysis_directory / f"{index:03d}-{name}.npz",
                )
            )
    manifest = _write_json(
        output / "manifest.json",
        {
            "artifact": "quantum-lattice-hamiltonian-bundle",
            "files": [str(item.relative_to(output)) for item in paths],
            "format": format,
            "schema_version": result.model.schema_version,
            "analyses": [
                str(item.relative_to(output)) for item in paths if item.parent.name == "analyses"
            ],
        },
    )
    paths.append(manifest)
    return tuple(paths)


def metadata_path(matrix_path: str | Path) -> Path:
    """Return the deterministic JSON sidecar path for an NPY matrix."""

    path = Path(matrix_path)
    return path.with_suffix(path.suffix + ".json")


def _load_external_matrix(path: Path) -> np.ndarray | sp.csr_matrix:
    if path.suffix.lower() == ".npy":
        return np.asarray(np.load(path, allow_pickle=False))
    if path.suffix.lower() != ".npz":
        raise ValueError("External matrix path must end in .npy or .npz.")

    with np.load(path, allow_pickle=False) as archive:
        keys = set(archive.files)
        if "matrix" in keys:
            return np.asarray(archive["matrix"])
        if {"data", "indices", "indptr", "shape"} <= keys:
            shape_values = np.asarray(archive["shape"]).reshape(-1)
            if shape_values.size != 2:
                raise ValueError("Sparse NPZ shape must contain exactly two dimensions.")
            shape = tuple(int(value) for value in shape_values)
            return sp.csr_matrix(
                (archive["data"], archive["indices"], archive["indptr"]),
                shape=shape,
            )
        if len(keys) == 1:
            return np.asarray(archive[next(iter(keys))])
    raise ValueError(
        "NPZ matrix must contain 'matrix', one unnamed array, or CSR "
        "'data', 'indices', 'indptr', and 'shape' arrays."
    )


def _validate_external_matrix(
    matrix: np.ndarray | sp.spmatrix,
    *,
    basis_dimension: int | None,
    require_hermitian: bool,
) -> None:
    if len(matrix.shape) != 2 or matrix.shape[0] < 1 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("External Hamiltonian matrix must be nonempty and square.")
    values = matrix.data if sp.issparse(matrix) else np.asarray(matrix)
    if not np.issubdtype(values.dtype, np.number):
        raise ValueError("External Hamiltonian matrix must contain numeric values.")
    if not np.all(np.isfinite(values)):
        raise ValueError("External Hamiltonian matrix must contain only finite values.")
    if basis_dimension is not None and basis_dimension != matrix.shape[0]:
        raise ValueError(
            "External matrix dimension must equal metadata basis_dimension when supplied."
        )
    if require_hermitian and not is_hermitian(matrix):
        raise ValueError(
            "External Hamiltonian matrix must be Hermitian; "
            "set require_hermitian=False to import a non-Hermitian matrix."
        )


def _object_mapping(value: object, name: str) -> dict[str, object]:
    if not isinstance(value, dict) or not all(isinstance(key, str) for key in value):
        raise ValueError(f"External matrix metadata '{name}' must be an object.")
    return dict(value)


def _string_mapping(value: object, name: str) -> dict[str, str]:
    mapping = _object_mapping(value, name)
    if not all(isinstance(item, str) for item in mapping.values()):
        raise ValueError(f"External matrix metadata '{name}' values must be strings.")
    return {key: str(item) for key, item in mapping.items()}


def _string_list(value: object, name: str) -> tuple[str, ...]:
    if not isinstance(value, list) or not all(isinstance(item, str) for item in value):
        raise ValueError(f"External matrix metadata '{name}' must be a list of strings.")
    return tuple(value)


def _metadata_records(value: object, name: str) -> list[dict[str, object]]:
    if not isinstance(value, list) or not all(isinstance(item, dict) for item in value):
        raise ValueError(f"External matrix metadata '{name}' must be a list of objects.")
    return value


def _write_json(path: Path, value: object) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    return path


def _prepare_bundle_directory(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for name in (
        "matrix.npy",
        "matrix.npy.json",
        "matrix.npz",
        "model.json",
        "metadata.json",
        "lattice.json",
        "manifest.json",
        "analyses",
    ):
        candidate = path / name
        if candidate.is_dir():
            for item in candidate.iterdir():
                item.unlink()
            candidate.rmdir()
        elif candidate.exists():
            candidate.unlink()


def _with_suffix(path: Path, output_format: str) -> Path:
    suffix = f".{output_format}"
    return path if path.suffix.lower() == suffix else path.with_suffix(suffix)


def _safe_analysis_name(value: str) -> str:
    normalized = "".join(character if character.isalnum() else "-" for character in value.lower())
    return normalized.strip("-") or "analysis"
