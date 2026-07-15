"""Optional ecosystem and lightweight text-format adapters."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from quantum_lattice_models.specs import LatticeSpec, ModelSpec


def export_lattice_dot(lattice: LatticeSpec, path: str | Path) -> Path:
    """Export a deterministic Graphviz DOT interaction graph."""

    output = Path(path)
    lines = ["graph lattice {"]
    for site in range(lattice.n_sites):
        label = lattice.site_labels[site] if lattice.site_labels else str(site)
        lines.append(f'  {site} [label="{_dot_escape(label)}"];')
    for index, bond in enumerate(lattice.bonds):
        label = "" if bond.value is None else f' [label="{bond.value}"]'
        lines.append(f"  {bond.source} -- {bond.target}{label}; // bond {index}")
    lines.append("}")
    output.write_text("\n".join(lines) + "\n")
    return output


def export_table_json(
    coordinates: dict[str, np.ndarray],
    values: dict[str, np.ndarray],
    path: str | Path,
) -> Path:
    """Export named coordinates and values as deterministic plot-data JSON."""

    output = Path(path)
    payload = {
        "coordinates": {name: np.asarray(value).tolist() for name, value in coordinates.items()},
        "values": {name: np.asarray(value).tolist() for name, value in values.items()},
    }
    output.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    return output


def export_model_yaml(model: ModelSpec, path: str | Path) -> Path:
    """Export a model specification through the optional PyYAML adapter."""

    yaml = _yaml()
    output = Path(path)
    output.write_text(yaml.safe_dump(model.to_dict(), sort_keys=True))
    return output


def import_model_yaml(path: str | Path) -> ModelSpec:
    """Load and validate a model specification from optional YAML."""

    yaml = _yaml()
    data = yaml.safe_load(Path(path).read_text())
    if not isinstance(data, dict):
        raise ValueError("YAML model specification must contain a mapping.")
    return ModelSpec.from_dict(data)


def export_lattice_xyz(lattice: LatticeSpec, path: str | Path) -> Path:
    """Export lattice coordinates to XYZ with labels as element-like tokens."""

    if not lattice.positions or len(lattice.positions[0]) not in (2, 3):
        raise ValueError("XYZ export requires 2D or 3D lattice positions.")
    output = Path(path)
    lines = [str(lattice.n_sites), "quantum-lattice-models coordinates"]
    for site, position in enumerate(lattice.positions):
        xyz = tuple(position) + (0.0,) * (3 - len(position))
        label = lattice.site_labels[site] if lattice.site_labels else "X"
        token = "".join(character for character in label if character.isalpha()) or "X"
        lines.append(f"{token} {xyz[0]:.16g} {xyz[1]:.16g} {xyz[2]:.16g}")
    output.write_text("\n".join(lines) + "\n")
    return output


def import_lattice_xyz(path: str | Path) -> LatticeSpec:
    """Import XYZ coordinates without inferring bonds or model semantics."""

    lines = Path(path).read_text().splitlines()
    if len(lines) < 2:
        raise ValueError("XYZ file must contain a count and comment line.")
    n_sites = int(lines[0].strip())
    records = lines[2:]
    if len(records) != n_sites:
        raise ValueError("XYZ coordinate count does not match its header.")
    labels = []
    positions = []
    for record in records:
        fields = record.split()
        if len(fields) < 4:
            raise ValueError("Each XYZ record requires a label and three coordinates.")
        labels.append(fields[0])
        positions.append(tuple(float(value) for value in fields[1:4]))
    return LatticeSpec(
        n_sites=n_sites,
        positions=tuple(positions),
        site_labels=tuple(labels),
        conventions={"xyz_semantics": "coordinates only; bonds and orbitals are unspecified"},
        provenance=({"operation": "import_xyz", "parameters": {"path": str(path)}},),
    )


def from_ase(atoms: Any) -> LatticeSpec:
    """Translate an ASE ``Atoms`` object into coordinate-only portable geometry."""

    positions = np.asarray(atoms.get_positions(), dtype=float)
    labels = tuple(str(symbol) for symbol in atoms.get_chemical_symbols())
    cell = np.asarray(atoms.cell.array, dtype=float).tolist()
    return LatticeSpec(
        n_sites=len(labels),
        positions=tuple(tuple(row) for row in positions),
        site_labels=labels,
        units={"position": "angstrom"},
        conventions={"structure_semantics": "coordinates only; no Hamiltonian inferred"},
        metadata={"cell": cell, "periodic_axes": [bool(value) for value in atoms.pbc]},
        provenance=({"operation": "import_ase"},),
    )


def to_openfermion(model: ModelSpec) -> Any:
    """Translate portable fermion interaction records to ``FermionOperator``."""

    try:
        from openfermion import FermionOperator
    except ImportError as exc:
        raise ImportError(
            "OpenFermion export requires the optional 'openfermion' package."
        ) from exc
    operator = FermionOperator()
    for term in model.interactions:
        factors = []
        for degree, label in zip(term.degrees, term.operators, strict=True):
            if label == "create":
                factors.append((degree, 1))
            elif label == "annihilate":
                factors.append((degree, 0))
            elif label == "number":
                factors.extend(((degree, 1), (degree, 0)))
            else:
                raise ValueError(f"OpenFermion export does not support operator {label!r}.")
        operator += FermionOperator(tuple(factors), term.coefficient)
    return operator


def to_qiskit_sparse_pauli(model: ModelSpec) -> Any:
    """Translate portable spin terms to Qiskit's ``SparsePauliOp``."""

    try:
        from qiskit.quantum_info import SparsePauliOp
    except ImportError as exc:
        raise ImportError("Qiskit export requires the optional 'qiskit' package.") from exc
    n_qubits = len(model.local_degrees)
    labels = []
    coefficients = []
    for term in model.interactions:
        label = ["I"] * n_qubits
        for degree, axis in zip(term.degrees, term.operators, strict=True):
            if axis not in {"I", "X", "Y", "Z"}:
                raise ValueError("Qiskit Pauli export requires spin Pauli operators.")
            label[n_qubits - 1 - degree] = axis
        labels.append("".join(label))
        coefficients.append(term.coefficient)
    return SparsePauliOp(labels, coeffs=coefficients)


def to_qutip_qobj(model: ModelSpec) -> Any:
    """Construct a QuTiP ``Qobj`` while retaining portable model metadata."""

    try:
        import qutip
    except ImportError as exc:
        raise ImportError("QuTiP export requires the optional 'qutip' package.") from exc
    matrix = model.hamiltonian(sparse=False)
    dimensions = [degree.local_dimension for degree in model.local_degrees]
    dims = (
        [dimensions, dimensions] if dimensions and np.prod(dimensions) == matrix.shape[0] else None
    )
    return qutip.Qobj(matrix, dims=dims)


def to_netket_graph(lattice: LatticeSpec) -> Any:
    """Translate finite portable geometry to a NetKet graph."""

    try:
        import netket
    except ImportError as exc:
        raise ImportError("NetKet export requires the optional 'netket' package.") from exc
    edges = sorted({tuple(sorted((bond.source, bond.target))) for bond in lattice.bonds})
    return netket.graph.Graph(edges=edges, n_nodes=lattice.n_sites)


def to_quspin_hamiltonian(model: ModelSpec, **basis_parameters: Any) -> Any:
    """Translate portable spin Pauli terms to a QuSpin Hamiltonian."""

    try:
        from quspin.basis import spin_basis_1d
        from quspin.operators import hamiltonian
    except ImportError as exc:
        raise ImportError("QuSpin export requires the optional 'quspin' package.") from exc
    if not model.local_degrees or any(degree.kind != "spin" for degree in model.local_degrees):
        raise ValueError("QuSpin export currently supports spin models only.")
    grouped: dict[str, list[list[object]]] = {}
    for term in model.interactions:
        operator = "".join(axis.lower() for axis in term.operators)
        grouped.setdefault(operator, []).append(
            [term.coefficient, *[int(index) for index in term.degrees]]
        )
    static = sorted(grouped.items())
    basis = spin_basis_1d(len(model.local_degrees), pauli=True, **basis_parameters)
    return hamiltonian(static, [], basis=basis, dtype=np.complex128)


def _yaml() -> Any:
    try:
        import yaml
    except ImportError as exc:
        raise ImportError("YAML support requires the optional 'yaml' extra.") from exc
    return yaml


def _dot_escape(value: str) -> str:
    return value.replace("\\", "\\\\").replace('"', '\\"')
