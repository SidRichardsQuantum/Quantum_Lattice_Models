from __future__ import annotations

import numpy as np

from quantum_lattice_models import (
    HamiltonianResult,
    LatticeSpec,
    ModelSpec,
    SpinInteraction,
    adapter_capability_report,
    create_graph_spin_spec,
    create_model_spec,
    describe_model,
    lint_model,
)
from quantum_lattice_models.lattice import Bond


def test_describe_model_returns_user_facing_physical_summary() -> None:
    model = create_graph_spin_spec(
        3,
        interactions=(SpinInteraction(0, 1, "X", "X", 1.0),),
        positions=((0.0, 0.0), (1.0, 0.0), (2.0, 0.0)),
        site_labels=("left", "middle", "right"),
    )

    summary = describe_model(model)

    assert summary.family == "graph_spin"
    assert summary.dimension == 8
    assert summary.lattice_sites == 3
    assert summary.local_degrees == 3
    assert summary.interactions == 1
    assert summary.reconstructable
    assert summary.to_dict()["matrix_shape"] is None


def test_lint_model_reports_import_metadata_gaps_and_matrix_warnings() -> None:
    lattice = LatticeSpec(
        n_sites=3,
        bonds=(Bond(0, 1), Bond(1, 2)),
        site_labels=("A", "A", "B"),
    )
    model = create_model_spec("custom_tight_binding", lattice=lattice)
    result = HamiltonianResult(
        matrix=np.array([[0.0, 1.0], [0.0, 0.0]]),
        model=ModelSpec(family="external_matrix", basis="site basis"),
        basis="site basis",
        representation="dense",
        metadata={},
    )

    model_report = lint_model(model)
    matrix_report = lint_model(result)

    assert model_report.valid
    assert "duplicate site labels detected" in model_report.warnings
    assert "lattice positions are missing" in model_report.warnings
    assert "boundary conditions are unspecified" in model_report.warnings
    assert "matrix is not Hermitian" in matrix_report.warnings
    assert (
        "record basis conventions and provenance for external matrices" in matrix_report.suggestions
    )


def test_adapter_capability_reports_preserved_and_lost_semantics() -> None:
    graph_spin = create_graph_spin_spec(
        2,
        interactions=(SpinInteraction(0, 1, "Z", "Z", 1.0),),
        positions=((0.0, 0.0), (1.0, 0.0)),
        site_labels=("A", "B"),
    )
    tight_binding = create_model_spec("tight_binding_chain", parameters={"n_sites": 3})

    qiskit = adapter_capability_report(graph_spin, "qiskit")
    csv = adapter_capability_report(graph_spin, "csv")
    netket = adapter_capability_report(tight_binding, "netket")

    assert qiskit.supported
    assert "interactions" in qiskit.preserved
    assert "geometry" in qiskit.lost
    assert csv.supported
    assert "geometry" in csv.preserved
    assert "interactions" in csv.lost
    assert netket.supported
    assert "bonds" in netket.preserved
