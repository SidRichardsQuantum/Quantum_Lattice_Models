from __future__ import annotations

import json
import runpy

import pytest

from quantum_lattice_models import (
    LatticeSpec,
    export_graphml,
    export_lattice_csv,
    from_networkx,
    import_graphml,
    import_lattice_csv,
    to_networkx,
)
from quantum_lattice_models.cli import main
from quantum_lattice_models.lattice import Bond


def _lattice() -> LatticeSpec:
    return LatticeSpec(
        n_sites=3,
        positions=((0.0, 0.0), (1.0, 0.0), (0.5, 0.75)),
        bonds=(Bond(0, 1), Bond(1, 2, 0.25 + 0.5j), Bond(1, 2, -0.1j)),
        site_labels=("left", "middle", "right"),
        orbital_labels=("s", "p", "s"),
        sublattice_labels=("A", "B", "A"),
        unit_cells=(0, 0, 1),
        boundary_conditions={"x": "open"},
        units={"position": "angstrom", "energy": "eV"},
        conventions={"edge_interpretation": "directed matrix elements"},
        references=("https://example.test/model",),
        provenance=({"operation": "created", "parameters": {"source": "test"}},),
        metadata={"name": "triangle", "tags": ("test", "complex")},
    )


def test_csv_lattice_round_trip_preserves_tables_and_metadata(tmp_path) -> None:
    original = _lattice()
    sites, bonds, metadata = export_lattice_csv(
        original,
        tmp_path / "sites.csv",
        tmp_path / "bonds.csv",
    )

    restored = import_lattice_csv(sites, bonds)

    assert metadata == tmp_path / "sites.csv.json"
    assert restored == original
    assert json.loads(metadata.read_text())["units"]["energy"] == "eV"


def test_networkx_and_graphml_round_trip(tmp_path) -> None:
    pytest.importorskip("networkx")
    original = _lattice()

    graph = to_networkx(original)
    restored_graph = from_networkx(graph)
    graphml = export_graphml(original, tmp_path / "lattice.graphml")
    restored_graphml = import_graphml(graphml)

    assert graph.is_directed()
    assert graph.is_multigraph()
    assert restored_graph == original
    assert restored_graphml == original


def test_networkx_weighted_chain_case_study() -> None:
    pytest.importorskip("networkx")
    namespace = runpy.run_path("case_studies/networkx_weighted_chain.py")
    assert namespace["model"].lattice.n_sites == 3
    assert namespace["report"].supported


def test_csv_validation_rejects_incomplete_coordinates(tmp_path) -> None:
    sites = tmp_path / "sites.csv"
    bonds = tmp_path / "bonds.csv"
    sites.write_text(
        "site,x,y,z,site_label,orbital_label,sublattice_label,unit_cell\n"
        "0,0,0,,,,,\n"
        "1,1,,,,,,\n"
    )
    bonds.write_text("source,target,has_value,value_real,value_imag\n")

    with pytest.raises(ValueError, match="x and y"):
        import_lattice_csv(sites, bonds)


def test_cli_csv_import_and_export_round_trip(tmp_path, capsys) -> None:
    original = _lattice()
    sites, bonds, metadata = export_lattice_csv(
        original,
        tmp_path / "input-sites.csv",
        tmp_path / "input-bonds.csv",
    )
    model_path = tmp_path / "model.json"
    assert (
        main(
            [
                "import",
                str(sites),
                "--format",
                "csv",
                "--bonds",
                str(bonds),
                "--metadata",
                str(metadata),
                "--output",
                str(model_path),
            ]
        )
        == 0
    )
    capsys.readouterr()

    output_sites = tmp_path / "output-sites.csv"
    output_bonds = tmp_path / "output-bonds.csv"
    assert (
        main(
            [
                "export-lattice",
                str(model_path),
                "--format",
                "csv",
                "--sites",
                str(output_sites),
                "--bonds",
                str(output_bonds),
            ]
        )
        == 0
    )
    capsys.readouterr()

    assert import_lattice_csv(output_sites, output_bonds) == original
