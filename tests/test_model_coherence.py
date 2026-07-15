from __future__ import annotations

import json

import numpy as np
import scipy.sparse as sp

from quantum_lattice_models import create_model_spec, list_models
from quantum_lattice_models.cli import main


def test_every_registered_name_has_normalized_physical_records() -> None:
    for name in list_models(include_aliases=True):
        model = create_model_spec(name)
        assert model.lattice is not None, name
        assert model.local_degrees, name
        assert model.basis_mappings, name
        assert model.interactions, name
        assert len(model.basis_mappings) == len(model.local_degrees), name


def test_all_default_single_particle_records_match_builders() -> None:
    for name in list_models(include_aliases=True):
        model = create_model_spec(name)
        if not all(degree.kind == "orbital" for degree in model.local_degrees):
            continue
        symbolic = np.zeros((len(model.local_degrees), len(model.local_degrees)), dtype=complex)
        for term in model.interactions:
            if term.operators == ("number",):
                symbolic[term.degrees[0], term.degrees[0]] += term.coefficient
            elif term.operators == ("create", "annihilate"):
                source, target = term.degrees
                symbolic[source, target] += term.coefficient
                if term.metadata.get("hermitian_conjugate"):
                    symbolic[target, source] += complex(term.coefficient).conjugate()
        built = model.hamiltonian(sparse=model.representation == "sparse")
        dense = built.toarray() if sp.issparse(built) else np.asarray(built)
        assert np.allclose(symbolic, dense), name


def test_intake_cli_commands_report_portable_semantics(tmp_path, capsys) -> None:
    path = create_model_spec("tight_binding_chain", parameters={"n_sites": 3}).save(
        tmp_path / "model.json"
    )

    assert main(["describe", str(path), "--json"]) == 0
    description = json.loads(capsys.readouterr().out)
    assert description["local_degrees"] == 3
    assert description["lattice_bonds"] == 2

    assert main(["lint", str(path), "--json"]) == 0
    lint = json.loads(capsys.readouterr().out)
    assert lint["valid"]

    assert main(["adapter-capabilities", str(path), "graphml", "--json"]) == 0
    capability = json.loads(capsys.readouterr().out)
    assert capability["supported"]
    assert "geometry" in capability["preserved"]
