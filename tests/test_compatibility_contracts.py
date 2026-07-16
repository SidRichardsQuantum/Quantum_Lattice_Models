from __future__ import annotations

import hashlib
import json
from pathlib import Path

import quantum_lattice_models as qlm
from quantum_lattice_models.cli import main


def test_public_api_0_1_snapshot() -> None:
    snapshot = json.loads(Path("tests/fixtures/public_api_0_1.json").read_text())
    names = "\n".join(sorted(qlm.__all__)) + "\n"

    assert len(qlm.__all__) == snapshot["count"]
    assert len(set(qlm.__all__)) == len(qlm.__all__)
    assert hashlib.sha256(names.encode()).hexdigest() == snapshot["sha256"]
    assert all(hasattr(qlm, name) for name in qlm.__all__)


def test_cli_describe_schema_1_0_golden_output(capsys) -> None:
    expected = json.loads(Path("tests/fixtures/cli_describe_schema_1_0.json").read_text())

    assert (
        main(
            [
                "describe",
                "tests/fixtures/model_schema_1_0.json",
                "--json",
            ]
        )
        == 0
    )
    assert json.loads(capsys.readouterr().out) == expected
