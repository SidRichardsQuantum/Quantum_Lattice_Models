from __future__ import annotations

from quantum_lattice_models import (
    _cli_output,
    _cli_parameters,
    _cli_sources,
    _physical_inference,
    specs,
)
from quantum_lattice_models.cli import main


def test_specification_and_cli_boundaries_preserve_stable_imports() -> None:
    assert specs.ModelSpec.__module__ == "quantum_lattice_models.specs"
    assert specs.create_model_spec.__module__ == "quantum_lattice_models.specs"
    assert callable(_physical_inference.infer_physical_system)
    assert callable(_cli_parameters.parameter_values)
    assert callable(_cli_sources.build_source)
    assert callable(_cli_output.print_json)
    assert callable(main)
