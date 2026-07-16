"""NetworkX intake case study with weighted directed-edge metadata."""

from __future__ import annotations

import networkx as nx

from quantum_lattice_models import (
    adapter_capability_report,
    create_model_spec,
    from_networkx,
)

graph = nx.MultiDiGraph()
graph.add_node(0, x=0.0, y=0.0, site_label="left")
graph.add_node(1, x=1.0, y=0.0, site_label="center")
graph.add_node(2, x=2.0, y=0.0, site_label="right")
graph.add_edge(0, 1, has_value=True, value_real=-1.0, value_imag=0.0)
graph.add_edge(1, 2, has_value=True, value_real=0.0, value_imag=0.25)

lattice = from_networkx(graph)
model = create_model_spec(
    "custom_tight_binding",
    lattice=lattice,
    parameters={"hopping": 1.0, "onsite": 0.0, "hermitian": True},
    provenance={"source": "NetworkX MultiDiGraph case study"},
)
report = adapter_capability_report(model, "graphml")

print(model.summary())
print(report.to_dict())
