PYTHON ?= .venv/bin/python

.PHONY: format format-check lint typecheck test check examples docs-assets docs-models docs-models-check notebook-summary

format:
	find src tests examples scripts case_studies -name '*.py' -print0 | xargs -0 -n 1 $(PYTHON) -m black

format-check:
	find src tests examples scripts case_studies -name '*.py' -print0 | xargs -0 -n 1 $(PYTHON) -m black --check --quiet

lint:
	$(PYTHON) -m ruff check src tests examples scripts case_studies

typecheck:
	$(PYTHON) -m mypy

test:
	$(PYTHON) -m pytest -q

check: format-check lint typecheck docs-models-check test

examples:
	$(PYTHON) examples/heisenberg_density_of_states.py
	$(PYTHON) examples/ssh_edge_state.py
	$(PYTHON) examples/hofstadter_butterfly.py
	$(PYTHON) examples/kagome_graph.py
	$(PYTHON) examples/haldane_matrix_structure.py
	$(PYTHON) examples/physical_interaction_graph.py
	$(PYTHON) examples/portable_model_bundle.py
	$(PYTHON) examples/external_matrix_import.py

docs-assets: examples
	mkdir -p docs/images
	cp images/*.png docs/images/

docs-models:
	PYTHONPATH=src $(PYTHON) scripts/build_model_docs.py

docs-models-check:
	PYTHONPATH=src $(PYTHON) scripts/build_model_docs.py --check

notebook-summary:
	$(PYTHON) scripts/summarize_notebooks.py
