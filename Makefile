PYTHON ?= .venv/bin/python

.PHONY: format lint test check examples docs-assets notebook-summary

format:
	find src tests examples -name '*.py' -print0 | xargs -0 -n 1 $(PYTHON) -m black

lint:
	$(PYTHON) -m ruff check src tests examples

test:
	$(PYTHON) -m pytest -q

check: lint test

examples:
	$(PYTHON) examples/ising_spectrum.py
	$(PYTHON) examples/heisenberg_density.py
	$(PYTHON) examples/ssh_edge_state.py
	$(PYTHON) examples/tight_binding_spectrum.py
	$(PYTHON) examples/rice_mele_spectrum.py
	$(PYTHON) examples/hofstadter_butterfly.py
	$(PYTHON) examples/bose_hubbard_spectrum.py
	$(PYTHON) examples/haldane_spectrum.py
	$(PYTHON) examples/kagome_graph.py
	$(PYTHON) examples/hamiltonian_matrix.py

docs-assets: examples
	mkdir -p docs/images
	cp images/*.png docs/images/

notebook-summary:
	$(PYTHON) scripts/summarize_notebooks.py
