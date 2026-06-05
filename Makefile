PYTHON ?= .venv/bin/python

.PHONY: format lint test check examples

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
