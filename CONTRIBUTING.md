# Contributing

Contributions are welcome when they preserve the package's core distinctions:
geometry is not a Hamiltonian, single-particle and many-body bases are not
interchangeable, and conversions must state convention loss explicitly.

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
make check
```

`make check` verifies formatting, linting, strict typing, generated model
documentation, and the test suite. Before submitting a change, also validate
the distributions:

```bash
python -m build
python -m twine check dist/*
```

## Model and importer contributions

Please include:

- an explicit basis and index ordering;
- geometry, units, boundary conditions, gauge/sign conventions, and provenance;
- the smallest representative source file or construction script;
- analytic limits or a trusted cross-package reference;
- dense/sparse or full/reduced block comparisons where relevant;
- a capability report identifying preserved and lost semantics; and
- documentation of computational scaling.

Importers should report inferred values and assumptions. They must not infer a
complete Hamiltonian from coordinates alone.

## Public API and schema changes

Read `DEPRECATIONS.md` and `SCHEMA.md` before modifying public records.
Intentional public API changes require updating the compatibility snapshot and
changelog. Existing schema `1.0` fixtures must remain readable unless an
explicit migration is introduced.

## Testing expectations

Prefer small deterministic systems. Scientific features should compare against
an analytic result, a full-space block, a gauge-equivalent construction, or an
independent package. Performance tests should use operation or storage budgets
where possible; wall-clock thresholds belong in the scheduled benchmark job.

## Feedback reports

The issue forms request concrete files, package versions, conventions, and
expected outputs. Reproducible samples are substantially easier to support than
format names without representative data.
