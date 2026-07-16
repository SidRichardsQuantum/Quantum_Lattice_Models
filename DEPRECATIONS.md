# Deprecation and Compatibility Policy

The package is still alpha software, but portable files and documented public
entry points are treated as compatibility-sensitive.

## Public Python API

- Symbols listed in `quantum_lattice_models.__all__` are public.
- Removing or changing the meaning of a public symbol requires a deprecation
  period of at least two minor releases where practical.
- Deprecated symbols remain importable during that period and emit
  `DeprecationWarning` with the replacement and planned removal release.
- Adding optional parameters is compatible when existing calls retain their
  behavior. Renaming parameters, changing defaults, or narrowing accepted
  values requires deprecation.
- Modules and names beginning with an underscore are internal and may change
  without notice.

The `tests/fixtures/public_api_0_1.json` snapshot detects accidental additions,
removals, duplicates, or renames. Intentional public changes must update the
snapshot and changelog together.

## Portable schemas

Schema compatibility follows `SCHEMA.md`. Existing schema `1.0` files remain
readable. New optional fields may be omitted when empty so canonical output for
older records remains stable.

Removing, reinterpreting, or making a field mandatory requires:

1. a new schema version;
2. an explicit migration;
3. golden fixtures for both the old and new representation; and
4. release notes explaining semantic changes or loss.

## Command-line interface

Documented commands, option names, exit behavior, and JSON field meanings are
compatibility-sensitive. Human-readable terminal formatting may improve, but
machine-readable JSON changes require a golden-fixture update and changelog
entry.

## Deprecation process

Each deprecation should identify:

- the deprecated API or behavior;
- the supported replacement;
- the first release that warns;
- the earliest removal release; and
- any model, basis, schema, or numerical convention affected.

Emergency removals are reserved for security, data corruption, or scientifically
incorrect behavior and must be called out prominently in release notes.
