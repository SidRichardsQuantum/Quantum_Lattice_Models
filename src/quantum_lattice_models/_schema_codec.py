"""Internal JSON codec and validation helpers for portable schemas."""

from __future__ import annotations

import json
from typing import Any, cast

import numpy as np

from quantum_lattice_models.lattice import Bond


def encode_value(value: object) -> object:
    if isinstance(value, complex):
        return {"__complex__": [value.real, value.imag]}
    if isinstance(value, np.generic):
        return encode_value(value.item())
    if isinstance(value, dict):
        return {str(key): encode_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [encode_value(item) for item in value]
    return value


def decode_value(value: object) -> Any:
    if isinstance(value, dict):
        if set(value) == {"__complex__"}:
            parts = value["__complex__"]
            if not isinstance(parts, list) or len(parts) != 2:
                raise ValueError("Complex values must contain [real, imaginary].")
            return complex(float(parts[0]), float(parts[1]))
        return {str(key): decode_value(item) for key, item in value.items()}
    if isinstance(value, list):
        return tuple(decode_value(item) for item in value)
    return value


def mapping(value: object) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise ValueError("Expected a JSON object.")
    return dict(value)


def record_list(value: object, field_name: str) -> list[dict[str, Any]]:
    if not isinstance(value, list) or not all(isinstance(item, dict) for item in value):
        raise ValueError(f"{field_name} must be a list of objects.")
    return value


def validate_fields(data: dict[str, object], allowed: set[str], field_name: str) -> None:
    unknown = sorted(set(data) - allowed)
    if unknown:
        raise ValueError(f"{field_name} contains unknown fields: {', '.join(unknown)}.")


def require_fields(data: dict[str, object], required: set[str], field_name: str) -> None:
    missing = sorted(required - set(data))
    if missing:
        raise ValueError(f"{field_name} is missing required fields: {', '.join(missing)}.")


def require_type(value: object, expected: type, field_name: str) -> None:
    if expected is int:
        valid = isinstance(value, int) and not isinstance(value, bool)
    else:
        valid = isinstance(value, expected)
    if not valid:
        raise ValueError(f"{field_name} must have type {expected.__name__}.")


def required_int(record: dict[str, Any], name: str, field_name: str) -> int:
    if name not in record:
        raise ValueError(f"{field_name} entries require {name!r}.")
    require_type(record[name], int, f"{field_name}.{name}")
    return cast(int, record[name])


def bond_from_record(record: dict[str, Any]) -> Bond:
    validate_fields(record, {"source", "target", "value"}, "lattice.bonds entry")
    return Bond(
        source=required_int(record, "source", "lattice.bonds"),
        target=required_int(record, "target", "lattice.bonds"),
        value=decode_value(record.get("value")),
    )


def validate_json_value(value: object, field_name: str) -> None:
    try:
        json.dumps(encode_value(value))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must contain portable JSON values.") from exc


def validate_string_mapping(value: object, field_name: str) -> None:
    if not isinstance(value, dict) or not all(
        isinstance(key, str) and isinstance(item, str) for key, item in value.items()
    ):
        raise ValueError(f"{field_name} must map strings to strings.")
