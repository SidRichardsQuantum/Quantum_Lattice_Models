"""Deterministic CLI rendering helpers."""

from __future__ import annotations

import json

import numpy as np


def print_json(value: object) -> None:
    print(json.dumps(value, indent=2, sort_keys=True, default=json_default))


def json_default(value: object) -> object:
    if isinstance(value, complex):
        return {"real": value.real, "imag": value.imag}
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, tuple):
        return list(value)
    return str(value)


def print_key_values(value: dict[str, object]) -> None:
    for key in sorted(value):
        item = value[key]
        if isinstance(item, (dict, list, tuple)):
            item = json.dumps(item, sort_keys=True, default=json_default)
        print(f"{key}\t{item}")
