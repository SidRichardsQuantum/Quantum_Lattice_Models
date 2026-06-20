"""Portable physical-system records shared by model specifications and plots."""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

LOCAL_DEGREE_KINDS = ("spin", "fermion", "boson", "orbital", "nambu")
BASIS_MAPPING_ROLES = ("tensor_factor", "single_particle_state", "mode")
_OPERATORS_BY_KIND = {
    "spin": {"I", "X", "Y", "Z"},
    "fermion": {"identity", "create", "annihilate", "number"},
    "boson": {"identity", "create", "annihilate", "number", "number_pair"},
    "orbital": {"identity", "create", "annihilate", "number"},
    "nambu": {"identity", "create", "annihilate", "number", "particle", "hole"},
}


@dataclass(frozen=True)
class LocalDegreeOfFreedom:
    """One indexed local physical degree of freedom."""

    index: int
    site: int
    kind: str
    local_dimension: int
    label: str = ""
    component: str | None = None
    orbital: str | None = None
    metadata: dict[str, object] = field(default_factory=dict)

    def validate(self) -> None:
        if not isinstance(self.index, int) or isinstance(self.index, bool) or self.index < 0:
            raise ValueError("local degree index must be a nonnegative integer.")
        if not isinstance(self.site, int) or isinstance(self.site, bool) or self.site < 0:
            raise ValueError("local degree site must be a nonnegative integer.")
        if self.kind not in LOCAL_DEGREE_KINDS:
            raise ValueError(f"local degree kind must be one of {LOCAL_DEGREE_KINDS!r}.")
        if (
            not isinstance(self.local_dimension, int)
            or isinstance(self.local_dimension, bool)
            or self.local_dimension < 1
        ):
            raise ValueError("local degree local_dimension must be a positive integer.")
        if not isinstance(self.label, str):
            raise ValueError("local degree label must be a string.")
        for name, value in (("component", self.component), ("orbital", self.orbital)):
            if value is not None and not isinstance(value, str):
                raise ValueError(f"local degree {name} must be a string or null.")
        _validate_portable_metadata(self.metadata, "local degree metadata")

    def to_dict(self) -> dict[str, object]:
        self.validate()
        return {
            "index": self.index,
            "site": self.site,
            "kind": self.kind,
            "local_dimension": self.local_dimension,
            "label": self.label,
            "component": self.component,
            "orbital": self.orbital,
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> LocalDegreeOfFreedom:
        allowed = {
            "index",
            "site",
            "kind",
            "local_dimension",
            "label",
            "component",
            "orbital",
            "metadata",
        }
        _validate_record(
            data,
            allowed,
            {"index", "site", "kind", "local_dimension"},
            "local degree",
        )
        record = cls(
            index=_integer(data["index"], "local degree index"),
            site=_integer(data["site"], "local degree site"),
            kind=str(data["kind"]),
            local_dimension=_integer(data["local_dimension"], "local degree local_dimension"),
            label=str(data.get("label", "")),
            component=None if data.get("component") is None else str(data["component"]),
            orbital=None if data.get("orbital") is None else str(data["orbital"]),
            metadata=_mapping(data.get("metadata", {}), "local degree metadata"),
        )
        record.validate()
        return record


@dataclass(frozen=True)
class BasisIndexMapping:
    """Map one local degree to a basis factor, mode, or single-particle index."""

    local_degree: int
    basis_index: int
    role: str
    label: str = ""

    def validate(self) -> None:
        for name, value in (
            ("local_degree", self.local_degree),
            ("basis_index", self.basis_index),
        ):
            if not isinstance(value, int) or isinstance(value, bool) or value < 0:
                raise ValueError(f"basis mapping {name} must be a nonnegative integer.")
        if self.role not in BASIS_MAPPING_ROLES:
            raise ValueError(f"basis mapping role must be one of {BASIS_MAPPING_ROLES!r}.")
        if not isinstance(self.label, str):
            raise ValueError("basis mapping label must be a string.")

    def to_dict(self) -> dict[str, object]:
        self.validate()
        return {
            "local_degree": self.local_degree,
            "basis_index": self.basis_index,
            "role": self.role,
            "label": self.label,
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> BasisIndexMapping:
        _validate_record(
            data,
            {"local_degree", "basis_index", "role", "label"},
            {"local_degree", "basis_index", "role"},
            "basis mapping",
        )
        record = cls(
            local_degree=_integer(data["local_degree"], "basis mapping local_degree"),
            basis_index=_integer(data["basis_index"], "basis mapping basis_index"),
            role=str(data["role"]),
            label=str(data.get("label", "")),
        )
        record.validate()
        return record


@dataclass(frozen=True)
class InteractionTerm:
    """Portable onsite or two-body operator term."""

    degrees: tuple[int, ...]
    operators: tuple[str, ...]
    coefficient: complex
    kind: str
    label: str = ""
    unit: str | None = None
    metadata: dict[str, object] = field(default_factory=dict)

    def validate(self, local_degrees: tuple[LocalDegreeOfFreedom, ...]) -> None:
        if len(self.degrees) not in (1, 2):
            raise ValueError("interaction terms must be onsite or two-body.")
        if len(self.operators) != len(self.degrees):
            raise ValueError("interaction operators must align with interaction degrees.")
        if len(set(self.degrees)) != len(self.degrees):
            raise ValueError("two-body interaction degrees must be distinct.")
        by_index = {degree.index: degree for degree in local_degrees}
        for index, operator in zip(self.degrees, self.operators, strict=True):
            if index not in by_index:
                raise ValueError("interaction degree index is not defined by local_degrees.")
            allowed = _OPERATORS_BY_KIND[by_index[index].kind]
            if operator not in allowed:
                raise ValueError(
                    f"operator {operator!r} is incompatible with "
                    f"local degree kind {by_index[index].kind!r}."
                )
        if not isinstance(self.coefficient, (int, float, complex, np.number)):
            raise ValueError("interaction coefficient must be numeric.")
        if not np.isfinite(complex(self.coefficient)):
            raise ValueError("interaction coefficient must be finite.")
        if not isinstance(self.kind, str) or not self.kind:
            raise ValueError("interaction kind must be a nonempty string.")
        if not isinstance(self.label, str):
            raise ValueError("interaction label must be a string.")
        if self.unit is not None and not isinstance(self.unit, str):
            raise ValueError("interaction unit must be a string or null.")
        _validate_portable_metadata(self.metadata, "interaction metadata")

    def to_dict(self) -> dict[str, object]:
        return {
            "degrees": list(self.degrees),
            "operators": list(self.operators),
            "coefficient": {
                "__complex__": [complex(self.coefficient).real, complex(self.coefficient).imag]
            },
            "kind": self.kind,
            "label": self.label,
            "unit": self.unit,
            "metadata": dict(self.metadata),
        }

    @classmethod
    def from_dict(cls, data: dict[str, object]) -> InteractionTerm:
        _validate_record(
            data,
            {"degrees", "operators", "coefficient", "kind", "label", "unit", "metadata"},
            {"degrees", "operators", "coefficient", "kind"},
            "interaction",
        )
        degrees = data["degrees"]
        operators = data["operators"]
        if not isinstance(degrees, list) or not isinstance(operators, list):
            raise ValueError("interaction degrees and operators must be lists.")
        record = cls(
            degrees=tuple(_integer(value, "interaction degree") for value in degrees),
            operators=tuple(str(value) for value in operators),
            coefficient=_complex_value(data["coefficient"]),
            kind=str(data["kind"]),
            label=str(data.get("label", "")),
            unit=None if data.get("unit") is None else str(data["unit"]),
            metadata=_mapping(data.get("metadata", {}), "interaction metadata"),
        )
        return record


def validate_physical_system(
    local_degrees: tuple[LocalDegreeOfFreedom, ...],
    basis_mappings: tuple[BasisIndexMapping, ...],
    interactions: tuple[InteractionTerm, ...],
    *,
    n_sites: int | None,
) -> None:
    """Validate cross-record indices, mappings, and operator compatibility."""

    for degree in local_degrees:
        degree.validate()
        if n_sites is not None and degree.site >= n_sites:
            raise ValueError("local degree site must be less than lattice.n_sites.")
    indices = [degree.index for degree in local_degrees]
    if indices != list(range(len(local_degrees))):
        raise ValueError("local degree indices must be contiguous and start at zero.")
    for mapping in basis_mappings:
        mapping.validate()
        if mapping.local_degree >= len(local_degrees):
            raise ValueError("basis mapping references an undefined local degree.")
    if local_degrees and {mapping.local_degree for mapping in basis_mappings} != set(indices):
        raise ValueError("every local degree must have exactly one basis mapping.")
    if len({mapping.local_degree for mapping in basis_mappings}) != len(basis_mappings):
        raise ValueError("each local degree may have at most one basis mapping.")
    if len({(mapping.role, mapping.basis_index) for mapping in basis_mappings}) != len(
        basis_mappings
    ):
        raise ValueError("basis indices must be unique within each mapping role.")
    for interaction in interactions:
        interaction.validate(local_degrees)


def _validate_record(
    data: dict[str, object],
    allowed: set[str],
    required: set[str],
    name: str,
) -> None:
    if not isinstance(data, dict):
        raise ValueError(f"{name} must be an object.")
    unknown = sorted(set(data) - allowed)
    if unknown:
        raise ValueError(f"{name} contains unknown fields: {', '.join(unknown)}.")
    missing = sorted(required - set(data))
    if missing:
        raise ValueError(f"{name} is missing required fields: {', '.join(missing)}.")


def _integer(value: object, name: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool):
        raise ValueError(f"{name} must be an integer.")
    return value


def _mapping(value: object, name: str) -> dict[str, object]:
    if not isinstance(value, dict) or not all(isinstance(key, str) for key in value):
        raise ValueError(f"{name} must be an object with string keys.")
    return dict(value)


def _complex_value(value: object) -> complex:
    if isinstance(value, (int, float, complex)) and not isinstance(value, bool):
        return complex(value)
    if isinstance(value, dict) and set(value) == {"__complex__"}:
        parts = value["__complex__"]
        if isinstance(parts, list) and len(parts) == 2:
            return complex(float(parts[0]), float(parts[1]))
    raise ValueError("interaction coefficient must be numeric or complex encoded.")


def _validate_portable_metadata(value: object, name: str) -> None:
    import json

    try:
        json.dumps(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must contain portable JSON values.") from exc
