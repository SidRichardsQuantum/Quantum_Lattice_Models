"""Model builders organized by physical domain.

The historical flat imports from :mod:`quantum_lattice_models.models` remain
available. New code may use the focused ``spin``, ``particles``,
``topological``, ``periodic``, and ``benchmarks`` submodules.
"""

from __future__ import annotations

from . import benchmarks as benchmarks
from . import particles as particles
from . import periodic as periodic
from . import spin as spin
from . import topological as topological
from .particles import *  # noqa: F403
from .spin import *  # noqa: F403
from .topological import *  # noqa: F403

__all__ = [*particles.__all__, *spin.__all__, *topological.__all__]
