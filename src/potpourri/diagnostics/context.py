# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""What a diagnostic check is handed, and what it is allowed to assume.

A check must work on a model that has not had `add_OPF()` called, on one
that failed to solve, and on one that solved fine. `DiagnosticContext`
carries the state each check needs and answers the capability questions
that decide which checks apply at all — asking for reactive-power
diagnostics on a DC model is a category error, not a finding.

Capabilities are detected from the components the model actually declares
rather than from its class name. A subclass, a mix-in or a future
formulation then behaves correctly without this module knowing about it.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from potpourri.diagnostics.mappings import IndexMap


@dataclass
class DiagnosticContext:
    """The model under diagnosis, with its mappings and capabilities.

    Attributes:
        model_obj: The potpourri model object, e.g. an `ACOPF`.
        net: Its preprocessed pandapower network (`model_obj.net`).
        model: The Pyomo `ConcreteModel`, or `None` if not built.
        imap: Index map back to the caller's pandapower objects.
        tol: Absolute tolerance for calling a constraint violated, in the
            model's per-unit system.
        near_bound_fraction: How close to a bound counts as "nearly
            binding", as a fraction of the bound's range.
        options: Free-form switches passed through from `diagnose()`.
    """

    model_obj: Any
    net: Any
    model: Any
    imap: IndexMap
    tol: float = 1e-6
    near_bound_fraction: float = 0.01
    options: dict[str, Any] = field(default_factory=dict)

    @classmethod
    def from_model(cls, model_obj, **kwargs) -> DiagnosticContext:
        """Build a context for one potpourri model.

        Args:
            model_obj: A constructed potpourri model object.
            **kwargs: Overrides for any other field, e.g. `tol`.

        Returns:
            A ready `DiagnosticContext`. Missing pieces stay `None`; no
            check may assume the Pyomo model exists.
        """
        return cls(
            model_obj=model_obj,
            net=getattr(model_obj, "net", None),
            model=getattr(model_obj, "model", None),
            imap=IndexMap.from_model(model_obj),
            **kwargs,
        )

    # --- capability detection -----------------------------------------

    def has(self, component: str) -> bool:
        """Whether the Pyomo model declares a component of this name."""
        return self.model is not None and hasattr(self.model, component)

    @property
    def base_mva(self) -> float:
        """The per-unit base, in MVA. Defaults to 1.0 if unavailable."""
        value = getattr(self.model_obj, "baseMVA", None)
        try:
            return float(value)
        except (TypeError, ValueError):
            return 1.0

    @property
    def is_multi_period(self) -> bool:
        """Whether the model carries a time index."""
        return self.has("T")

    @property
    def time_steps(self) -> list[int]:
        """The model's time steps, or `[None]` for single-period models."""
        if not self.is_multi_period:
            return [None]
        return [int(t) for t in self.model.T]

    @property
    def has_reactive(self) -> bool:
        """Whether the formulation models reactive power.

        Detected from the presence of a reactive nodal balance rather than
        from the class name, so DC and any future active-only formulation
        are both excluded without being enumerated here.
        """
        return self.has("KCL_reactive")

    @property
    def has_voltage_magnitude(self) -> bool:
        """Whether bus voltage magnitude is a variable of the model."""
        return self.has("v")

    @property
    def is_opf(self) -> bool:
        """Whether `add_OPF()` has run and attached operating limits."""
        return self.has("Vmax") or self.has("SLmax") or self.has("v_pyo")

    @property
    def has_solution(self) -> bool:
        """Whether a solver result has been loaded onto the model object."""
        return getattr(self.model_obj, "results", None) is not None

    @property
    def formulation(self) -> str:
        """A short label for the formulation, for the report header."""
        name = type(self.model_obj).__name__
        kind = "AC" if self.has_reactive else "DC"
        if "LPAC" in name:
            kind = "LPAC"
        suffix = " multi-period" if self.is_multi_period else ""
        return f"{kind}{suffix} ({name})"

    def value(self, obj, default=None):
        """Read a Pyomo value without raising on an uninitialised one.

        Args:
            obj: A Pyomo variable data, parameter or expression.
            default: Returned when the value cannot be obtained.

        Returns:
            The float value, or `default` when it is `None`, not yet
            initialised, or raises on evaluation. Diagnostics run on
            half-solved models, where all three happen.
        """
        import pyomo.environ as pyo

        try:
            raw = pyo.value(obj, exception=False)
        except (ValueError, TypeError, ZeroDivisionError):
            return default
        if raw is None:
            return default
        try:
            out = float(raw)
        except (TypeError, ValueError):
            return default
        return default if out != out else out
