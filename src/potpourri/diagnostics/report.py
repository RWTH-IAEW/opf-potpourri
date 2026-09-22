# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""The types a diagnostic run produces.

A finding is a `DiagnosticIssue`: a typed record, not a sentence. It carries
a stable machine-readable `code`, the pandapower object it concerns, the
Pyomo component behind it, and the numbers — so a caller can filter on
`code == "LINE_THERMAL_OVERLOAD"` without parsing English, and a human
still gets a readable line when the report is printed.

Severities are `INFO`, `WARNING` and `ERROR`. There is deliberately no
`CRITICAL`: nothing here is more urgent than "this is why your OPF did not
solve", and an unused level invites inconsistent use.

The wording rule for every message in this package: **state what was
measured, not what caused what.** "These constraints form the identified
conflicting set" is supportable; "line 8 caused the solver to fail" is not,
unless the analysis actually established it.
"""

from __future__ import annotations

import enum
from dataclasses import dataclass, field
from typing import Any

from potpourri.diagnostics.mappings import ElementRef


class DiagnosticSeverity(enum.IntEnum):
    """How much attention a finding deserves.

    Ordered, so `max(...)` and `>=` filtering work: `INFO < WARNING < ERROR`.
    """

    INFO = 10
    WARNING = 20
    ERROR = 30

    def __str__(self) -> str:
        """The name, for printing."""
        return self.name


class DiagnosticCategory(enum.Enum):
    """What area a finding belongs to.

    Used for grouping in the printed report and for filtering
    programmatically. The value is the string that appears in serialised
    output.
    """

    NETWORK = "NETWORK"
    DATA = "DATA"
    STRUCTURE = "STRUCTURE"
    BOUNDS = "BOUNDS"
    INITIALIZATION = "INITIALIZATION"
    INFEASIBILITY = "INFEASIBILITY"
    NUMERICAL = "NUMERICAL"
    SOLVER = "SOLVER"
    POWER_BALANCE = "POWER_BALANCE"
    VOLTAGE = "VOLTAGE"
    THERMAL = "THERMAL"
    GENERATOR_CAPABILITY = "GENERATOR_CAPABILITY"
    STORAGE = "STORAGE"
    RESULT_MAPPING = "RESULT_MAPPING"
    CROSS_CHECK = "CROSS_CHECK"

    def __str__(self) -> str:
        """The category name, for printing."""
        return self.value


@dataclass
class DiagnosticIssue:
    """One finding, in both engineering and Pyomo terms.

    Attributes:
        severity: `INFO`, `WARNING` or `ERROR`.
        category: Which `DiagnosticCategory` this belongs to.
        code: Stable identifier, e.g. `"BUS_VOLTAGE_HIGH"`. Filter on this
            rather than on `message`, which is prose and may be reworded.
        message: One sentence for a human, in engineering terms.
        recommendation: What to look at next, where there is something
            useful to say. Never a claim of causality.
        element: The pandapower object concerned, in the caller's numbering.
        pyomo_component: Name of the Pyomo component, e.g.
            `"line_lim_from"`.
        pyomo_index: Index into that component.
        value: The measured quantity.
        lower_bound: Applicable lower limit, if any.
        upper_bound: Applicable upper limit, if any.
        violation: How far outside the limits `value` is, in `unit`.
        relative_violation: `violation` scaled by the limit, where that is
            meaningful. A fraction, not a percentage.
        unit: Unit of `value`, `violation` and the bounds.
        time_step: Time index for multi-period models; `None` otherwise.
        context: Anything else worth keeping, serialised as-is. Must hold
            only JSON-friendly values.
    """

    severity: DiagnosticSeverity
    category: DiagnosticCategory
    code: str
    message: str
    recommendation: str | None = None
    element: ElementRef | None = None
    pyomo_component: str | None = None
    pyomo_index: Any = None
    value: float | None = None
    lower_bound: float | None = None
    upper_bound: float | None = None
    violation: float | None = None
    relative_violation: float | None = None
    unit: str | None = None
    time_step: int | None = None
    context: dict[str, Any] = field(default_factory=dict)

    # --- rendering ----------------------------------------------------

    def headline(self) -> str:
        """The first line: severity, category and the object concerned."""
        where = f" {self.element}" if self.element is not None else ""
        return f"[{self.severity}] {self.category}{where}"

    def pyomo_label(self) -> str | None:
        """`line_lim_from[23]`, or `None` if no component is attached."""
        if self.pyomo_component is None:
            return None
        if self.pyomo_index is None:
            return self.pyomo_component
        index = self.pyomo_index
        if isinstance(index, tuple):
            inner = ",".join(str(i) for i in index)
        else:
            inner = str(index)
        return f"{self.pyomo_component}[{inner}]"

    def __str__(self) -> str:
        """A short multi-line rendering, numbers included where present."""
        lines = [self.headline(), f"  {self.message}"]
        if self.time_step is not None:
            lines.append(f"  time step: {self.time_step}")
        if self.value is not None:
            unit = f" {self.unit}" if self.unit else ""
            lines.append(f"  value: {self.value:.6g}{unit}")
        bounds = []
        if self.lower_bound is not None:
            bounds.append(f"min {self.lower_bound:.6g}")
        if self.upper_bound is not None:
            bounds.append(f"max {self.upper_bound:.6g}")
        if bounds:
            lines.append("  limits: " + ", ".join(bounds))
        if self.violation is not None:
            rel = (
                f" ({self.relative_violation * 100:.2f} %)"
                if self.relative_violation is not None
                else ""
            )
            unit = f" {self.unit}" if self.unit else ""
            lines.append(f"  violation: {self.violation:.6g}{unit}{rel}")
        label = self.pyomo_label()
        if label:
            lines.append(f"  pyomo: {label}")
        if self.recommendation:
            lines.append(f"  -> {self.recommendation}")
        return "\n".join(lines)

    def to_dict(self) -> dict[str, Any]:
        """A JSON-serialisable view.

        No Pyomo objects survive this: the component is reduced to its
        name and index, both of which are plain data.

        Returns:
            A flat dict with the element split into `element_type`,
            `element_index` and `element_name`.
        """
        index = self.pyomo_index
        if isinstance(index, tuple):
            index = list(index)
        return {
            "severity": str(self.severity),
            "category": str(self.category),
            "code": self.code,
            "message": self.message,
            "recommendation": self.recommendation,
            "element_type": self.element.table if self.element else None,
            "element_index": self.element.index if self.element else None,
            "element_name": self.element.name if self.element else None,
            "time": self.time_step,
            "value": self.value,
            "lower_bound": self.lower_bound,
            "upper_bound": self.upper_bound,
            "violation": self.violation,
            "relative_violation": self.relative_violation,
            "unit": self.unit,
            "pyomo_component": self.pyomo_component,
            "pyomo_index": index,
            "context": dict(self.context),
        }


@dataclass
class DiagnosticReport:
    """Everything one `diagnose()` call found.

    Attributes:
        issues: Every finding, in the order the checks produced them.
        summary: Headline facts about the model, solver and network, for
            the top of the printed report.
        skipped: Checks that could not run, mapped to why. A diagnostic
            that cannot run says so rather than staying silent.
    """

    issues: list[DiagnosticIssue] = field(default_factory=list)
    summary: dict[str, Any] = field(default_factory=dict)
    skipped: dict[str, str] = field(default_factory=dict)

    # --- collection ---------------------------------------------------

    def add(self, issue: DiagnosticIssue) -> DiagnosticIssue:
        """Record one finding and hand it back."""
        self.issues.append(issue)
        return issue

    def skip(self, check: str, reason: str) -> None:
        """Record that a check could not run, and why."""
        self.skipped[check] = reason

    def extend(self, other: DiagnosticReport) -> None:
        """Fold another report's findings into this one."""
        self.issues.extend(other.issues)
        self.summary.update(other.summary)
        self.skipped.update(other.skipped)

    # --- filtered views -----------------------------------------------

    @property
    def errors(self) -> list[DiagnosticIssue]:
        """Findings at `ERROR`."""
        return self.by_severity(DiagnosticSeverity.ERROR)

    @property
    def warnings(self) -> list[DiagnosticIssue]:
        """Findings at `WARNING`."""
        return self.by_severity(DiagnosticSeverity.WARNING)

    @property
    def info(self) -> list[DiagnosticIssue]:
        """Findings at `INFO`."""
        return self.by_severity(DiagnosticSeverity.INFO)

    def by_severity(
        self, severity: DiagnosticSeverity
    ) -> list[DiagnosticIssue]:
        """Findings at exactly this severity."""
        return [i for i in self.issues if i.severity == severity]

    def by_category(
        self, category: DiagnosticCategory
    ) -> list[DiagnosticIssue]:
        """Findings in this category."""
        return [i for i in self.issues if i.category == category]

    def by_code(self, code: str) -> list[DiagnosticIssue]:
        """Findings carrying this code."""
        return [i for i in self.issues if i.code == code]

    @property
    def ok(self) -> bool:
        """Whether nothing was reported at `ERROR`."""
        return not self.errors

    # --- export -------------------------------------------------------

    def to_dict(self) -> dict[str, Any]:
        """A JSON-serialisable view of the whole report."""
        return {
            "summary": dict(self.summary),
            "skipped": dict(self.skipped),
            "issues": [i.to_dict() for i in self.issues],
        }

    def to_dataframe(self):
        """The findings as a pandas DataFrame, one row per issue.

        Returns:
            A DataFrame with one column per `DiagnosticIssue.to_dict` key
            except `context`, which is dropped because a nested dict does
            not belong in a flat table. Empty with the right columns when
            there are no findings, so downstream code can rely on them.
        """
        import pandas as pd

        columns = [
            "severity",
            "category",
            "code",
            "element_type",
            "element_index",
            "element_name",
            "time",
            "value",
            "lower_bound",
            "upper_bound",
            "violation",
            "relative_violation",
            "unit",
            "pyomo_component",
            "pyomo_index",
            "message",
            "recommendation",
        ]
        rows = []
        for issue in self.issues:
            row = issue.to_dict()
            row.pop("context", None)
            rows.append(row)
        return pd.DataFrame(rows, columns=columns)

    # --- printing -----------------------------------------------------

    def __str__(self) -> str:
        """The terminal report: summary, then findings grouped by category."""
        width = 62
        out = ["potpourri OPF diagnostics", "=" * width, ""]

        if self.summary:
            for section, entries in _grouped_summary(self.summary).items():
                out.append(section)
                for key, value in entries:
                    dots = "." * max(3, 20 - len(key))
                    out.append(f"  {key} {dots} {value}")
                out.append("")

        counts = {
            sev: len(self.by_severity(sev)) for sev in DiagnosticSeverity
        }
        out.append(
            "Findings: "
            + ", ".join(f"{n} {sev.name.lower()}" for sev, n in counts.items())
        )
        out.append("")

        for category in DiagnosticCategory:
            found = self.by_category(category)
            if not found:
                continue
            out.append(f"{category}")
            for issue in sorted(found, key=lambda i: -int(i.severity)):
                for line in str(issue).splitlines():
                    out.append("  " + line)
                out.append("")

        if self.skipped:
            out.append("Not run")
            for check, reason in self.skipped.items():
                out.append(f"  {check}: {reason}")
            out.append("")

        if not self.issues:
            out.append("No issues found.")
            out.append("")
        return "\n".join(out).rstrip() + "\n"


def _grouped_summary(summary: dict[str, Any]) -> dict[str, list]:
    """Split flat `section.key` summary entries into printable sections.

    Args:
        summary: Summary mapping, whose keys may be `"model.buses"` style.

    Returns:
        An ordered mapping of section title to `(key, value)` pairs, with
        unprefixed entries collected under `"General"`.
    """
    sections: dict[str, list] = {}
    for key, value in summary.items():
        section, _, leaf = key.partition(".")
        if not leaf:
            section, leaf = "General", key
        sections.setdefault(section.capitalize(), []).append((leaf, value))
    return sections
