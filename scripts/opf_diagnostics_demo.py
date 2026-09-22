# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Diagnosing an OPF that will not solve, and one that solved oddly.

``opf.diagnose()`` answers the question a solver cannot: *why*. A solver
says ``infeasible``; it has no idea that ``line_lim_from[23]`` is the
from-side current limit of the cable feeding your worst-supplied street.
The diagnostics report in terms of the pandapower objects you built, and
keep the Pyomo identifier alongside so you can get back to the equation.

This script walks five deliberately broken four-bus networks plus one
healthy one, and prints what the diagnostics make of each:

  1. an island with load and no source
  2. ``min_p_mw`` above ``max_p_mw`` on the PV units
  3. a slack voltage setpoint outside the OPF's own voltage band
  4. demand beyond anything the sources can supply
  5. a line whose current rating is zero
  6. a healthy network, solved — which must produce **no** errors and
     **no** warnings, or nobody will read the reports that matter

Cases 1-5 are all caught *before* a solver runs, by ``level="basic"``.
That is the point: a contradiction in the data does not need an
optimiser to find it.

The last section solves the healthy case and shows the other half of the
feature — which limits are shaping the answer, and whether an
independent pandapower power flow reproduces it.

Network: ``simple_four_bus_system`` (small and deterministic)
Solver : IPOPT (only for the final, healthy case)

Author: Steffen Kortmann (2026)
"""

from __future__ import annotations

import copy
import logging
import warnings

import pandapower as pp

from potpourri.models.ACOPF_base import ACOPF

warnings.filterwarnings("ignore")
# Same as battery_multi_period_opf.py: Pyomo logs a deprecation notice
# about an implicit Param domain while the model is built, which is noise
# at the top of a worked example.
logging.getLogger("pyomo.core").setLevel(logging.ERROR)

# ── Configuration ────────────────────────────────────────────────────────────

#: Voltage band the OPF enforces on every bus.
VM_MIN = 0.95
VM_MAX = 1.05

#: External-grid limits, wide enough not to bind in the healthy case.
EXT_GRID_LIMIT_MW = 100.0

#: Diagnostic level for the broken cases. "basic" runs no solver and no
#: power flow, which is the whole point: these faults are found in the data.
BROKEN_LEVEL = "basic"

#: Level for the solved case, which adds plausibility checks and the
#: independent pandapower cross-check.
SOLVED_LEVEL = "standard"

SOLVER = "ipopt"


# ── Network preparation ──────────────────────────────────────────────────────


def healthy_net():
    """A four-bus network carrying the limits an OPF needs.

    Returns:
        A pandapower network with a voltage band on every bus and finite
        active and reactive limits on the external grid.
    """
    net = pp.networks.simple_four_bus_system()
    net.bus["min_vm_pu"] = VM_MIN
    net.bus["max_vm_pu"] = VM_MAX
    net.ext_grid["min_p_mw"] = -EXT_GRID_LIMIT_MW
    net.ext_grid["max_p_mw"] = EXT_GRID_LIMIT_MW
    net.ext_grid["min_q_mvar"] = -EXT_GRID_LIMIT_MW
    net.ext_grid["max_q_mvar"] = EXT_GRID_LIMIT_MW
    return net


def break_island(net):
    """Strand a load on a bus nothing supplies."""
    bus = pp.create_bus(net, vn_kv=0.4, name="Stranded bus")
    pp.create_load(net, bus, p_mw=0.02, name="Stranded load")
    return net


def break_bounds(net):
    """Set a lower active-power limit above the upper one."""
    net.sgen["min_p_mw"] = 5.0
    net.sgen["max_p_mw"] = 1.0
    return net


def break_voltage_setpoint(net):
    """Hold the slack outside the band the OPF enforces at its bus."""
    net.ext_grid["vm_pu"] = 1.20
    return net


def break_adequacy(net):
    """Add demand that no source could ever cover."""
    net.ext_grid["max_p_mw"] = 0.001
    pp.create_load(net, 3, p_mw=5.0, name="Oversized load")
    return net


def break_thermal_rating(net):
    """Give a line a zero current rating, forbidding any flow."""
    net.line.loc[net.line.index[0], "max_i_ka"] = 0.0
    return net


#: Each entry is a label and the function that introduces one known fault.
BROKEN_CASES = [
    ("island with load and no source", break_island),
    ("min_p_mw above max_p_mw", break_bounds),
    ("slack setpoint outside the voltage band", break_voltage_setpoint),
    ("demand beyond any possible supply", break_adequacy),
    ("line with a zero current rating", break_thermal_rating),
]


# ── Reporting ────────────────────────────────────────────────────────────────


def show(label, report, severities=("ERROR", "WARNING")):
    """Print the findings of one report at the given severities.

    Args:
        label: Heading for this case.
        report: The `DiagnosticReport` to print.
        severities: Severity names to include.

    Returns:
        None. Everything is printed.
    """
    print("=" * 76)
    print(label.upper())
    print("=" * 76)
    shown = 0
    for issue in report.issues:
        if issue.severity.name not in severities:
            continue
        shown += 1
        print(f"  {issue.headline()}")
        print(f"    {issue.message}")
        if issue.recommendation:
            print(f"    -> {issue.recommendation}")
        label_pyomo = issue.pyomo_label()
        if label_pyomo:
            print(f"    pyomo: {label_pyomo}")
    if not shown:
        print("  nothing at this severity")
    print(
        f"  [{len(report.errors)} errors, {len(report.warnings)} warnings, "
        f"{len(report.info)} info]"
    )
    print()


def diagnose_broken_cases():
    """Run each deliberately broken network through the diagnostics.

    Returns:
        None. Everything is printed.
    """
    for label, break_it in BROKEN_CASES:
        net = break_it(healthy_net())
        opf = ACOPF(copy.deepcopy(net))
        opf.add_OPF()
        opf.add_voltage_deviation_objective()
        # No solve: every fault below is found in the data itself.
        show(label, opf.diagnose(level=BROKEN_LEVEL))


def diagnose_solved_case():
    """Solve the healthy network and inspect the answer.

    Shows the other half of the feature: which limits shape the result,
    whether the power balance closes, and whether an independent
    pandapower power flow reproduces the OPF state.

    Returns:
        None. Everything is printed.
    """
    opf = ACOPF(healthy_net())
    opf.add_OPF()
    opf.add_voltage_deviation_objective()
    opf.solve(solver=SOLVER, print_solver_output=False)

    report = opf.diagnose(level=SOLVED_LEVEL)
    show("healthy network, solved", report)

    print("Limits shaping the solution")
    binding = report.by_code("SOLUTION_CONSTRAINT_BINDING")
    if not binding:
        print("  none are close to binding")
    for issue in binding:
        print(f"  {issue.element}: {issue.message}")
    print()

    print("Independent pandapower cross-check")
    for issue in report.by_code("REPLAY_AGREES") + report.by_code(
        "REPLAY_MISMATCH"
    ):
        print(f"  {issue.message}")
    print()

    print("System totals")
    for key, value in report.summary.items():
        if key.startswith("balance."):
            print(f"  {key.split('.', 1)[1]:<20} {value}")
    print()

    print("Machine-readable view")
    frame = report.to_dataframe()
    print(f"  {len(frame)} findings as a DataFrame")
    print(f"  columns: {', '.join(list(frame.columns)[:6])}, ...")
    print("  filter on `code`, never on `message` — codes are stable")
    print()


def main():
    """Run every case and print the reports.

    Returns:
        None.
    """
    print()
    print("potpourri OPF diagnostics — worked examples")
    print()
    print(
        "Cases 1-5 are caught before any solver runs. A contradiction in\n"
        "the data does not need an optimiser to find it.\n"
    )
    diagnose_broken_cases()
    diagnose_solved_case()
    print(
        "Full guide: docs/user-guide/diagnostics.md\n"
        "Levels: basic (no solver), standard (default), deep (relaxation)\n"
    )


if __name__ == "__main__":
    main()
