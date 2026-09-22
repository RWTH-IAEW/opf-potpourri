# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""A diagnostic framework has to be tested on broken networks.

Each network below has one known fault, and the test asserts that the
report names it — with the right pandapower object, not just the right
category. Two of them exist specifically to catch the failure mode that
would make the whole feature untrustworthy:

* `test_mapping_survives_non_contiguous_indices` uses a network whose bus,
  line and load indices have holes in them, plus a bus-bus switch that
  fuses two buses away. potpourri renumbers buses internally, so a report
  that names `net.bus[2]` when the user's network has buses 10, 20, 30 is
  worse than no report at all.
* `test_healthy_network_is_quiet` checks the other direction: a sound
  model must not produce a wall of warnings, or nobody will read the ones
  that matter.
"""

import copy

import pandapower as pp
import pytest

from potpourri.diagnostics import (
    DiagnosticCategory,
    DiagnosticReport,
    DiagnosticSeverity,
)
from potpourri.diagnostics.mappings import IndexMap
from potpourri.models.ACOPF_base import ACOPF
from potpourri.models.DCOPF import DCOPF


def healthy_net():
    """A four-bus network with the limits an OPF needs."""
    net = pp.networks.simple_four_bus_system()
    net.bus["min_vm_pu"] = 0.95
    net.bus["max_vm_pu"] = 1.05
    net.ext_grid["min_p_mw"] = -100.0
    net.ext_grid["max_p_mw"] = 100.0
    net.ext_grid["min_q_mvar"] = -100.0
    net.ext_grid["max_q_mvar"] = 100.0
    return net


def build(net, objective=True):
    """Construct an ACOPF on a copy of `net`, without solving it."""
    opf = ACOPF(copy.deepcopy(net))
    opf.add_OPF()
    if objective:
        opf.add_voltage_deviation_objective()
    return opf


def codes(report):
    """The set of issue codes a report carries."""
    return {issue.code for issue in report.issues}


# --- the report object ------------------------------------------------


def test_report_is_machine_readable():
    report = build(healthy_net()).diagnose(level="basic")
    assert isinstance(report, DiagnosticReport)

    payload = report.to_dict()
    assert set(payload) == {"summary", "skipped", "issues"}

    import json

    json.dumps(payload)  # must not contain Pyomo objects

    frame = report.to_dataframe()
    for column in ("severity", "category", "code", "pyomo_component"):
        assert column in frame.columns


def test_severity_is_ordered():
    assert DiagnosticSeverity.INFO < DiagnosticSeverity.WARNING
    assert DiagnosticSeverity.WARNING < DiagnosticSeverity.ERROR


def test_printing_a_report_mentions_the_formulation():
    text = str(build(healthy_net()).diagnose(level="basic"))
    assert "potpourri OPF diagnostics" in text
    assert "AC" in text


# --- index mapping, the part everything else depends on ---------------


def gappy_net():
    """Non-contiguous indices, a fused bus pair and an impedance branch."""
    net = pp.create_empty_network()
    buses = {
        i: pp.create_bus(net, vn_kv=20.0, index=i, name=f"Bus {i}")
        for i in (10, 20, 30, 40, 50)
    }
    pp.create_ext_grid(net, buses[10])
    pp.create_line_from_parameters(
        net,
        buses[10],
        buses[20],
        1.0,
        0.1,
        0.1,
        0.0,
        0.4,
        index=3,
        name="Cable A",
    )
    pp.create_line_from_parameters(
        net,
        buses[20],
        buses[30],
        1.0,
        0.1,
        0.1,
        0.0,
        0.4,
        index=7,
        name="Cable B",
    )
    pp.create_switch(
        net, bus=buses[30], element=buses[40], et="b", closed=True, z_ohm=0.0
    )
    pp.create_impedance(
        net,
        buses[30],
        buses[50],
        rft_pu=0.01,
        xft_pu=0.05,
        sn_mva=1.0,
        index=2,
        name="Imp Z",
    )
    pp.create_load(net, buses[50], p_mw=0.05, index=9, name="Load L")
    return net


def test_mapping_survives_non_contiguous_indices():
    """The report must name the caller's objects, not internal ones."""
    net = gappy_net()
    opf = ACOPF(copy.deepcopy(net))
    imap = IndexMap.from_model(opf)

    # potpourri renumbers buses to 0..n-1 internally; the caller's are not
    assert list(opf.net.bus.index) != list(net.bus.index)

    resolved = {
        imap.bus_ref(b).index
        for b in opf.model.B
        if imap.bus_ref(b).is_resolved
    }
    assert resolved <= set(net.bus.index), (
        "a bus was named that the caller's network does not contain"
    )
    assert 10 in resolved and 50 in resolved

    # names travel with the bus through the renumbering
    names = {
        imap.bus_ref(b).index: imap.bus_ref(b).name
        for b in opf.model.B
        if imap.bus_ref(b).is_resolved
    }
    assert names[10] == "Bus 10"

    # model.L spans lines and impedance rows; both must resolve correctly
    line_refs = {int(i): imap.line_ref(i) for i in opf.model.L}
    tables = {(ref.table, ref.index) for ref in line_refs.values()}
    assert ("line", 3) in tables
    assert ("line", 7) in tables
    assert ("impedance", 2) in tables


def test_merged_buses_are_recorded_not_hidden():
    """A bus fused away by a switch must not vanish silently."""
    opf = ACOPF(gappy_net())
    imap = IndexMap.from_model(opf)
    merged = {m for group in imap.merged_buses.values() for m in group}
    assert 40 in merged or 30 in merged


def test_unknown_element_is_reported_as_unmapped():
    imap = IndexMap.from_model(ACOPF(healthy_net()))
    ref = imap.element_ref("load", 9999)
    assert not ref.is_resolved
    assert str(ref) == "<unmapped>"


# --- broken cases -----------------------------------------------------


def test_disconnected_island_is_found():
    net = healthy_net()
    bus = pp.create_bus(net, vn_kv=0.4, name="Island bus")
    pp.create_load(net, bus, p_mw=0.02, name="Island load")

    report = build(net).diagnose(level="basic")
    assert "NET_ISLAND_WITHOUT_SOURCE" in codes(report)
    issue = report.by_code("NET_ISLAND_WITHOUT_SOURCE")[0]
    assert issue.severity is DiagnosticSeverity.ERROR
    assert issue.element is not None and issue.element.table == "bus"


def test_contradictory_bounds_are_found():
    net = healthy_net()
    net.sgen["min_p_mw"] = 5.0
    net.sgen["max_p_mw"] = 1.0

    report = build(net).diagnose(level="basic")
    found = report.by_code("BOUND_MIN_ABOVE_MAX")
    assert found, "min_p_mw > max_p_mw was not reported"
    assert all(i.severity is DiagnosticSeverity.ERROR for i in found)
    assert {i.element.table for i in found} == {"sgen"}
    assert found[0].lower_bound == 5.0 and found[0].upper_bound == 1.0


def test_voltage_setpoint_outside_band_is_found():
    net = healthy_net()
    net.ext_grid["vm_pu"] = 1.20

    report = build(net).diagnose(level="basic")
    issue = report.by_code("VOLTAGE_SETPOINT_OUTSIDE_BAND")[0]
    assert issue.severity is DiagnosticSeverity.ERROR
    assert issue.value == pytest.approx(1.20)
    assert issue.upper_bound == pytest.approx(1.05)
    assert issue.category is DiagnosticCategory.VOLTAGE


def test_impossible_demand_is_found():
    net = healthy_net()
    net.ext_grid["max_p_mw"] = 0.001
    pp.create_load(net, 3, p_mw=5.0, name="Huge load")

    report = build(net).diagnose(level="basic")
    issue = report.by_code("ADEQUACY_INSUFFICIENT_GENERATION")[0]
    assert issue.severity is DiagnosticSeverity.ERROR
    assert issue.unit == "MW"
    assert "Losses are not included" in issue.recommendation


def test_nonpositive_thermal_rating_is_found():
    net = healthy_net()
    net.line.loc[net.line.index[0], "max_i_ka"] = 0.0

    report = build(net).diagnose(level="basic")
    issue = report.by_code("THERMAL_RATING_NONPOSITIVE")[0]
    assert issue.element.table == "line"
    assert issue.element.index == int(net.line.index[0])


def test_missing_objective_is_found():
    report = build(healthy_net(), objective=False).diagnose(level="basic")
    assert "MODEL_NO_OBJECTIVE" in codes(report)


def test_healthy_network_is_quiet():
    """No errors and no warnings on a sound model, or nobody will read it."""
    report = build(healthy_net()).diagnose(level="basic")
    assert report.errors == [], [i.message for i in report.errors]
    assert report.warnings == [], [i.message for i in report.warnings]
    assert report.ok


# --- state the diagnostics must survive -------------------------------


def test_runs_before_add_opf():
    """A model with no OPF layer must still be diagnosable."""
    opf = ACOPF(healthy_net())
    report = opf.diagnose(level="basic")
    assert isinstance(report, DiagnosticReport)
    assert "MODEL_NO_OBJECTIVE" not in codes(report)


def test_checks_that_cannot_run_say_so():
    report = build(healthy_net()).diagnose(level="basic")
    assert "solver_result" in report.skipped
    assert "solve" in report.skipped["solver_result"]


def test_diagnose_does_not_modify_the_model():
    opf = build(healthy_net())
    before = len(list(opf.model.component_objects()))
    opf.diagnose(level="basic")
    assert len(list(opf.model.component_objects())) == before


def test_invalid_level_is_rejected():
    with pytest.raises(ValueError, match="level must be one of"):
        build(healthy_net()).diagnose(level="exhaustive")


# --- formulation awareness --------------------------------------------


def test_dc_model_is_not_given_reactive_diagnostics():
    """Reactive checks on a DC model would be a category error."""
    dc = DCOPF(healthy_net())
    dc.add_OPF()
    from potpourri.diagnostics.context import DiagnosticContext

    ctx = DiagnosticContext.from_model(dc)
    assert not ctx.has_reactive
    assert "DC" in ctx.formulation


def test_ac_model_reports_reactive_capability():
    from potpourri.diagnostics.context import DiagnosticContext

    ctx = DiagnosticContext.from_model(build(healthy_net()))
    assert ctx.has_reactive
    assert ctx.has_voltage_magnitude


# --- Pyomo identifiers stay usable ------------------------------------


def test_findings_keep_the_pyomo_component():
    """The trail from report to constraint must not need guesswork."""
    net = healthy_net()
    net.ext_grid["vm_pu"] = 1.20
    opf = build(net)
    report = opf.diagnose(level="basic")

    labelled = [i for i in report.issues if i.pyomo_component]
    assert labelled, "no finding carried a Pyomo component"
    for issue in labelled:
        assert hasattr(opf.model, issue.pyomo_component), issue.pyomo_component
        assert issue.pyomo_label()


def test_metadata_registry_matches_the_model():
    """Every registered constraint family that exists must be indexable."""
    from potpourri.diagnostics.metadata import CONSTRAINTS

    opf = build(healthy_net())
    present = [name for name in CONSTRAINTS if hasattr(opf.model, name)]
    assert len(present) >= 5, "the registry does not describe this model"
    for name in present:
        meta = CONSTRAINTS[name]
        assert meta.description
        assert isinstance(meta.category, DiagnosticCategory)


# --- solved models ----------------------------------------------------


@pytest.mark.integration
def test_solved_model_reports_binding_limits_and_agrees_with_pandapower():
    opf = build(healthy_net())
    opf.solve(solver="ipopt", print_solver_output=False)

    report = opf.diagnose(level="standard")
    assert "SOLVER_OPTIMAL" in codes(report)
    assert report.by_code("SOLUTION_CONSTRAINT_VIOLATED") == []

    replay = report.by_code("REPLAY_AGREES")
    assert replay, "the pandapower cross-check did not agree"
    assert replay[0].value < 1e-6

    assert "balance.residual" in report.summary
    assert report.errors == [], [i.message for i in report.errors]


@pytest.mark.integration
def test_explain_bus_lists_what_is_attached():
    from potpourri.diagnostics.context import DiagnosticContext
    from potpourri.diagnostics.solution import explain_bus

    opf = build(healthy_net())
    opf.solve(solver="ipopt", print_solver_output=False)
    ctx = DiagnosticContext.from_model(opf)

    detail = explain_bus(ctx, int(opf.net.bus.index[1]))
    assert "vm_pu" in detail
    assert detail["branches"], "no branches were reported for a connected bus"
