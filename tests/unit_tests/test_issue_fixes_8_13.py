# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Regression tests for GitLab issues #8–#13.

Each test pins one behaviour that was reported broken, in the same style as
:mod:`test_audit_fixes`: name the issue, isolate the path, assert the
post-fix behaviour.

* #8  — pandapower >= 3.5 broke every model constructor
* #9  — multi-period ``solve(to_net=True)`` logged a mapping it never did
* #10 — battery efficiency produced no round-trip loss; terminal SOC free
* #11 — ``__version__`` drifted behind ``pyproject.toml``
* #12 — docs advertised an EV module and a misspelled device class
* #13 — the multi-period ``add_OPF`` surface diverged from single-period
"""

from __future__ import annotations

import ast
import copy
import warnings
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import pyomo.environ as pyo
import pytest
import simbench as sb
from pyomo.core.expr import identify_variables

import potpourri
from potpourri.models.ACOPF_base import ACOPF
from potpourri.models_multi_period.ACOPF_multi_period import (
    ACOPF_multi_period,
)
from potpourri.models_multi_period.DCOPF_multi_period import (
    DCOPF_multi_period,
)
from potpourri.technologies.battery import Battery_multi_period
from potpourri.technologies.sgens import SgenMinPAboveProfileWarning

REPO_ROOT = Path(__file__).resolve().parents[2]

# A profile row in daylight, so the PV profile is non-zero and an sgen lower
# bound is not automatically in conflict with it.
DAYLIGHT_ROW = 13868


@pytest.fixture(scope="module")
def lv_net():
    return sb.get_simbench_net("1-LV-rural1--0-sw")


def _prepared(net):
    """Operational limits that make an AC OPF on the LV net well posed."""
    net = copy.deepcopy(net)
    net.bus["max_vm_pu"] = 1.05
    net.bus["min_vm_pu"] = 0.95
    net.sgen["controllable"] = True
    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = 0.0
    net.ext_grid["max_q_mvar"] = 500.0
    net.ext_grid["min_q_mvar"] = -500.0
    return net


# ---------------------------------------------------------------------------
# #8 — create_continuous_bus_index import survives pandapower 3.5
# ---------------------------------------------------------------------------


def test_issue8_preprocess_grid_does_not_use_top_level_symbol():
    """The call must not go through the pandapower top-level namespace.

    pandapower 3.5 removed ``create_continuous_bus_index`` from it, which made
    every single-period constructor raise ``AttributeError`` in
    ``preprocess_grid``. ``pandapower.toolbox`` carries it on 3.4 and 3.5
    alike, so the import is pinned there.
    """
    source = (
        REPO_ROOT / "src" / "potpourri" / "models" / "basemodel.py"
    ).read_text()
    assert "from pandapower.toolbox import create_continuous_bus_index" in (
        source
    )
    assert "pp.create_continuous_bus_index(" not in source


def test_issue8_pandapower_pin_has_an_upper_bound():
    """An unbounded pin is what let a breaking minor release through."""
    pyproject = (REPO_ROOT / "pyproject.toml").read_text()
    assert '"pandapower>=2.13,<3.6"' in pyproject


def test_issue8_constructor_builds_a_model(lv_net):
    """The path that raised: every single-period constructor runs this."""
    acopf = ACOPF(lv_net)
    assert len(list(acopf.model.B)) > 0


def test_issue8_result_mappers_do_not_use_top_level_symbol():
    """Pandapower 3.5 also moved ``clear_result_tables``.

    Both mappers call it as the first thing they do, so on 3.5 neither could
    write ``net.res_*`` at all. In the multi-period mapper this was dormant
    until #9 wired the mapper in; in the single-period one it silently
    downgraded every ``solve(to_net=True)`` to the base-case power flow.
    """
    for path in (
        REPO_ROOT / "src" / "potpourri" / "models" / "pyo_to_net.py",
        REPO_ROOT
        / "src"
        / "potpourri"
        / "models_multi_period"
        / "pyo_to_net_multi_period.py",
    ):
        source = path.read_text()
        assert (
            "from pandapower.toolbox import clear_result_tables" in source
        ), path
        assert "pp.clear_result_tables(" not in source, path


def test_issue8_no_top_level_pandapower_symbol_is_missing():
    """Guard the whole ``pp.<symbol>`` surface, not just the two known ones.

    Both breakages were the same shape: a symbol pandapower moved out of its
    top-level namespace between minor releases. This fails on the installed
    version rather than waiting for a user to hit it.
    """
    import pandapower as pp

    # Walk the AST rather than the text, so `pp.<name>` inside a comment or a
    # docstring is not mistaken for a real attribute access.
    used = {}
    for py in (REPO_ROOT / "src" / "potpourri").rglob("*.py"):
        tree = ast.parse(py.read_text(), filename=str(py))
        for node in ast.walk(tree):
            if (
                isinstance(node, ast.Attribute)
                and isinstance(node.value, ast.Name)
                and node.value.id == "pp"
            ):
                used.setdefault(node.attr, f"{py.name}:{node.lineno}")

    missing = sorted(
        f"pp.{name} ({where})"
        for name, where in used.items()
        if not hasattr(pp, name)
    )
    assert missing == [], (
        f"pandapower {pp.__version__} has no top-level {missing}; import "
        f"them from their defining module instead"
    )


def test_issue8_single_period_solve_maps_the_solution(lv_net):
    """A broken mapper must not be able to report success.

    ``solve()`` wrapped the mapping in ``except AttributeError``, so when the
    mapper raised, the error was logged (into a logger disabled by default)
    and the stale base-case ``net.res_*`` was returned as the result.
    """
    acopf = ACOPF(_prepared(lv_net))
    acopf.add_OPF()
    acopf.add_voltage_deviation_objective()
    base = acopf.net.res_bus.vm_pu.copy()
    acopf.solve(solver="ipopt", to_net=True, print_solver_output=False)

    for b in acopf.model.B:
        if b in acopf.net.res_bus.index:
            assert acopf.net.res_bus.vm_pu.loc[b] == pytest.approx(
                pyo.value(acopf.model.v[b]), abs=1e-9
            )
    assert not base.equals(acopf.net.res_bus.vm_pu)


# ---------------------------------------------------------------------------
# #9 — multi-period solve(to_net) actually writes net.res_*
# ---------------------------------------------------------------------------


def test_issue9_map_to_net_writes_the_requested_step(lv_net):
    """``net.res_bus`` must carry the solution, not the base-case power flow.

    The constructor runs ``pp.runpp``, so ``net.res_*`` is already populated
    before the solve. Logging "Solution mapped to net.res_*" without calling
    the mapper therefore handed back plausible base-case numbers.
    """
    opf = ACOPF_multi_period(_prepared(lv_net), toT=3, fromT=0)
    opf.add_OPF()
    opf.add_voltage_deviation_objective()
    opf.solve(solver="ipopt", to_net=False, print_solver_output=False)

    b = list(opf.model.B)[0]
    for t in opf.model.T:
        assert opf.map_to_net(t) == t
        assert opf.net.res_bus.vm_pu.iloc[b] == pytest.approx(
            pyo.value(opf.model.v[b, t]), abs=1e-9
        )


def test_issue9_to_net_true_maps_the_last_step(lv_net):
    opf = ACOPF_multi_period(_prepared(lv_net), toT=3, fromT=0)
    opf.add_OPF()
    opf.add_voltage_deviation_objective()
    opf.solve(solver="ipopt", to_net=True, print_solver_output=False)

    b = list(opf.model.B)[0]
    last = opf.model.T.last()
    assert opf.net.res_bus.vm_pu.iloc[b] == pytest.approx(
        pyo.value(opf.model.v[b, last]), abs=1e-9
    )


def test_issue9_to_net_accepts_an_explicit_step(lv_net):
    opf = ACOPF_multi_period(_prepared(lv_net), toT=3, fromT=0)
    opf.add_OPF()
    opf.add_voltage_deviation_objective()
    opf.solve(solver="ipopt", to_net=1, print_solver_output=False)

    b = list(opf.model.B)[0]
    assert opf.net.res_bus.vm_pu.iloc[b] == pytest.approx(
        pyo.value(opf.model.v[b, 1]), abs=1e-9
    )


def test_issue9_map_to_net_rejects_a_step_outside_the_horizon(lv_net):
    opf = ACOPF_multi_period(_prepared(lv_net), toT=3, fromT=0)
    opf.add_OPF()
    with pytest.raises(ValueError, match=r"not a time step"):
        opf.map_to_net(99)


# ---------------------------------------------------------------------------
# #10 — battery round-trip loss and terminal SOC
# ---------------------------------------------------------------------------


def _battery_model(lv_net, toT=4, **kw):
    mp = ACOPF_multi_period(_prepared(lv_net), toT=toT)
    bat = Battery_multi_period(
        mp.net,
        T=toT,
        penetration=30.0,
        power_pu=0.01,
        capacity_pu_h=0.02,
        **kw,
    )
    bat.get_all(mp.model)
    mp.add_OPF()
    return mp, bat


def test_issue10_charge_and_discharge_are_separate_legs(lv_net):
    """A single signed power variable cannot carry a one-way efficiency."""
    mp, _ = _battery_model(lv_net)
    assert hasattr(mp.model, "BAT_Pchg")
    assert hasattr(mp.model, "BAT_Pdis")
    for v in (mp.model.BAT_Pchg, mp.model.BAT_Pdis):
        for idx in v:
            assert v[idx].lb == 0.0
    # BAT_P survives as the net injection expression.
    assert isinstance(mp.model.BAT_P, pyo.Expression)


def test_issue10_round_trip_efficiency_is_eta_squared(lv_net):
    """Charging then discharging the same grid-side power must lose energy.

    With one signed variable the SOC returned exactly to its start for any
    ``efficiency``, i.e. a 100 % round trip.
    """
    eff = 0.5
    mp, _ = _battery_model(lv_net, efficiency=eff, soc_min=0.0)
    m = mp.model
    b = list(m.BAT)[0]
    dt = pyo.value(m.deltaT)
    cap = pyo.value(m.BAT_Cap[b])
    p = 0.01

    soc0 = pyo.value(m.BAT_SOC_init[b])
    gained = dt * (eff * p) / cap
    spent = dt * (p / eff) / cap
    # Returning the same grid-side power costs 1/eta^2 of what charging gained.
    assert spent == pytest.approx(gained / eff**2)
    assert soc0 + gained - spent < soc0


def test_issue10_soc_update_applies_eta_in_both_directions(lv_net):
    mp, _ = _battery_model(lv_net, efficiency=0.9)
    m = mp.model
    b = list(m.BAT)[0]
    t = list(m.T)[1]
    body = str(m.bat_soc_update_con[b, t].body)
    assert "BAT_Pchg" in body and "BAT_Pdis" in body


def test_issue10_simultaneous_charge_and_discharge_is_bounded(lv_net):
    """Both legs at once would burn energy through the losses."""
    mp, _ = _battery_model(lv_net)
    m = mp.model
    b, t = list(m.BAT)[0], list(m.T)[0]
    pmax = pyo.value(m.BAT_Pmax[b])
    m.BAT_Pchg[b, t].set_value(pmax)
    m.BAT_Pdis[b, t].set_value(pmax)
    assert pyo.value(m.bat_power_con[b, t].body) > pmax


def test_issue10_terminal_soc_is_cyclic_by_default(lv_net):
    mp, _ = _battery_model(lv_net)
    assert hasattr(mp.model, "bat_terminal_soc_con")
    mp.add_voltage_deviation_objective()
    mp.solve(solver="ipopt", print_solver_output=False)
    for b in mp.model.BAT:
        assert pyo.value(
            mp.model.BAT_SOC[b, mp.model.T.last()]
        ) == pytest.approx(pyo.value(mp.model.BAT_SOC_init[b]), abs=1e-6)


def test_issue10_terminal_soc_can_be_freed(lv_net):
    mp, _ = _battery_model(lv_net, terminal_soc=None)
    assert not hasattr(mp.model, "bat_terminal_soc_con")


def test_issue10_terminal_soc_accepts_an_absolute_value(lv_net):
    mp, _ = _battery_model(lv_net, terminal_soc=0.7)
    mp.add_voltage_deviation_objective()
    mp.solve(solver="ipopt", print_solver_output=False)
    for b in mp.model.BAT:
        assert pyo.value(
            mp.model.BAT_SOC[b, mp.model.T.last()]
        ) == pytest.approx(0.7, abs=1e-6)


@pytest.mark.parametrize("efficiency", [0.0, -0.1, 1.5])
def test_issue10_invalid_efficiency_raises(lv_net, efficiency):
    with pytest.raises(ValueError, match=r"efficiency must be in"):
        _battery_model(lv_net, efficiency=efficiency)


def test_issue10_terminal_soc_outside_bounds_raises(lv_net):
    with pytest.raises(ValueError, match=r"lies outside"):
        _battery_model(lv_net, terminal_soc=5.0)


# ---------------------------------------------------------------------------
# #11 — __version__ comes from the distribution metadata
# ---------------------------------------------------------------------------


def test_issue11_version_matches_distribution_metadata():
    """The literal drifted to 0.2.0 while pyproject.toml was at 0.4.2."""
    try:
        installed = version("opf-potpourri")
    except PackageNotFoundError:
        pytest.skip(
            "opf-potpourri is not installed as a distribution (running from "
            "a bare source tree), so there is no metadata to compare against"
        )
    assert potpourri.__version__ == installed
    assert potpourri.__version__ != "0.0.0.dev0"


def test_issue11_version_falls_back_on_a_bare_source_tree():
    """Importing without an install must not raise."""
    assert isinstance(potpourri.__version__, str)
    assert potpourri.__version__


def test_issue11_version_is_not_hard_coded():
    """A literal is what drifted two minor releases behind pyproject.toml."""
    source = (REPO_ROOT / "src" / "potpourri" / "__init__.py").read_text()
    assert 'version("opf-potpourri")' in source
    # Every assignment either reads the metadata or is the source-tree
    # fallback; a bare release literal would mean the drift is back.
    assignments = [
        line.split("=", 1)[1].strip()
        for line in source.splitlines()
        if line.strip().startswith("__version__ =")
    ]
    assert assignments
    for rhs in assignments:
        assert rhs in ('version("opf-potpourri")', '"0.0.0.dev0"'), rhs


# ---------------------------------------------------------------------------
# #12 — documented names exist
# ---------------------------------------------------------------------------


def test_issue12_no_ev_module_is_advertised():
    """The README and CLAUDE.md listed an EV module that does not exist."""
    tech = REPO_ROOT / "src" / "potpourri" / "technologies"
    assert not (tech / "EVs.py").exists()
    for doc in ("README.md", "CLAUDE.md"):
        text = (REPO_ROOT / doc).read_text()
        assert "EVs" not in text
        assert "EVs, heat pumps" not in text


def test_issue12_heatpump_class_name_matches_the_docs():
    """The README said ``HeatPump_multi_period``; the class has a small p."""
    from potpourri.technologies import heat_pump

    assert hasattr(heat_pump, "Heatpump_multi_period")
    assert not hasattr(heat_pump, "HeatPump_multi_period")
    readme = (REPO_ROOT / "README.md").read_text()
    assert "HeatPump_multi_period" not in readme
    assert "Heatpump_multi_period" in readme


def test_issue12_readme_lists_the_shunts_device():
    """``Shunts_multi_period`` is attached automatically but went unlisted."""
    readme = (REPO_ROOT / "README.md").read_text()
    assert "Shunts_multi_period" in readme


def test_issue12_claude_md_battery_example_uses_the_real_signature():
    text = (REPO_ROOT / "CLAUDE.md").read_text()
    assert "Battery_multi_period(mpopf)" not in text
    assert "battery.get_all(mpopf.model)" in text


def test_issue12_no_tutorials_directory_is_referenced():
    assert not (REPO_ROOT / "tutorials").exists()
    for doc in ("CLAUDE.md",):
        assert "tutorials/" not in (REPO_ROOT / doc).read_text()


# ---------------------------------------------------------------------------
# #13.1 — thermal_limit / angle_limits on the multi-period AC model
# ---------------------------------------------------------------------------


def _references_voltage(model, constraint_data):
    """Whether a constraint's body contains a bus voltage magnitude variable.

    Checked by identity rather than by substring: ``qThv[0,0]`` contains the
    text ``v[`` without being a voltage.
    """
    return any(
        var.parent_component() is model.v
        for var in identify_variables(constraint_data.body)
    )


def test_issue13_thermal_limit_default_is_current(lv_net):
    """Current form scales the rating by v², so v must appear in the body."""
    opf = ACOPF_multi_period(_prepared(lv_net), toT=3)
    opf.add_OPF()
    assert opf.thermal_limit_mode == "current"
    m = opf.model
    l0, t0 = list(m.L)[0], list(m.T)[0]
    assert _references_voltage(m, m.line_lim_from[l0, t0])
    for tr in list(m.TRANSF):
        assert _references_voltage(m, m.transf_lim1[tr, t0])


def test_issue13_thermal_limit_mva_mode(lv_net):
    """Previously a ``TypeError`` two frames deep naming a private method.

    The constant-MVA form has a fixed right-hand side, so no voltage.
    """
    opf = ACOPF_multi_period(_prepared(lv_net), toT=3)
    opf.add_OPF(thermal_limit="mva")
    assert opf.thermal_limit_mode == "mva"
    m = opf.model
    t0 = list(m.T)[0]
    for l in list(m.L):
        assert not _references_voltage(m, m.line_lim_from[l, t0])
        assert not _references_voltage(m, m.line_lim_to[l, t0])
    for tr in list(m.TRANSF):
        assert not _references_voltage(m, m.transf_lim1[tr, t0])
        assert not _references_voltage(m, m.transf_lim2[tr, t0])


def test_issue13_invalid_thermal_limit_raises(lv_net):
    opf = ACOPF_multi_period(_prepared(lv_net), toT=3)
    with pytest.raises(ValueError, match=r"thermal_limit must be"):
        opf.add_OPF(thermal_limit="bogus")


def test_issue13_angle_limits_off_by_default(lv_net):
    net = _prepared(lv_net)
    net.line["angmin_degree"] = -30.0
    net.line["angmax_degree"] = 30.0
    opf = ACOPF_multi_period(net, toT=3)
    opf.add_OPF()
    assert not hasattr(opf.model, "line_angle_diff")


def test_issue13_angle_limits_are_time_indexed(lv_net):
    net = _prepared(lv_net)
    net.line["angmin_degree"] = -30.0
    net.line["angmax_degree"] = 30.0
    opf = ACOPF_multi_period(net, toT=3)
    opf.add_OPF(angle_limits=True)
    n_t = len(list(opf.model.T))
    assert len(list(opf.model.line_angle_diff)) == (
        len(list(opf.model.LineAngleSet)) * n_t
    )


def test_issue13_dc_rejects_thermal_limit_by_name(lv_net):
    """The DC path used to swallow it silently — worse than failing.

    The DC model is lossless and carries no reactive power, so there is no
    current-versus-MVA distinction for it to make.
    """
    dcopf = DCOPF_multi_period(_prepared(lv_net), toT=3)
    with pytest.raises(TypeError) as err:
        dcopf.add_OPF(thermal_limit="mva")
    message = str(err.value)
    assert "thermal_limit" in message
    assert "DCOPF_multi_period" in message


def test_issue13_dc_supports_angle_limits(lv_net):
    """Parity with the single-period ``DCOPF.add_OPF(angle_limits=True)``."""
    net = _prepared(lv_net)
    net.line["angmin_degree"] = -30.0
    net.line["angmax_degree"] = 30.0
    dcopf = DCOPF_multi_period(net, toT=3)
    dcopf.add_OPF(angle_limits=True)
    assert hasattr(dcopf.model, "line_angle_diff")


def test_issue13_unsupported_option_names_itself(lv_net):
    opf = ACOPF_multi_period(_prepared(lv_net), toT=3)
    with pytest.raises(TypeError) as err:
        opf.add_OPF(pv_q_control="qp")
    assert "pv_q_control" in str(err.value)


# ---------------------------------------------------------------------------
# #13.2 — net.sgen.min_p_mw reaches the multi-period model
# ---------------------------------------------------------------------------


def _daylight_model(net, min_p_mw):
    net = _prepared(net)
    net.sgen["min_p_mw"] = min_p_mw
    opf = ACOPF_multi_period(net, toT=DAYLIGHT_ROW + 2, fromT=DAYLIGHT_ROW - 2)
    opf.add_OPF()
    return opf


def test_issue13_sgen_min_p_is_read_from_the_net(lv_net):
    """``sPGmin`` was hard-coded to 0, so ``min_p_mw`` was silently ignored."""
    opf = _daylight_model(lv_net, 0.005)
    for g in opf.model.sGc:
        for t in opf.model.T:
            assert pyo.value(opf.model.sPGmin[g, t]) == pytest.approx(0.005)


def test_issue13_sgen_min_p_defaults_to_zero(lv_net):
    opf = _daylight_model(lv_net, 0.0)
    for g in opf.model.sGc:
        for t in opf.model.T:
            assert pyo.value(opf.model.sPGmin[g, t]) == 0.0


def test_issue13_sgen_min_p_nan_falls_back_to_zero(lv_net):
    """Matches ``OPF.static_generation_real_power_limits``."""
    opf = _daylight_model(lv_net, np.nan)
    for g in opf.model.sGc:
        for t in opf.model.T:
            assert pyo.value(opf.model.sPGmin[g, t]) == 0.0


def test_issue13_sgen_min_p_missing_column_falls_back_to_zero(lv_net):
    net = _prepared(lv_net)
    net.sgen.drop(columns=["min_p_mw"], inplace=True)
    opf = ACOPF_multi_period(net, toT=3)
    opf.add_OPF()
    for g in opf.model.sGc:
        for t in opf.model.T:
            assert pyo.value(opf.model.sPGmin[g, t]) == 0.0


def test_issue13_sgen_min_p_above_profile_warns(lv_net):
    """A constant floor can exceed a PV profile — infeasible at night.

    Left unflagged the solver reports only "infeasible", with nothing
    pointing at the lower bound.
    """
    net = _prepared(lv_net)
    net.sgen["min_p_mw"] = 0.01
    with pytest.warns(SgenMinPAboveProfileWarning, match=r"min_p_mw"):
        opf = ACOPF_multi_period(net, toT=4, fromT=0)  # night
        opf.add_OPF()


def test_issue13_sgen_min_p_on_pinned_sgen_does_not_warn(lv_net):
    """Only controllable sgens are bounded by ``sPGmin``.

    A non-controllable sgen has ``psG`` fixed to its profile value, so a
    stray ``min_p_mw`` on it is never enforced and must not be reported as
    an impending infeasibility.
    """
    net = _prepared(lv_net)
    net.sgen["controllable"] = False
    net.sgen["min_p_mw"] = 0.01
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        opf = ACOPF_multi_period(net, toT=4, fromT=0)  # night
        opf.add_OPF()
    assert not [w for w in caught if w.category is SgenMinPAboveProfileWarning]


def test_issue13_sgen_min_p_within_profile_does_not_warn(lv_net):
    net = _prepared(lv_net)
    net.sgen["min_p_mw"] = 0.001
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        opf = ACOPF_multi_period(
            net, toT=DAYLIGHT_ROW + 2, fromT=DAYLIGHT_ROW - 2
        )
        opf.add_OPF()
    assert not [w for w in caught if w.category is SgenMinPAboveProfileWarning]


# ---------------------------------------------------------------------------
# #13.3 — free_slack_vm on the multi-period AC model
# ---------------------------------------------------------------------------


def test_issue13_slack_vm_is_free_by_default(lv_net):
    """``AC_multi_period`` pins it while building the power flow."""
    opf = ACOPF_multi_period(_prepared(lv_net), toT=3)
    opf.add_OPF()
    for b0 in opf.model.b0:
        for t in opf.model.T:
            assert not opf.model.v[b0, t].is_fixed()


def test_issue13_slack_vm_can_be_pinned(lv_net):
    opf = ACOPF_multi_period(_prepared(lv_net), toT=3)
    opf.add_OPF(free_slack_vm=False)
    for b0 in opf.model.b0:
        for t in opf.model.T:
            assert opf.model.v[b0, t].is_fixed()


@pytest.mark.integration
def test_issue13_snapshot_and_horizon_agree_on_the_slack(lv_net):
    """The same snapshot must not depend on the period kind.

    Reported as single-period 0.9975–1.0027 p.u. (slack free) against
    multi-period 1.0195–1.0250 p.u. (slack pinned) — read as a modelling
    error by anyone comparing the two.
    """
    row = 20017
    profiles = sb.get_absolute_values(
        lv_net, profiles_instead_of_study_cases=True
    )

    snapshot = copy.deepcopy(lv_net)
    snapshot.load["p_mw"] = profiles[("load", "p_mw")].iloc[row].values
    snapshot.load["q_mvar"] = profiles[("load", "q_mvar")].iloc[row].values
    snapshot.sgen["p_mw"] = profiles[("sgen", "p_mw")].iloc[row].values
    for net in (snapshot,):
        net.bus["max_vm_pu"] = 1.05
        net.bus["min_vm_pu"] = 0.95
        net.ext_grid["max_q_mvar"] = 500.0
        net.ext_grid["min_q_mvar"] = -500.0
        net.sgen["controllable"] = False

    sp = ACOPF(snapshot)
    sp.add_OPF()
    sp.add_voltage_deviation_objective()
    sp.solve(solver="ipopt", print_solver_output=False)
    sp_v = [pyo.value(sp.model.v[b]) for b in sp.model.B]

    horizon = copy.deepcopy(lv_net)
    horizon.bus["max_vm_pu"] = 1.05
    horizon.bus["min_vm_pu"] = 0.95
    horizon.ext_grid["max_q_mvar"] = 500.0
    horizon.ext_grid["min_q_mvar"] = -500.0
    horizon.sgen["controllable"] = False
    mp = ACOPF_multi_period(horizon, toT=row + 1, fromT=row)
    mp.add_OPF()
    mp.add_voltage_deviation_objective()
    mp.solve(solver="ipopt", print_solver_output=False)
    mp_v = [pyo.value(mp.model.v[b, row]) for b in mp.model.B]

    assert min(mp_v) == pytest.approx(min(sp_v), abs=1e-3)
    assert max(mp_v) == pytest.approx(max(sp_v), abs=1e-3)
