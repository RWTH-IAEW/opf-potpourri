# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Multi-period hosting capacity: it has to build, and be time-indexed.

`HC_ACOPF_multi_period` could not run at all before this suite existed: the
constructor failed on any network without a `wind_hc` column, `add_OPF()`
built no hosting-capacity layer unless an optional column happened to be
present, `add_loss_obj()` reached for a name nothing created, and the whole
capability layer was single-period code indexing `psG[w]` against a model
whose `psG` is indexed `(g, t)`.

Each test below pins one of those, so none of them can come back quietly.
"""

import math

import pytest

pyo = pytest.importorskip("pyomo.environ")

HORIZON = 3
CASE = "1-LV-rural1--0-sw"


@pytest.fixture(scope="module")
def plain_net():
    """A SimBench network with nothing hosting-capacity specific on it."""
    sb = pytest.importorskip("simbench")
    return sb.get_simbench_net(CASE)


def build(net, horizon=HORIZON, windpot=None):
    import copy

    from potpourri.models_multi_period.HC_ACOPF_multi_period import (
        HC_ACOPF_multi_period,
    )

    net = copy.deepcopy(net)
    if windpot is not None:
        net.bus["windpot_p_mw"] = windpot
    model = HC_ACOPF_multi_period(net, horizon)
    model.add_OPF()
    return model


def test_constructs_from_a_plain_network(plain_net):
    """The candidates get a profile, so the base class can resolve one.

    Before, `pp.create_sgens(...)` left 14 sgens with no SimBench profile
    and construction died in `get_absolute_values` with "These profiles are
    set to be applied but are missing in the profiles data".
    """
    hc = build(plain_net)
    assert len(hc.model.WIND_HC) > 0


def test_builds_the_hc_layer_without_windpot(plain_net):
    """Wind potential is an optional cap, not a precondition.

    Keying the whole layer on `net.bus.windpot_p_mw` meant `add_OPF()`
    returned a plain ACOPF — no candidates, no selection variable, no
    objective — and reported success.
    """
    model = build(plain_net).model
    assert hasattr(model, "y"), "no selection variable: the HC layer is absent"
    assert hasattr(model, "SW2")
    active = [
        o.name
        for o in model.component_data_objects(pyo.Objective, active=True)
    ]
    assert active == ["obj"]


def test_sizing_is_per_candidate_and_dispatch_is_per_step(plain_net):
    """A plant is built once; only its dispatch moves over the horizon."""
    model = build(plain_net).model
    n_cand, n_t = len(model.WIND_HC), len(model.T)
    assert n_t == HORIZON

    # decided once
    assert len(model.y) == n_cand
    assert len(model.SW2) == n_cand
    assert len(model.hc_size_upper) == n_cand
    assert len(model.hc_size_lower) == n_cand

    # decided per step
    for name in ("SW_max", "QW_min", "QW_max", "QU_min_hc", "QU_max_hc"):
        comp = model.component(name)
        assert comp is not None, f"{name} is missing"
        assert len(comp) == n_cand * n_t, f"{name} is not time-indexed"


def test_objective_covers_the_whole_horizon(plain_net):
    """A single-step objective would let a candidate look free."""
    model = build(plain_net).model
    text = str(model.obj.expr)
    for t in range(HORIZON):
        assert f",{t}]" in text, f"time step {t} is missing from the objective"


def test_windpot_caps_every_step(plain_net):
    """`PW_max` appears only with the column, and then binds per step."""
    without = build(plain_net).model
    assert without.component("PW_max") is None

    with_cap = build(plain_net, windpot=1.5).model
    pw = with_cap.component("PW_max")
    assert pw is not None
    assert len(pw) == len(with_cap.WIND_HC) * len(with_cap.T)


def test_add_loss_obj_swaps_the_objective(plain_net):
    """It used to raise AttributeError on `obj_hc`, which never existed."""
    hc = build(plain_net)
    hc.add_loss_obj()
    active = [
        o.name
        for o in hc.model.component_data_objects(pyo.Objective, active=True)
    ]
    assert active == ["OBJ_with_loss"]
    text = str(hc.model.OBJ_with_loss.expr)
    for t in range(HORIZON):
        assert f",{t}]" in text


def test_hosting_capacity_needs_a_solve(plain_net):
    hc = build(plain_net)
    with pytest.raises(RuntimeError, match="solve"):
        hc.hosting_capacity_mva()


@pytest.mark.integration
def test_solves_and_respects_the_installed_rating(plain_net):
    """Dispatch stays inside the rating that was chosen, at every step."""
    hc = build(plain_net, windpot=1.5)
    # the binary makes this a MINLP; relax it so IPOPT can take it
    for w in hc.model.WIND_HC:
        hc.model.y[w].domain = pyo.UnitInterval
    hc.solve(solver="ipopt", print_solver_output=False)

    model = hc.model
    for w in model.WIND_HC:
        rating = pyo.value(model.SW2[w])
        for t in model.T:
            s2 = (
                pyo.value(model.psG[w, t]) ** 2
                + pyo.value(model.qsG[w, t]) ** 2
            )
            assert s2 <= rating + 1e-6

    capacity = hc.hosting_capacity_mva()
    assert set(capacity) == set(model.WIND_HC)
    assert all(math.isfinite(v) and v >= 0 for v in capacity.values())
