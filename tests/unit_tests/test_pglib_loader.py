"""The PGLib loader must hand potpourri the network MATPOWER describes.

Four conversions went wrong before and each is pinned here with a small
MATPOWER case written to ``tmp_path`` (no submodule needed):

* dropping out-of-service generators without renumbering ``poly_cost``
  attached every later cost curve to the wrong unit (case200_activ: +14 %);
* rebalancing the initial dispatch zeroed the sgens ``from_mpc`` creates
  from negative demand (case240_pserc lost 4.6 GW of injection: +4.9 %);
* those sgens were marked controllable although they are fixed injections;
* every transformer tap was encoded on the high-voltage side although
  MATPOWER's ``TAP`` acts on the from bus (case162_ieee_dtc, case300_ieee:
  AC-OPF infeasible within the voltage band).
"""

import numpy as np
import pandapower as pp
import pytest

from potpourri.benchmarks.pglib import load_pglib_case

# bus 1: slack (type 3), 110 kV; bus 2: PV bus with two generators, one of
# them out of service; bus 3: negative demand (a fixed 20 MW injection);
# bus 4: 20 kV behind a transformer whose MATPOWER from bus is the LV side.
TINY_CASE = """function mpc = tiny
mpc.version = '2';
mpc.baseMVA = 100;
mpc.bus = [
\t1\t3\t0\t0\t0\t0\t1\t1\t0\t110\t1\t1.1\t0.9;
\t2\t2\t50\t10\t0\t0\t1\t1\t0\t110\t1\t1.1\t0.9;
\t3\t1\t-20\t0\t0\t0\t1\t1\t0\t110\t1\t1.1\t0.9;
\t4\t1\t10\t2\t0\t0\t1\t1\t0\t20\t1\t1.1\t0.9;
];
mpc.gen = [
\t1\t0\t0\t100\t-100\t1\t100\t1\t200\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0;
\t2\t30\t0\t50\t-50\t1\t100\t1\t80\t5\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0;
\t2\t10\t0\t50\t-50\t1\t100\t0\t40\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0;
];
mpc.branch = [
\t1\t2\t0.01\t0.05\t0.02\t250\t250\t250\t0\t0\t1\t-30\t30;
\t2\t3\t0.01\t0.05\t0.02\t250\t250\t250\t0\t0\t1\t-30\t30;
\t4\t2\t0.002\t0.04\t0\t100\t100\t100\t1.05\t0\t1\t-30\t30;
];
mpc.gencost = [
\t2\t0\t0\t3\t0.01\t20\t100;
\t2\t0\t0\t3\t0.02\t30\t200;
\t2\t0\t0\t3\t0.03\t40\t300;
];
"""


@pytest.fixture
def tiny_case(tmp_path):
    path = tmp_path / "pglib_opf_tiny.m"
    path.write_text(TINY_CASE)
    return path


def test_out_of_service_generator_takes_its_cost_row_along(tiny_case):
    net = load_pglib_case(tiny_case)
    assert len(net.gen) == 1
    gen_costs = net.poly_cost[net.poly_cost.et == "gen"]
    assert list(gen_costs.element) == [0]
    # the surviving generator keeps its own curve, not the dropped unit's
    assert gen_costs.cp1_eur_per_mw.iloc[0] == pytest.approx(30.0)
    assert 40.0 not in net.poly_cost.cp1_eur_per_mw.values
    assert net.poly_cost.cp0_eur.sum() == pytest.approx(300.0)


def test_negative_demand_stays_a_fixed_injection(tiny_case):
    net = load_pglib_case(tiny_case)
    sgen = net.sgen[net.sgen.bus == 2]
    assert len(sgen) == 1
    assert sgen.p_mw.iloc[0] == pytest.approx(
        20.0
    )  # not zeroed by the rebalance
    assert not bool(sgen.controllable.iloc[0])
    # the real generator is dispatchable and was rebalanced towards the load
    assert bool(net.gen.controllable.iloc[0])
    assert 5.0 <= net.gen.p_mw.iloc[0] <= 80.0


def test_transformer_tap_sits_on_matpower_from_bus(tiny_case):
    net = load_pglib_case(tiny_case)
    assert len(net.trafo) == 1
    trafo = net.trafo.iloc[0]
    assert trafo.lv_bus == 3  # MATPOWER from bus 4 is the 20 kV side
    assert trafo.tap_side == "lv"
    # opting out keeps pandapower's encoding
    raw = load_pglib_case(tiny_case, align_tap_sides=False)
    assert raw.trafo.iloc[0].tap_side == "hv"


def test_loaded_case_runs_a_power_flow_that_matches_matpower_admittance(
    tiny_case,
):
    """With the tap on the right side the diagonal of Ybus follows MATPOWER:
    the LV bus sees ys / TAP², the HV bus sees ys."""
    from pandapower.pypower.makeYbus import makeYbus

    net = load_pglib_case(tiny_case)
    pp.rundcpp(net)
    ppc = net._ppc
    Y = makeYbus(ppc["baseMVA"], ppc["bus"], ppc["branch"])[0].toarray()
    lk = net._pd2ppc_lookups["bus"]
    # series admittance of the transformer branch in p.u. on the 100 MVA base
    ys = 1 / complex(0.002, 0.04)
    # bus 4 (index 3) only carries the transformer; its diagonal is ys / TAP²
    assert Y[lk[3], lk[3]] == pytest.approx(ys / 1.05**2, rel=2e-3)
    # the off-diagonal is -ys / TAP either way
    assert Y[lk[3], lk[1]] == pytest.approx(-ys / 1.05, rel=2e-3)
    assert np.isfinite(Y).all()


def test_baseline_parser_keeps_rows_powermodels_found_infeasible(tmp_path):
    """Most SAD rows carry ``inf.`` in the DC column; their AC reference must
    survive parsing and the DC value must come back as ``inf``."""
    import math

    from potpourri.benchmarks.pglib import parse_baseline_md

    md = tmp_path / "BASELINE.md"
    md.write_text(
        "## Typical Operating Conditions (TYP)\n"
        "| pglib_opf_case5_pjm | 5 | 6 | 1.7480e+04 | 1.7552e+04 | 0.1 | 0.1 | <1 | <1 |\n"
        "## Small Angle Difference Conditions (SAD)\n"
        "| pglib_opf_case5_pjm__sad | 5 | 6 | inf. | 2.6109e+04 | 0.99 | 3.62 | <1 | <1 |\n"
    )
    parsed = parse_baseline_md(md)
    assert parsed["TYP"]["pglib_opf_case5_pjm"] == {
        "dc": 17480.0,
        "ac": 17552.0,
    }
    sad = parsed["SAD"]["pglib_opf_case5_pjm__sad"]
    assert math.isinf(sad["dc"]) and sad["ac"] == 26109.0


# bus 5 (20 kV) hangs off bus 3 through a branch with a nominal ratio between
# different voltage levels: pandapower turns that into an ``impedance``
# element, the third table a MATPOWER branch can land in.
IMPEDANCE_CASE = TINY_CASE.replace(
    "\t4\t1\t10\t2\t0\t0\t1\t1\t0\t20\t1\t1.1\t0.9;\n];",
    "\t4\t1\t10\t2\t0\t0\t1\t1\t0\t20\t1\t1.1\t0.9;\n"
    "\t5\t1\t5\t1\t0\t0\t1\t1\t0\t20\t1\t1.1\t0.9;\n];",
).replace(
    "\t4\t2\t0.002\t0.04\t0\t100\t100\t100\t1.05\t0\t1\t-30\t30;\n];",
    "\t4\t2\t0.002\t0.04\t0\t100\t100\t100\t1.05\t0\t1\t-30\t30;\n"
    "\t3\t5\t0.01\t0.05\t0\t100\t100\t100\t0\t0\t1\t-2\t2;\n];",
)


@pytest.fixture
def impedance_case(tmp_path):
    path = tmp_path / "pglib_opf_tiny_imp.m"
    path.write_text(IMPEDANCE_CASE)
    return path


def test_angle_limits_reach_impedance_branches(impedance_case):
    from potpourri.models.ACOPF_base import ACOPF
    from potpourri.models.DCOPF import DCOPF

    net = load_pglib_case(impedance_case)
    assert len(net.impedance) == 1
    assert list(net.impedance.angmin_degree) == [-2.0]
    assert list(net.impedance.angmax_degree) == [2.0]
    synthetic = len(net.line)  # impedance rows follow the lines in model.L
    for builder, kwargs in (
        (DCOPF, dict(angle_limits=True)),
        (ACOPF, dict(thermal_limit="mva", angle_limits=True)),
    ):
        model = builder(net)
        model.add_OPF(**kwargs)
        assert synthetic in model.model.LineAngleSet
        lo, body, hi = model.model.line_angle_diff[
            synthetic
        ].to_bounded_expression()
        assert lo == pytest.approx(np.deg2rad(-2.0))
        assert hi == pytest.approx(np.deg2rad(2.0))


# the MATPOWER slack bus 1 carries no generator (both units sit at bus 2), as
# in the RTE cases; pandapower then creates no ext_grid
NO_SLACK_GEN_CASE = TINY_CASE.replace(
    "\t1\t0\t0\t100\t-100\t1\t100\t1\t200\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0;\n",
    "",
).replace(
    "\t2\t0\t0\t3\t0.01\t20\t100;\n",
    "",
)


def test_slack_bus_without_generator_gets_a_zero_capacity_ext_grid(tmp_path):
    from potpourri.models.DCOPF import DCOPF

    path = tmp_path / "pglib_opf_tiny_noslack.m"
    path.write_text(NO_SLACK_GEN_CASE)
    net = load_pglib_case(path)
    assert len(net.ext_grid) == 1
    eg = net.ext_grid.iloc[0]
    assert eg.bus == 0 and eg.max_p_mw == 0.0 and eg.min_p_mw == 0.0
    assert eg.max_q_mvar == 0.0 and eg.min_q_mvar == 0.0
    # the model builds (a reference bus exists) and the reference does not
    # appear in the cost objective
    model = DCOPF(net)
    model.add_OPF(angle_limits=True)
    assert 0 in model.model.b0
    assert (net.poly_cost.et == "ext_grid").sum() == 0
