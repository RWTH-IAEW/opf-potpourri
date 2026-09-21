# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Basemodel needs a power flow only for the ppc tables and a starting point.

When Newton-Raphson diverges from the flat start it must fall back to a DC
power flow instead of refusing to build the model.
"""

import numpy as np
import pandapower as pp
import pytest
from loguru import logger

from potpourri.models.ACOPF_base import ACOPF
from potpourri.models.basemodel import Basemodel


def _feeder(load_mw):
    net = pp.create_empty_network(sn_mva=1.0)
    b0 = pp.create_bus(net, 0.4, max_vm_pu=1.1, min_vm_pu=0.9)
    b1 = pp.create_bus(net, 0.4, max_vm_pu=1.1, min_vm_pu=0.9)
    pp.create_ext_grid(net, b0, vm_pu=1.0)
    pp.create_line(net, b0, b1, 1.0, "NAYY 4x150 SE")
    pp.create_load(net, b1, load_mw, 0.2 * load_mw)
    return net


def test_diverging_feeder_is_the_premise():
    net = _feeder(5.0)  # far beyond what a 1 km 0.4 kV cable can carry
    with pytest.raises(pp.LoadflowNotConverged):
        pp.runpp(net, voltage_depend_loads=False)


def test_basemodel_falls_back_to_dc_power_flow():
    messages = []
    # the package silences its logger for library use; listen in for this test
    logger.enable("potpourri")
    sink = logger.add(lambda m: messages.append(m), level="WARNING")
    try:
        model = Basemodel(_feeder(5.0))
    finally:
        logger.remove(sink)
        logger.disable("potpourri")
    assert any("DC power flow" in m for m in messages)
    # tables are there and the start is the DC one: flat magnitudes, angles set
    assert model.net._ppc["branch"].shape[0] == 1
    assert np.allclose(model.bus_data.v_m.values, 1.0)
    assert model.bus_data.v_a_rad.abs().max() > 0.0
    # the start stays inside the (-π, π) bounds of the angle variables
    assert model.bus_data.v_a_rad.abs().max() <= np.pi


def test_converging_feeder_keeps_the_ac_start():
    model = Basemodel(_feeder(0.05))
    assert (
        model.bus_data.v_m.min() < 1.0
    )  # a real voltage drop, not the flat fallback


def test_acopf_builds_on_the_fallback():
    acopf = ACOPF(_feeder(5.0))
    acopf.add_OPF(thermal_limit="mva")
    assert len(acopf.model.B) == 2


def test_isolated_bus_gets_no_degenerate_kcl():
    """A fully out-of-service bus contributes no constraint.

    A bus whose every branch is out of service has a balance without any
    variable. Its sums are numpy floats, so the equality is a numpy bool, not
    a Python bool; the KCL rules must skip it instead of handing Pyomo a
    constant (case78484_epigrids has such buses).
    """
    import pandapower as pp

    from potpourri.models.ACOPF_base import ACOPF
    from potpourri.models.DCOPF import DCOPF

    net = pp.create_empty_network(sn_mva=100.0)
    b0 = pp.create_bus(net, 110.0)
    b1 = pp.create_bus(net, 110.0)
    b2 = pp.create_bus(net, 110.0)
    pp.create_ext_grid(net, b0)
    pp.create_line(net, b0, b1, 10.0, "149-AL1/24-ST1A 110.0")
    pp.create_line(
        net, b1, b2, 10.0, "149-AL1/24-ST1A 110.0", in_service=False
    )
    pp.create_load(net, b1, 20.0, 5.0)
    pp.create_shunt(
        net, b2, q_mvar=0.5, p_mw=0.0
    )  # only a shunt: no variable in the balance
    for builder, kwargs in ((DCOPF, {}), (ACOPF, dict(thermal_limit="mva"))):
        model = builder(net)
        model.add_OPF(**kwargs)
        assert len(model.model.B) == 3


def _feeder_with_zero_reactance():
    net = pp.create_empty_network(sn_mva=100.0)
    b0 = pp.create_bus(net, 110.0)
    b1 = pp.create_bus(net, 110.0)
    b2 = pp.create_bus(net, 110.0)
    pp.create_ext_grid(net, b0)
    pp.create_line(net, b0, b1, 10.0, "149-AL1/24-ST1A 110.0")
    # a purely resistive tie, as PGLib case1803_snem ships two of
    pp.create_line_from_parameters(
        net,
        b1,
        b2,
        1.0,
        r_ohm_per_km=0.5,
        x_ohm_per_km=0.0,
        c_nf_per_km=0.0,
        max_i_ka=1.0,
    )
    pp.create_load(net, b2, 5.0, 1.0)
    return net


def test_zero_reactance_breaks_every_pandapower_power_flow():
    net = _feeder_with_zero_reactance()
    with pytest.raises(FloatingPointError):
        pp.runpp(net, voltage_depend_loads=False)
    with pytest.raises(FloatingPointError):
        pp.rundcpp(net)


def test_basemodel_builds_tables_without_a_power_flow():
    messages = []
    logger.enable("potpourri")
    sink = logger.add(lambda m: messages.append(m), level="WARNING")
    try:
        model = Basemodel(_feeder_with_zero_reactance())
    finally:
        logger.remove(sink)
        logger.disable("potpourri")
    assert any("No power flow could be run" in m for m in messages)
    assert model.net._ppc["branch"].shape[0] == 2
    assert np.allclose(model.bus_data.v_m.values, 1.0)
    assert np.allclose(model.bus_data.v_a_rad.values, 0.0)
    assert list(model.net._ppc["internal"]["ref_gens"]) == [0]
    acopf = ACOPF(_feeder_with_zero_reactance())
    acopf.add_OPF(thermal_limit="mva")
    assert len(acopf.model.B) == 3
