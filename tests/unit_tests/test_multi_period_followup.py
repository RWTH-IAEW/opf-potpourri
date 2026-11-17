# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Heat balances, dispatch bounds, warm starts and result indices over time."""

import copy

import numpy as np
import pyomo.environ as pyo
import pytest
import simbench as sb

from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period
from potpourri.models_multi_period.DCOPF_multi_period import DCOPF_multi_period
from potpourri.technologies.demand import Demand_multi_period
from potpourri.technologies.heat_pump import Heatpump_multi_period


@pytest.fixture(scope="module")
def raw_net():
    """An LV network with profiles, kept unchanged by each constructor."""
    return sb.get_simbench_net("1-LV-rural1--0-sw")


def test_heat_capacity_uses_mwh_per_kelvin(raw_net):
    opf = ACOPF_multi_period(raw_net, fromT=4, toT=7)
    hp = Heatpump_multi_period(opf.net, penetration=50.0)
    capacity_j_per_k = (
        500.0 * 1.2 * 1005.0 + 4 * (500.0 / 3) * 0.2 * 2400 * 840
    )
    assert hp.heat_cap == pytest.approx(capacity_j_per_k / 3_600_000_000)


@pytest.mark.parametrize("base_mva", [1.0, 10.0])
@pytest.mark.parametrize("steps", [1, 3])
def test_heat_balance_includes_every_interval(raw_net, base_mva, steps):
    net = copy.deepcopy(raw_net)
    net.sn_mva = base_mva
    opf = ACOPF_multi_period(net, fromT=4, toT=4 + steps)
    hp = Heatpump_multi_period(
        opf.net,
        penetration=50.0,
        power_max_pu=0.004 / base_mva,
        cop=3.0,
        qloss_max=0.002 / base_mva,
        initial_temp_c=18.0,
    )
    hp.heat_cap = 0.05  # MWh/K, independent of the electrical system base.
    hp.get_all(opf.model)
    m = opf.model
    assert len(m.hp_temp_update_con) == len(m.HP) * steps
    for h in m.HP:
        temperature = 18.0
        for t in m.T:
            m.hp_p[h, t].set_value(0.004 / base_mva)
            loss_mw = pyo.value(m.Qloss[h, t]) * base_mva
            temperature += 0.25 * (0.004 * 3.0 - loss_mw) / 0.05
            m.temp[h, t].set_value(temperature)
            assert pyo.value(m.hp_temp_update_con[h, t].body) == pytest.approx(
                0.0, abs=1e-12
            )


def test_heatpump_default_initial_temperature_is_mid_band(raw_net):
    opf = ACOPF_multi_period(raw_net, toT=1)
    hp = Heatpump_multi_period(
        opf.net, penetration=50.0, temp_min_c=18.0, temp_max_c=22.0
    )
    hp.get_all(opf.model)
    for h in opf.model.HP:
        assert pyo.value(opf.model.HP_TempInit[h]) == 20.0


@pytest.mark.parametrize("initial_temp_c", [14.0, 21.0, np.nan])
def test_heatpump_rejects_invalid_initial_temperature(raw_net, initial_temp_c):
    opf = ACOPF_multi_period(raw_net, toT=1)
    with pytest.raises(ValueError, match="initial_temp_c"):
        Heatpump_multi_period(
            opf.net, penetration=50.0, initial_temp_c=initial_temp_c
        )


@pytest.mark.parametrize(
    "model_class",
    [ACOPF_multi_period, DCOPF_multi_period],
)
def test_transformer_results_preserve_network_index(raw_net, model_class):
    net = copy.deepcopy(raw_net)
    assert list(net.trafo.index) == [0]
    net.trafo.index = [5]
    net.switch.loc[net.switch.et == "t", "element"] = 5
    opf = model_class(net, fromT=4, toT=6)
    assert opf.warm_start_from_pf() == 2
    expected = {}
    for t in opf.model.T:
        # Distinct flows identify the requested time step and each winding.
        expected[t] = (0.01 + t * 0.001, -0.009 - t * 0.001)
        opf.model.pThv[0, t].set_value(expected[t][0] / net.sn_mva)
        opf.model.pTlv[0, t].set_value(expected[t][1] / net.sn_mva)
    for t in opf.model.T:
        opf.map_to_net(t=t)
        assert list(opf.net.res_trafo.index) == [5]
        assert opf.net.res_trafo.at[5, "p_hv_mw"] == pytest.approx(
            expected[t][0]
        )
        assert opf.net.res_trafo.at[5, "p_lv_mw"] == pytest.approx(
            expected[t][1]
        )


@pytest.mark.parametrize("element", ["sgen", "load"])
@pytest.mark.parametrize(
    "model_class",
    [ACOPF_multi_period, DCOPF_multi_period],
)
def test_out_of_service_injections_are_excluded(raw_net, element, model_class):
    net = copy.deepcopy(raw_net)
    index = net[element].index[0]
    net[element].at[index, "in_service"] = False
    net[element]["controllable"] = True
    opf = model_class(net, fromT=4, toT=6)
    opf.add_OPF()
    assert opf.warm_start_from_pf() == 2
    opf.map_to_net(t=4)
    result = opf.net["res_" + element]
    assert result.at[index, "p_mw"] == 0.0
    if "AC" in opf.model.name:
        assert result.at[index, "q_mvar"] == 0.0


@pytest.mark.parametrize("controllable", [False, True])
def test_warm_start_preserves_fixed_dispatch(raw_net, controllable):
    net = copy.deepcopy(raw_net)
    net.sgen["controllable"] = controllable
    opf = ACOPF_multi_period(net, fromT=13868, toT=13870)
    opf.add_OPF()
    m = opf.model
    expected = {index: pyo.value(m.PsG[index]) for index in m.psG}
    assert max(expected.values()) > 0.0
    assert opf.warm_start_from_pf(curtailment=0.5) == 2
    for index, profile in expected.items():
        assert m.psG[index].fixed == (not controllable)
        assert pyo.value(m.psG[index]) == pytest.approx(
            0.5 * profile if controllable else profile
        )


@pytest.mark.parametrize("base_mva", [1.0, 10.0])
def test_controllable_load_uses_network_bounds(raw_net, base_mva):
    net = copy.deepcopy(raw_net)
    net.sn_mva = base_mva
    net.load["controllable"] = True
    net.load["min_p_mw"] = 0.002
    net.load["max_p_mw"] = 0.008
    net.load["min_q_mvar"] = -0.003
    net.load["max_q_mvar"] = 0.004
    opf = ACOPF_multi_period(net, fromT=4, toT=6)
    opf.add_OPF()
    m = opf.model
    for d in m.Dc:
        for t in m.T:
            for name, expected in (
                ("PDmin", 0.002),
                ("PDmax", 0.008),
                ("QDmin", -0.003),
                ("QDmax", 0.004),
            ):
                assert pyo.value(getattr(m, name)[d, t]) == pytest.approx(
                    expected / base_mva
                )


def test_load_bound_fallbacks_and_explicit_overrides(raw_net):
    net = copy.deepcopy(raw_net)
    net.sn_mva = 10.0
    net.load["controllable"] = True
    for column in ("min_p_mw", "max_p_mw", "min_q_mvar", "max_q_mvar"):
        net.load[column] = np.nan
    opf = ACOPF_multi_period(net, fromT=4, toT=6)
    opf.add_OPF()
    m = opf.model
    for d in m.Dc:
        for t in m.T:
            assert pyo.value(m.PDmin[d, t]) == 0.0
            assert pyo.value(m.PDmax[d, t]) == pyo.value(m.PD[d, t])
            assert pyo.value(m.QDmin[d, t]) == -abs(pyo.value(m.QD[d, t]))
            assert pyo.value(m.QDmax[d, t]) == abs(pyo.value(m.QD[d, t]))
    demand = next(
        f for f in opf.flexibilities if isinstance(f, Demand_multi_period)
    )
    demand.get_demand_real_power_data(m, max_p_mw=0.008, min_p_mw=0.002)
    demand.get_demand_reactive_data(m, max_q_mvar=0.004, min_q_mvar=-0.003)
    for data, physical_value in (
        (demand.PDmax_data_dict, 0.008),
        (demand.PDmin_data_dict, 0.002),
        (demand.QDmax_data_dict, 0.004),
        (demand.QDmin_data_dict, -0.003),
    ):
        assert all(
            value == pytest.approx(physical_value / 10.0)
            for value in data.values()
        )
