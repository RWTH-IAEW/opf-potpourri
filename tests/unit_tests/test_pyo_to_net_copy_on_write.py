# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""The result mappers must keep writing under pandas copy-on-write.

pandas 3 makes copy-on-write the default, and under it a chained assignment
such as `net.res_bus.vm_pu[net.gen.bus] = ...` updates a temporary copy and
leaves the table untouched, without an error. These tests switch
copy-on-write on under pandas 2 so that a chained write fails here instead
of silently after the upgrade.

The networks deliberately give generators an index that differs from their
bus index: a `.loc` write fed a Series aligns on the Series' index, so it
only looks right when the two happen to coincide.
"""

from __future__ import annotations

from types import SimpleNamespace

import pandapower as pp
import pandas as pd
import pytest

from potpourri.models import pyo_to_net
from potpourri.models_multi_period import pyo_to_net_multi_period


@pytest.fixture(autouse=True)
def copy_on_write():
    """Run each test with pandas 3's copy-on-write semantics."""
    with pd.option_context("mode.copy_on_write", True):
        yield


class _Var:
    """Stand-in for an indexed Pyomo variable: `get_values()` and `[i].value`."""

    def __init__(self, values):
        self._values = values

    def get_values(self):
        return dict(self._values)

    def __getitem__(self, index):
        return SimpleNamespace(value=self._values[index])


class _Stop(Exception):
    """Raised to end `pyo_sol_to_net_res` right after the HC sgen writes."""


def _stop(net):
    raise _Stop


def _dc_net():
    """Five buses; the gens sit on buses 3 and 1, the slack on bus 4."""
    net = pp.create_empty_network()
    for _ in range(5):
        pp.create_bus(net, vn_kv=20.0)
    for b in range(4):
        pp.create_line(net, b, b + 1, 1.0, "NA2XS2Y 1x95 RM/25 12/20 kV")
    pp.create_ext_grid(net, 4, vm_pu=1.01)
    pp.create_gen(net, 3, p_mw=0.5, vm_pu=1.02)
    pp.create_gen(net, 1, p_mw=0.5, vm_pu=1.05)
    pp.create_load(net, 0, p_mw=1.0)
    pp.create_load(net, 2, p_mw=0.5)
    pp.rundcpp(net)
    return net


def _expected_vm(net):
    vm = pd.Series(1.0, index=net.bus.index)
    vm[3], vm[1], vm[4] = 1.02, 1.05, 1.01
    return vm


def _hc_net():
    """Three sgens, of which 7 and 3 are hosting-capacity candidates."""
    net = pp.create_empty_network()
    b = pp.create_bus(net, vn_kv=20.0)
    for idx in (5, 7, 3):
        pp.create_sgen(net, b, p_mw=0.0, q_mvar=0.0, index=idx)
    return net


def test_dc_bus_voltages_single_period():
    net = _dc_net()
    ppc = net._pd2ppc_lookups["bus"][net.bus.index.values]
    model = SimpleNamespace(delta=_Var({b: 0.0 for b in ppc}))

    pyo_to_net._bus_voltage_results_to_net(net, model)

    pd.testing.assert_series_equal(
        net.res_bus.vm_pu, _expected_vm(net), check_names=False
    )


def test_dc_bus_voltages_multi_period():
    net = _dc_net()
    ppc = net._pd2ppc_lookups["bus"][net.bus.index.values]
    model = SimpleNamespace(
        name="DCOPF_multi_period", delta=_Var({(b, 0): 0.0 for b in ppc})
    )

    pyo_to_net_multi_period._bus_voltage_results_to_net(net, model, 0)

    pd.testing.assert_series_equal(
        net.res_bus.vm_pu, _expected_vm(net), check_names=False
    )


def test_hc_candidate_dispatch_single_period(monkeypatch):
    net = _hc_net()
    model = SimpleNamespace(
        WIND_HC=[7, 3],
        baseMVA=SimpleNamespace(value=10.0),
        y=_Var({7: 1.0, 3: 0.0}),
        psG=_Var({7: 0.2, 3: 0.4}),
        qsG=_Var({7: 0.1, 3: 0.3}),
    )
    monkeypatch.setattr(pyo_to_net, "clear_result_tables", _stop)

    with pytest.raises(_Stop):
        pyo_to_net.pyo_sol_to_net_res(net, model)

    assert net.sgen.p_mw.to_dict() == {5: 0.0, 7: 2.0, 3: 0.0}
    assert net.sgen.q_mvar.to_dict() == {5: 0.0, 7: 1.0, 3: 0.0}


def test_hc_candidate_dispatch_multi_period(monkeypatch):
    net = _hc_net()
    model = SimpleNamespace(
        name="HC_ACOPF_multi_period",
        WIND_HC=[7, 3],
        baseMVA=SimpleNamespace(value=10.0),
        y=_Var({7: 1.0, 3: 1.0}),
        psG=_Var({(7, 2): 0.2, (3, 2): 0.4}),
        qsG=_Var({(7, 2): 0.1, (3, 2): 0.3}),
    )
    monkeypatch.setattr(pyo_to_net_multi_period, "clear_result_tables", _stop)

    with pytest.raises(_Stop):
        pyo_to_net_multi_period.pyo_sol_to_net_res(net, model, 2)

    assert net.sgen.p_mw.to_dict() == {5: 0.0, 7: 2.0, 3: 4.0}
    assert net.sgen.q_mvar.to_dict() == {5: 0.0, 7: 1.0, 3: 3.0}
