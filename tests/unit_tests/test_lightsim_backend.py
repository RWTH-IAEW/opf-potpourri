# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Default AC power-flow backend, result parity and compatibility fallback."""

import copy

import numpy as np
import pandapower as pp

from potpourri.models.ACOPF_base import ACOPF
from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period


def test_default_backend_matches_pandapower_reference(four_bus_net):
    reference = copy.deepcopy(four_bus_net)
    pp.runpp(reference, voltage_depend_loads=False, lightsim2grid=False)
    opf = ACOPF(four_bus_net)
    assert opf.net._options["lightsim2grid"] is True
    np.testing.assert_allclose(
        opf.net.res_bus[["vm_pu", "va_degree"]],
        reference.res_bus[["vm_pu", "va_degree"]],
        atol=1e-8,
        rtol=1e-8,
    )
    np.testing.assert_allclose(
        opf.net.res_line[["p_from_mw", "q_from_mvar"]],
        reference.res_line[["p_from_mw", "q_from_mvar"]],
        atol=1e-8,
        rtol=1e-8,
    )


def test_multiple_slacks_fall_back_to_pandapower(four_bus_net):
    net = copy.deepcopy(four_bus_net)
    pp.create_ext_grid(net, int(net.load.bus.iloc[-1]), vm_pu=1.0)
    opf = ACOPF(net)
    assert opf.net.converged
    assert opf.net._options["lightsim2grid"] is False


def test_multi_period_base_and_every_warm_start_use_lightsim(
    lv_rural_net, monkeypatch
):
    opf = ACOPF_multi_period(lv_rural_net, fromT=4, toT=7)
    assert opf.net._options["lightsim2grid"] is True
    original = pp.runpp
    backends = []

    def record_backend(net, **kwargs):
        original(net, **kwargs)
        backends.append(net._options["lightsim2grid"])

    monkeypatch.setattr(pp, "runpp", record_backend)
    assert opf.warm_start_from_pf() == 3
    assert backends == [True, True, True]
