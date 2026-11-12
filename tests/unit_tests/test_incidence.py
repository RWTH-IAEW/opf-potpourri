# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Bus incidence lists replace the per-bus element scans of the balances.

The nodal balance used to scan every element for every bus. The lists in
`potpourri.models.incidence` must hold exactly the elements those scans
selected, in the order they visited them, or the balance would differ term
for term. The tests recompute the scans the slow way and compare.
"""

import pandapower as pp
import pandapower.networks as pn
import pytest
import simbench as sb

from potpourri.models.ACOPF_base import ACOPF
from potpourri.models.incidence import bus_incidence
from potpourri.models_multi_period.ACOPF_multi_period import (
    ACOPF_multi_period,
)


def _scan(model, b, storage):
    """The element lists at bus `b`, found by scanning every element."""
    return {
        "sG": [g for g in model.sG if (g, b) in model.sGbs],
        "G": [g for g in model.G if (g, b) in model.Gbs],
        "STOR": (
            [s for s in model.STOR if model.STOR_bus[s] == b]
            if storage
            else []
        ),
        "D": [d for d in model.D if (b, d) in model.Dbs],
        "L_from": [line for line in model.L if model.A[line, 1] == b],
        "L_to": [line for line in model.L if model.A[line, 2] == b],
        "T_hv": [t for t in model.TRANSF if model.AT[t, 1] == b],
        "T_lv": [t for t in model.TRANSF if model.AT[t, 2] == b],
        "SHUNT": [s for s in model.SHUNT if (b, s) in model.SHUNTbs],
    }


def _assert_matches_scan(model, storage):
    inc = bus_incidence(model, storage=storage)
    for b in model.B:
        expected = _scan(model, b, storage)
        for name, elements in expected.items():
            assert list(getattr(inc, name).get(b, ())) == elements, (
                name,
                b,
            )


def _case9_with_storage_and_shunts():
    net = pn.case9()
    for bus, p in ((4, 10.0), (4, 5.0), (6, 5.0)):
        pp.create_storage(
            net,
            bus,
            p_mw=p,
            max_e_mwh=20.0,
            q_mvar=1.0,
            sn_mva=10.0,
            efficiency_percent=95.0,
        )
    net.shunt.drop(net.shunt.index, inplace=True)
    pp.create_shunt(net, 4, q_mvar=-5.0, p_mw=0.5)
    pp.create_shunt(net, 4, q_mvar=-2.0, p_mw=0.2)
    pp.create_shunt(net, 7, q_mvar=3.0, p_mw=0.0)
    return net


def test_single_period_lists_match_the_scans_with_storage_and_shunts():
    model = ACOPF(_case9_with_storage_and_shunts()).model
    assert len(model.STOR) == 3 and len(model.SHUNT) == 3
    _assert_matches_scan(model, storage=True)


def test_lists_match_the_scans_on_a_grid_with_auxiliary_buses():
    """The MV grid has switch nodes, so ppc buses outnumber `net.bus`."""
    net = sb.get_simbench_net("1-MV-rural--0-sw")
    model = ACOPF(net).model
    assert len(model.B) > len(net.bus)
    _assert_matches_scan(model, storage=True)


def test_multi_period_lists_match_the_scans():
    net = sb.get_simbench_net("1-LV-rural1--0-sw")
    model = ACOPF_multi_period(net, toT=2).model
    _assert_matches_scan(model, storage=False)


@pytest.mark.parametrize("name", ["KCL_real", "KCL_reactive"])
def test_rule_called_outside_a_build_gives_the_built_constraint(name):
    """Without the cache a rule falls back to fresh lists, same terms."""
    opf = ACOPF(_case9_with_storage_and_shunts())
    rule = getattr(opf, f"_{name.lower()}_rule")
    constraint = getattr(opf.model, name)
    for b in opf.model.B:
        direct = rule(opf.model, b)
        assert str(direct) == str(constraint[b].expr)
