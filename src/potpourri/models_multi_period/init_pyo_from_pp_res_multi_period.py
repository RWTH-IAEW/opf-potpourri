# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Warm-start a multi-period model from per-step power flows.

Initialises the Pyomo variables from a pandapower power flow
solved at each time step.

A cold-started multi-period AC OPF begins with the bus state of the base
power flow and every branch flow at zero, which violates Kirchhoff's laws at
every bus by the full nodal injection (and a flat start, ``v = 1`` with every
angle at zero, is worse still). On a nonconvex problem IPOPT can fail to
recover from that: it reports a locally infeasible point on models that are
demonstrably feasible. Seeding the whole state from a power flow — voltages,
angles, branch flows and generation together, so the starting point is
*consistent* — is what avoids it. The seeded operating point does not have to
be near the optimum; it has to satisfy the power flow.

The bus state is read from pandapower's internal ppc bus table rather than
from ``net.res_bus``, because the model is indexed over ppc buses and the ppc
can hold buses that ``net.bus`` does not: pandapower reconnects the open end
of a line with an open switch to an auxiliary bus so the line charging stays
in the power flow. Those buses are in service and solved, but ``res_bus`` has
no row for them. Seeding through ``res_bus`` left them at the flat start,
about 150° away from their neighbours behind an HV/MV Dyn5 transformer, which
put a residual of the order of $10^3$ p.u. into the branch-flow equations of
short lines and sent IPOPT to a locally infeasible point on every SimBench MV
network.
"""

import copy
from math import pi

import numpy as np
import pandapower as pp
from loguru import logger
from pandapower.pypower.idx_bus import VA, VM

DEG_TO_RAD = pi / 180.0


def _seed_bus_state(model, scratch, bus_lookup, t, set_value):
    """Seed ``v`` and ``delta`` at step ``t`` for every ppc bus of the model.

    Reads the solved ``VM`` / ``VA`` columns of ``scratch._ppc["bus"]``, which
    pandapower fills for every in-service bus including the auxiliary ones
    that have no ``net.bus`` row (see the module docstring). Isolated buses
    carry NaN there and keep their default initial value.

    The ppc numbering of ``scratch`` is checked against the model's
    ``bus_lookup`` first. The two agree whenever ``scratch`` is a copy of the
    model's own network, which is the only way this module is called; should
    they ever differ, the seed falls back to ``net.res_bus`` through
    ``bus_lookup`` and leaves the auxiliary buses alone.

    Args:
        model: Multi-period Pyomo ConcreteModel.
        scratch: Network on which the power flow for step ``t`` just ran.
        bus_lookup: The model's pandapower-bus to ppc-bus map.
        t: Time step whose variables are seeded.
        set_value: Callback ``(component, index, value)`` that assigns an
            initial value when the variable exists.
    """
    ppc_bus = scratch._ppc["bus"]
    pd_buses = scratch.bus.index.values
    scratch_lookup = scratch._pd2ppc_lookups["bus"]
    same_numbering = len(scratch_lookup) > pd_buses.max() and np.array_equal(
        scratch_lookup[pd_buses], bus_lookup[pd_buses]
    )
    if same_numbering:
        for ppc_index in model.B:
            if ppc_index >= ppc_bus.shape[0]:
                continue
            v_m, v_a = ppc_bus[ppc_index, VM], ppc_bus[ppc_index, VA]
            if np.isfinite(v_m) and np.isfinite(v_a):
                set_value("v", (ppc_index, t), v_m)
                set_value("delta", (ppc_index, t), v_a * DEG_TO_RAD)
        return

    logger.warning(
        "Warm start at t={}: the power flow's ppc bus numbering differs from "
        "the model's; seeding only the buses net.res_bus can reach.",
        t,
    )
    for bus in pd_buses:
        ppc_index = int(bus_lookup[bus])
        set_value("v", (ppc_index, t), scratch.res_bus.vm_pu[bus])
        set_value(
            "delta",
            (ppc_index, t),
            scratch.res_bus.va_degree[bus] * DEG_TO_RAD,
        )


def init_pyo_from_pp_res_multi_period(net, model, bus_lookup, curtailment=1.0):
    """Warm-start every state variable of a multi-period model.

    Runs one pandapower power flow per time step, with the loads and static
    generators set to that step's profile values, and copies the result into
    the variables indexed by that step.

    Args:
        net: pandapower network carrying ``net.profiles`` (the model's own
            ``self.net``). Not modified — the power flows run on a copy.
        model: Multi-period Pyomo ConcreteModel.
        bus_lookup: pandapower-bus to ppc-bus map (``self.bus_lookup``), since
            the model is indexed over ppc bus numbers.
        curtailment: Factor applied to the static-generation profile in the
            seeding power flow. ``1.0`` seeds the uncurtailed state. The
            optimum is insensitive to this — what matters is that the seed
            satisfies the power flow — so it is exposed only for the rare case
            where the uncurtailed state will not converge in pandapower.

    Returns:
        Number of time steps successfully seeded.

    Variables that do not exist on the given model kind are skipped, so this
    works for the AC, LPAC and DC formulations alike.
    """
    base = model.baseMVA.value
    scratch = copy.deepcopy(net)
    profiles = net.profiles
    seeded = 0

    def _set(component, index, value):
        """Assign an initial value when the variable exists on this model."""
        var = getattr(model, component, None)
        if var is None or index not in var:
            return
        var[index].set_value(float(value))

    for t in model.T:
        for element, columns in (
            ("load", ("p_mw", "q_mvar")),
            ("sgen", ("p_mw", "q_mvar")),
        ):
            scale = curtailment if element == "sgen" else 1.0
            for column in columns:
                key = (element, column)
                if key in profiles and t in profiles[key].index:
                    scratch[element][column] = (
                        profiles[key].loc[t].values * scale
                    )

        try:
            pp.runpp(scratch, voltage_depend_loads=False)
        except (pp.LoadflowNotConverged, KeyError) as err:
            logger.warning(
                "Warm-start power flow did not converge at t={}: {}. "
                "Leaving that step at its default initial values.",
                t,
                err,
            )
            continue

        _seed_bus_state(model, scratch, bus_lookup, t, _set)

        for line in scratch.line.index:
            _set("pLfrom", (line, t), scratch.res_line.p_from_mw[line] / base)
            _set("pLto", (line, t), scratch.res_line.p_to_mw[line] / base)
            _set(
                "qLfrom", (line, t), scratch.res_line.q_from_mvar[line] / base
            )
            _set("qLto", (line, t), scratch.res_line.q_to_mvar[line] / base)

        # Transformers are indexed positionally in the model.
        for position, trafo in enumerate(scratch.trafo.index):
            _set(
                "pThv", (position, t), scratch.res_trafo.p_hv_mw[trafo] / base
            )
            _set(
                "pTlv", (position, t), scratch.res_trafo.p_lv_mw[trafo] / base
            )
            _set(
                "qThv",
                (position, t),
                scratch.res_trafo.q_hv_mvar[trafo] / base,
            )
            _set(
                "qTlv",
                (position, t),
                scratch.res_trafo.q_lv_mvar[trafo] / base,
            )

        for sgen in model.sG:
            if sgen < len(scratch.res_sgen):
                _set("psG", (sgen, t), scratch.res_sgen.p_mw.iloc[sgen] / base)
                _set(
                    "qsG",
                    (sgen, t),
                    scratch.res_sgen.q_mvar.iloc[sgen] / base,
                )

        for position, ext_grid in enumerate(scratch.ext_grid.index):
            _set(
                "pG",
                (position, t),
                scratch.res_ext_grid.p_mw[ext_grid] / base,
            )
            _set(
                "qG",
                (position, t),
                scratch.res_ext_grid.q_mvar[ext_grid] / base,
            )

        seeded += 1

    logger.debug("Warm-started {} of {} time steps", seeded, len(model.T))
    return seeded
