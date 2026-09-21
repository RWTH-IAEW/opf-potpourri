# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Warm-start a multi-period model from per-step power flows.

Initialises the Pyomo variables from a pandapower power flow
solved at each time step.

A cold-started multi-period AC OPF begins with ``v = 1``, every angle at zero
and every branch flow at zero, which violates Kirchhoff's laws at every bus by
the full nodal injection. On a nonconvex problem IPOPT can fail to recover from
that: it reports a locally infeasible point on models that are demonstrably
feasible. Seeding the whole state from a power flow — voltages, angles, branch
flows and generation together, so the starting point is *consistent* — is what
avoids it. The seeded operating point does not have to be near the optimum; it
has to satisfy the power flow.
"""

import copy
from math import pi

import pandapower as pp
from loguru import logger

DEG_TO_RAD = pi / 180.0


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

        for bus in scratch.bus.index:
            ppc_bus = int(bus_lookup[bus])
            _set("v", (ppc_bus, t), scratch.res_bus.vm_pu[bus])
            _set(
                "delta",
                (ppc_bus, t),
                scratch.res_bus.va_degree[bus] * DEG_TO_RAD,
            )

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
