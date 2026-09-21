# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Generator cost objective for OPF models.

Reads polynomial cost coefficients from ``net.poly_cost`` (the table that
pandapower's MATPOWER importer fills with ``mpc.gencost``) and builds a Pyomo
objective expression of the form

    Σ_e  c2_e · P_e^2 + c1_e · P_e + c0_e

where ``P_e`` is the active-power dispatch of generation element ``e`` in MW.
Internal Pyomo variables are per-unit on ``baseMVA``; this helper handles the
conversion.

Used by ACOPF (quadratic, ``c2`` allowed) and DCOPF (linear by default;
quadratic terms can be enabled when the solver supports QP).
"""

from __future__ import annotations

import pyomo.environ as pyo


def _ext_grid_or_gen_index(net, et: str, element: int) -> int:
    """Row index in generation_data for a poly_cost entry.

    Map a (et, element) pair from net.poly_cost to the row index used by
    Basemodel.generation_data / model.G.

    Pandapower stacks ext_grid and gen elements into a single _ppc['gen']
    array; the starting offset for each element type is recorded in
    ``net._gen_order``. ``Basemodel.generation_data`` indexes match those
    _ppc rows.
    """
    if et not in ("ext_grid", "gen"):
        raise ValueError(f"_ext_grid_or_gen_index called with et={et!r}")
    f, _ = net._gen_order[et]
    return int(f) + int(element)


def build_poly_cost_expression(opf_model):
    """Build a Pyomo expression for the total polynomial generation cost.

    Args:
        opf_model: A Basemodel-derived object with ``.model`` and ``.net``
            populated. The model must already expose ``pG`` (over ``G``) and
            ``psG`` (over ``sG``) variables, i.e. ``create_model`` and (for
            OPF) ``add_OPF`` have been called.

    Returns:
        A Pyomo expression representing the total $/h cost. Returns 0 when
        ``net.poly_cost`` is empty or missing.
    """
    net = opf_model.net
    model = opf_model.model
    baseMVA = float(opf_model.baseMVA)

    if "poly_cost" not in net or net.poly_cost.empty:
        return pyo.Expression(expr=0.0)

    G_idx = set(model.G)
    sG_idx = set(model.sG)
    expr = 0.0

    for _, row in net.poly_cost.iterrows():
        et = row["et"]
        element = int(row["element"])
        c0 = float(row.get("cp0_eur", 0.0))
        c1 = float(row.get("cp1_eur_per_mw", 0.0))
        c2 = float(row.get("cp2_eur_per_mw2", 0.0))

        if et in ("ext_grid", "gen"):
            g_idx = _ext_grid_or_gen_index(net, et, element)
            if g_idx not in G_idx:
                continue
            p_mw = model.pG[g_idx] * baseMVA
        elif et == "sgen":
            if element not in sG_idx:
                continue
            p_mw = model.psG[element] * baseMVA
        else:
            # storage etc. — PGLib doesn't use them; ignore silently
            continue

        if c2:
            expr = expr + c2 * p_mw * p_mw
        if c1:
            expr = expr + c1 * p_mw
        if c0:
            expr = expr + c0

    return expr


def add_poly_cost_objective(opf_model, allow_quadratic: bool = True):
    """Attach a polynomial generation-cost objective to the Pyomo model.

    Args:
        opf_model: A Basemodel-derived object whose ``.model`` already has the
            OPF sets/variables defined.
        allow_quadratic: When False, raises ``ValueError`` if any cost row has
            a non-zero ``cp2`` term. Set ``False`` for LP DC-OPF runs with
            solvers that don't support QP (e.g. GLPK, CBC).

    Raises:
        ValueError: If ``net.pwl_cost`` contains rows. Piecewise-linear cost
            curves (MATPOWER ``gencost`` type 1) are not yet supported; we
            refuse rather than silently producing the wrong objective.

    Returns:
        The created ``pyo.Objective`` component.
    """
    net = opf_model.net
    model = opf_model.model

    if "pwl_cost" in net and not net.pwl_cost.empty:
        raise ValueError(
            "Piecewise-linear generator cost curves (MATPOWER gencost type "
            "1, pandapower net.pwl_cost) are not supported by "
            "add_poly_cost_objective. Convert them to polynomial form or "
            "implement a PWL objective handler before solving."
        )

    if not allow_quadratic and "poly_cost" in net and not net.poly_cost.empty:
        has_quad = (net.poly_cost["cp2_eur_per_mw2"].abs() > 0).any()
        if has_quad:
            raise ValueError(
                "Network contains quadratic cost terms (cp2 != 0) but "
                "allow_quadratic=False. Use an LP-capable formulation or "
                "switch to a QP solver."
            )

    expr = build_poly_cost_expression(opf_model)
    model.obj_poly_cost = pyo.Objective(expr=expr, sense=pyo.minimize)
    return model.obj_poly_cost
