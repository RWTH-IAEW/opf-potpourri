# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Custom objective function: per-voltage-level weighted voltage deviation.

potpourri does not restrict you to the built-in objectives (such as
``add_voltage_deviation_objective``).  After calling ``add_OPF()`` the model's
Pyomo ``ConcreteModel`` is fully accessible via ``ac.model``, so you can attach
any ``pyo.Objective`` directly — for example one that assigns different weights
per voltage level:

    min  C_lv * Σ (v_b - 1)²   (0.4 kV LV buses)
       + C_mv * Σ (v_b - 1)²   (10 kV MV buses)

This script solves the 1-MVLV-urban-5.303-0-no_sw network three times with
different (C_mv, C_lv) weight pairs:

  1. Equal weights  (C_mv=1, C_lv=1)  → same as add_voltage_deviation_objective
  2. LV-priority   (C_mv=1, C_lv=10) → optimizer pushes LV buses harder toward 1 p.u.
  3. MV-priority   (C_mv=10, C_lv=1) → optimizer pushes MV buses harder toward 1 p.u.

The hPV scenario (high load + near-rated PV output) produces a natural
voltage gradient: LV buses sit below 1 p.u. (load dominates) while MV buses
sit slightly above 1 p.u.  With a reactive-power trade-off across the MV/LV
transformer, the optimizer can push one level closer to 1 p.u. only at the
cost of the other — the weight signal makes this choice explicit.

Network: 1-MVLV-urban-5.303-0-no_sw  (246 buses: 110 kV / 10 kV / 0.4 kV)
Solver : IPOPT

Author: Steffen Kortmann (2026)
"""

from __future__ import annotations

import copy
import warnings

import pandas as pd
import pyomo.environ as pyo
import simbench as sb

from potpourri.models.ACOPF_base import ACOPF

warnings.filterwarnings("ignore")

NET_NAME = "1-MVLV-urban-5.303-0-no_sw"
LOADCASE = (
    "hPV"  # high load + near-rated PV → LV below 1 p.u., MV above 1 p.u.
)

VM_MAX = 1.06
VM_MIN = 0.94

SCENARIOS: list[tuple[str, dict[float, float]]] = [
    # label                                weights_by_kv
    ("equal weights  (C_mv=1,  C_lv=1) ", {}),
    ("LV-priority    (C_mv=1,  C_lv=10)", {0.4: 10.0}),
    ("MV-priority    (C_mv=10, C_lv=1) ", {10.0: 10.0}),
]


# ── Network preparation ───────────────────────────────────────────────────────


def configure_for_opf(net) -> None:
    """Configure OPF bounds and drop elements that break preprocess_grid."""
    f = net.loadcases.loc[LOADCASE]
    net.load.p_mw *= f["pload"]
    net.load.q_mvar *= f["qload"]
    net.ext_grid.vm_pu = f["Slack_vm"]

    # sgen Q bounds from rated apparent power (independent of current dispatch)
    sn = net.sgen["sn_mva"].clip(lower=0.01)
    net.sgen["max_q_mvar"] = 0.4 * sn
    net.sgen["min_q_mvar"] = -0.4 * sn
    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = 0.0
    net.sgen["controllable"] = True

    net.bus["max_vm_pu"] = VM_MAX
    net.bus["min_vm_pu"] = VM_MIN
    net.line["max_loading_percent"] = 100.0
    net.trafo["max_loading_percent"] = 100.0

    net.ext_grid["max_p_mw"] = 1.0e4
    net.ext_grid["min_p_mw"] = -1.0e4
    net.ext_grid["max_q_mvar"] = 1.0e4
    net.ext_grid["min_q_mvar"] = -1.0e4

    # Measurements carry bus references that pandapower's create_continuous_bus_index
    # remaps; if the bus no longer exists after topology cleanup the remap fails.
    if len(net.measurement):
        net.measurement.drop(net.measurement.index, inplace=True)

    # Any remaining switches (the no_sw variant has very few)
    if len(net.switch):
        net.switch.drop(net.switch.index, inplace=True)


# ── Custom objective ──────────────────────────────────────────────────────────


def add_weighted_voltage_deviation_objective(
    ac: ACOPF,
    weights_by_kv: dict[float, float],
    default_weight: float = 1.0,
) -> pyo.Objective:
    """Attach a per-voltage-level weighted voltage-deviation objective.

    This is a custom objective defined *outside* the ACOPF class, attached
    directly to ``ac.model`` — exactly as you would with any Pyomo component.

    The minimum viable form is three lines::

        ac.model.my_obj = pyo.Objective(
            expr=sum(w[b] * (ac.model.v[b] - 1.0)**2 for b in ac.model.B),
            sense=pyo.minimize,
        )

    This function wraps that pattern with a ``net.bus.vn_kv`` weight lookup.

    Args:
        ac: A fully constructed ACOPF instance (after ``add_OPF()``).
        weights_by_kv: Mapping ``{nominal_kv: weight}``.  Buses at unlisted
            voltage levels receive ``default_weight``.
        default_weight: Fallback weight (default 1.0).

    Returns:
        The newly created ``pyo.Objective``, also stored as
        ``ac.model.obj_weighted_v_dev``.
    """
    # Build a reverse map: Pyomo internal bus index → pandapower bus index.
    #
    # ac.bus_lookup  maps  pp_bus_idx  →  ppc (Pyomo) internal bus idx
    # ac.model.B     is the set of ppc bus indices
    pyomo_to_pp: dict[int, int] = {
        int(ac.bus_lookup[pp_idx]): pp_idx for pp_idx in ac.net.bus.index
    }

    # Pre-compute per-bus weights once (avoids repeated lookups inside Pyomo)
    w: dict[int, float] = {}
    for b in ac.model.B:
        pp_idx = pyomo_to_pp.get(b)
        if pp_idx is None:
            w[b] = default_weight
        else:
            kv = float(ac.net.bus.at[pp_idx, "vn_kv"])
            w[b] = weights_by_kv.get(kv, default_weight)

    @ac.model.Objective(sense=pyo.minimize)
    def obj_weighted_v_dev(model: pyo.ConcreteModel) -> pyo.Expression:
        # Non-slack buses: penalise deviation from 1 p.u.
        non_slack = sum(
            w[b] * (model.v[b] - 1.0) ** 2 for b in model.B - model.b0
        )
        # Slack bus: penalise deviation from the load-flow setpoint v_b0
        # so the solver does not push the reference voltage to its limit.
        slack = sum(w[b] * (model.v[b] - model.v_b0[b]) ** 2 for b in model.b0)
        return non_slack + slack

    return ac.model.obj_weighted_v_dev


# ── Result helper ─────────────────────────────────────────────────────────────


def voltage_stats_by_level(ac: ACOPF) -> pd.DataFrame:
    """Return mean / min / max vm_pu grouped by nominal voltage level (kV)."""
    vm = ac.net.res_bus["vm_pu"].copy()
    kv = ac.net.bus["vn_kv"]
    return (
        pd.DataFrame({"vm_pu": vm, "vn_kv": kv})
        .groupby("vn_kv")["vm_pu"]
        .agg(["mean", "min", "max"])
        .rename(columns={"mean": "mean_pu", "min": "min_pu", "max": "max_pu"})
    )


# ── Main ──────────────────────────────────────────────────────────────────────


def main() -> None:
    """Run the analysis this script demonstrates.

    Configuration comes from the module-level constants above, not from the
    command line. Edit those, or import and call this function, to change what
    is run.

    Returns:
        None. Results are printed, and written to the paths named in the
        configuration block where the script produces files.
    """
    print(f"Network  : {NET_NAME}")
    print(f"Loadcase : {LOADCASE}  (high load + near-rated PV)")
    print(f"Bounds   : Vmin={VM_MIN}  Vmax={VM_MAX} p.u.\n")

    net_template = sb.get_simbench_net(NET_NAME)
    configure_for_opf(net_template)

    kv_levels = sorted(net_template.bus.vn_kv.unique())
    print(f"Voltage levels  : {[float(k) for k in kv_levels]} kV")
    print(f"Buses           : {len(net_template.bus)}")
    print(
        f"Controllable sgen: {net_template.sgen.controllable.sum()} ({len(net_template.sgen)} total)"
    )
    print(f"Load total      : {net_template.load.p_mw.sum():.1f} MW")
    print(f"sGen total      : {net_template.sgen.p_mw.sum():.1f} MW\n")

    results = []

    for label, weights_by_kv in SCENARIOS:
        print(f"{'─' * 60}")
        print(f"Scenario: {label}")

        net = copy.deepcopy(net_template)
        ac = ACOPF(net)
        ac.add_OPF(thermal_limit="current", free_slack_vm=True)

        # ── Key pattern: attach any pyo.Objective to ac.model ────────────
        #
        # The minimum form (uniform weights) is just:
        #
        #   ac.model.obj = pyo.Objective(
        #       expr=sum((ac.model.v[b] - 1.0)**2 for b in ac.model.B),
        #       sense=pyo.minimize,
        #   )
        #
        # add_weighted_voltage_deviation_objective() builds the same
        # expression but looks up w[b] from net.bus.vn_kv.
        add_weighted_voltage_deviation_objective(ac, weights_by_kv)

        res = ac.solve(solver="ipopt", print_solver_output=False)
        ok = res is not None and pyo.check_optimal_termination(res)

        tag = (
            "✓ optimal"
            if ok
            else f"✗ {res.solver.termination_condition if res else 'error'}"
        )
        print(f"  Status : {tag}")

        if ok:
            stats = voltage_stats_by_level(ac)
            print("  Voltage stats by voltage level:")
            print(stats.to_string(float_format=lambda x: f"{x:.4f}"))
            results.append({"label": label, "stats": stats})

        print()

    # ── Comparison table ──────────────────────────────────────────────────────
    if len(results) == len(SCENARIOS):
        print("=" * 60)
        print("Comparison: mean vm_pu per voltage level across scenarios")
        print("=" * 60)
        rows = []
        for r in results:
            row = {"Scenario": r["label"].split("(")[0].strip()}
            for kv, stat in r["stats"].iterrows():
                row[f"{float(kv):.1f} kV mean"] = stat["mean_pu"]
            rows.append(row)
        df = pd.DataFrame(rows).set_index("Scenario")
        print(df.to_string(float_format=lambda x: f"{x:.4f}"))


if __name__ == "__main__":
    main()
