"""DC Optimal Power Flow example for potpourri.

Demonstrates the linearised DC power flow formulation in two steps:

1. Pure DC power flow — compare angle results against pandapower's built-in
   DC solver (``pp.rundcpp``).
2. DC OPF — maximise local renewable generation (minimise power drawn from
   the external grid) subject to a line-loading constraint.  Solved with
   GLPK as a linear programme.

The DC approximation ignores voltage magnitudes and reactive power, reducing
the power flow equations to a linear system in bus angle differences.  It is
orders of magnitude faster than AC OPF and useful for rapid operational
pre-screening.

Network: 1-LV-rural1--0-sw (15 buses, summer midday snapshot).

Institut für Elektrische Anlagen und Netze, Digitalisierung und
Energiewirtschaft (IAEW)
(c) 2023, Steffen Kortmann
"""

import copy
import warnings

import numpy as np
import pandapower as pp
import pyomo.environ as pe
import simbench as sb

from potpourri.models.DC import DC
from potpourri.models.DCOPF import DCOPF

warnings.filterwarnings("ignore")

SOLVER = "glpk"
NET_NAME = "1-LV-rural1--0-sw"
PROFILE_IDX = 672  # midday in summer (~slot 7 × 96)


if __name__ == "__main__":
    # ── 1. Load network and apply a single-snapshot profile ───────────────
    net = sb.get_simbench_net(NET_NAME)
    profiles = sb.get_absolute_values(
        net, profiles_instead_of_study_cases=True
    )

    net.sgen["p_mw"] = profiles[("sgen", "p_mw")].iloc[PROFILE_IDX]
    net.load["p_mw"] = profiles[("load", "p_mw")].iloc[PROFILE_IDX]
    net.load["q_mvar"] = profiles[("load", "q_mvar")].iloc[PROFILE_IDX]

    print(f"Network  : {NET_NAME}")
    print(f"Snapshot : profile index {PROFILE_IDX}  (summer midday)")
    print(
        f"PV total : {net.sgen['p_mw'].sum():.4f} MW  "
        f"| Load : {net.load['p_mw'].sum():.5f} MW\n"
    )

    # ── 2. Pure DC power flow — compare against pandapower ────────────────
    pp.rundcpp(net)  # pandapower reference

    pf = DC(net)
    pf.solve(solver=SOLVER, print_solver_output=False)

    print("Bus voltage angles — pandapower vs potpourri DC:")
    print(f"  {'Bus':>5}  {'pp (°)':>10}  {'pyo (°)':>10}  {'|diff|':>10}")
    for b in pf.model.B:
        pp_va = net.res_bus.va_degree.iloc[b]
        pyo_va = pe.value(pf.model.delta[b]) * 180 / np.pi
        print(
            f"  {b:>5}  {pp_va:>10.4f}  {pyo_va:>10.4f}  "
            f"{abs(pp_va - pyo_va):>10.6f}"
        )

    # ── 3. DC OPF — minimise external grid import ─────────────────────────
    net_opf = copy.deepcopy(net)
    net_opf.sgen["controllable"] = True
    net_opf.sgen["max_p_mw"] = net_opf.sgen["p_mw"]
    net_opf.sgen["min_p_mw"] = 0.0  # allow curtailment to zero
    net_opf.ext_grid["max_p_mw"] = 10_000.0
    net_opf.ext_grid["min_p_mw"] = -10_000.0
    net_opf.line["max_loading_percent"] = 80.0

    dcopf = DCOPF(net_opf)
    dcopf.add_OPF()

    dcopf.model.obj = pe.Objective(
        expr=sum(dcopf.model.pG[g] for g in dcopf.model.G),
        sense=pe.minimize,
    )
    dcopf.solve(solver=SOLVER, print_solver_output=False)

    base = dcopf.model.baseMVA

    print("\n== External grid dispatch ==")
    for g in dcopf.model.G:
        p_mw = pe.value(dcopf.model.pG[g]) * base
        print(f"  Generator {g}: {p_mw:+.3f} MW  (+ = import, − = export)")

    print("\n== Static generator dispatch ==")
    for g in dcopf.model.sG:
        p_mw = pe.value(dcopf.model.psG[g]) * base
        print(f"  sgen {g}: {p_mw:.4f} MW")

    print("\n== Line loading — five most loaded ==")
    loadings = {}
    for line in dcopf.model.L:
        p_max = abs(
            max(
                pe.value(dcopf.model.pLfrom[line]),
                abs(pe.value(dcopf.model.pLto[line])),
            )
        )
        smax = pe.value(dcopf.model.SLmax[line]) or 1e-9
        loadings[line] = p_max / smax * 100
    for line, pct in sorted(loadings.items(), key=lambda x: -x[1])[:5]:
        print(f"  Line {line:>3}: {pct:.1f} %")

    print(
        "\nKey takeaway: The DC OPF maximises local PV use as a linear "
        "programme.  It runs in milliseconds but ignores voltage magnitudes "
        "and reactive power — use AC OPF for final verification."
    )
