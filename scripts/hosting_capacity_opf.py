"""Hosting capacity analysis with HC_ACOPF for potpourri.

The hosting capacity (HC) of a distribution grid is the maximum amount of
distributed generation that can be connected without violating operational
limits.  ``HC_ACOPF`` extends the AC OPF with:

  - Binary variables ``y[w] ∈ {0, 1}`` — whether wind generator w is active.
  - Apparent power envelope per generator: ``SWmin·y ≤ S² ≤ SWmax·y``.
  - VDE-AR-N 4105 Q-P and Q-U droop curves enforced by grid-code constraints.
  - Default objective: maximise total wind generation minus network losses.

If the network's sgen table has no ``wind_hc`` column, HC_ACOPF automatically
places one candidate wind generator at every non-slack bus.

Workflow:
  1. Run a baseline AC power flow.
  2. Solve the HC optimisation and identify which buses can host wind power.
  3. Sweep the ``eps`` parameter to trade off wind generation against losses.
  4. Sweep the ``SWmin`` apparent-power floor to control minimum turbine size.

Network: 1-LV-rural1--0-sw (high-load winter snapshot, PV × 1).

Institut für Elektrische Anlagen und Netze, Digitalisierung und
Energiewirtschaft (IAEW)
(c) 2023, Steffen Kortmann
"""

import copy
import warnings

import pandapower as pp
import pyomo.environ as pe
import simbench as sb

from potpourri.models.HC_ACOPF import HC_ACOPF

warnings.filterwarnings("ignore")

SOLVER = (
    "mindtpy"  # HC uses MindtPy (MINLP); replace with gurobi for real runs
)
# For a production run use: SOLVER = "gurobi_direct" or mindtpy + ipopt/glpk
NET_NAME = "1-LV-rural1--0-sw"
PROFILE_IDX = 1190  # high-load winter evening


if __name__ == "__main__":
    # ── 1. Load network, apply snapshot ──────────────────────────────────
    net = sb.get_simbench_net(NET_NAME)
    profiles = sb.get_absolute_values(
        net, profiles_instead_of_study_cases=True
    )

    net.sgen["p_mw"] = profiles[("sgen", "p_mw")].iloc[PROFILE_IDX]
    net.load["p_mw"] = profiles[("load", "p_mw")].iloc[PROFILE_IDX]
    net.load["q_mvar"] = profiles[("load", "q_mvar")].iloc[PROFILE_IDX]

    net.bus["max_vm_pu"] = 1.05
    net.bus["min_vm_pu"] = 0.95
    net.line["max_loading_percent"] = 80.0
    net.ext_grid["max_q_mvar"] = 500.0
    net.ext_grid["min_q_mvar"] = -500.0

    pp.runpp(net, voltage_depend_loads=False)
    print(f"Network  : {NET_NAME}  |  snapshot {PROFILE_IDX}")
    print(
        f"Baseline : max loading {net.res_line.loading_percent.max():.1f}%"
        f"  |  V ∈ [{net.res_bus.vm_pu.min():.3f}, "
        f"{net.res_bus.vm_pu.max():.3f}] p.u.\n"
    )

    # ── 2. Solve HC-OPF ──────────────────────────────────────────────────
    hc = HC_ACOPF(net)
    hc.add_OPF()
    hc.solve(solver=SOLVER, print_solver_output=False)

    baseMVA = hc.model.baseMVA
    print(f"Objective (wind − losses): {pe.value(hc.model.obj):.4f} p.u.")

    active, inactive = [], []
    for w in hc.model.WIND_HC:
        y = pe.value(hc.model.y[w])
        bus = hc.net.sgen.bus.iloc[w]
        p_mw = pe.value(hc.model.psG[w]) * baseMVA
        q_mvar = pe.value(hc.model.qsG[w]) * baseMVA
        if y and y > 0.5:
            active.append((w, bus, p_mw, q_mvar))
        else:
            inactive.append((w, bus))

    print(f"\nActive wind sites ({len(active)}):")
    print(f"  {'sgen':>5}  {'bus':>5}  {'P (MW)':>8}  {'Q (Mvar)':>10}")
    for w, bus, p, q in active:
        print(f"  {w:>5}  {bus:>5}  {p:>8.3f}  {q:>10.3f}")
    print(
        f"\nTotal hosted wind capacity: {sum(p for _, _, p, _ in active):.2f} MW"
    )

    # ── 3. Minimum turbine size sweep (SWmin) ─────────────────────────────
    print("\nSWmin sweep (minimum generator size):")
    print(f"  {'SWmin (MW)':>12}  {'Active sites':>14}  {'Total (MW)':>12}")
    for swmin_pu in [0.0, 0.01, 0.02, 0.05]:
        hc_s = HC_ACOPF(copy.deepcopy(net))
        hc_s._calc_opf_parameters(SWmin=swmin_pu * baseMVA)
        hc_s.add_OPF()
        hc_s.solve(solver=SOLVER, print_solver_output=False)

        n = sum(
            1
            for w in hc_s.model.WIND_HC
            if (v := pe.value(hc_s.model.y[w])) and v > 0.5
        )
        tot = sum(
            pe.value(hc_s.model.psG[w]) * baseMVA
            for w in hc_s.model.WIND_HC
            if (v := pe.value(hc_s.model.y[w])) and v > 0.5
        )
        print(f"  {swmin_pu * baseMVA:>12.2f}  {n:>14}  {tot:>12.2f}")

    # ── 4. Wind-vs-loss trade-off sweep (eps) ─────────────────────────────
    print("\neps sweep (wind−loss trade-off):")
    print(
        f"  {'eps':>6}  {'Active sites':>14}  {'Total (MW)':>12}  "
        f"{'Losses (MW)':>13}"
    )
    for eps in [1.0, 0.8, 0.5]:
        hc_w = HC_ACOPF(copy.deepcopy(net))
        hc_w.add_OPF()
        hc_w.add_loss_obj()
        hc_w.model.eps.set_value(eps)
        hc_w.solve(solver=SOLVER, print_solver_output=False)

        n = sum(
            1
            for w in hc_w.model.WIND_HC
            if (v := pe.value(hc_w.model.y[w])) and v > 0.5
        )
        tot = sum(
            pe.value(hc_w.model.psG[w]) * baseMVA
            for w in hc_w.model.WIND_HC
            if (v := pe.value(hc_w.model.y[w])) and v > 0.5
        )
        losses = sum(
            (pe.value(hc_w.model.pLfrom[l]) + pe.value(hc_w.model.pLto[l]))
            * baseMVA
            for l in hc_w.model.L
        )
        print(f"  {eps:>6.1f}  {n:>14}  {tot:>12.2f}  {losses:>13.4f}")

    print(
        "\nKey takeaway: HC_ACOPF finds which buses can host wind generation "
        "while satisfying AC power flow physics and grid-code Q constraints.  "
        "The eps and SWmin parameters enable sensitivity analysis without "
        "rebuilding the model."
    )
