"""Minimal AC power flow: potpourri vs pandapower.

End-to-end pipeline:

  pandapower network
    → AC Pyomo model (Basemodel + AC mixin)
    → IPOPT solve
    → pyo_to_net writes net.res_*
    → compare against pandapower Newton-Raphson

Network: 1-LV-rural1--0-sw  (15 buses, 1 ext grid, 4 PV sgens, ~0.4 kV)

Key mappings
------------
Basemodel.__init__   extracts buses/lines/loads into Pyomo sets and params
AC.__init__          adds line/trafo admittances and reactive-power variables
solve(solver="ipopt")  calls SolverFactory("ipopt"), then pyo_to_net writes
                       results back to ac.net.res_bus and ac.net.res_line

Institut für Elektrische Anlagen und Netze, Digitalisierung und
Energiewirtschaft (IAEW)
(c) 2023, Steffen Kortmann
"""

import warnings

import pandas as pd
import pandapower as pp
import simbench as sb

from potpourri.models.AC import AC

warnings.filterwarnings("ignore")


def _bus_comparison(pp_net, ac_net):
    """Side-by-side bus voltage table with absolute deviations."""
    df = pd.DataFrame(
        {
            "pp_vm_pu": pp_net.res_bus.vm_pu,
            "ac_vm_pu": ac_net.res_bus.vm_pu,
            "Δvm_pu": (ac_net.res_bus.vm_pu - pp_net.res_bus.vm_pu).abs(),
            "pp_va_deg": pp_net.res_bus.va_degree,
            "ac_va_deg": ac_net.res_bus.va_degree,
            "Δva_deg": (
                ac_net.res_bus.va_degree - pp_net.res_bus.va_degree
            ).abs(),
        }
    )
    return df


def _line_comparison(pp_net, ac_net):
    """Side-by-side line flow table with absolute deviations.

    pl_mw  = total active power loss  (p_from_mw + p_to_mw)
    ql_mvar = total reactive loss       (q_from_mvar + q_to_mvar)
    """
    df = pd.DataFrame(
        {
            "pp_pl_mw": pp_net.res_line.pl_mw,
            "ac_pl_mw": ac_net.res_line.pl_mw,
            "Δpl_mw": (ac_net.res_line.pl_mw - pp_net.res_line.pl_mw).abs(),
            "pp_ql_mvar": pp_net.res_line.ql_mvar,
            "ac_ql_mvar": ac_net.res_line.ql_mvar,
            "Δql_mvar": (
                ac_net.res_line.ql_mvar - pp_net.res_line.ql_mvar
            ).abs(),
        }
    )
    return df


def _summary(bus_df, line_df):
    """Compact MAE / max absolute error for each compared quantity."""
    rows = [
        {
            "quantity": col,
            "MAE": df[col].mean(),
            "max_AE": df[col].max(),
        }
        for df, col in [
            (bus_df, "Δvm_pu"),
            (bus_df, "Δva_deg"),
            (line_df, "Δpl_mw"),
            (line_df, "Δql_mvar"),
        ]
    ]
    return pd.DataFrame(rows).set_index("quantity")


if __name__ == "__main__":
    NET_NAME = "1-LV-rural1--0-sw"

    # --- load network ---
    # AC.__init__ deep-copies the net internally, so the original `net` stays
    # unmodified and can be used for the pandapower reference solve below.
    net = sb.get_simbench_net(NET_NAME)

    print(f"Network: {NET_NAME}")
    print(
        f"  buses: {len(net.bus)}  |  "
        f"lines: {len(net.line)}  |  "
        f"ext_grids: {len(net.ext_grid)}  |  "
        f"sgens: {len(net.sgen)}  |  "
        f"loads: {len(net.load)}"
    )

    # --- potpourri AC model ---
    ac = AC(net)

    # The bare AC model leaves sgen reactive power (qsG) as a free variable.
    # Without reactive bounds (added by add_OPF()), the power flow equations
    # are under-determined when sgens are present: IPOPT can find a feasible
    # but physically different solution.  Fixing qsG to the nominal net values
    # (typically 0 for PV) makes the problem equivalent to a fixed-PQ power
    # flow and ensures agreement with pandapower.
    for g in ac.model.sG:
        ac.model.qsG[g].fix(ac.static_generation_data["q"].get(g, 0.0))

    result = ac.solve(solver="ipopt", print_solver_output=False)

    converged = (
        result is not None
        and result.solver.termination_condition.value == "optimal"
    )
    print(f"\nIPOPT termination: {result.solver.termination_condition.value}")

    if not converged:
        print("Solver did not converge – aborting comparison.")
        raise SystemExit(1)

    # --- pandapower reference (Newton-Raphson) ---
    pp.runpp(net, voltage_depend_loads=False)

    # --- comparison ---
    pd.set_option("display.float_format", "{:.6f}".format)
    pd.set_option("display.width", 120)

    bus_df = _bus_comparison(net, ac.net)
    line_df = _line_comparison(net, ac.net)
    summary = _summary(bus_df, line_df)

    print("\nBus voltages (pandapower vs potpourri AC):")
    print(bus_df.to_string())

    print("\nLine losses (pandapower vs potpourri AC):")
    print(line_df.to_string())

    print("\nError summary:")
    print(summary.to_string(float_format="{:.2e}".format))
