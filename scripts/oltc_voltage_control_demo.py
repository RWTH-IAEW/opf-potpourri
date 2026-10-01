# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""On-load tap changers as OPF decisions on a SimBench LV network.

A sunny-day window of ``1-LV-rural1--0-sw`` (15 buses, a 160 kVA 20/0.4 kV
regulated distribution transformer with a ±2 × 2.5 % tap changer, four PV
units). The PV is scaled up so that it pushes the feeder past the voltage
limit at the neutral tap, the PV is curtailable, and the DSO objective is

    minimise  curtailed energy  +  W_V · Σ (v − 1)²  +  c · tap operations

over the horizon. Four ways to treat the tap are compared:

  A. fixed       -- the tap stays where the network data puts it (neutral);
  B. continuous  -- the tap position is a real variable (a relaxation);
  C. discrete    -- integer positions, obtained by relax → round → fix →
                    re-solve with IPOPT (``solve_oltc_round_and_fix``), with
                    at most one position per 15 min and a limited number of
                    operations;
  D. minlp       -- the same discrete model solved globally with Gurobi 12+
                    (optional, slow on a horizon; off by default);
  E. controller  -- pandapower's local ``DiscreteTapControl`` holding the
                    0.4 kV busbar in a voltage band, one power flow per step.
                    Not an OPF: it curtails nothing and knows nothing about
                    the feeder ends, so it is the reference the OPF is
                    judged against, not a competitor on the same objective.

The script prints objective, curtailed energy, tap operations and voltage
extremes per case and writes four panels (min/max voltage, tap position,
transformer loading, curtailment over time) to ``results/``.

SimBench delivers the tap data but leaves ``tap_changer_type`` None, which
pandapower 3.x reads as "no tap changer"; the script sets it to ``"Ratio"``
first -- the one piece of data preparation the feature needs.

Why an LV network: it is the smaller example, and the regulated 20/0.4 kV
distribution transformer with ±2 positions shows the effect of each tap
step clearly. (The MV rural network was the first choice; when this script
was written the multi-period AC model did not converge on the SimBench MV
networks, which came down to an incomplete warm start that has since been
fixed -- see ``docs/research/dso_controllable_equipment.md``, Section 11.)
The single-period OLTC tests cover 110/20 kV units with ±9 positions.
"""

import copy
import os
import time
import warnings

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandapower as pp  # noqa: E402
import pandas as pd  # noqa: E402
import pyomo.environ as pyo  # noqa: E402
import simbench as sb  # noqa: E402
from pandapower.control import DiscreteTapControl  # noqa: E402

from potpourri.models_multi_period.ACOPF_multi_period import (  # noqa: E402
    ACOPF_multi_period,
)

warnings.filterwarnings("ignore")

# ── Configuration ────────────────────────────────────────────────────────────
NET_NAME = "1-LV-rural1--0-sw"
DAY = 146  # day of the SimBench year (0-based); the sunniest noon
FIRST_STEP = 40  # 15-min step of that day to start at (40 = 10:00)
N_STEPS = 24  # horizon length (24 × 15 min = 6 h)
V_MIN, V_MAX = 0.95, 1.03  # voltage band for every bus (LV share of ±10 %)
PV_SCALING = (
    2.0  # installed PV × 2: 1.037 p.u. at the neutral tap, 92 % loading
)
W_VOLTAGE = 0.01  # small weight of Σ (v − 1)²: a tie-breaker, not a goal
TAP_COST = 0.002  # cost of one tap operation, in objective units
MAX_CHANGE_PER_STEP = 1  # positions per 15 min
MAX_OPERATIONS = 4  # operations over the horizon (the range is ±2)
SOLVER = "ipopt"
USE_GUROBI = False  # case D: global MINLP with gurobi_direct_minlp
GUROBI_TIME_LIMIT = 600  # seconds
CONTROLLER_BAND = (0.99, 1.02)  # DiscreteTapControl band at the 0.4 kV busbar
RESULTS_DIR = "results"
PRINT_SOLVER_OUTPUT = False
# ─────────────────────────────────────────────────────────────────────────────


def prepare_network():
    """Load the SimBench network and make it an OLTC/curtailment case.

    Returns:
        The prepared pandapower network (with SimBench profiles attached).
    """
    net = sb.get_simbench_net(NET_NAME)
    net.trafo["tap_changer_type"] = "Ratio"  # pandapower 3.x needs this
    net.bus["min_vm_pu"] = V_MIN
    net.bus["max_vm_pu"] = V_MAX
    net.sgen["controllable"] = True
    net.sgen["min_p_mw"] = 0.0
    # SimBench profiles are per unit of p_mw, so scaling the installed PV
    # power scales the whole infeed profile.
    pv = net.sgen["type"].str.contains("PV", case=False, na=False)
    net.sgen.loc[pv, "p_mw"] *= PV_SCALING
    # The multi-period model enforces the inverter rating sn_mva when the
    # column exists; scale it with the installed power, or the rating caps
    # the infeed and looks like curtailment.
    if "sn_mva" in net.sgen.columns:
        net.sgen.loc[pv, "sn_mva"] *= PV_SCALING
    return net


def dso_objective(mp):
    """Attach the DSO objective: curtailed energy plus voltage deviation.

    Args:
        mp: The multi-period model object (after `add_OPF`).
    """
    model = mp.model
    dt = float(pyo.value(model.deltaT))

    @model.Objective(sense=pyo.minimize)
    def obj_dso(model):
        """Curtailed sgen energy (p.u.·h) plus weighted voltage deviation.

        Args:
            model: The Pyomo model.

        Returns:
            A Pyomo expression to minimise.
        """
        curtailed = sum(
            (model.sPGmax[g, t] - model.psG[g, t]) * dt
            for g in model.sGc
            for t in model.T
        )
        deviation = sum(
            (model.v[b, t] - 1.0) ** 2 for b in model.Bpd for t in model.T
        )
        return curtailed + W_VOLTAGE * deviation


def build(net, from_t, to_t, mode=None):
    """Build the multi-period model for one case.

    Args:
        net: Prepared network.
        from_t: First profile step.
        to_t: Last profile step (exclusive).
        mode: `None` (fixed tap), `"continuous"` or `"discrete"`.

    Returns:
        The model object, with objective attached.
    """
    mp = ACOPF_multi_period(net, fromT=from_t, toT=to_t)
    # The 20 kV side of an LV feeder is set by the upstream network, not by
    # this OPF: pin the slack magnitude at the SimBench value so that only
    # the tap and the PV can act on the LV voltages.
    mp.add_OPF(free_slack_vm=False)
    if mode is not None:
        mp.enable_oltc(
            mode=mode,
            max_change_per_step=MAX_CHANGE_PER_STEP,
            max_operations=MAX_OPERATIONS,
        )
    dso_objective(mp)
    if mode is not None:
        mp.penalize_tap_movement(cost=TAP_COST)
    return mp


def summarise(mp, label):
    """Collect per-step series and scalar indicators of a solved model.

    Args:
        mp: Solved model object.
        label: Case name.

    Returns:
        A dict with the time series (DataFrames) and the indicators.
    """
    model = mp.model
    steps = list(model.T)
    base = float(pyo.value(model.baseMVA))
    dt = float(pyo.value(model.deltaT))
    v = pd.DataFrame(
        {b: [pyo.value(model.v[b, t]) for t in steps] for b in model.Bpd},
        index=steps,
    )
    available = sum(
        pyo.value(model.sPGmax[g, t]) for g in model.sGc for t in steps
    )
    dispatched = sum(
        pyo.value(model.psG[g, t]) for g in model.sGc for t in steps
    )
    curtail_t = [
        sum(
            (pyo.value(model.sPGmax[g, t]) - pyo.value(model.psG[g, t]))
            for g in model.sGc
        )
        * base
        for t in steps
    ]
    loading = {}
    taps = {}
    for tr in model.TRANSF:
        s_hv = [
            math_hypot(
                pyo.value(model.pThv[tr, t]), pyo.value(model.qThv[tr, t])
            )
            * base
            for t in steps
        ]
        loading[tr] = np.array(s_hv) / float(mp.net.trafo.sn_mva.at[tr]) * 100
        if hasattr(model, "TRANSF_OLTC") and tr in model.TRANSF_OLTC:
            taps[tr] = [
                pyo.value(model.trafo_tap_position[tr, t]) for t in steps
            ]
        else:
            taps[tr] = [float(mp.net.trafo.tap_pos.at[tr])] * len(steps)
    operations = (
        mp.tap_operations().to_dict() if hasattr(model, "TRANSF_OLTC") else {}
    )
    # Line loading through the result mapper, one step at a time (net.res_*
    # has no time dimension).
    line_loading = []
    for t in steps:
        mp.map_to_net(t)
        line_loading.append(float(mp.net.res_line.loading_percent.max()))
    return {
        "label": label,
        "steps": steps,
        "vmin": v.min(axis=1).values,
        "vmax": v.max(axis=1).values,
        "taps": pd.DataFrame(taps, index=steps),
        "loading": pd.DataFrame(loading, index=steps),
        "line_max": max(line_loading),
        "curtail_mw": np.array(curtail_t),
        "objective": float(pyo.value(model.obj_dso)),
        "curtailed_mwh": max(0.0, (available - dispatched) * base * dt),
        "operations": operations,
    }


def math_hypot(p, q):
    """Apparent power from active and reactive power.

    Args:
        p: Active power.
        q: Reactive power.

    Returns:
        sqrt(p² + q²).
    """
    return float(np.hypot(p, q))


def controller_case(net, from_t, to_t):
    """Step-by-step pandapower power flows with a local tap controller.

    Args:
        net: Prepared network (profiles attached).
        from_t: First profile step.
        to_t: Last profile step (exclusive).

    Returns:
        A summary dict in the same shape as `summarise`, plus the number
        of steps on which the OPF voltage band was violated.
    """
    sim = copy.deepcopy(net)
    profiles = sb.get_absolute_values(
        sim, profiles_instead_of_study_cases=True
    )
    for tr in sim.trafo.index:
        DiscreteTapControl(
            sim,
            element_index=int(tr),
            vm_lower_pu=CONTROLLER_BAND[0],
            vm_upper_pu=CONTROLLER_BAND[1],
            side="lv",
        )
    steps = list(range(from_t, to_t))
    vmin, vmax, taps, loading, violations = [], [], {}, {}, 0
    line_loading = []
    for t in steps:
        sim.load["p_mw"] = profiles[("load", "p_mw")].loc[t].values
        sim.load["q_mvar"] = profiles[("load", "q_mvar")].loc[t].values
        sim.sgen["p_mw"] = profiles[("sgen", "p_mw")].loc[t].values
        pp.runpp(sim, run_control=True)
        line_loading.append(float(sim.res_line.loading_percent.max()))
        vmin.append(sim.res_bus.vm_pu.min())
        vmax.append(sim.res_bus.vm_pu.max())
        violations += int(vmax[-1] > V_MAX + 1e-6 or vmin[-1] < V_MIN - 1e-6)
        for tr in sim.trafo.index:
            taps.setdefault(tr, []).append(float(sim.trafo.tap_pos.at[tr]))
            loading.setdefault(tr, []).append(
                float(sim.res_trafo.loading_percent.at[tr])
            )
    taps = pd.DataFrame(taps, index=steps)
    operations = {
        tr: int(np.abs(np.diff([0.0] + taps[tr].tolist())).sum())
        for tr in taps.columns
    }
    return {
        "label": "E controller",
        "steps": steps,
        "vmin": np.array(vmin),
        "vmax": np.array(vmax),
        "taps": taps,
        "loading": pd.DataFrame(loading, index=steps),
        "line_max": max(line_loading),
        "curtail_mw": np.zeros(len(steps)),
        "objective": float("nan"),
        "curtailed_mwh": 0.0,
        "operations": operations,
        "violations": violations,
    }


def plot(cases, path):
    """Four panels over time for every case.

    Args:
        cases: List of summary dicts.
        path: Output PNG path.
    """
    fig, axes = plt.subplots(2, 2, figsize=(13, 8), sharex=True)
    hours = None
    for case in cases:
        steps = np.array(case["steps"])
        hours = (steps % 96) / 4.0
        axes[0, 0].plot(hours, case["vmax"], label=f"{case['label']} max")
        axes[0, 0].plot(
            hours, case["vmin"], "--", label=f"{case['label']} min"
        )
        tr0 = case["taps"].columns[0]
        axes[0, 1].step(
            hours, case["taps"][tr0].values, where="post", label=case["label"]
        )
        axes[1, 0].plot(
            hours, case["loading"][tr0].values, label=case["label"]
        )
        axes[1, 1].plot(hours, case["curtail_mw"], label=case["label"])
    axes[0, 0].axhline(V_MAX, color="k", lw=0.8, ls=":")
    axes[0, 0].axhline(V_MIN, color="k", lw=0.8, ls=":")
    axes[0, 0].set_ylabel("bus voltage [p.u.]")
    axes[0, 0].set_title("network voltage extremes")
    axes[0, 1].set_ylabel("tap position")
    axes[0, 1].set_title("tap position, transformer 0")
    axes[1, 0].set_ylabel("loading [%]")
    axes[1, 0].set_title("transformer 0 loading")
    axes[1, 0].set_xlabel("hour of day")
    axes[1, 1].set_ylabel("curtailed PV [MW]")
    axes[1, 1].set_title("total curtailment")
    axes[1, 1].set_xlabel("hour of day")
    for ax in axes.flat:
        ax.grid(alpha=0.3)
    axes[0, 0].legend(fontsize=7, ncol=2)
    axes[0, 1].legend(fontsize=8)
    fig.suptitle(
        f"{NET_NAME}: on-load tap changer as an OPF decision "
        f"(day {DAY}, {N_STEPS} × 15 min)"
    )
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    """Run the four cases, print the comparison and write the figure."""
    os.makedirs(RESULTS_DIR, exist_ok=True)
    net = prepare_network()
    from_t = DAY * 96 + FIRST_STEP
    to_t = from_t + N_STEPS
    cases = []

    for label, mode in (("A fixed", None), ("B continuous", "continuous")):
        tic = time.time()
        mp = build(net, from_t, to_t, mode)
        results = mp.solve(
            solver=SOLVER, print_solver_output=PRINT_SOLVER_OUTPUT
        )
        print(
            f"{label}: {results.solver.termination_condition} in "
            f"{time.time() - tic:.1f} s"
        )
        cases.append(summarise(mp, label))

    tic = time.time()
    mp = build(net, from_t, to_t, "discrete")
    results = mp.solve_oltc_round_and_fix(
        solver=SOLVER, print_solver_output=PRINT_SOLVER_OUTPUT
    )
    print(
        f"C discrete (round & fix): {results.solver.termination_condition} in "
        f"{time.time() - tic:.1f} s; relaxed objective "
        f"{mp.rounding_info['objective_relaxed']:.5f}"
    )
    cases.append(summarise(mp, "C discrete"))
    schedule = mp.tap_schedule()

    if USE_GUROBI:
        tic = time.time()
        mg = build(net, from_t, to_t, "discrete")
        results = mg.solve(
            solver="gurobi_direct_minlp",
            time_limit=GUROBI_TIME_LIMIT,
            print_solver_output=PRINT_SOLVER_OUTPUT,
        )
        print(
            f"D minlp (Gurobi): {results.solver.termination_condition} in "
            f"{time.time() - tic:.1f} s"
        )
        cases.append(summarise(mg, "D minlp"))

    controller = controller_case(net, from_t, to_t)
    cases.append(controller)

    print(
        "\ncase             objective  curtailed [MWh]  tap ops  "
        "v_min   v_max   trafo max [%]  line max [%]"
    )
    for case in cases:
        ops = ", ".join(
            f"{k}:{v:g}" if isinstance(v, float) else f"{k}:{v}"
            for k, v in case["operations"].items()
        )
        print(
            f"{case['label']:<16} {case['objective']:>9.5f}  "
            f"{case['curtailed_mwh']:>14.3f}  {ops:<8} "
            f"{case['vmin'].min():.4f}  {case['vmax'].max():.4f}  "
            f"{case['loading'].max().max():>13.1f}  "
            f"{case['line_max']:>12.1f}"
        )
    print(
        f"\nE controller: {controller['violations']} of {N_STEPS} steps "
        f"outside the [{V_MIN}, {V_MAX}] band (it curtails nothing)"
    )
    print("\nC discrete tap schedule:")
    print(schedule.T.to_string())

    path = os.path.join(RESULTS_DIR, "oltc_voltage_control_demo.png")
    plot(cases, path)
    print(f"\nfigure written to {path}")


if __name__ == "__main__":
    main()
