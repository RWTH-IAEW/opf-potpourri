# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Time-series AC OPF via pandapower's run_timeseries.

Demonstrates how to integrate potpourri's single-period AC OPF into
pandapower's time-series framework by passing a custom run function to
run_timeseries().

At every time step pandapower's controllers write the simbench profile values
into net.load / net.sgen; then run_acopf() builds a fresh ACOPF model, solves
it with IPOPT, and copies the Pyomo solution back to net.res_* so the
OutputWriter can record it.

Workflow:
  1. Load a simbench LV network and fetch absolute 15-min profiles.
  2. Create ConstControl objects for load P/Q and sgen P from DFData sources.
  3. Wire an OutputWriter to log vm_pu, line loading, ext-grid P, and sgen P.
  4. Define run_acopf() as the custom run function.
  5. Call run_timeseries(net, run=run_acopf) for one full day (96 steps).
  6. Print a summary of the logged results.

Network: 1-LV-rural1--0-sw (400 V rural, diverse load / PV mix).

Institut für Elektrische Anlagen und Netze, Digitalisierung und
Energiewirtschaft (IAEW)
Author: Steffen Kortmann (2023)
"""

import os
import warnings

import numpy as np
import pandapower as pp
import simbench as sb
from pandapower.control import ConstControl
from pandapower.timeseries import OutputWriter, run_timeseries
from pandapower.timeseries.data_sources.frame_data import DFData

from potpourri.models.ACOPF_base import ACOPF

warnings.filterwarnings("ignore")

# ── Configuration ─────────────────────────────────────────────────────────────
NET_NAME = "1-LV-rural1--0-sw"
SOLVER = "ipopt"
RESULTS_DIR = "results"
N_TS = 96  # one full day: 96 × 15 min = 24 h
PV_SCALE = 3.0  # scale PV output to create realistic voltage stress
PF_MIN = 0.95  # PV inverter minimum power factor (sets Q capability)
VM_MAX = 1.06  # upper voltage bound [p.u.]
VM_MIN = 0.95  # lower voltage bound [p.u.]
# ──────────────────────────────────────────────────────────────────────────────

_TAN_PHI = np.sqrt(1 - PF_MIN**2) / PF_MIN


# ── custom run function ────────────────────────────────────────────────────────


def run_acopf(net, **kwargs):
    """Build and solve a single-period AC OPF for the current time step.

    Called by run_timeseries at every step after the controllers have updated
    net.load and net.sgen from the profile.  The Pyomo results are written
    back to net.res_* so the OutputWriter can log them.

    Raises RuntimeError on solver failure so run_timeseries can register the
    step as failed (handled gracefully when continue_on_divergence=True).
    """
    # Active-power curtailment ceiling follows the profile value for this step.
    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = 0.0
    # Q capability from cos(phi)_min envelope at the current P set-point.
    net.sgen["max_q_mvar"] = _TAN_PHI * net.sgen["p_mw"]
    net.sgen["min_q_mvar"] = -_TAN_PHI * net.sgen["p_mw"]

    ac = ACOPF(net)
    ac.add_OPF()
    ac.add_voltage_deviation_objective()
    result = ac.solve(solver=SOLVER, print_solver_output=False)

    term = (
        result.solver.termination_condition.value
        if result is not None
        else "no_result"
    )
    if term != "optimal":
        raise RuntimeError(f"AC OPF did not converge (termination: {term})")

    # Write Pyomo results back to the original net for the OutputWriter.
    # Use .update() (partial overwrite) rather than full replacement so that
    # bus indices merged by preprocess_grid (bus-to-bus switch collapsing)
    # are preserved in net.res_bus with their initialised values.
    for table in (
        "res_bus",
        "res_line",
        "res_sgen",
        "res_load",
        "res_ext_grid",
        "res_trafo",
    ):
        src = getattr(ac.net, table, None)
        if src is None or src.empty:
            continue
        dst = getattr(net, table, None)
        if dst is None or dst.empty:
            setattr(net, table, src.copy())
        else:
            dst.update(src)


# ── main ───────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # ── 1. Load simbench network and profiles ─────────────────────────────
    print(f"Network  : {NET_NAME}")
    print("Loading simbench profiles ...", end="  ", flush=True)
    net = sb.get_simbench_net(NET_NAME)
    profiles = sb.get_absolute_values(
        net, profiles_instead_of_study_cases=True
    )
    n_slots = len(profiles[("load", "p_mw")])
    print(f"OK  ({n_slots} slots × 15 min = {n_slots // 96} days)")

    # Static OPF limit columns required by ACOPF
    net.bus["max_vm_pu"] = VM_MAX
    net.bus["min_vm_pu"] = VM_MIN
    net.line["max_loading_percent"] = 80.0
    net.sgen["controllable"] = True
    net.sgen["max_p_mw"] = net.sgen["p_mw"]
    net.sgen["min_p_mw"] = 0.0
    net.ext_grid["max_q_mvar"] = 500.0
    net.ext_grid["min_q_mvar"] = -500.0

    # Pre-populate net.res_* with all original bus/line indices so the
    # OutputWriter can find them even after preprocess_grid merges some buses.
    pp.runpp(net, voltage_depend_loads=False)

    # ── 2. Create DFData sources from the first N_TS profile rows ─────────
    # Profile DataFrames: rows = time steps (int index), cols = element indices.
    time_steps = range(N_TS)

    load_p = profiles[("load", "p_mw")].iloc[:N_TS].reset_index(drop=True)
    load_q = profiles[("load", "q_mvar")].iloc[:N_TS].reset_index(drop=True)
    sgen_p = (
        profiles[("sgen", "p_mw")].iloc[:N_TS].reset_index(drop=True)
        * PV_SCALE
    )

    ds_load_p = DFData(load_p)
    ds_load_q = DFData(load_q)
    ds_sgen_p = DFData(sgen_p)

    # ── 3. Wire ConstControl objects ──────────────────────────────────────
    load_idx = net.load.index.tolist()
    sgen_idx = net.sgen.index.tolist()

    ConstControl(
        net,
        element="load",
        variable="p_mw",
        element_index=load_idx,
        data_source=ds_load_p,
        profile_name=load_idx,
    )
    ConstControl(
        net,
        element="load",
        variable="q_mvar",
        element_index=load_idx,
        data_source=ds_load_q,
        profile_name=load_idx,
    )
    ConstControl(
        net,
        element="sgen",
        variable="p_mw",
        element_index=sgen_idx,
        data_source=ds_sgen_p,
        profile_name=sgen_idx,
    )

    # ── 4. Configure OutputWriter ─────────────────────────────────────────
    os.makedirs(RESULTS_DIR, exist_ok=True)
    ow = OutputWriter(
        net,
        time_steps=time_steps,
        output_path=RESULTS_DIR,
        output_file_type=".csv",
    )
    ow.log_variable("res_bus", "vm_pu")
    ow.log_variable("res_line", "loading_percent")
    ow.log_variable("res_ext_grid", "p_mw")
    ow.log_variable("res_sgen", "p_mw")

    # ── 5. Run time series with AC OPF at each step ───────────────────────
    print(
        f"\nRunning {N_TS} × single-period AC OPF "
        f"(solver={SOLVER}, Vmax={VM_MAX} p.u., PV×{PV_SCALE}) ...\n",
        flush=True,
    )
    run_timeseries(
        net,
        time_steps=time_steps,
        run=run_acopf,
        continue_on_divergence=True,
        verbose=True,
    )

    # ── 6. Print summary of logged results ────────────────────────────────
    res_vm = ow.output["res_bus.vm_pu"]
    res_ll = ow.output["res_line.loading_percent"]
    res_eg = ow.output["res_ext_grid.p_mw"]
    res_sg = ow.output["res_sgen.p_mw"]

    n_failed = ow.output["Parameters"]["powerflow_failed"].sum()
    n_ok = N_TS - n_failed

    print("\n" + "─" * 60)
    print(f"RESULTS  ({n_ok}/{N_TS} time steps converged)")
    print("─" * 60)
    print(
        f"  Voltage:      max = {res_vm.max().max():.4f} p.u. "
        f"  min = {res_vm.min().min():.4f} p.u."
    )
    print(f"  Line loading: max = {res_ll.max().max():.2f} %")
    print(
        f"  Ext-grid P  : [{res_eg.min().min():.4f}, "
        f"{res_eg.max().max():.4f}] MW  (neg = export)"
    )
    print(
        f"  Sgen P total: max = {res_sg.sum(axis=1).max():.4f} MW "
        f"(after curtailment)"
    )
    print("─" * 60)

    print(
        "\nCSV results written to:",
        os.path.abspath(RESULTS_DIR),
    )
    print(
        "\nKey takeaway: run_timeseries() with run=run_acopf solves one "
        "independent AC OPF per time step.  Unlike ACOPF_multi_period, "
        "there is no coupling across steps — no ramp constraints, no storage "
        "state of charge.  Use ACOPF_multi_period for coupled optimisation."
    )
