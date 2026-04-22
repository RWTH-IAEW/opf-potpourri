"""Compare NEOS solver backends on SimBench networks.

Solves six SimBench benchmark networks with the ``AC`` power flow model via
several NEOS NLP solvers and records both solve time and result accuracy
relative to pandapower Newton-Raphson.

Useful for identifying which NEOS solvers converge reliably and how solve
time scales with network size.

Workflow:
  For each network in ``NETS`` × each solver in ``SOLVERS``:
    1. Build an ``AC`` model and call ``solve(solver="neos", neos_opt=...)``.
    2. Run a reference power flow with ``pp.runpp``.
    3. Record elapsed time, MAE for vm_pu / va_degree / pl_mw / ql_mvar.
  Results are saved to ``results/solver_comparison_results.csv``.

Requirements:
  - A NEOS account email in the NEOS_EMAIL environment variable::

        export NEOS_EMAIL=your@email.com

  - Internet access (NEOS submissions are sent to neos-server.org).
"""

import os
import pathlib
import time
import warnings

import pandas as pd
import pandapower as pp
import simbench as sb
from tqdm import tqdm

from potpourri.models.AC import AC

warnings.filterwarnings("ignore")

SOLVERS = ["ipopt", "knitro", "bonmin", "couenne", "conopt", "snopt"]

NETS = [
    "1-HV-mixed--0-sw",
    "1-HV-urban--0-sw",
    "1-MV-rural--0-sw",
    "1-MV-urban--0-sw",
    "1-LV-urban6--0-sw",
    "1-LV-rural1--0-sw",
]

RESULT_KEYS = {
    "res_bus": ["vm_pu", "va_degree"],
    "res_line": ["pl_mw", "ql_mvar"],
}

RESULTS_DIR = pathlib.Path("results")


if __name__ == "__main__":
    os.environ.setdefault("NEOS_EMAIL", "your@email.com")

    records = []

    for net_name in tqdm(NETS, desc="Networks"):
        net = sb.get_simbench_net(net_name)
        pp.runpp(net, voltage_depend_loads=False)  # reference power flow

        for solver in tqdm(SOLVERS, desc="Solvers", leave=False):
            try:
                ac = AC(net)
                t0 = time.perf_counter()
                result = ac.solve(
                    solver="neos",
                    neos_opt=solver,
                    print_solver_output=False,
                )
                elapsed = time.perf_counter() - t0

                converged = (
                    result is not None
                    and result.solver.termination_condition.value == "optimal"
                )

                row = {
                    "network": net_name,
                    "solver": solver,
                    "time_s": round(elapsed, 2),
                    "converged": converged,
                }
                for table, keys in RESULT_KEYS.items():
                    for key in keys:
                        mae = (
                            (ac.net[table][key] - net[table][key]).abs().mean()
                        )
                        row[f"mae_{key}"] = mae

                records.append(row)

            except Exception as exc:
                print(f"  [{net_name}] {solver}: {exc}")
                records.append(
                    {
                        "network": net_name,
                        "solver": solver,
                        "time_s": None,
                        "converged": False,
                    }
                )

    df = pd.DataFrame(records)
    RESULTS_DIR.mkdir(exist_ok=True)
    out = RESULTS_DIR / "solver_comparison_results.csv"
    df.to_csv(out, index=False)

    print("\nSolver comparison results:")
    print(df.to_string(index=False))
    print(f"\nSaved to {out}")
