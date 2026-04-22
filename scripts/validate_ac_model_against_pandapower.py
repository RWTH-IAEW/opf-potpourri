"""Validation of the potpourri AC power flow model against pandapower.

Solves the AC power flow on six SimBench benchmark networks with the
potpourri Pyomo model (IPOPT) and compares the results against pandapower's
built-in Newton-Raphson solver.

Workflow:
  For each network in ``NETS``:
    1. Load the network via SimBench.
    2. Build a potpourri ``AC`` model and solve with IPOPT.
    3. Run a reference power flow with pandapower (``pp.runpp``).
    4. Compare ``res_bus.vm_pu``, ``res_bus.va_degree``, ``res_line.pl_mw``,
       and ``res_line.ql_mvar`` using MAE, max absolute error, RMSE, and std.
    5. Print a detailed table and compact pivot tables.

Expected behaviour:
  The bare AC model leaves sgen reactive power as a free variable.  The
  optimizer may therefore find a reactive dispatch that differs from
  pandapower's PQ assumption.  Small voltage magnitude differences (< 1e-3)
  are expected; active power losses should match closely.
"""

import warnings

import pandas as pd
import pandapower as pp
import simbench as sb
from tqdm import tqdm

from potpourri.models.AC import AC

warnings.filterwarnings("ignore")


DELTA_KEYS = {
    "res_bus": ["vm_pu", "va_degree"],
    "res_line": ["pl_mw", "ql_mvar"],
}

NETS = [
    "1-HV-mixed--0-sw",
    "1-HV-urban--0-sw",
    "1-MV-rural--0-sw",
    "1-MV-urban--0-sw",
    "1-LV-urban6--0-sw",
    "1-LV-rural1--0-sw",
]


def calculate_error_metrics(
    reference: pd.Series, candidate: pd.Series
) -> dict:
    """Compute element-wise error statistics between two result series.

    Args:
        reference: Ground-truth values (e.g. pandapower result).
        candidate: Values to evaluate (e.g. POTPOURRI AC result).

    Returns:
        A dict with keys ``mean_abs``, ``max_abs``, ``rmse``, ``std_abs``,
        each holding a scalar float error metric.
    """
    delta = candidate - reference
    abs_delta = delta.abs()

    return {
        "mean_abs": abs_delta.mean(),
        "max_abs": abs_delta.max(),
        "rmse": (delta.pow(2).mean()) ** 0.5,
        "std_abs": abs_delta.std(),
    }


def compare_results(pp_net, ac_net, delta_keys: dict) -> list[dict]:
    """Compare selected result tables and columns of two solved networks.

    For every (table, column) pair in *delta_keys*, the function calculates
    error metrics of *ac_net* relative to *pp_net* and records the bus/line
    index where the largest deviation occurs.

    Args:
        pp_net: pandapower network after ``pp.runpp``, used as reference.
        ac_net: pandapower network after POTPOURRI ``AC.solve``.
        delta_keys: Mapping ``{result_table: [column, ...]}`` describing which
            result fields to compare.

    Returns:
        A list of dicts, one per (table, column) pair, each containing:
        ``element``, ``variable``, ``mean_abs``, ``max_abs``, ``rmse``,
        ``std_abs``, ``worst_index``.
    """
    records = []

    for res_element, keys in delta_keys.items():
        for key in keys:
            metrics = calculate_error_metrics(
                reference=pp_net[res_element][key],
                candidate=ac_net[res_element][key],
            )

            max_idx = (
                (ac_net[res_element][key] - pp_net[res_element][key])
                .abs()
                .idxmax()
            )

            records.append(
                {
                    "element": res_element,
                    "variable": key,
                    "mean_abs": metrics["mean_abs"],
                    "max_abs": metrics["max_abs"],
                    "rmse": metrics["rmse"],
                    "std_abs": metrics["std_abs"],
                    "worst_index": max_idx,
                }
            )

    return records


if __name__ == "__main__":
    all_results = []

    for net_name in tqdm(NETS, desc="Comparing networks"):
        net = sb.get_simbench_net(net_name)

        # Solve with POTPOURRI AC model via NEOS cloud solver or ipopt.
        ac = AC(net)
        # results = ac.solve(solver="neos")
        results = ac.solve(solver="ipopt")
        converged = (
            results is not None
            and results.solver.termination_condition.value == "optimal"
        )

        # Reference: pandapower Newton-Raphson power flow (voltage-independent loads).
        pp.runpp(net, voltage_depend_loads=False)

        vm_tol = 1e-3
        vm_close = (
            ac.net.res_bus.vm_pu - net.res_bus.vm_pu
        ).abs().max() < vm_tol

        comparison = compare_results(
            pp_net=net, ac_net=ac.net, delta_keys=DELTA_KEYS
        )

        for row in comparison:
            row["net"] = net_name
            row["solver_converged"] = converged
            row["vm_pu_close"] = vm_close
            all_results.append(row)

    results_df = pd.DataFrame(all_results)

    results_df = results_df[
        [
            "net",
            "element",
            "variable",
            "mean_abs",
            "max_abs",
            "rmse",
            "std_abs",
            "worst_index",
            "solver_converged",
            "vm_pu_close",
        ]
    ]

    print("\nDetailed comparison:")
    print(results_df)

    print("\nCompact pivot table (mean absolute error):")
    pivot_mean = results_df.pivot(
        index="net", columns="variable", values="mean_abs"
    )
    print(pivot_mean)

    print("\nCompact pivot table (max absolute error):")
    pivot_max = results_df.pivot(
        index="net", columns="variable", values="max_abs"
    )
    print(pivot_max)

    print("\nNetworks ranked by total mean absolute error:")
    ranking = (
        results_df.groupby("net")["mean_abs"]
        .sum()
        .sort_values()
        .rename("total_mean_abs_error")
    )
    print(ranking)

    # Optional CSV / LaTeX exports
    # results_df.to_csv("detailed_comparison.csv", index=False)
    # pivot_mean.to_csv("comparison_mean_abs.csv")
    # pivot_max.to_csv("comparison_max_abs.csv")
    # latex_table = pivot_mean.to_latex(float_format="{:.2e}".format)
    # print(latex_table)
