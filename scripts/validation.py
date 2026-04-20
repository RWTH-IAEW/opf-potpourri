import os
import warnings

import pandas as pd
import pandapower as pp
import simbench as sb
from tqdm import tqdm

from potpourri.models.AC import AC

# Set a valid email address registered with NEOS (https://neos-server.org)
os.environ["NEOS_EMAIL"] = os.environ.get(
    "NEOS_EMAIL", "your-email@example.com"
)
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
    """Calculate comparison metrics between two result series."""
    delta = candidate - reference
    abs_delta = delta.abs()

    return {
        "mean_abs": abs_delta.mean(),
        "max_abs": abs_delta.max(),
        "rmse": (delta.pow(2).mean()) ** 0.5,
        "std_abs": abs_delta.std(),
    }


def compare_results(pp_net, ac_net, delta_keys: dict) -> list[dict]:
    """Compare selected result tables and columns of two solved networks."""
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

        ac = AC(net)
        results = ac.solve(solver="neos")
        converged = (
            results is not None
            and results.solver.termination_condition.value == "optimal"
        )

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

    # Optional exports
    # results_df.to_csv("detailed_comparison.csv", index=False)
    # pivot_mean.to_csv("comparison_mean_abs.csv")
    # pivot_max.to_csv("comparison_max_abs.csv")
    # latex_table = pivot_mean.to_latex(float_format="{:.2e}".format)
    # print(latex_table)
