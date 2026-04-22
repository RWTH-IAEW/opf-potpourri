"""Solver performance benchmark using perfplot.

Measures AC power flow solve time as a function of network size across six
SimBench benchmark networks and multiple NLP solvers (via NEOS).  Uses
`perfplot <https://github.com/nschloe/perfplot>`_ to produce a scaling plot.

Requirements:
  - ``perfplot`` optional dependency::

        pip install potpourri[performance-test]

  - NEOS email::

        export NEOS_EMAIL=your@email.com

  - Internet access.

Output is saved to ``results/perf_comparison.png``.
"""

import os
import pathlib
import warnings

import simbench as sb

from potpourri.models.AC import AC

warnings.filterwarnings("ignore")

SOLVERS = ["ipopt", "knitro", "bonmin", "couenne", "conopt", "snopt"]

NETS = [
    "1-HV-mixed--2-sw",
    "1-HV-urban--2-sw",
    "1-MV-rural--2-sw",
    "1-MV-urban--2-sw",
    "1-LV-urban6--2-sw",
    "1-LV-rural1--2-sw",
]

RESULTS_DIR = pathlib.Path("results")


def _load_networks():
    return [(name, sb.get_simbench_net(name)) for name in NETS]


def setup(n, net_data):
    """Return the network whose bus count is closest to n."""
    return min(net_data, key=lambda x: abs(x[1]["bus"].shape[0] - n))[1]


def _make_kernel(solver):
    def kernel(net):
        ac = AC(net)
        ac.solve(solver="neos", neos_opt=solver, print_solver_output=False)
        return ac

    kernel.__name__ = solver.upper()
    return kernel


if __name__ == "__main__":
    try:
        import perfplot
    except ImportError as exc:
        raise ImportError(
            "perfplot is required: pip install potpourri[performance-test]"
        ) from exc

    os.environ.setdefault("NEOS_EMAIL", "your@email.com")

    net_data = _load_networks()
    bus_counts = sorted({net["bus"].shape[0] for _, net in net_data})

    out = perfplot.bench(
        setup=lambda n: setup(n, net_data),
        kernels=[_make_kernel(s) for s in SOLVERS],
        labels=[s.upper() for s in SOLVERS],
        n_range=bus_counts,
        xlabel="Number of buses",
        equality_check=None,
        target_time_per_measurement=100.0,
    )

    out.show(time_unit="s")

    RESULTS_DIR.mkdir(exist_ok=True)
    out.save(
        RESULTS_DIR / "perf_comparison.png",
        transparent=True,
        bbox_inches="tight",
    )
    print(f"Plot saved to {RESULTS_DIR / 'perf_comparison.png'}")
