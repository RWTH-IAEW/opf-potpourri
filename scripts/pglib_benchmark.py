"""Benchmark potpourri AC-OPF and DC-OPF against the PGLib-OPF reference.

PGLib-OPF (https://github.com/power-grid-lib/pglib-opf) is the IEEE PES Power
Grid Library benchmark suite for optimal power flow. Each case ships as a
MATPOWER ``.m`` file with a published reference objective value (DC and AC,
solved by PowerModels.jl + IPOPT) in ``BASELINE.md``.

This script:

1. Loads a PGLib case via :func:`potpourri.benchmarks.load_pglib_case`
   (which attaches MATPOWER ``ANGMIN``/``ANGMAX`` to ``net.line`` /
   ``net.trafo`` so phase-angle limits can be enforced).
2. Builds the potpourri DC- and AC-OPF with PGLib-compatible flags:
   * ``thermal_limit='mva'`` — constant-MVA branch limit (matches MATPOWER /
     PowerModels' ``constraint_thermal_limit_*``).
   * ``free_slack_vm=True`` — slack-bus voltage magnitude bounded by
     ``[Vmin, Vmax]`` rather than pinned.
   * ``angle_limits=True`` — branch phase-angle-difference constraints.
3. Wires the polynomial generator cost from ``net.poly_cost`` as the
   objective via :func:`add_poly_cost_objective`.
4. Reports the objective vs. the PGLib reference.
"""

from __future__ import annotations

import os
import time
import warnings
from pathlib import Path
from typing import Iterable

import pandas as pd
import pyomo.environ as pyo

from potpourri.benchmarks import (
    PGLIB_BASELINE_TYP,
    load_pglib_case,
)
from potpourri.models.ACOPF_base import ACOPF
from potpourri.models.DCOPF import DCOPF
from potpourri.models.cost_objective import add_poly_cost_objective

warnings.filterwarnings("ignore")

# ── Configuration ─────────────────────────────────────────────────────────────
SOLVER = "ipopt"
MAX_BUSES = 300  # skip cases with more buses than this
RUN_DC = True  # include DC-OPF column
RUN_AC = True  # include AC-OPF column
CASES = None  # None → all cases within MAX_BUSES; list of names to override
RESULTS_DIR = Path(__file__).parent / "results"
# ──────────────────────────────────────────────────────────────────────────────


# Cases to skip even when within the size budget, with a short reason.
SKIP_CASES: dict[str, str] = {}


def _select_cases(
    requested: Iterable[str] | None, max_buses: int
) -> list[tuple[str, int]]:
    """Choose which cases to benchmark.

    Returns a list of ``(case_name_without_prefix, n_nodes)`` tuples.
    """
    if requested:
        names = []
        for c in requested:
            n = c if c.startswith("pglib_opf_") else f"pglib_opf_{c}"
            if n not in PGLIB_BASELINE_TYP:
                print(f"  ! '{c}' not in baseline; skipping")
                continue
            names.append(n)
    else:
        names = list(PGLIB_BASELINE_TYP.keys())

    out: list[tuple[str, int]] = []
    for n in names:
        nodes = _case_node_count(n)
        if nodes is None or nodes > max_buses:
            continue
        out.append((n.replace("pglib_opf_", ""), nodes))
    out.sort(key=lambda t: t[1])
    return out


def _case_node_count(case_name: str) -> int | None:
    """Return the node count from the baseline metadata, parsed lazily."""
    return _NODE_COUNTS.get(case_name)


def _build_node_counts() -> dict[str, int]:
    """Read BASELINE.md once to extract node counts per case."""
    import re
    from potpourri.benchmarks.pglib import PGLIB_ROOT

    out: dict[str, int] = {}
    md = PGLIB_ROOT / "BASELINE.md"
    if not md.is_file():
        return out
    in_typ = False
    row_re = re.compile(r"^\|\s*(pglib_opf_\S+?)\s*\|\s*(\d+)\s*\|\s*\d+\s*\|")
    with md.open() as fh:
        for line in fh:
            if line.startswith("## "):
                in_typ = "Typical" in line
                continue
            if not in_typ:
                continue
            m = row_re.match(line)
            if m:
                out[m.group(1)] = int(m.group(2))
    return out


_NODE_COUNTS = _build_node_counts()


def run_dcopf(case_name: str) -> dict:
    """Solve DC-OPF for ``case_name`` and return result summary."""
    net = load_pglib_case(case_name)
    dcopf = DCOPF(net)
    dcopf.add_OPF(angle_limits=True)

    # PGLib DC-OPF uses linear costs only (c2 dropped in pure LP).
    # We still accept quadratic terms but warn — IPOPT can handle them as QP.
    add_poly_cost_objective(dcopf, allow_quadratic=True)

    t0 = time.perf_counter()
    res = dcopf.solve(solver=SOLVER, print_solver_output=False)
    elapsed = time.perf_counter() - t0

    ok = res is not None and pyo.check_optimal_termination(res)
    obj = pyo.value(dcopf.model.obj_poly_cost) if ok else float("nan")
    return {
        "obj": obj,
        "time_s": elapsed,
        "ok": ok,
        "termination": str(res.solver.termination_condition),
    }


def run_acopf(case_name: str) -> dict:
    """Solve AC-OPF for ``case_name`` and return result summary.

    Uses PGLib-compatible defaults:
      * ``thermal_limit='mva'`` (constant-MVA branch limit per MATPOWER)
      * ``free_slack_vm=True``  (slack voltage magnitude free in [Vmin, Vmax])
      * ``fix_hv_buses=False``  (no 110 kV pinning)
      * ``angle_limits=True``   (uses ANGMIN/ANGMAX from the .m file)
    """
    net = load_pglib_case(case_name)
    acopf = ACOPF(net)
    acopf.add_OPF(
        thermal_limit="mva",
        free_slack_vm=True,
        fix_hv_buses=False,
        angle_limits=True,
    )
    add_poly_cost_objective(acopf, allow_quadratic=True)

    t0 = time.perf_counter()
    res = acopf.solve(solver=SOLVER, print_solver_output=False)
    elapsed = time.perf_counter() - t0

    ok = res is not None and pyo.check_optimal_termination(res)
    obj = pyo.value(acopf.model.obj_poly_cost) if ok else float("nan")
    return {
        "obj": obj,
        "time_s": elapsed,
        "ok": ok,
        "termination": str(res.solver.termination_condition),
    }


def main():
    cases = _select_cases(CASES, MAX_BUSES)
    if not cases:
        print("No cases selected.")
        return

    print(f"Running {len(cases)} case(s) with solver={SOLVER!r}:")
    for c, n in cases:
        print(f"  - {c}  ({n} buses)")
    print()

    rows = []
    for case, nodes in cases:
        full_name = f"pglib_opf_{case}"
        baseline = PGLIB_BASELINE_TYP[full_name]
        ref_dc, ref_ac = baseline["dc"], baseline["ac"]

        row = {
            "case": case,
            "nodes": nodes,
            "ref_dc": ref_dc,
            "ref_ac": ref_ac,
        }

        if RUN_DC:
            try:
                dc = run_dcopf(case)
                row["dc_obj"] = dc["obj"]
                row["dc_gap_%"] = (dc["obj"] - ref_dc) / ref_dc * 100
                row["dc_t"] = dc["time_s"]
                row["dc_ok"] = dc["ok"]
            except Exception as e:
                row["dc_obj"] = float("nan")
                row["dc_gap_%"] = float("nan")
                row["dc_t"] = float("nan")
                row["dc_ok"] = False
                row["dc_err"] = str(e)[:60]

        if RUN_AC:
            try:
                ac = run_acopf(case)
                row["ac_obj"] = ac["obj"]
                row["ac_gap_%"] = (ac["obj"] - ref_ac) / ref_ac * 100
                row["ac_t"] = ac["time_s"]
                row["ac_ok"] = ac["ok"]
            except Exception as e:
                row["ac_obj"] = float("nan")
                row["ac_gap_%"] = float("nan")
                row["ac_t"] = float("nan")
                row["ac_ok"] = False
                row["ac_err"] = str(e)[:60]

        rows.append(row)

        msg = [f"{case:25s} ({nodes:5d} buses)"]
        if "dc_obj" in row:
            tag = "✓" if row["dc_ok"] else "✗"
            msg.append(
                f"DC: {row['dc_obj']:11.2f} (ref {ref_dc:11.2f}, "
                f"{row['dc_gap_%']:+6.2f}%, {row['dc_t']:5.2f}s) {tag}"
            )
        if "ac_obj" in row:
            tag = "✓" if row["ac_ok"] else "✗"
            msg.append(
                f"AC: {row['ac_obj']:11.2f} (ref {ref_ac:11.2f}, "
                f"{row['ac_gap_%']:+6.2f}%, {row['ac_t']:5.2f}s) {tag}"
            )
        print("  ".join(msg))

    df = pd.DataFrame(rows)

    print()
    print("=" * 78)
    print("Summary (sorted by node count):")
    print(df.to_string(index=False, float_format=lambda v: f"{v:.3f}"))

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_csv = RESULTS_DIR / "pglib_benchmark.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv}")

    out_md = RESULTS_DIR / "pglib_benchmark.md"
    _write_markdown_table(df, out_md, ac=RUN_AC, dc=RUN_DC)
    print(f"Wrote {out_md}")


def _fmt_obj(val: float) -> str:
    """Format an objective value in the BASELINE.md scientific style."""
    if val != val:  # NaN
        return "—"
    return f"{val:.4e}"


def _fmt_time(t: float) -> str:
    if t != t:
        return "—"
    if t < 1:
        return "<1"
    return f"{t:.0f}"


def _fmt_gap(g: float) -> str:
    if g != g:
        return "—"
    return f"{g:+.2f}"


def _write_markdown_table(
    df: pd.DataFrame, path: os.PathLike, ac: bool, dc: bool
) -> None:
    """Write a results table in the same shape as PGLib's BASELINE.md."""
    headers = ["**Case Name**", "**Nodes**"]
    if dc:
        headers += ["**DC ($/h)**", "**DC gap (%)**", "**DC Time (s)**"]
    if ac:
        headers += ["**AC ($/h)**", "**AC gap (%)**", "**AC Time (s)**"]

    lines = ["# potpourri results vs PGLib-OPF baseline", ""]
    lines.append(
        "Solver: IPOPT (potpourri Pyomo model) — PGLib-OPF v23.07 reference "
        "values from upstream `BASELINE.md`."
    )
    lines.append("")
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join("---" for _ in headers) + " |")
    for _, row in df.iterrows():
        cells = [f"pglib_opf_{row['case']}", str(int(row["nodes"]))]
        if dc:
            cells += [
                _fmt_obj(row.get("dc_obj", float("nan"))),
                _fmt_gap(row.get("dc_gap_%", float("nan"))),
                _fmt_time(row.get("dc_t", float("nan"))),
            ]
        if ac:
            cells += [
                _fmt_obj(row.get("ac_obj", float("nan"))),
                _fmt_gap(row.get("ac_gap_%", float("nan"))),
                _fmt_time(row.get("ac_t", float("nan"))),
            ]
        lines.append("| " + " | ".join(cells) + " |")

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
