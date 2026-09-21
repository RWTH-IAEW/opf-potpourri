# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Benchmark potpourri AC-OPF and DC-OPF against the PGLib-OPF reference.

PGLib-OPF (https://github.com/power-grid-lib/pglib-opf) is the IEEE PES Power
Grid Library benchmark suite for optimal power flow. Each case ships as a
MATPOWER ``.m`` file with a published reference objective value (DC and AC,
solved by PowerModels.jl + IPOPT) in ``BASELINE.md``, for three operating
conditions: Typical (TYP), Congested (API, loads scaled up so thermal limits
bind) and Small Angle Difference (SAD, tight phase-angle limits).

This script:

1. Loads a PGLib case via :func:`potpourri.benchmarks.load_pglib_case`
   (which repairs what pandapower's MATPOWER import gets wrong for OPF use
   and attaches ``ANGMIN``/``ANGMAX`` so phase-angle limits can be enforced).
2. Builds the potpourri DC- and AC-OPF with PGLib-compatible flags:
   * ``thermal_limit='mva'`` — constant-MVA branch limit (matches MATPOWER /
     PowerModels' ``constraint_thermal_limit_*``).
   * ``free_slack_vm=True`` — slack-bus voltage magnitude bounded by
     ``[Vmin, Vmax]`` rather than pinned.
   * ``angle_limits=True`` — branch phase-angle-difference constraints.
3. Wires the polynomial generator cost from ``net.poly_cost`` as the
   objective via :func:`add_poly_cost_objective`.
4. Reports the objective vs. the PGLib reference, per group.

Cases run in parallel worker processes (one IPOPT each); the largest cases
are dispatched first so they overlap with the many small ones. A full run
over all three groups and all sizes takes hours and tens of gigabytes for
the 20 000+ bus cases; set ``MAX_BUSES`` for a quick pass.
"""

from __future__ import annotations

import math
import os
import re
import time
import traceback
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable

import pandas as pd

warnings.filterwarnings("ignore")

# IPOPT's linear solver and the BLAS underneath it are multi-threaded by
# default; one solve then takes about ten cores, and ``N_WORKERS`` solves
# oversubscribe the host several times over, which slows every one of them
# down. Pin them to a single thread each — the solver runs as a child
# process and inherits this environment — so that ``N_WORKERS`` alone
# decides how much of the machine the run uses. Export the variables to
# override.
for _threads in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
):
    os.environ.setdefault(_threads, "1")

# ── Configuration ─────────────────────────────────────────────────────────────
SOLVER = "ipopt"
GROUPS = ("typ", "api", "sad")  # PGLib operating conditions to run
MAX_BUSES = None  # None → every case; an int skips cases with more buses
RUN_DC = True  # include DC-OPF column
RUN_AC = True  # include AC-OPF column
CASES = None  # None → all cases of each group; list of bare names to override
N_WORKERS = 8  # parallel worker processes, one solver thread each
TIME_LIMIT_S = 3600  # IPOPT wall-time limit per solve
RETRY_STARTS = (
    "midrange",
)  # retry starts for a failed AC solve; "dc" also available
DC_CONVENTION = "powermodels"  # DC linearisation convention, see DCOPF
RESULTS_DIR = Path(__file__).parent / "results"
# ──────────────────────────────────────────────────────────────────────────────

_GROUP_SUFFIX = {"typ": "", "api": "__api", "sad": "__sad"}
_GROUP_TITLE = {
    "typ": "Typical Operating Conditions (TYP)",
    "api": "Congested Operating Conditions (API)",
    "sad": "Small Angle Difference Conditions (SAD)",
}


def _baselines():
    from potpourri.benchmarks import (
        PGLIB_BASELINE_API,
        PGLIB_BASELINE_SAD,
        PGLIB_BASELINE_TYP,
    )

    return {
        "typ": PGLIB_BASELINE_TYP,
        "api": PGLIB_BASELINE_API,
        "sad": PGLIB_BASELINE_SAD,
    }


def _build_node_counts() -> dict[str, int]:
    """Read BASELINE.md once to extract node counts per case (all groups)."""
    from potpourri.benchmarks.pglib import PGLIB_ROOT

    out: dict[str, int] = {}
    md = PGLIB_ROOT / "BASELINE.md"
    if not md.is_file():
        return out
    row_re = re.compile(r"^\|\s*(pglib_opf_\S+?)\s*\|\s*(\d+)\s*\|\s*\d+\s*\|")
    with md.open() as fh:
        for line in fh:
            m = row_re.match(line)
            if m:
                out[m.group(1)] = int(m.group(2))
    return out


def _select_cases(
    group: str, requested: Iterable[str] | None, max_buses: int | None
) -> list[tuple[str, int]]:
    """Choose which cases of ``group`` to benchmark.

    Returns ``(full_case_name, n_nodes)`` tuples sorted by descending size so
    the long-running cases enter the pool first.
    """
    baseline = _baselines()[group]
    nodes = _build_node_counts()
    suffix = _GROUP_SUFFIX[group]
    if requested:
        names = []
        for c in requested:
            n = c if c.startswith("pglib_opf_") else f"pglib_opf_{c}"
            if not n.endswith(suffix):
                n = f"{n}{suffix}"
            if n not in baseline:
                print(
                    f"  ! '{n}' not in the {group.upper()} baseline; skipping"
                )
                continue
            names.append(n)
    else:
        names = list(baseline)

    out = []
    for n in names:
        size = nodes.get(n)
        if size is None or (max_buses is not None and size > max_buses):
            continue
        out.append((n, size))
    out.sort(key=lambda t: -t[1])
    return out


def _midrange_start(model) -> None:
    """Put the model on a start that ignores the network's own setpoint.

    Flat voltages and every dispatchable unit at the middle of its range. The
    shipped setpoint of the congested (API) files is far from anything
    feasible, and IPOPT then converges to a locally infeasible point although
    the problem has a solution. Nothing about the problem changes here, only
    the point the solver starts from.
    """
    import pyomo.environ as pyo

    m = model.model
    for b in m.B:
        if hasattr(m, "v") and not m.v[b].fixed:
            m.v[b].set_value(1.0)
        if not m.delta[b].fixed:
            m.delta[b].set_value(0.0)
    for g in m.G:
        if not m.pG[g].fixed:
            lo = pyo.value(m.PGmin[g])
            hi = pyo.value(m.PGmax[g])
            m.pG[g].set_value((lo + hi) / 2)
    for s in getattr(m, "sGc", ()):
        if not m.psG[s].fixed:
            lo = pyo.value(m.sPGmin[s])
            hi = pyo.value(m.sPGmax[s])
            m.psG[s].set_value((lo + hi) / 2)


def _dc_start(model, case_name: str) -> bool:
    """Put the model on the DC-OPF solution of the same case.

    The DC-OPF is a linear program: its own start does not matter, it solves
    in seconds, and its dispatch and bus angles are a much better guess for
    the AC-OPF than a setpoint the case file never meant as one. Returns
    False when the DC-OPF itself does not solve.
    """
    import pyomo.environ as pyo

    from potpourri.benchmarks import load_pglib_case
    from potpourri.models.cost_objective import add_poly_cost_objective
    from potpourri.models.DCOPF import DCOPF

    dc = DCOPF(load_pglib_case(case_name), dc_convention=DC_CONVENTION)
    dc.add_OPF(angle_limits=True)
    add_poly_cost_objective(dc, allow_quadratic=True)
    res = dc.solve(
        solver=SOLVER,
        print_solver_output=False,
        time_limit=TIME_LIMIT_S,
        to_net=False,
    )
    if not pyo.check_optimal_termination(res):
        return False

    m, d = model.model, dc.model
    for b in m.B:
        if hasattr(m, "v") and not m.v[b].fixed:
            m.v[b].set_value(1.0)
        if not m.delta[b].fixed and b in d.delta:
            m.delta[b].set_value(pyo.value(d.delta[b]))
    for g in m.G:
        if not m.pG[g].fixed and g in d.pG:
            m.pG[g].set_value(pyo.value(d.pG[g]))
    for s in m.sG:
        if not m.psG[s].fixed and s in d.psG:
            m.psG[s].set_value(pyo.value(d.psG[s]))
    return True


_STARTS = {"midrange": _midrange_start}


def _solve(
    builder,
    case_name: str,
    opf_kwargs: dict,
    model_kwargs: dict | None = None,
    retries: tuple[str, ...] = (),
) -> dict:
    """Build and solve one model; return objective, status and timings.

    A solve that does not reach optimality is repeated from each start named
    in ``retries`` until one succeeds; ``start`` says which one produced the
    reported result.
    """
    import pyomo.environ as pyo

    from potpourri.benchmarks import load_pglib_case
    from potpourri.models.cost_objective import add_poly_cost_objective

    t0 = time.perf_counter()
    net = load_pglib_case(case_name)
    model = builder(net, **(model_kwargs or {}))
    model.add_OPF(**opf_kwargs)
    add_poly_cost_objective(model, allow_quadratic=True)
    t_build = time.perf_counter() - t0

    def _run():
        return model.solve(
            solver=SOLVER, print_solver_output=False, time_limit=TIME_LIMIT_S
        )

    t0 = time.perf_counter()
    res = _run()
    ok = res is not None and pyo.check_optimal_termination(res)
    start = "setpoint"
    for name in retries:
        if ok:
            break
        if name == "dc":
            if not _dc_start(model, case_name):
                continue
        else:
            _STARTS[name](model)
        res = _run()
        ok = res is not None and pyo.check_optimal_termination(res)
        start = name
    t_solve = time.perf_counter() - t0

    term = str(res.solver.termination_condition) if res is not None else "none"
    return {
        "obj": pyo.value(model.model.obj_poly_cost) if ok else float("nan"),
        "ok": ok,
        "termination": term,
        "start": start,
        "t_build": t_build,
        "t_solve": t_solve,
    }


def run_dcopf(case_name: str) -> dict:
    """Solve DC-OPF for ``case_name`` (full PGLib name) and return a summary."""
    from potpourri.models.DCOPF import DCOPF

    # The PGLib DC references come from PowerModels' DCPPowerModel: branch
    # susceptance -x/(r²+x²) and no transformer phase shift in the DC flow;
    # potpourri's default is MATPOWER's -1/x with the shift. Using the
    # PowerModels convention here makes the DC column a like-for-like
    # comparison (it closed gaps of up to 2.8 % on the API cases and
    # reproduces the SAD infeasibilities).
    # No retries: the DC-OPF is a linear program, its optimum does not depend
    # on the starting point.
    return _solve(
        DCOPF,
        case_name,
        dict(angle_limits=True),
        model_kwargs=dict(dc_convention=DC_CONVENTION),
    )


def run_acopf(case_name: str) -> dict:
    """Solve AC-OPF for ``case_name`` with the PGLib-compatible flags."""
    from potpourri.models.ACOPF_base import ACOPF

    return _solve(
        ACOPF,
        case_name,
        dict(
            thermal_limit="mva",
            free_slack_vm=True,
            fix_hv_buses=False,
            angle_limits=True,
        ),
        retries=RETRY_STARTS,
    )


def _gap(obj: float, ref: float) -> float:
    if not math.isfinite(ref) or ref == 0.0 or not math.isfinite(obj):
        return float("nan")
    return (obj - ref) / ref * 100


def run_case(group: str, full_name: str, nodes: int) -> dict:
    """Worker entry point: one case, DC and/or AC. Never raises."""
    warnings.filterwarnings("ignore")
    ref = _baselines()[group][full_name]
    row = {
        "group": group,
        "case": full_name.replace("pglib_opf_", ""),
        "nodes": nodes,
        "ref_dc": ref["dc"],
        "ref_ac": ref["ac"],
    }
    for tag, runner, flag in (
        ("dc", run_dcopf, RUN_DC),
        ("ac", run_acopf, RUN_AC),
    ):
        if not flag:
            continue
        try:
            r = runner(full_name)
            row[f"{tag}_obj"] = r["obj"]
            row[f"{tag}_ok"] = r["ok"]
            row[f"{tag}_termination"] = r["termination"]
            row[f"{tag}_start"] = r["start"]
            row[f"{tag}_t_build"] = r["t_build"]
            row[f"{tag}_t_solve"] = r["t_solve"]
        except Exception as e:  # noqa: BLE001 — one bad case must not end the run
            row[f"{tag}_obj"] = float("nan")
            row[f"{tag}_ok"] = False
            row[f"{tag}_termination"] = "exception"
            row[f"{tag}_start"] = "none"
            row[f"{tag}_t_build"] = float("nan")
            row[f"{tag}_t_solve"] = float("nan")
            row[f"{tag}_err"] = f"{type(e).__name__}: {str(e)[:120]}"
            traceback.print_exc()
        row[f"{tag}_gap_%"] = _gap(row[f"{tag}_obj"], row[f"ref_{tag}"])
        # PowerModels infeasible and potpourri infeasible: agreement, not failure
        row[f"{tag}_agree"] = bool(
            row[f"{tag}_ok"]
            if math.isfinite(row[f"ref_{tag}"])
            else not row[f"{tag}_ok"]
        )
    return row


def _describe(row: dict, tag: str) -> str:
    ref = row[f"ref_{tag}"]
    if f"{tag}_obj" not in row:
        return ""
    mark = "✓" if row[f"{tag}_agree"] else "✗"
    if not math.isfinite(ref):
        return f"{tag.upper()}: ref inf, {row[f'{tag}_termination']} {mark}"
    return (
        f"{tag.upper()}: {row[f'{tag}_obj']:12.2f} (ref {ref:12.2f}, "
        f"{row[f'{tag}_gap_%']:+7.2f}%, build {row[f'{tag}_t_build']:6.1f}s "
        f"solve {row[f'{tag}_t_solve']:7.1f}s) {mark}"
    )


def main():
    jobs = []
    for group in GROUPS:
        for full_name, nodes in _select_cases(group, CASES, MAX_BUSES):
            jobs.append((group, full_name, nodes))
    if not jobs:
        print("No cases selected.")
        return
    print(
        f"Running {len(jobs)} case(s) over groups {GROUPS} with solver={SOLVER!r}, "
        f"{N_WORKERS} workers, IPOPT time limit {TIME_LIMIT_S} s",
        flush=True,
    )

    rows = []
    t_start = time.perf_counter()
    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        futures = {pool.submit(run_case, *job): job for job in jobs}
        for done, fut in enumerate(as_completed(futures), start=1):
            group, full_name, nodes = futures[fut]
            row = fut.result()
            rows.append(row)
            print(
                f"[{done:3d}/{len(jobs)} {time.perf_counter() - t_start:7.0f}s] "
                f"{group.upper()} {row['case']:28s} ({nodes:6d} buses)  "
                f"{_describe(row, 'dc')}  {_describe(row, 'ac')}",
                flush=True,
            )

    df = pd.DataFrame(rows).sort_values(["group", "nodes", "case"])
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    out_csv = RESULTS_DIR / "pglib_benchmark.csv"
    df.to_csv(out_csv, index=False)
    print(f"\nWrote {out_csv}")
    for group in GROUPS:
        part = df[df.group == group]
        if part.empty:
            continue
        out_md = RESULTS_DIR / f"pglib_benchmark_{group}.md"
        _write_markdown_table(part, out_md, group, ac=RUN_AC, dc=RUN_DC)
        print(f"Wrote {out_md}")
    _print_summary(df)


def _print_summary(df: pd.DataFrame) -> None:
    print("\nSummary per group:")
    for group, part in df.groupby("group"):
        line = [f"  {group.upper():4s} {len(part):3d} cases"]
        for tag in ("dc", "ac"):
            if f"{tag}_obj" not in part:
                continue
            finite_ref = part[f"ref_{tag}"].apply(math.isfinite)
            solved = part[f"{tag}_ok"] & finite_ref
            gaps = part.loc[solved, f"{tag}_gap_%"].abs()
            line.append(
                f"{tag.upper()}: agree {int(part[f'{tag}_agree'].sum())}, "
                f"solved {int(solved.sum())}/{int(finite_ref.sum())}, "
                f"|gap| median {gaps.median():.3f} % max {gaps.max():.2f} %"
            )
        print(" | ".join(line))


def _fmt_obj(val: float) -> str:
    if val != val:  # NaN
        return "—"
    return f"{val:.4e}"


def _fmt_ref(val: float) -> str:
    return "inf." if math.isinf(val) else f"{val:.4e}"


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
    df: pd.DataFrame, path: os.PathLike, group: str, ac: bool, dc: bool
) -> None:
    """Write a results table in the same shape as PGLib's BASELINE.md."""
    headers = ["**Case Name**", "**Nodes**"]
    for tag, flag in (("DC", dc), ("AC", ac)):
        if flag:
            headers += [
                f"**{tag} ($/h)**",
                f"**{tag} ref ($/h)**",
                f"**{tag} gap (%)**",
                f"**{tag} status**",
                f"**{tag} build (s)**",
                f"**{tag} solve (s)**",
            ]

    lines = [
        f"# potpourri results vs PGLib-OPF baseline — {_GROUP_TITLE[group]}",
        "",
    ]
    lines.append(
        f"Solver: {SOLVER} (potpourri Pyomo model) — PGLib-OPF v23.07 reference "
        "values from upstream `BASELINE.md`; `inf.` marks a reference problem "
        "PowerModels.jl found infeasible."
    )
    lines.append("")
    lines.append("| " + " | ".join(headers) + " |")
    lines.append("| " + " | ".join("---" for _ in headers) + " |")
    for _, row in df.iterrows():
        cells = [f"pglib_opf_{row['case']}", str(int(row["nodes"]))]
        for tag, flag in (("dc", dc), ("ac", ac)):
            if not flag:
                continue
            cells += [
                _fmt_obj(row.get(f"{tag}_obj", float("nan"))),
                _fmt_ref(row[f"ref_{tag}"]),
                _fmt_gap(row.get(f"{tag}_gap_%", float("nan"))),
                str(row.get(f"{tag}_termination", "")),
                _fmt_time(row.get(f"{tag}_t_build", float("nan"))),
                _fmt_time(row.get(f"{tag}_t_solve", float("nan"))),
            ]
        lines.append("| " + " | ".join(cells) + " |")

    with open(path, "w") as fh:
        fh.write("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
