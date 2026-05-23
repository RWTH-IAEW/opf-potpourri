"""PGLib-OPF benchmark loader.

PGLib-OPF (https://github.com/power-grid-lib/pglib-opf) is the IEEE PES Power
Grid Library benchmark suite for optimal power flow. Each case is a MATPOWER
``.m`` file. This module provides:

* ``load_pglib_case`` -- parse a PGLib ``.m`` file into a pandapower network
  configured for OPF (controllable generators, generation/voltage limits,
  polynomial cost coefficients in ``net.poly_cost``).
* ``parse_baseline_md`` -- read the reference DC- and AC-OPF objective values
  published in ``BASELINE.md``.
* ``PGLIB_BASELINE_{TYP,API,SAD}`` -- curated dicts of the published reference
  objective values for the three benchmark groups (Typical Operations,
  Congested Operations, Small Angle Difference) for quick comparison.

The MATPOWER → pandapower conversion is delegated to pandapower's built-in
``from_mpc`` (which uses ``matpowercaseframes`` to parse the ``.m`` source).
On top of that we (a) ensure generators are flagged ``controllable`` so the
OPF actually optimises them, and (b) keep the per-MW polynomial cost
coefficients in ``net.poly_cost`` so a cost objective can be wired up.
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pandapower as pp
from pandapower.converter.matpower import from_mpc

# Repository root -> benchmarks/pglib-opf/
PGLIB_ROOT = Path(__file__).resolve().parents[3] / "benchmarks" / "pglib-opf"


def _resolve_case_path(case: str | Path) -> Path:
    """Resolve a case name or path to an absolute ``.m`` file path."""
    p = Path(case)
    if p.is_file():
        return p

    name = str(case)
    if not name.startswith("pglib_opf_"):
        name = f"pglib_opf_{name}"
    if not name.endswith(".m"):
        name = f"{name}.m"

    candidate = PGLIB_ROOT / name
    if candidate.is_file():
        return candidate

    api_candidate = PGLIB_ROOT / "api" / name
    if api_candidate.is_file():
        return api_candidate
    sad_candidate = PGLIB_ROOT / "sad" / name
    if sad_candidate.is_file():
        return sad_candidate

    raise FileNotFoundError(
        f"PGLib case '{case}' not found under {PGLIB_ROOT}. "
        "Check the case name or pass an absolute path."
    )


def load_pglib_case(
    case: str | Path,
    f_hz: int = 60,
    make_controllable: bool = True,
    rebalance_initial_dispatch: bool = True,
    attach_angle_limits: bool = True,
) -> pp.pandapowerNet:
    """Load a PGLib-OPF benchmark case as a pandapower network.

    Args:
        case: Either a case name (``"case5_pjm"``, ``"pglib_opf_case5_pjm"``,
            with or without ``.m`` suffix) or a full path to the ``.m`` file.
        f_hz: System frequency (Hz). PGLib cases are MATPOWER cases without
            an explicit frequency; 60 Hz is the conventional default.
        make_controllable: Flag all generators/sgens/ext_grids as
            ``controllable=True`` so the OPF optimises them.
        rebalance_initial_dispatch: Scale generator ``p_mw`` so the initial
            power flow run inside ``Basemodel.__init__`` converges. Several
            PGLib cases ship a ``mpc.gen.Pg`` setpoint that is far from the
            load total (it's only meant as a flat-start guess), which makes
            Newton-Raphson diverge.

    Returns:
        A pandapower network with ``net.poly_cost`` populated (one row per
        generation element, polynomial coefficients in ``cp{0,1,2}_eur``).
    """
    path = _resolve_case_path(case)
    net = from_mpc(str(path), f_hz=f_hz, validate_conversion=False)
    net.name = path.stem

    if make_controllable:
        if not net.gen.empty:
            net.gen["controllable"] = True
        if not net.sgen.empty:
            net.sgen["controllable"] = True
        if not net.ext_grid.empty:
            net.ext_grid["controllable"] = True

    # Some PGLib cases ship lines without an active-power limit but with a
    # current rating of 0. pandapower's converter then sets max_i_ka=0, which
    # makes the OPF infeasible. Replace 0 with a large value (essentially
    # unrated).
    if "max_i_ka" in net.line.columns:
        unlimited = net.line["max_i_ka"].fillna(0) <= 0
        if unlimited.any():
            net.line.loc[unlimited, "max_i_ka"] = 1e6

    # potpourri's OPF._calc_opf_parameters() uses net._gen_order slices
    # (which only count in-service gens) against the full net.gen/sgen
    # DataFrames, so a mix of in-service and out-of-service rows triggers
    # a shape mismatch. Drop out-of-service generation entirely.
    for table_name in ("gen", "sgen"):
        table = net[table_name]
        if not table.empty and "in_service" in table.columns:
            in_service = table.in_service.astype(bool)
            if (~in_service).any():
                net[table_name] = table.loc[in_service].reset_index(drop=True)

    if attach_angle_limits:
        _attach_branch_angle_limits(net, path)

    if rebalance_initial_dispatch:
        _rebalance_initial_dispatch(net)

    return net


def _attach_branch_angle_limits(net: pp.pandapowerNet, mpc_path: Path) -> None:
    """Read MATPOWER ``ANGMIN``/``ANGMAX`` from the source ``.m`` file and
    write them to ``net.line.angmin_degree`` / ``net.line.angmax_degree``
    (and the transformer equivalent).

    pandapower's ``from_mpc`` doesn't preserve these, but they are needed
    for PGLib-OPF benchmark fidelity. We match each MATPOWER branch to a
    pandapower line or trafo by (from_bus, to_bus), since pandapower's
    line/trafo split is not purely a function of the TAP column.
    """
    try:
        from matpowercaseframes import CaseFrames
    except ImportError:
        return

    cf = CaseFrames(str(mpc_path))
    branch = cf.branch
    if "ANGMIN" not in branch.columns or "ANGMAX" not in branch.columns:
        return

    # MATPOWER bus ids are 1-based; pandapower's from_mpc made them 0-based.
    fr = branch["F_BUS"].astype(int).values - 1
    to = branch["T_BUS"].astype(int).values - 1
    amin = branch["ANGMIN"].astype(float).values
    amax = branch["ANGMAX"].astype(float).values

    line_amin = np.full(len(net.line), -360.0, dtype=float)
    line_amax = np.full(len(net.line), 360.0, dtype=float)
    trafo_amin = np.full(len(net.trafo), -360.0, dtype=float)
    trafo_amax = np.full(len(net.trafo), 360.0, dtype=float)

    line_from = net.line["from_bus"].astype(int).values
    line_to = net.line["to_bus"].astype(int).values
    line_used = np.zeros(len(net.line), dtype=bool)
    if not net.trafo.empty:
        trafo_hv = net.trafo["hv_bus"].astype(int).values
        trafo_lv = net.trafo["lv_bus"].astype(int).values
        trafo_used = np.zeros(len(net.trafo), dtype=bool)

    for i in range(len(branch)):
        # Try lines first (lookups against both directions)
        line_hit = np.where(
            ~line_used
            & (
                ((line_from == fr[i]) & (line_to == to[i]))
                | ((line_from == to[i]) & (line_to == fr[i]))
            )
        )[0]
        if len(line_hit):
            j = line_hit[0]
            line_amin[j] = amin[i]
            line_amax[j] = amax[i]
            line_used[j] = True
            continue
        if not net.trafo.empty:
            t_hit = np.where(
                ~trafo_used
                & (
                    ((trafo_hv == fr[i]) & (trafo_lv == to[i]))
                    | ((trafo_hv == to[i]) & (trafo_lv == fr[i]))
                )
            )[0]
            if len(t_hit):
                k = t_hit[0]
                trafo_amin[k] = amin[i]
                trafo_amax[k] = amax[i]
                trafo_used[k] = True

    net.line["angmin_degree"] = line_amin
    net.line["angmax_degree"] = line_amax
    if not net.trafo.empty:
        net.trafo["angmin_degree"] = trafo_amin
        net.trafo["angmax_degree"] = trafo_amax


def _rebalance_initial_dispatch(net: pp.pandapowerNet) -> None:
    """Set a balanced initial dispatch so a flat-start power flow converges.

    For each in-service controllable gen/sgen, set ``p_mw`` to ``max_p_mw``
    scaled by ``total_load / total_max_p``. Gens with ``max_p_mw == 0`` (e.g.
    PGLib synchronous condensers) keep ``p_mw = 0``. Tables without
    ``max_p_mw`` / ``min_p_mw`` columns are skipped — they're unbounded so
    their current ``p_mw`` is already valid.
    """
    total_load = float(net.load.p_mw.sum()) if not net.load.empty else 0.0
    if total_load <= 0:
        return

    def _max(table) -> float:
        if table.empty or "max_p_mw" not in table.columns:
            return 0.0
        return float(table["max_p_mw"].fillna(0).clip(lower=0).sum())

    total_max = _max(net.gen) + _max(net.sgen)
    if total_max <= 0:
        return

    # Aim a hair below 1.0 so the slack absorbs losses rather than spilling.
    scale = min(1.0, 0.9 * total_load / total_max)

    for table in (net.gen, net.sgen):
        if table.empty or "max_p_mw" not in table.columns:
            continue
        max_p = table["max_p_mw"].fillna(0).clip(lower=0)
        if "min_p_mw" in table.columns:
            min_p = table["min_p_mw"].fillna(0)
        else:
            min_p = 0.0
        table["p_mw"] = (max_p * scale).clip(lower=min_p, upper=max_p)


def list_available_cases() -> list[str]:
    """Return the sorted list of bare case names available in PGLIB_ROOT."""
    if not PGLIB_ROOT.is_dir():
        return []
    cases = []
    for p in PGLIB_ROOT.glob("pglib_opf_case*.m"):
        cases.append(p.stem.replace("pglib_opf_", ""))
    return sorted(cases)


_ROW_RE = re.compile(
    r"^\|\s*(pglib_opf_\S+?)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*"
    r"([0-9.eE+\-]+)\s*\|\s*([0-9.eE+\-]+)\s*\|"
)


def parse_baseline_md(
    baseline_path: str | Path | None = None,
) -> dict[str, dict[str, dict[str, float]]]:
    """Parse PGLib's ``BASELINE.md`` to extract DC/AC reference values.

    Returns a nested dict ``{group: {case_name: {"dc": $/h, "ac": $/h}}}``
    with ``group`` in ``{"TYP", "API", "SAD"}``.

    Args:
        baseline_path: Path to ``BASELINE.md`` (defaults to the one inside
            ``PGLIB_ROOT``).
    """
    if baseline_path is None:
        baseline_path = PGLIB_ROOT / "BASELINE.md"
    baseline_path = Path(baseline_path)

    out: dict[str, dict[str, dict[str, float]]] = {
        "TYP": {},
        "API": {},
        "SAD": {},
    }
    current = None
    with baseline_path.open() as fh:
        for line in fh:
            stripped = line.strip()
            if stripped.startswith("## "):
                if "Typical" in stripped:
                    current = "TYP"
                elif "Congested" in stripped or "(API)" in stripped:
                    current = "API"
                elif "Small Angle" in stripped or "(SAD)" in stripped:
                    current = "SAD"
                else:
                    current = None
                continue

            if current is None:
                continue
            m = _ROW_RE.match(line)
            if not m:
                continue
            case_name, _, _, dc_val, ac_val = m.groups()
            out[current][case_name] = {
                "dc": float(dc_val),
                "ac": float(ac_val),
            }
    return out


def _parse_or_empty() -> dict[str, dict[str, dict[str, float]]]:
    try:
        return parse_baseline_md()
    except FileNotFoundError:
        return {"TYP": {}, "API": {}, "SAD": {}}


_BASELINE = _parse_or_empty()
PGLIB_BASELINE_TYP: dict[str, dict[str, float]] = _BASELINE["TYP"]
PGLIB_BASELINE_API: dict[str, dict[str, float]] = _BASELINE["API"]
PGLIB_BASELINE_SAD: dict[str, dict[str, float]] = _BASELINE["SAD"]


def baseline_as_dataframe() -> pd.DataFrame:
    """Return the full baseline table as a tidy DataFrame."""
    rows = []
    for group, cases in _BASELINE.items():
        for case_name, vals in cases.items():
            rows.append(
                {
                    "group": group,
                    "case": case_name,
                    "dc_ref": vals["dc"],
                    "ac_ref": vals["ac"],
                }
            )
    return pd.DataFrame(rows)
