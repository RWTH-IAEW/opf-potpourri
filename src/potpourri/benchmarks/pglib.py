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
    align_tap_sides: bool = True,
) -> pp.pandapowerNet:
    """Load a PGLib-OPF benchmark case as a pandapower network.

    Args:
        case: Either a case name (``"case5_pjm"``, ``"pglib_opf_case5_pjm"``,
            with or without ``.m`` suffix) or a full path to the ``.m`` file.
        f_hz: System frequency (Hz). PGLib cases are MATPOWER cases without
            an explicit frequency; 60 Hz is the conventional default.
        make_controllable: Flag all generators and ext_grids as
            ``controllable=True`` so the OPF optimises them, and likewise the
            sgens that represent MATPOWER generators (those with a cost row
            or active-power limits). sgens that ``from_mpc`` created from
            buses with negative demand stay fixed injections.
        rebalance_initial_dispatch: Scale generator ``p_mw`` so the initial
            power flow run inside ``Basemodel.__init__`` converges. Several
            PGLib cases ship a ``mpc.gen.Pg`` setpoint that is far from the
            load total (it's only meant as a flat-start guess), which makes
            Newton-Raphson diverge.
        attach_angle_limits: Copy MATPOWER ``ANGMIN``/``ANGMAX`` onto
            ``net.line``/``net.trafo`` (needs ``matpowercaseframes``).
        align_tap_sides: Put each transformer tap on the side MATPOWER puts
            it. MATPOWER's ``TAP`` acts on the from bus; ``from_mpc`` always
            encodes it on the high-voltage side, which distorts the
            admittance matrix wherever the from bus is the low-voltage one
            (18 of the 31 transformers in ``case162_ieee_dtc``; the AC-OPF
            was infeasible within ``[Vmin, Vmax]`` before). Needs
            ``matpowercaseframes``.

    Returns:
        A pandapower network with ``net.poly_cost`` populated (one row per
        generation element, polynomial coefficients in ``cp{0,1,2}_eur``).
    """
    path = _resolve_case_path(case)
    net = from_mpc(str(path), f_hz=f_hz, validate_conversion=False)
    net.name = path.stem

    # Reference bus first: the units it re-creates take part in the
    # in-service filter below like every other generator.
    _normalise_reference_bus(net, path)

    # potpourri's OPF._calc_opf_parameters() uses net._gen_order slices
    # (which only count in-service gens) against the full net.gen/sgen
    # DataFrames, so a mix of in-service and out-of-service rows triggers
    # a shape mismatch. Drop out-of-service generation entirely — and carry
    # the cost tables along: ``net.poly_cost.element`` addresses rows of
    # ``net.gen``/``net.sgen`` by index, so dropping rows without renumbering
    # the cost table attaches every later generator's cost curve to the
    # wrong unit and keeps charging the constant terms of units that no
    # longer exist (case200_activ: 11 of 49 gens are out of service; the
    # objective was 14 % above the PGLib reference before this remap).
    for table_name in ("gen", "sgen"):
        table = net[table_name]
        if not table.empty and "in_service" in table.columns:
            in_service = table.in_service.astype(bool)
            if (~in_service).any():
                _drop_generation_rows(net, table_name, in_service)

    if make_controllable:
        if not net.gen.empty:
            net.gen["controllable"] = True
        if not net.sgen.empty:
            # pandapower's converter turns buses with negative demand into
            # sgens (case162_ieee_dtc has nine, case300_ieee eight). Those
            # carry neither active-power limits nor a cost row and are fixed
            # injections in the PGLib formulation, not dispatchable units;
            # leaving them controllable hands the OPF free generation with
            # zero cost and a P bound at the shipped setpoint.
            net.sgen["controllable"] = _dispatchable_sgens(net)
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

    if attach_angle_limits:
        _attach_branch_angle_limits(net, path)

    if align_tap_sides:
        _align_transformer_tap_sides(net, path)

    if rebalance_initial_dispatch:
        _rebalance_initial_dispatch(net)

    return net


def _normalise_reference_bus(net: pp.pandapowerNet, mpc_path: Path) -> int:
    """Guarantee an in-service reference bus without inventing capacity.

    ``from_mpc`` turns the *first* generator at the MATPOWER type-3 bus into
    ``net.ext_grid`` and keeps every further unit there as an ``sgen`` with
    its own cost row, which is a faithful picture. Two situations break it:

    * the slack bus carries no generator at all (the eight RTE cases, where
      MATPOWER and PowerModels use the bus purely as the angle reference):
      no external grid exists and pandapower cannot run a power flow;
    * the first unit is out of service (``case2746wop_k``): the external grid
      is created out of service, so again there is no reference bus, while
      its capacity and cost row belong to a unit that is not running.

    Both are repaired the same way: an in-service external grid with zero
    active and reactive capacity and no cost row, i.e. an angle reference
    only. Everything else stays as ``from_mpc`` made it.

    Returns the number of external grids created or reduced.
    """
    changed = 0
    eg = net.ext_grid
    if not eg.empty:
        off = ~eg["in_service"].astype(bool)
        if off.any():
            for col in ("min_p_mw", "max_p_mw", "min_q_mvar", "max_q_mvar"):
                eg.loc[off, col] = 0.0
            eg.loc[off, "in_service"] = True
            if "poly_cost" in net and not net.poly_cost.empty:
                drop = (net.poly_cost["et"] == "ext_grid") & net.poly_cost[
                    "element"
                ].isin(eg.index[off])
                net.poly_cost = net.poly_cost.loc[~drop].reset_index(drop=True)
            changed += int(off.sum())
        return changed

    try:
        from matpowercaseframes import CaseFrames
    except ImportError:
        return 0
    bus = CaseFrames(str(mpc_path)).bus
    slack_ids = bus.loc[bus["BUS_TYPE"].astype(int) == 3, "BUS_I"].astype(int)
    for bus_id in slack_ids:
        idx = int(bus_id) - 1  # from_mpc numbers pandapower buses as id - 1
        if idx not in net.bus.index:
            continue
        pp.create_ext_grid(
            net,
            bus=idx,
            vm_pu=1.0,
            min_p_mw=0.0,
            max_p_mw=0.0,
            min_q_mvar=0.0,
            max_q_mvar=0.0,
            name="angle reference (MATPOWER slack bus without generator)",
        )
        changed += 1
    return changed


def _drop_generation_rows(
    net: pp.pandapowerNet, table_name: str, keep: pd.Series
) -> None:
    """Drop the rows of ``net[table_name]`` where ``keep`` is False and
    renumber the surviving rows *and* their ``poly_cost``/``pwl_cost``
    entries consistently.

    The cost tables address generators by ``(et, element)``; ``element`` is
    the row index of the generator table. After ``reset_index`` every row
    behind a dropped one moves up, so the cost rows must move with them and
    the cost rows of dropped generators must go.
    """
    table = net[table_name]
    new_index = {old: new for new, old in enumerate(table.index[keep])}
    for cost_table in ("poly_cost", "pwl_cost"):
        if cost_table not in net or net[cost_table].empty:
            continue
        cost = net[cost_table]
        is_table = cost["et"] == table_name
        cost = cost.loc[~is_table | cost["element"].isin(new_index)].copy()
        is_table = cost["et"] == table_name
        cost.loc[is_table, "element"] = cost.loc[is_table, "element"].map(
            new_index
        )
        net[cost_table] = cost.reset_index(drop=True)
    net[table_name] = table.loc[keep].reset_index(drop=True)


def _dispatchable_sgens(net: pp.pandapowerNet) -> np.ndarray:
    """Boolean mask of ``net.sgen`` rows that are generators in the MATPOWER
    sense: they have a cost row or explicit active-power limits.

    ``from_mpc`` creates sgens for two unrelated things: generators sitting
    at PQ buses (with ``min_p_mw``/``max_p_mw`` and a ``gencost`` row) and
    buses with negative demand (neither). Only the former are dispatchable.
    """
    sgen = net.sgen
    costed = np.zeros(len(sgen), dtype=bool)
    if "poly_cost" in net and not net.poly_cost.empty:
        elements = net.poly_cost.loc[net.poly_cost["et"] == "sgen", "element"]
        costed = sgen.index.isin(elements.astype(int))
    has_limits = np.zeros(len(sgen), dtype=bool)
    if "max_p_mw" in sgen.columns:
        has_limits = sgen["max_p_mw"].notna().to_numpy()
    return costed | has_limits


def _attach_branch_angle_limits(net: pp.pandapowerNet, mpc_path: Path) -> None:
    """Read MATPOWER ``ANGMIN``/``ANGMAX`` from the source ``.m`` file and
    write them to ``angmin_degree`` / ``angmax_degree`` on ``net.line``,
    ``net.trafo`` and ``net.impedance``.

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

    # Every MATPOWER branch ends up in one of three pandapower tables: lines
    # (same voltage level, nominal ratio), transformers (off-nominal ratio or
    # phase shift) or impedance elements (different voltage levels without a
    # tap). All three carry the angle limit; leaving the impedance rows out
    # let case60_c__sad settle 3.7 % below the PGLib reference with six
    # angle limits violated.
    tables = [
        ("line", net.line, "from_bus", "to_bus"),
        ("trafo", net.trafo, "hv_bus", "lv_bus"),
        ("impedance", net.get("impedance"), "from_bus", "to_bus"),
    ]
    tables = [
        (n, t, a, b) for n, t, a, b in tables if t is not None and not t.empty
    ]
    ends = {
        n: (t[a].astype(int).values, t[b].astype(int).values)
        for n, t, a, b in tables
    }
    used = {n: np.zeros(len(t), dtype=bool) for n, t, _, _ in tables}
    lim = {
        n: (
            np.full(len(t), -360.0, dtype=float),
            np.full(len(t), 360.0, dtype=float),
        )
        for n, t, _, _ in tables
    }

    for i in range(len(branch)):
        for name, _, _, _ in tables:
            a, b = ends[name]
            hit = np.where(
                ~used[name]
                & (
                    ((a == fr[i]) & (b == to[i]))
                    | ((a == to[i]) & (b == fr[i]))
                )
            )[0]
            if len(hit):
                j = hit[0]
                lim[name][0][j] = amin[i]
                lim[name][1][j] = amax[i]
                used[name][j] = True
                break

    for name, table, _, _ in tables:
        table["angmin_degree"] = lim[name][0]
        table["angmax_degree"] = lim[name][1]


def _align_transformer_tap_sides(net: pp.pandapowerNet, mpc_path: Path) -> int:
    """Move the tap of every transformer whose MATPOWER from bus is the
    pandapower ``lv_bus`` to ``tap_side="lv"``.

    MATPOWER models an off-nominal ratio ``TAP`` as an ideal transformer at
    the from bus (``Yff = ys / TAP²``, ``Ytt = ys``). pandapower's converter
    keeps the ratio but always writes ``tap_side="hv"``; when the from bus is
    the low-voltage side that divides the wrong diagonal of the admittance
    matrix by ``TAP²``. With the tap on the right side the bus admittance
    matrix built from the pandapower network matches MATPOWER's to the
    rounding of the transformer parameters (1e-2 p.u. on case162/case300
    instead of 5–12 p.u.).

    Returns the number of transformers re-encoded. Transformers are matched
    to MATPOWER branches by their bus pair, like the angle limits.
    """
    if net.trafo.empty:
        return 0
    try:
        from matpowercaseframes import CaseFrames
    except ImportError:
        return 0

    branch = CaseFrames(str(mpc_path)).branch
    if "TAP" not in branch.columns:
        return 0
    # MATPOWER bus ids are 1-based; pandapower's from_mpc made them 0-based.
    fr = branch["F_BUS"].astype(int).values - 1
    to = branch["T_BUS"].astype(int).values - 1
    tap = branch["TAP"].astype(float).fillna(0.0).values

    hv = net.trafo["hv_bus"].astype(int).values
    lv = net.trafo["lv_bus"].astype(int).values
    used = np.zeros(len(net.trafo), dtype=bool)
    moved = 0
    for i in range(len(branch)):
        hit = np.where(
            ~used
            & (
                ((hv == fr[i]) & (lv == to[i]))
                | ((hv == to[i]) & (lv == fr[i]))
            )
        )[0]
        if not len(hit):
            continue
        k = hit[0]
        used[k] = True
        if tap[i] in (0.0, 1.0) or hv[k] == fr[i]:
            continue  # nominal ratio, or the tap already sits on the from bus
        net.trafo.iat[k, net.trafo.columns.get_loc("tap_side")] = "lv"
        moved += 1
    return moved


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
        # Only rows with an active-power limit are generators. sgens that
        # ``from_mpc`` made from negative demand have no limit and must keep
        # their setpoint (case240_pserc carries 4.6 GW of such injections;
        # zeroing them raised the objective by 4.9 %).
        dispatchable = table["max_p_mw"].notna()
        if not dispatchable.any():
            continue
        max_p = table.loc[dispatchable, "max_p_mw"].clip(lower=0)
        if "min_p_mw" in table.columns:
            min_p = table.loc[dispatchable, "min_p_mw"].fillna(0)
        else:
            min_p = 0.0
        table.loc[dispatchable, "p_mw"] = (max_p * scale).clip(
            lower=min_p, upper=max_p
        )


def list_available_cases() -> list[str]:
    """Return the sorted list of bare case names available in PGLIB_ROOT."""
    if not PGLIB_ROOT.is_dir():
        return []
    cases = []
    for p in PGLIB_ROOT.glob("pglib_opf_case*.m"):
        cases.append(p.stem.replace("pglib_opf_", ""))
    return sorted(cases)


# DC/AC cells hold a number or "inf." (PowerModels found the problem
# infeasible; 45 of the 66 SAD cases have no DC solution under the tight
# angle limits). Both must parse, or the AC reference of those rows is lost.
_ROW_RE = re.compile(
    r"^\|\s*(pglib_opf_\S+?)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*"
    r"([0-9.eE+\-]+|inf\.?)\s*\|\s*([0-9.eE+\-]+|inf\.?)\s*\|"
)


def _baseline_value(cell: str) -> float:
    """``"inf."`` means PowerModels reported the problem infeasible."""
    return float("inf") if cell.startswith("inf") else float(cell)


def parse_baseline_md(
    baseline_path: str | Path | None = None,
) -> dict[str, dict[str, dict[str, float]]]:
    """Parse PGLib's ``BASELINE.md`` to extract DC/AC reference values.

    Returns a nested dict ``{group: {case_name: {"dc": $/h, "ac": $/h}}}``
    with ``group`` in ``{"TYP", "API", "SAD"}``. A value of ``inf`` records
    that PowerModels.jl found that problem infeasible (``"inf."`` in the
    table; most SAD cases have no DC-OPF solution).

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
                "dc": _baseline_value(dc_val),
                "ac": _baseline_value(ac_val),
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
