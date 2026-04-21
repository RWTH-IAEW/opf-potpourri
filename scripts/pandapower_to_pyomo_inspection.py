"""Pandapower → Pyomo model inspection for potpourri.

Demonstrates how each pandapower table maps to Pyomo sets, parameters,
variables, and constraints in the AC power flow model.

Data flow::

  net.bus      → Sets  B, b0, bPV
  net.ext_grid → Sets  G, eG          + Params PG, v_b0, delta_b0
  net.gen      → Sets  G, gG          + Params PG
  net.sgen     → Set   sG             + Params PsG, QsG
  net.load     → Set   D              + Params PD, QD
  net.line     → Set   L              + Params Gii/Bii/Gik/Bik (admittances)
  net.trafo    → Set   TRANSF         + Params GiiT/BiiT/… shift
  net.shunt    → Set   SHUNT          + Params GB, BB

All power quantities are stored in p.u. on net.sn_mva (baseMVA).

Institut für Elektrische Anlagen und Netze, Digitalisierung und
Energiewirtschaft (IAEW)
(c) 2023, Steffen Kortmann
"""

import warnings

import pyomo.environ as pyo
import simbench as sb

from potpourri.models.AC import AC

warnings.filterwarnings("ignore")

NET_NAME = "1-LV-rural1--0-sw"

# ── helpers ──────────────────────────────────────────────────────────────────


def _set_len(s):
    return len(list(s))


def _param_len(p):
    try:
        return len(list(p))
    except TypeError:
        return 1  # scalar Param


def _var_stats(v):
    data = list(v.values())
    total = len(data)
    fixed = sum(1 for vd in data if vd.is_fixed())
    return total, total - fixed, fixed


def _con_count(c):
    return sum(1 for cd in c.values() if cd.active)


def _hr(width=64):
    print("─" * width)


def _section(title):
    print(f"\n{'━' * 64}")
    print(f"  {title}")
    print(f"{'━' * 64}")


def _table(rows, headers, col_widths):
    fmt = "  " + "  ".join(f"{{:<{w}}}" for w in col_widths)
    print(fmt.format(*headers))
    print("  " + "  ".join("─" * w for w in col_widths))
    for row in rows:
        print(fmt.format(*[str(x) for x in row]))


# ── main ─────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # ── 1. Pandapower network ─────────────────────────────────────────────────
    _section("1. PANDAPOWER NETWORK  ·  " + NET_NAME)

    net = sb.get_simbench_net(NET_NAME)

    pp_rows = [
        ("bus", len(net.bus), net.bus.in_service.sum()),
        ("ext_grid", len(net.ext_grid), net.ext_grid.in_service.sum()),
        ("line", len(net.line), net.line.in_service.sum()),
        ("trafo", len(net.trafo), net.trafo.in_service.sum()),
        ("sgen", len(net.sgen), net.sgen.in_service.sum()),
        ("load", len(net.load), net.load.in_service.sum()),
        (
            "shunt",
            len(net.shunt),
            int(net.shunt.in_service.sum()) if len(net.shunt) else 0,
        ),
    ]
    _table(pp_rows, ["table", "rows", "in-service"], [10, 8, 12])
    print(f"\n  baseMVA : {net.sn_mva} MVA")
    extgrid_bus = net.ext_grid.bus.iloc[0]
    print(
        f"  vn_kv   : {net.bus.vn_kv[extgrid_bus]:.3g} kV  (ext_grid bus {extgrid_bus})"
    )

    # ── 2. Build AC Pyomo model ───────────────────────────────────────────────
    _section("2. PYOMO SETS  (how pandapower tables become index sets)")

    ac = AC(net)
    m = ac.model

    set_rows = [
        ("B", _set_len(m.B), "all buses in network"),
        ("b0", _set_len(m.b0), "slack / ext-grid buses  (bus type 3)"),
        ("bPV", _set_len(m.bPV), "voltage-controlled buses (bus type 2)"),
        ("G", _set_len(m.G), "all generators (ext_grid + gen)"),
        ("eG", _set_len(m.eG), "└─ ext-grid / slack generators"),
        ("gG", _set_len(m.gG), "└─ non-slack synchronous generators"),
        ("sG", _set_len(m.sG), "static generators  (net.sgen)"),
        ("D", _set_len(m.D), "loads  (net.load, in-service)"),
        ("L", _set_len(m.L), "lines  (net.line, in-service)"),
        (
            "TRANSF",
            _set_len(m.TRANSF),
            "transformers  (net.trafo, in-service)",
        ),
        ("SHUNT", _set_len(m.SHUNT), "shunts  (net.shunt, in-service)"),
        ("", "", ""),
        ("sGbs", _set_len(m.sGbs), "sgen → bus incidence pairs  (sG × B)"),
        ("Gbs", _set_len(m.Gbs), "gen  → bus incidence pairs  (G  × B)"),
        ("Dbs", _set_len(m.Dbs), "load → bus incidence pairs  (B  × D)"),
    ]
    _table(set_rows, ["set", "size", "meaning"], [8, 6, 44])

    # ── 3. Parameters ─────────────────────────────────────────────────────────
    _section("3. PYOMO PARAMETERS  (read-only data extracted from net.*)")

    param_rows = [
        # topology
        (
            "A",
            _param_len(m.A),
            "L×{1,2}",
            "bus-line incidence (from/to bus index per line)",
        ),
        (
            "AT",
            _param_len(m.AT),
            "T×{1,2}",
            "bus-trafo incidence (HV/LV bus index per trafo)",
        ),
        ("baseMVA", 1, "scalar", f"system base = {ac.baseMVA} MVA"),
        # injections
        ("PsG", _param_len(m.PsG), "sG", "sgen real power  [p.u. on baseMVA]"),
        ("QsG", _param_len(m.QsG), "sG", "sgen reactive power  [p.u.]"),
        ("PG", _param_len(m.PG), "G", "generator real power  [p.u.]"),
        ("PD", _param_len(m.PD), "D", "load real power  [p.u.]"),
        ("QD", _param_len(m.QD), "D", "load reactive power  [p.u.]"),
        # voltage references
        (
            "v_b0",
            _param_len(m.v_b0),
            "b0",
            "slack bus voltage magnitude  [p.u.]",
        ),
        (
            "delta_b0",
            _param_len(m.delta_b0),
            "b0",
            "slack bus voltage angle  [rad]",
        ),
        # line admittances
        (
            "Gii/Bii",
            _param_len(m.Gii),
            "L",
            "line self-conductance / self-susceptance  [p.u.]",
        ),
        (
            "Gik/Bik",
            _param_len(m.Gik),
            "L",
            "line mutual conductance / susceptance  [p.u.]",
        ),
        # transformer admittances
        (
            "GiiT/BiiT",
            _param_len(m.GiiT),
            "TRANSF",
            "trafo self-conductance / self-susceptance  [p.u.]",
        ),
        (
            "GikT/BikT",
            _param_len(m.GikT),
            "TRANSF",
            "trafo mutual conductance / susceptance  [p.u.]",
        ),
        (
            "shift",
            _param_len(m.shift),
            "TRANSF",
            "transformer phase-shift angle  [rad]",
        ),
    ]
    _table(
        param_rows, ["param", "size", "index", "description"], [12, 6, 8, 48]
    )

    # show one concrete row: PsG values to illustrate p.u. conversion
    if _set_len(m.sG) > 0:
        print("\n  PsG sample (sgen real power in p.u.):")
        for g in m.sG:
            p_mw = ac.static_generation_data.p[g] * ac.baseMVA
            print(
                f"    sgen {g}: PsG = {pyo.value(m.PsG[g]):.5f} p.u. "
                f"= {p_mw:.4f} MW  (net.sgen.p_mw × scaling / baseMVA)"
            )

    # ── 4. Variables ──────────────────────────────────────────────────────────
    _section("4. PYOMO VARIABLES  (unknowns solved by IPOPT)")

    var_rows = []
    var_groups = [
        ("─── bus state", None, None),
        ("v", "voltage magnitude [p.u.]  bounds (0, 2)", True),
        ("delta", "voltage angle [rad]  bounds (−π, π)", True),
        ("─── generator output", None, None),
        ("pG", "real power: ext_grid + gen  [p.u.]", True),
        ("qG", "reactive power: ext_grid + gen  [p.u.]", True),
        ("psG", "sgen real power  [p.u.]  → fixed to PsG", True),
        ("qsG", "sgen reactive power [p.u.]  → free*", True),
        ("─── load consumption", None, None),
        ("pD", "real power demand [p.u.]  → fixed to PD", True),
        ("qD", "reactive power demand [p.u.]  → fixed to QD", True),
        ("─── line flows", None, None),
        ("pLfrom", "real power from-end  [p.u.]", True),
        ("pLto", "real power to-end  [p.u.]", True),
        ("qLfrom", "reactive power from-end  [p.u.]", True),
        ("qLto", "reactive power to-end  [p.u.]", True),
        ("─── transformer flows", None, None),
        ("pThv", "real power HV-end  [p.u.]", True),
        ("pTlv", "real power LV-end  [p.u.]", True),
        ("qThv", "reactive power HV-end  [p.u.]", True),
        ("qTlv", "reactive power LV-end  [p.u.]", True),
        ("Tap", "tap ratio  → fixed to pp.runpp value", True),
    ]

    print(
        f"\n  {'variable':<12}  {'total':>5}  {'free':>5}  {'fixed':>5}  description"
    )
    print(f"  {'─' * 12}  {'─' * 5}  {'─' * 5}  {'─' * 5}  {'─' * 44}")

    total_all = total_free = total_fixed = 0
    for entry in var_groups:
        name, desc, is_var = entry
        if not is_var:
            print(f"\n  {name}")
            continue
        obj = m.component(name)
        if obj is None:
            continue
        tot, free, fixed = _var_stats(obj)
        total_all += tot
        total_free += free
        total_fixed += fixed
        print(f"  {name:<12}  {tot:>5}  {free:>5}  {fixed:>5}  {desc}")

    print(
        f"\n  {'TOTAL':<12}  {total_all:>5}  {total_free:>5}  {total_fixed:>5}"
    )

    # ── 5. Constraints ────────────────────────────────────────────────────────
    _section(
        "5. PYOMO CONSTRAINTS  (equality equations from Kirchhoff's laws)"
    )

    con_groups = [
        ("─── KCL: power balance at each bus", None),
        ("KCL_real", "P-balance at every bus  (one per bus in B)"),
        ("KCL_reactive", "Q-balance at every bus  (one per bus in B)"),
        ("─── KVL: line power equations", None),
        ("KVL_real_from", "P from-end = f(v, δ, G, B)  (one per line)"),
        ("KVL_real_to", "P to-end   = f(v, δ, G, B)  (one per line)"),
        ("KVL_reactive_from", "Q from-end = f(v, δ, G, B)  (one per line)"),
        ("KVL_reactive_to", "Q to-end   = f(v, δ, G, B)  (one per line)"),
        ("─── KVL: transformer power equations", None),
        (
            "KVL_real_fromTransf",
            "P HV-end = f(v, δ, tap, shift)  (one per trafo)",
        ),
        (
            "KVL_real_toTransf",
            "P LV-end = f(v, δ, tap, shift)  (one per trafo)",
        ),
        (
            "KVL_reactive_fromTransf",
            "Q HV-end = f(v, δ, tap, shift)  (one per trafo)",
        ),
        (
            "KVL_reactive_toTransf",
            "Q LV-end = f(v, δ, tap, shift)  (one per trafo)",
        ),
        ("─── voltage setpoints", None),
        ("v_bPV_setpoint", "v[b] = v_bPV[b]  for PV buses"),
    ]

    print(f"\n  {'constraint':<28}  {'active':>6}  description")
    print(f"  {'─' * 28}  {'─' * 6}  {'─' * 44}")

    total_constraints = 0
    for entry in con_groups:
        name, desc = entry
        if desc is None:
            print(f"\n  {name}")
            continue
        obj = m.component(name)
        if obj is None:
            count = 0
        else:
            count = _con_count(obj)
        total_constraints += count
        print(f"  {name:<28}  {count:>6}  {desc}")

    print(f"\n  {'TOTAL':<28}  {total_constraints:>6}")

    # ── 6. Size summary ───────────────────────────────────────────────────────
    _section("6. PANDAPOWER → PYOMO : SIZE SUMMARY")

    n_b = len(net.bus)
    n_l = len(net.line)
    n_t = len(net.trafo)
    n_sg = len(net.sgen)
    n_ld = len(net.load)
    n_eg = len(net.ext_grid)

    rows = [
        ("buses", n_b, f"{2 * n_b} vars   (v[b], δ[b] per bus)"),
        ("ext_grids", n_eg, f"{2 * n_eg} vars   (pG[g], qG[g] per ext_grid)"),
        ("lines", n_l, f"{4 * n_l} vars   (pLfrom/to, qLfrom/to per line)"),
        (
            "transformers",
            n_t,
            f"{4 * n_t} vars   (pThv/Tlv, qThv/Tlv per trafo) + {n_t} fixed Tap",
        ),
        ("sgens", n_sg, f"{2 * n_sg} vars   (psG fixed, qsG free*)"),
        ("loads", n_ld, f"{2 * n_ld} vars   (pD fixed, qD fixed)"),
    ]
    _table(rows, ["element", "pp rows", "→ pyomo"], [14, 9, 50])

    dof = total_free - total_constraints
    print(f"""
  Total variables    : {total_all}
  Free variables     : {total_free}
  Active constraints : {total_constraints}
  Net DOF            : {total_free} − {total_constraints} = {dof}

  * qsG (sgen reactive power) is free in the bare AC model.
    For a well-determined power flow (DOF = 0), fix qsG to the
    nominal net.sgen.q_mvar values before calling solve(), or use
    ACOPF.add_OPF() which adds reactive-power bounds as constraints.
""")
