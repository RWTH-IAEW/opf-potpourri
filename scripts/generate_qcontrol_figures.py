"""Generate the reactive-power-control figures used in the user guide.

Every curve is computed from
:mod:`potpourri.technologies.q_control` rather than transcribed, so the figures
cannot drift from the implementation: change a grid code and the figures change
with it.

One figure per constraint group, each showing the operating area that group
permits and — just as importantly — what it leaves unbounded.

Rendered with rwthplots (RWTH corporate style, LaTeX typesetting, IEEE column
width).  Install with ``pip install rwthplots``; a LaTeX installation is
required for ``usetex``.

Colour choices follow the encoding job rather than taste:

* A single constraint is one hue (RWTH blue) — region fill plus boundary.
* The three ``var_q`` variants are *ordinal*, not categorical: they are ordered
  by envelope width, so they use one hue at three tints (RWTH blue 100/75/50),
  which validates as a monotone ramp with a light end clearing 2:1 on white.
* Where two different constraints appear together, RWTH blue and red are used;
  that pair measures ΔE 33.7 normal-vision and 20.1 under simulated CVD.
* Line style varies alongside colour throughout, so identity never rests on
  hue alone.
"""

import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import rwthplots  # noqa: E402

from potpourri.technologies.q_control import (  # noqa: E402
    VDE_AR_N_4120,
    compute_q_curves,
)

# --- configuration -------------------------------------------------------
OUT_DIR = (
    pathlib.Path(__file__).resolve().parents[1]
    / "docs"
    / "assets"
    / "q-control"
)
STYLES = ("rwth-latex", "color.standard", "size.ieee-column")
FORMATS = ("svg",)

GRID_CODE = VDE_AR_N_4120
DEADBAND = (0.98, 1.02)  # operator parameterisation, not a normative value
COS_PHI_MIN = 0.90  # cos(phi) cone / cos(phi)(P) target
COS_PHI_FIXED = 0.95  # fixed cos(phi) mode
S_INV_PU = 1.10  # inverter rating over Pn for the S^2 circle

BLUE, BLUE_75, BLUE_50 = "#00549F", "#407FB7", "#8EBAE5"
RED = "#CC071E"
FILL_ALPHA = 0.18
GRID_KW = dict(color="#E1E0D9", linewidth=0.5, zorder=0)

# Keep the SVG output reproducible: matplotlib otherwise stamps a
# <dc:date> and derives element ids from a random salt, so regenerating
# unchanged figures would produce a diff every time.
matplotlib.rcParams["svg.hashsalt"] = "potpourri-q-control"
# -------------------------------------------------------------------------

CURVES = compute_q_curves(GRID_CODE)
TAN_MIN = float(np.tan(np.arccos(COS_PHI_MIN)))
TAN_FIX = float(np.tan(np.arccos(COS_PHI_FIXED)))


def _axes(ax, xlabel, ylabel, title):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, fontsize=8)
    ax.grid(True, **GRID_KW)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    # rwth-latex draws ticks on all four sides; with the top and right spines
    # hidden those ticks float free, so switch them off too.
    ax.tick_params(top=False, right=False, which="both")


def _save(fig, name):
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for fmt in FORMATS:
        path = OUT_DIR / f"{name}.{fmt}"
        fig.savefig(
            path,
            format=fmt,
            bbox_inches="tight",
            pad_inches=0.02,
            metadata={"Date": None} if fmt == "svg" else None,
        )
        print(f"  wrote {path.relative_to(OUT_DIR.parents[2])}", flush=True)
    plt.close(fig)


def fig_qp():
    """Q(P) capability area for the three var_q variants.

    Bounds are drawn as line pairs rather than filled bands: three overlapping
    translucent fills composite into darker tones where they intersect, which
    destroys the ordinal reading the tints are there to carry.

    The bound ramps between the two active-power breakpoints and then holds.
    Before 0.4.1 the shelf was missing and the line ran on to +3.5 Pn at
    rated output, seven times the limit the standard sets.
    """
    p = np.linspace(0.0, 1.0, 400)
    fig, ax = plt.subplots()
    p_ref = GRID_CODE.qp_p_low
    for v, (c, ls) in enumerate(
        zip((BLUE, BLUE_75, BLUE_50), ("-", "--", ":"))
    ):
        lo, hi = GRID_CODE.pq_area.q_flexibility(p, v)
        ax.plot(p, hi, ls, color=c, label=rf"$\mathrm{{var\_q}}={v}$")
        ax.plot(p, lo, ls, color=c)
        q_lo, q_hi = GRID_CODE.pq_area.q_flexibility(p_ref, v)
        ax.plot([p_ref, p_ref], [q_lo, q_hi], "o", ms=2.5, color=c)

    ax.axvline(p_ref, color="#898781", lw=0.6, ls=(0, (1, 2)))
    ax.axhline(0.0, color="#C3C2B7", lw=0.6)
    ax.annotate(
        r"the bound holds beyond $0.2\,P_n$",
        xy=(0.34, 0.56),
        fontsize=5.5,
        color="#52514E",
    )
    ax.set_xlim(0, 1)
    ax.set_ylim(-0.52, 0.62)
    _axes(ax, r"$P/P_n$", r"$Q/P_n$", r"Q(P) capability area")
    ax.legend(loc="lower left", fontsize=6, frameon=False)
    _save(fig, "qp-characteristic")


def fig_qu():
    """Q(U) capability area, over the grid code's own voltage breakpoints.

    The area is a hexagon: the band is pinned to one limit below V1, opens
    out across V1-V2, spans the full range over the plateau, closes again
    across V3-V4 and is pinned above V4.  Before 0.4.1 both bounds were
    single unclipped lines, so the band was 3.4x too wide at every voltage
    and never saturated.
    """
    pts = np.sort(GRID_CODE.qv_area.x_points)
    v_pu = np.linspace(pts[0] - 0.02, pts[-1] + 0.02, 600)
    fig, ax = plt.subplots()
    for var, (c, ls) in enumerate(
        zip((BLUE, BLUE_75, BLUE_50), ("-", "--", ":"))
    ):
        lo, hi = GRID_CODE.qv_area.q_flexibility(v_pu, var)
        ax.plot(v_pu, hi, ls, color=c, label=rf"$\mathrm{{var\_q}}={var}$")
        ax.plot(v_pu, lo, ls, color=c)

    for i, x in enumerate(pts, start=1):
        ax.axvline(x, color="#898781", lw=0.5, ls=(0, (1, 2)))
        ax.annotate(
            rf"$V_{i}$",
            xy=(x, 0.60),
            fontsize=5.5,
            color="#52514E",
            ha="center",
        )
    ax.axhline(0.0, color="#C3C2B7", lw=0.6)
    ax.set_xlim(v_pu[0], v_pu[-1])
    ax.set_ylim(-0.52, 0.68)
    _axes(ax, r"$v$ [p.u.]", r"$Q/P_n$", r"Q(U) capability area")
    ax.legend(loc="lower left", fontsize=6, frameon=False)
    _save(fig, "qu-droop")


def fig_qu_deadband():
    """Q(U) characteristic with a dead band, against the area it sits in.

    The area *bounds* Q and leaves the optimiser free inside it; the
    characteristic *assigns* Q, so a dead band -- a voltage span over which
    Q is held at zero -- becomes expressible.  That pinch makes the feasible
    set non-convex, which is why this mode needs a MIP-capable solver.
    """
    curve = GRID_CODE.deadband_curve(deadband=DEADBAND)
    pts = np.sort(GRID_CODE.qv_area.x_points)
    v_pu = np.linspace(pts[0] - 0.02, pts[-1] + 0.02, 600)
    fig, ax = plt.subplots()

    lo, hi = GRID_CODE.qv_area.q_flexibility(v_pu, 0)
    ax.fill_between(v_pu, lo, hi, color=BLUE, alpha=FILL_ALPHA, linewidth=0)
    ax.plot(v_pu, hi, "-", color=BLUE, lw=0.8, label=r"capability area")
    ax.plot(v_pu, lo, "-", color=BLUE, lw=0.8)
    ax.plot(
        v_pu,
        curve.step(v_pu, 0),
        "-",
        color=RED,
        lw=1.2,
        label=r"characteristic",
    )

    db_lo, db_hi = curve.deadband(0)
    ax.axvspan(db_lo, db_hi, color="#898781", alpha=0.12, linewidth=0)
    ax.annotate(
        r"dead band",
        xy=((db_lo + db_hi) / 2, 0.13),
        fontsize=5.5,
        color="#52514E",
        ha="center",
    )
    ax.axhline(0.0, color="#C3C2B7", lw=0.6)
    ax.set_xlim(v_pu[0], v_pu[-1])
    ax.set_ylim(-0.52, 0.68)
    _axes(
        ax,
        r"$v$ [p.u.]",
        r"$Q/P_n$",
        r"Q(U) characteristic with dead band",
    )
    ax.legend(loc="lower left", fontsize=6, frameon=False)
    _save(fig, "qu-deadband")


def fig_s2():
    """Inverter apparent-power circle."""
    fig, ax = plt.subplots()
    th = np.linspace(-np.pi / 2, np.pi / 2, 400)
    ax.fill(
        S_INV_PU * np.cos(th),
        S_INV_PU * np.sin(th),
        color=BLUE,
        alpha=FILL_ALPHA,
        linewidth=0,
    )
    ax.plot(S_INV_PU * np.cos(th), S_INV_PU * np.sin(th), "-", color=BLUE)
    ax.annotate(
        rf"$S_{{\mathrm{{inv}}}}={S_INV_PU:.2f}\,P_n$",
        xy=(0.55, 0.95),
        fontsize=6,
        color=BLUE,
    )
    ax.set_xlim(0, 1.3)
    ax.set_ylim(-1.3, 1.3)
    _axes(ax, r"$P/P_n$", r"$Q/P_n$", r"Inverter $S^2$ circle")
    _save(fig, "inverter-s2-circle")


def fig_cone():
    """cos(phi) cone."""
    p = np.linspace(0.0, 1.2, 200)
    fig, ax = plt.subplots()
    ax.fill_between(
        p,
        -TAN_MIN * p,
        TAN_MIN * p,
        color=BLUE,
        alpha=FILL_ALPHA,
        linewidth=0,
    )
    ax.plot(p, TAN_MIN * p, "-", color=BLUE)
    ax.plot(p, -TAN_MIN * p, "-", color=BLUE)
    ax.annotate(
        rf"$|Q|\leq P\tan(\arccos {COS_PHI_MIN:.2f})$",
        xy=(0.12, 0.92),
        fontsize=6,
        color=BLUE,
    )
    ax.set_xlim(0, 1.2)
    ax.set_ylim(-1.3, 1.3)
    _axes(ax, r"$P/P_n$", r"$Q/P_n$", r"$\cos\varphi$ cone")
    _save(fig, "cos-phi-cone")


def fig_operating_region():
    """The PV operating region: P >= 0, S^2 circle and cone together."""
    fig, ax = plt.subplots()
    th = np.linspace(-np.pi / 2, np.pi / 2, 600)
    cx, cy = S_INV_PU * np.cos(th), S_INV_PU * np.sin(th)

    p = np.linspace(0, S_INV_PU, 600)
    cone_hi, cone_lo = TAN_MIN * p, -TAN_MIN * p
    circ = np.sqrt(np.maximum(S_INV_PU**2 - p**2, 0.0))
    hi, lo = np.minimum(cone_hi, circ), np.maximum(cone_lo, -circ)
    feasible = hi >= lo
    ax.fill_between(
        p[feasible],
        lo[feasible],
        hi[feasible],
        color=BLUE,
        alpha=FILL_ALPHA,
        linewidth=0,
    )
    ax.plot(cx, cy, "--", color=RED, label=r"$S^2$ circle")
    ax.plot(p, cone_hi, "-", color=BLUE, label=r"$\cos\varphi$ cone")
    ax.plot(p, cone_lo, "-", color=BLUE)
    p_cross = S_INV_PU * COS_PHI_MIN
    ax.plot([p_cross], [TAN_MIN * p_cross], "o", ms=3, color="#000000")
    ax.annotate(
        r"$P_{\mathrm{cross}}=S_{\mathrm{inv}}\cos\varphi$",
        xy=(p_cross, TAN_MIN * p_cross),
        xytext=(0.20, 0.85),
        fontsize=6,
        color="#0B0B0B",
        arrowprops=dict(arrowstyle="-", lw=0.5, color="#898781"),
    )
    ax.set_xlim(0, 1.3)
    ax.set_ylim(-1.3, 1.3)
    _axes(ax, r"$P/P_n$", r"$Q/P_n$", r"PV operating region")
    ax.legend(loc="lower right", fontsize=6, frameon=False)
    _save(fig, "pv-operating-region")


def fig_pu_curtail():
    """P(U) active-power curtailment."""
    vc, vmax = GRID_CODE.vpu_v_curtail, GRID_CODE.vpu_v_max
    v_pu = np.linspace(1.00, 1.12, 400)
    p_lim = np.clip((vmax - v_pu) / (vmax - vc), 0.0, 1.0)
    fig, ax = plt.subplots()
    ax.fill_between(v_pu, 0, p_lim, color=BLUE, alpha=FILL_ALPHA, linewidth=0)
    ax.plot(v_pu, p_lim, "-", color=BLUE)
    for x, lab in ((vc, r"$V_{\mathrm{curtail}}$"), (vmax, r"$V_{\max}$")):
        ax.axvline(x, color="#898781", lw=0.6, ls=(0, (1, 2)))
        ax.annotate(
            lab, xy=(x, 1.03), fontsize=6, color="#52514E", ha="center"
        )
    ax.set_xlim(1.00, 1.12)
    ax.set_ylim(0, 1.12)
    _axes(ax, r"$v$ [p.u.]", r"$P/P_n$", r"P(U) curtailment")
    _save(fig, "pu-curtailment")


def fig_fixed_cos_phi():
    """Fixed cos(phi): an equality, so the operating area collapses."""
    p = np.linspace(0.0, 1.0, 200)
    fig, ax = plt.subplots()
    ax.fill_between(
        p,
        -TAN_MIN * p,
        TAN_MIN * p,
        color=BLUE,
        alpha=0.08,
        linewidth=0,
    )
    ax.plot(p, TAN_MIN * p, ":", color=BLUE, lw=0.8)
    ax.plot(p, -TAN_MIN * p, ":", color=BLUE, lw=0.8)
    ax.plot(
        p,
        TAN_FIX * p,
        "-",
        color=RED,
        label=rf"$Q=P\tan(\arccos {COS_PHI_FIXED:.2f})$",
    )
    ax.annotate(
        r"$\cos\varphi$ cone, for comparison",
        xy=(0.05, 0.50),
        fontsize=5.5,
        color=BLUE,
    )
    ax.annotate(
        r"an equality: a line, not an area",
        xy=(0.05, 0.40),
        fontsize=5.5,
        color="#52514E",
    )
    ax.set_xlim(0, 1.0)
    ax.set_ylim(-0.7, 0.7)
    _axes(ax, r"$P/P_n$", r"$Q/P_n$", r"Fixed $\cos\varphi$")
    ax.legend(loc="lower left", fontsize=6, frameon=False)
    _save(fig, "fixed-cos-phi")


def fig_cpp():
    """cos(phi)(P) profile: a quadratic equality through the P-Q plane."""
    pt = GRID_CODE.cpp_p_threshold_pu
    p = np.linspace(0.0, 1.0, 400)
    q = TAN_MIN * p * (p - pt) / (1.0 - pt)
    fig, ax = plt.subplots()
    ax.fill_between(
        p,
        -TAN_MIN * p,
        TAN_MIN * p,
        color=BLUE,
        alpha=0.08,
        linewidth=0,
    )
    ax.plot(p, TAN_MIN * p, ":", color=BLUE, lw=0.8)
    ax.plot(p, -TAN_MIN * p, ":", color=BLUE, lw=0.8)
    ax.plot(p, q, "-", color=RED, label=r"$\cos\varphi(P)$ profile")
    ax.annotate(
        r"$\cos\varphi$ cone, for comparison",
        xy=(0.30, -0.45),
        fontsize=5.5,
        color=BLUE,
    )
    ax.axvline(pt, color="#898781", lw=0.6, ls=(0, (1, 2)))
    ax.annotate(
        rf"$P_t={pt:.1f}\,P_n$",
        xy=(pt + 0.02, -0.55),
        fontsize=6,
        color="#52514E",
    )
    ax.set_xlim(0, 1.0)
    ax.set_ylim(-0.7, 0.7)
    _axes(ax, r"$P/P_n$", r"$Q/P_n$", r"$\cos\varphi(P)$ profile")
    ax.legend(loc="upper left", fontsize=6, frameon=False)
    _save(fig, "cos-phi-p-profile")


FIGURES = (
    fig_qp,
    fig_qu,
    fig_qu_deadband,
    fig_s2,
    fig_cone,
    fig_operating_region,
    fig_pu_curtail,
    fig_fixed_cos_phi,
    fig_cpp,
)


def main():
    print(
        f"grid code: {GRID_CODE.title} ({GRID_CODE.voltage_level})", flush=True
    )
    with rwthplots.context(*STYLES):
        for make in FIGURES:
            make()


if __name__ == "__main__":
    main()
