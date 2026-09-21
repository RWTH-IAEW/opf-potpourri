# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Grid-code reactive-power capability parameters, aligned with pandapower.

The German technical connection rules (TAR) are represented as selectable
:class:`GridCode` parameter sets, so a study can target the rule that applies
to its voltage level:

* :data:`VDE_AR_N_4105` — low voltage, 2 variants.
* :data:`VDE_AR_N_4110` — medium voltage, 1 variant.
* :data:`VDE_AR_N_4120` — high voltage (110 kV), 3 variants.

Select one by short name (``"4105"``, ``"4110"``, ``"4120"``) or by passing
the :class:`GridCode` itself.

Alignment with pandapower
-------------------------
Every capability area here reproduces the corresponding class in
:mod:`pandapower.control.controller.DERController` to machine precision --
:class:`~pandapower.control.controller.DERController.PQVAreas.PQArea4105`,
``PQArea4110``, ``PQArea4120``, ``QVArea4105``, ``QVArea4110`` and
``QVArea4120``.  pandapower models these as shapely polygons (LV/MV) or as
explicit branch logic (HV); both collapse to the same piecewise-linear
envelope sampled with :func:`numpy.interp`, which is the form the Pyomo
constraints need.  ``tests/unit_tests/test_q_control.py`` asserts the
equivalence against pandapower's own ``q_flexibility()``.

Two deliberate differences:

* **Q is bounded, not assigned.** pandapower's controllers *set* Q during a
  time-series simulation; the OPF treats Q as a decision variable and these
  areas as its feasible region.  :class:`QVCurve` covers the assigning case.
* **VDE-AR-N 4130 (EHV) is not included.**  It needs ``vn_kv``-dependent
  breakpoints for 380/220 kV, and potpourri targets distribution grids.  Use
  pandapower's ``PQVArea4130*`` directly if you need it.

Consumed by :class:`~potpourri.technologies.pv.PV_multi_period`,
:class:`~potpourri.technologies.sgens.Sgens_multi_period` and the
single-period :class:`~potpourri.models.ACOPF_base.ACOPF`.
"""

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd


class ProvisionalGridCodeWarning(UserWarning):
    """Raised when a grid code whose parameters are placeholders is used."""


class SgenTypeOverlapWarning(UserWarning):
    """Raised when an sgen is claimed by both the PV and the wind Q path."""


class QuCurveOutsidePqAreaWarning(UserWarning):
    """Raised when a Q(U) characteristic conflicts with the Q(P) area.

    Both are imposed on the same reactive power, and they express different
    things: the Q(P) area *bounds* Q from active power, a Q(U)
    characteristic *assigns* it from voltage.  Where the assigned value lies
    outside the bound the model is infeasible, and the solver says only
    "infeasible".  See :func:`warn_if_curve_leaves_pq_area`.
    """


class EnvelopeRangeWarning(UserWarning):
    """Raised when the operating range leaves an envelope's exact span.

    The Pyomo constraints represent each capability bound as a set of linear
    inequalities, which reproduces the envelope exactly only between its
    outermost breakpoints.  Beyond them the bound is replaced by its concave
    majorant (upper) or convex minorant (lower), so the model permits
    slightly *more* reactive power than the grid code strictly requires.
    That direction is deliberate: the exact area is non-convex there, and
    the alternative — extrapolating the end segment — makes the two bounds
    cross and the model infeasible.  See :meth:`Envelope.exact_range`.
    """


# sgen ``type`` values treated as wind by the wind Q-control path.  SimBench
# spells wind differently per voltage level — its RES dataset uses "Wind" in
# HV, "Wind_MV" in MV and "wind onshore"/"wind offshore" in EHV — so matching
# only "Wind" reaches nothing on any SimBench MV or EHV grid.  Matching is
# exact and case-sensitive.  Lives here rather than in a model module so the
# single-period and multi-period wind paths cannot drift apart.
DEFAULT_WIND_SGEN_TYPES = (
    "Wind",
    "Wind_MV",
    "wind onshore",
    "wind offshore",
)

# Range of P/Pn a generator can be dispatched over: curtailed to zero, up to
# the installed power.  Passed as ``x_range`` when building the Q(P) pieces
# so they stay valid below the standard's first active-power breakpoint.
# Without it the two bounds extrapolate past each other and no reactive power
# is feasible at all -- for VDE-AR-N 4120 that happened below P = 0.061 Pn.
DEFAULT_P_RANGE_PU = (0.0, 1.0)

# Fallback bus voltage band when net.bus carries no limits, matching
# pandapower's own default.
DEFAULT_V_RANGE_PU = (0.9, 1.1)


def _cross(o, a, b):
    """Z-component of (a-o) x (b-o); sign gives the turn direction."""
    return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])


def _hull(points, upper):
    """Vertices of the concave majorant (``upper``) or convex minorant.

    A capability bound is only representable as a min/max of affine pieces
    where it is concave/convex.  Taking the hull first makes the piece form
    valid over the whole requested range: it can only *widen* the feasible
    band, never narrow it below what the grid code allows.
    """
    out = []
    for p in points:
        while len(out) >= 2 and (
            _cross(out[-2], out[-1], p) >= 0
            if upper
            else _cross(out[-2], out[-1], p) <= 0
        ):
            out.pop()
        out.append(p)
    return out


def _affine_pieces(x_points, y_points, tol=1e-12):
    """Return ``(slope, intercept)`` for each distinct linear segment.

    Consecutive segments with the same slope are collapsed, so a bound made
    of a flat shelf and a ramp yields two pieces rather than one per
    breakpoint interval.
    """
    pieces = []
    for i in range(len(x_points) - 1):
        x0, x1 = float(x_points[i]), float(x_points[i + 1])
        y0, y1 = float(y_points[i]), float(y_points[i + 1])
        if abs(x1 - x0) < tol:
            continue  # vertical step: carries no affine piece
        m = (y1 - y0) / (x1 - x0)
        b = y0 - m * x0
        if pieces:
            m_prev, b_prev = pieces[-1]
            if abs(m_prev - m) < tol and abs(b_prev - b) < tol:
                continue
        pieces.append((m, b))
    return pieces


@dataclass(frozen=True, eq=False)
class Envelope:
    """Piecewise-linear lower and upper Q bounds over a driving quantity.

    The driving quantity is active power (``P/Pn``) for a PQ area and bus
    voltage (p.u.) for a QV area.  Evaluation follows :func:`numpy.interp`,
    which holds the end values outside the breakpoint span — exactly the
    saturation behaviour of pandapower's ``q_flexibility()``.

    Attributes:
        x_points: Breakpoints, shape ``(n_points,)``, ascending.
        q_min: Lower bound per variant, shape ``(n_variants, n_points)``.
        q_max: Upper bound per variant, same shape.  All Q values are
            normalised to the installed active power ``Pn``.
    """

    x_points: np.ndarray
    q_min: np.ndarray
    q_max: np.ndarray

    def __post_init__(self):
        """Validate and normalise the Q(P) envelope on construction.

        Coerces the point arrays to float, promotes a single-variant envelope
        to 2-D so every consumer can index by variant, and rejects a
        non-monotonic abscissa, which would make the piecewise interpolation
        ambiguous.

        Raises:
            ValueError: If the points are not strictly increasing or the bound
                arrays disagree in shape.
        """
        x = np.asarray(self.x_points, dtype=float)
        lo = np.atleast_2d(np.asarray(self.q_min, dtype=float))
        hi = np.atleast_2d(np.asarray(self.q_max, dtype=float))
        if np.any(np.diff(x) < 0):
            raise ValueError("x_points must be ascending")
        if lo.shape != hi.shape:
            raise ValueError(f"q_min {lo.shape} and q_max {hi.shape} differ")
        if lo.shape[1] != x.shape[0]:
            raise ValueError(
                f"q bounds have {lo.shape[1]} points but x_points has "
                f"{x.shape[0]}"
            )
        if np.any(lo > hi + 1e-12):
            raise ValueError("q_min exceeds q_max at some breakpoint")
        object.__setattr__(self, "x_points", x)
        object.__setattr__(self, "q_min", lo)
        object.__setattr__(self, "q_max", hi)

    @property
    def n_variants(self) -> int:
        """Number of variants this envelope defines."""
        return int(self.q_min.shape[0])

    def q_flexibility(self, x, variant: int = 0):
        """Return ``(q_min, q_max)`` at ``x``, pandapower-compatible.

        Args:
            x: Driving quantity — ``P/Pn`` or voltage [p.u.]; scalar or array.
            variant: Variant index, i.e. the value of ``net.sgen.var_q``.

        Returns:
            Tuple of lower and upper Q/Pn bounds, matching the shape of ``x``.
        """
        self._check_variant(variant)
        return (
            np.interp(x, self.x_points, self.q_min[variant]),
            np.interp(x, self.x_points, self.q_max[variant]),
        )

    def _bound_points(self, variant, upper, x_range):
        """Breakpoints of one bound, padded to cover ``x_range``.

        Padding holds the end value, matching how the envelope itself
        saturates.  Without it the outermost affine piece extrapolates, and
        for a bound whose ramp starts above the operating minimum the two
        extrapolations cross: the lower bound rises past the upper one and
        no Q is feasible at all.  That is what made every VDE-AR-N 4120
        model infeasible below P = 0.061 Pn before 0.4.1.
        """
        y = (self.q_max if upper else self.q_min)[variant]
        xs = list(map(float, self.x_points))
        ys = list(map(float, y))
        if x_range is not None:
            lo, hi = float(x_range[0]), float(x_range[1])
            if lo < xs[0] - 1e-12:
                xs.insert(0, lo)
                ys.insert(0, ys[0])
            if hi > xs[-1] + 1e-12:
                xs.append(hi)
                ys.append(ys[-1])
        return list(zip(xs, ys))

    def upper_pieces(self, variant: int = 0, x_range=None):
        """Affine pieces whose pointwise **minimum** is the upper bound.

        Args:
            variant: Variant index (``net.sgen.var_q``).
            x_range: ``(min, max)`` the driving quantity can take.  When
                given, the bound is padded to span it and reduced to its
                concave majorant, so the piece form stays valid — and the
                model stays feasible — outside the breakpoints.  Where the
                bound is already concave over the range this changes
                nothing.

        Returns:
            List of ``(slope, intercept)`` pairs.
        """
        self._check_variant(variant)
        pts = self._bound_points(variant, True, x_range)
        if x_range is not None:
            pts = _hull(pts, upper=True)
        return _affine_pieces([p[0] for p in pts], [p[1] for p in pts])

    def lower_pieces(self, variant: int = 0, x_range=None):
        """Affine pieces whose pointwise **maximum** is the lower bound.

        The mirror of :meth:`upper_pieces`; ``x_range`` reduces the bound to
        its convex minorant.
        """
        self._check_variant(variant)
        pts = self._bound_points(variant, False, x_range)
        if x_range is not None:
            pts = _hull(pts, upper=False)
        return _affine_pieces([p[0] for p in pts], [p[1] for p in pts])

    def max_pieces(self, x_range=None) -> int:
        """Largest piece count over both bounds and all variants."""
        return max(
            max(
                len(self.upper_pieces(v, x_range)),
                len(self.lower_pieces(v, x_range)),
            )
            for v in range(self.n_variants)
        )

    def exact_range(self):
        """Span over which the linear-inequality form is exact.

        ``Q <= min(pieces)`` reproduces a *concave* upper bound and
        ``Q >= max(pieces)`` a *convex* lower bound.  Both hold between the
        outermost breakpoints.  Beyond them the shape is non-convex, so
        :meth:`upper_pieces` and :meth:`lower_pieces` fall back to the hull
        when given an ``x_range`` that reaches outside — exact here, a
        bounded relaxation there.
        """
        return float(self.x_points[0]), float(self.x_points[-1])

    def is_linearisable(self, variant: int = 0) -> bool:
        """Whether min/max-of-affines reproduces this variant exactly.

        True when the upper bound is concave (non-increasing slopes) and the
        lower bound convex (non-decreasing slopes) across the breakpoint
        span.  Every grid code shipped here satisfies both; a hand-built
        :class:`Envelope` need not.
        """
        upper = [m for m, _ in self.upper_pieces(variant)]
        lower = [m for m, _ in self.lower_pieces(variant)]
        return all(np.diff(upper) <= 1e-12) and all(np.diff(lower) >= -1e-12)

    def _check_variant(self, variant):
        """Reject an out-of-range grid-code variant.

        Returns:
            A Pyomo expression.
        """
        # Reject fractional values rather than truncating them: silently
        # rounding var_q would select a neighbouring capability column.
        if not float(variant).is_integer():
            raise IndexError(
                f"variant {variant!r} is not a whole number; var_q selects "
                f"a capability column and must be an integer"
            )
        if not 0 <= int(variant) < self.n_variants:
            raise IndexError(
                f"variant {variant} out of range: this envelope defines "
                f"{self.n_variants} variant(s), so var_q must be in "
                f"0..{self.n_variants - 1}"
            )


@dataclass(frozen=True, eq=False)
class QVCurve:
    r"""Piecewise-linear Q(U) **setpoint** characteristic, with dead band.

    Mirrors ``QVCurve`` in
    :mod:`pandapower.control.controller.DERController.DERBasics`.
    Unlike an :class:`Envelope`, this *assigns* Q as a function of voltage
    rather than bounding it, so a dead band — a voltage span over which Q is
    held at zero — is expressible::

        Q/Pn
          |
      q_max +-------\
          |          \
        0 +           \________            <- dead band, Q = 0
          |                     \
      q_min|                      \-------
          +---+------+--------+------+----> V [p.u.]
             V1     V2       V3     V4

    The graph of a non-constant curve is not a convex set, so the OPF cannot
    express this with plain inequalities.  It is built with
    :class:`pyomo.core.base.piecewise.Piecewise` and therefore needs a
    MIP-capable solver (MindtPy, Gurobi) rather than IPOPT alone.

    Attributes:
        v_points_pu: Voltage breakpoints [p.u.], shape ``(n_points,)``.
        q_points_pu: Q/Pn at each breakpoint, shape
            ``(n_variants, n_points)``.
    """

    v_points_pu: np.ndarray
    q_points_pu: np.ndarray

    def __post_init__(self):
        """Validate and normalise the Q(U) envelope on construction.

        As the Q(P) validator, for a characteristic in voltage: the voltage
        points must be strictly increasing.

        Raises:
            ValueError: If the voltage points are not strictly increasing or
                the bound arrays disagree in shape.
        """
        v = np.asarray(self.v_points_pu, dtype=float)
        q = np.atleast_2d(np.asarray(self.q_points_pu, dtype=float))
        if np.any(np.diff(v) <= 0):
            raise ValueError("v_points_pu must be strictly ascending")
        if q.shape[1] != v.shape[0]:
            raise ValueError(
                f"q_points_pu has {q.shape[1]} points but v_points_pu has "
                f"{v.shape[0]}"
            )
        object.__setattr__(self, "v_points_pu", v)
        object.__setattr__(self, "q_points_pu", q)

    @property
    def n_variants(self) -> int:
        """Number of variants this curve defines."""
        return int(self.q_points_pu.shape[0])

    def step(self, vm_pu, variant: int = 0):
        """Return Q/Pn at ``vm_pu`` — the pandapower ``QVCurve.step`` API."""
        if not 0 <= int(variant) < self.n_variants:
            raise IndexError(
                f"variant {variant} out of range: this curve defines "
                f"{self.n_variants} variant(s)"
            )
        return np.interp(vm_pu, self.v_points_pu, self.q_points_pu[variant])

    def padded(self, v_min, v_max) -> "QVCurve":
        """Extend the curve flat so it spans at least ``[v_min, v_max]``.

        The SOS2 formulation needs its input variable bounded inside the
        breakpoint span, so a curve narrower than the bus voltage limits
        would make the model infeasible rather than saturating.  Padding
        holds the end values, which is what the characteristic does anyway
        outside its outermost breakpoints.

        Args:
            v_min: Lowest voltage the model can reach [p.u.].
            v_max: Highest voltage the model can reach [p.u.].

        Returns:
            A curve covering the range; ``self`` if it already does.
        """
        v = self.v_points_pu
        q = self.q_points_pu
        lo_needed = float(v_min) < v[0] - 1e-12
        hi_needed = float(v_max) > v[-1] + 1e-12
        if not (lo_needed or hi_needed):
            return self
        v_new = list(v)
        q_new = [list(row) for row in q]
        if lo_needed:
            v_new.insert(0, float(v_min))
            for row in q_new:
                row.insert(0, row[0])
        if hi_needed:
            v_new.append(float(v_max))
            for row in q_new:
                row.append(row[-1])
        return QVCurve(np.asarray(v_new), np.asarray(q_new))

    def deadband(self, variant: int = 0, tol: float = 1e-12):
        """Return the ``(v_low, v_high)`` zero-Q span, or ``None``."""
        q = self.q_points_pu[variant]
        zero = np.flatnonzero(np.abs(q) <= tol)
        if zero.size < 2:
            return None
        return (
            float(self.v_points_pu[zero[0]]),
            float(self.v_points_pu[zero[-1]]),
        )


def deadband_qv_curve(v_points_pu, q_max, q_min=None) -> QVCurve:
    """Build a dead-band Q(U) characteristic from four breakpoints.

    Args:
        v_points_pu: ``[V1, V2, V3, V4]`` — full capacitive support below
            ``V1``, ramp to zero over ``V1..V2``, dead band over ``V2..V3``,
            ramp to full inductive over ``V3..V4``.
        q_max: Capacitive (overexcited) limit Q/Pn, scalar or per variant.
        q_min: Inductive limit Q/Pn.  Defaults to ``-q_max``.

    Returns:
        The :class:`QVCurve`.
    """
    v = np.asarray(v_points_pu, dtype=float)
    if v.shape != (4,):
        raise ValueError(f"expected 4 voltage breakpoints, got {v.shape[0]}")
    hi = np.atleast_1d(np.asarray(q_max, dtype=float))
    lo = (
        -hi if q_min is None else np.atleast_1d(np.asarray(q_min, dtype=float))
    )
    if hi.shape != lo.shape:
        raise ValueError(
            "q_max and q_min must have the same number of variants"
        )
    zero = np.zeros_like(hi)
    return QVCurve(v, np.column_stack([hi, zero, zero, lo]))


@dataclass(frozen=True, eq=False)
class GridCode:
    """Reactive-power capability parameters for one grid code.

    Attributes:
        name: Short identifier, e.g. ``"4110"``.
        title: Full designation, e.g. ``"VDE-AR-N 4110"``.
        voltage_level: Voltage level the rule applies to.
        pq_area: Q bounds as a function of ``P/Pn``.
        qv_area: Q bounds as a function of bus voltage [p.u.].
        vpu_v_curtail: Voltage above which P(U) curtailment begins [p.u.].
        vpu_v_max: Voltage at which P(U) curtailment reaches zero [p.u.].
        cpp_p_threshold_pu: P/Pn below which the cos(phi)(P) profile
            requires no reactive power.
        provisional: ``True`` when the values are placeholders rather than
            normative figures taken from the standard.
        provisional_note: Explanation surfaced in the warning when
            ``provisional`` is set.
    """

    name: str
    title: str
    voltage_level: str
    pq_area: Envelope
    qv_area: Envelope
    vpu_v_curtail: float
    vpu_v_max: float
    cpp_p_threshold_pu: float
    provisional: bool = False
    provisional_note: str = ""

    def __post_init__(self):
        """Check that the grid code's two envelopes agree.

        A code carries a Q(P) and a Q(U) area, and both must describe the same
        set of operating variants -- otherwise `var_q` would select different
        behaviour in each.

        Raises:
            ValueError: If the two areas declare a different number of
                variants.
        """
        if self.pq_area.n_variants != self.qv_area.n_variants:
            raise ValueError(
                f"{self.title}: PQ area has {self.pq_area.n_variants} "
                f"variants but QV area has {self.qv_area.n_variants}"
            )

    @property
    def n_variants(self) -> int:
        """Number of ``var_q`` variants this grid code defines."""
        return self.qv_area.n_variants

    def q_flexibility(self, p_pu, vm_pu, variant: int = 0):
        """Return ``(q_min, q_max)`` from both areas intersected.

        The pandapower ``BasePQVArea.q_flexibility`` semantics: take the
        tighter of the PQ and QV bounds.

        Args:
            p_pu: Active power P/Pn.
            vm_pu: Bus voltage [p.u.].
            variant: Variant index (``net.sgen.var_q``).

        Returns:
            Tuple of lower and upper Q/Pn bounds.
        """
        pq_lo, pq_hi = self.pq_area.q_flexibility(p_pu, variant)
        qv_lo, qv_hi = self.qv_area.q_flexibility(vm_pu, variant)
        return np.maximum(pq_lo, qv_lo), np.minimum(pq_hi, qv_hi)

    def deadband_curve(self, deadband=None, variant=None) -> QVCurve:
        """Build a dead-band Q(U) characteristic from this code's own limits.

        The result reuses the QV area's voltage breakpoints and reactive
        limits but passes through zero in the middle, turning the capability
        *area* into a control *characteristic*.

        Note:
            The dead band is a network-operator parameterisation, not a
            normative value.  The default reproduces the plateau of this
            code's QV area; pass ``deadband`` to set it explicitly.

        Args:
            deadband: ``(v_low, v_high)`` of the zero-Q span [p.u.].
                Defaults to this code's QV plateau.
            variant: Restrict to a single variant, or ``None`` for all.

        Returns:
            The :class:`QVCurve`.
        """
        v = self.qv_area.x_points
        if v.shape[0] != 4:
            raise ValueError(
                f"{self.title}: deadband_curve needs a 4-breakpoint QV area, "
                f"this one has {v.shape[0]}"
            )
        v_pts = np.array([v[0], v[1], v[2], v[3]], dtype=float)
        if deadband is not None:
            lo_db, hi_db = float(deadband[0]), float(deadband[1])
            if not v[0] <= lo_db <= hi_db <= v[3]:
                raise ValueError(
                    f"deadband {deadband} must lie inside "
                    f"[{v[0]:.4f}, {v[3]:.4f}] and be ordered"
                )
            v_pts = np.array([v[0], lo_db, hi_db, v[3]], dtype=float)
        variants = (
            range(self.n_variants) if variant is None else [int(variant)]
        )
        q_hi = np.array([self.qv_area.q_max[i].max() for i in variants])
        q_lo = np.array([self.qv_area.q_min[i].min() for i in variants])
        return deadband_qv_curve(v_pts, q_hi, q_lo)

    # --- deprecated shape-specific accessors -----------------------------
    # Kept so callers written against the pre-envelope registry keep working.
    @property
    def vqu_v_points(self) -> np.ndarray:
        """Q(U) voltage breakpoints as ``[[V1, V2], [V3, V4]]``."""
        return self.qv_area.x_points.reshape(2, 2)

    @property
    def vqu_q_max(self) -> np.ndarray:
        """``[[capacitive at low V], [inductive at high V]]`` per variant."""
        return np.vstack([self.qv_area.q_max[:, 0], self.qv_area.q_min[:, -1]])

    @property
    def qp_p_high(self) -> float:
        """Lower P/Pn breakpoint of the Q(P) characteristic."""
        return float(self.pq_area.x_points[0])

    @property
    def qp_p_low(self) -> float:
        """Upper P/Pn breakpoint of the Q(P) characteristic."""
        return float(self.pq_area.x_points[1])


def _vde_qv_area(v_points, q_min, q_max, low_v_q_max=None, high_v_q_min=None):
    """Build the "hold, ramp, hold" QV envelope shared by the VDE rules.

    Args:
        v_points: ``[V1, V2, V3, V4]``.
        q_min: Inductive limit per variant (negative).
        q_max: Capacitive limit per variant (positive).
        low_v_q_max: Lower bound at ``V1``.  Defaults to ``q_max`` (the 4120
            shape, where the area pinches to a point below ``V1``); 4105 and
            4110 pass zeros instead.
        high_v_q_min: Upper bound at ``V4``.  Defaults to ``q_min``.
    """
    lo = np.atleast_1d(np.asarray(q_min, dtype=float))
    hi = np.atleast_1d(np.asarray(q_max, dtype=float))
    at_v1 = hi if low_v_q_max is None else np.atleast_1d(low_v_q_max)
    at_v4 = lo if high_v_q_min is None else np.atleast_1d(high_v_q_min)
    return Envelope(
        x_points=np.asarray(v_points, dtype=float),
        q_min=np.column_stack([at_v1, lo, lo, lo]),
        q_max=np.column_stack([hi, hi, hi, at_v4]),
    )


# cos(phi) = 0.95 and 0.90 expressed as tan(phi) — the reactive limits the
# VDE rules are written in terms of.
_Q_COSPHI_095 = 0.328684
_Q_COSPHI_090 = 0.484322


# --- VDE-AR-N 4105 (low voltage) -----------------------------------------
# Variant 1 applies to plants with S_E,max <= 4.6 kVA (cos phi 0.95), variant
# 2 above that (cos phi 0.90).  Matches pandapower PQArea4105 / QVArea4105.
_Q4105 = np.array([_Q_COSPHI_095, _Q_COSPHI_090])

VDE_AR_N_4105 = GridCode(
    name="4105",
    title="VDE-AR-N 4105",
    voltage_level="low voltage",
    # Q ramps linearly from zero at P = 0 to the full limit at P = Pn.
    pq_area=Envelope(
        x_points=np.array([0.0, 1.0]),
        q_min=np.column_stack([np.zeros(2), -_Q4105]),
        q_max=np.column_stack([np.zeros(2), _Q4105]),
    ),
    qv_area=_vde_qv_area(
        [0.90, 0.95, 1.05, 1.10],
        q_min=-_Q4105,
        q_max=_Q4105,
        low_v_q_max=np.zeros(2),
        high_v_q_min=np.zeros(2),
    ),
    vpu_v_curtail=1.06,  # VDE-AR-N 4105 section 8.5
    vpu_v_max=1.10,
    cpp_p_threshold_pu=0.2,
)


# --- VDE-AR-N 4110 (medium voltage) --------------------------------------
# Single variant at cos phi 0.90, with a 0.05 p.u. active-power threshold
# below which only a token reactive range is required.  Matches pandapower
# PQArea4110 / QVArea4110.
_Q4110 = np.array([_Q_COSPHI_090])

VDE_AR_N_4110 = GridCode(
    name="4110",
    title="VDE-AR-N 4110",
    voltage_level="medium voltage",
    pq_area=Envelope(
        x_points=np.array([0.05, 1.0]),
        q_min=np.array([[-0.01961505, -_Q_COSPHI_090]]),
        q_max=np.array([[0.01961505, _Q_COSPHI_090]]),
    ),
    qv_area=_vde_qv_area(
        [0.90, 0.95, 1.05, 1.10],
        q_min=-_Q4110,
        q_max=_Q4110,
        low_v_q_max=np.zeros(1),
        high_v_q_min=np.zeros(1),
    ),
    vpu_v_curtail=1.06,
    vpu_v_max=1.10,
    cpp_p_threshold_pu=0.2,
)


# --- VDE-AR-N 4120 (high voltage, 110 kV) --------------------------------
# Three variants, distinguished by how the capability is split between
# capacitive and inductive.  The voltage breakpoints are 96 / 103 / 120 /
# 127 kV on the 110 kV base, which is why the span reaches 0.87-1.15 p.u.
# rather than the 0.9-1.1 of the LV and MV rules.  Matches pandapower
# PQArea4120(version=2015) / QVArea4120.
_Q4120_MIN = np.array([-0.227902, -0.328684, -0.410775])
_Q4120_MAX = np.array([0.484322, 0.410775, 0.328684])

VDE_AR_N_4120 = GridCode(
    name="4120",
    title="VDE-AR-N 4120",
    voltage_level="high voltage",
    pq_area=Envelope(
        x_points=np.array([0.1, 0.2, 1.0]),
        q_min=np.column_stack([np.full(3, -0.1), _Q4120_MIN, _Q4120_MIN]),
        q_max=np.column_stack([np.full(3, 0.1), _Q4120_MAX, _Q4120_MAX]),
    ),
    qv_area=_vde_qv_area(
        np.array([96, 103, 120, 127]) / 110.0,
        q_min=_Q4120_MIN,
        q_max=_Q4120_MAX,
    ),
    vpu_v_curtail=1.06,
    vpu_v_max=1.10,
    cpp_p_threshold_pu=0.2,
)


GRID_CODES = {
    VDE_AR_N_4105.name: VDE_AR_N_4105,
    VDE_AR_N_4110.name: VDE_AR_N_4110,
    VDE_AR_N_4120.name: VDE_AR_N_4120,
}

# 4120 is the default because it is what the pre-0.4.1 constants encoded:
# the breakpoints and variant values shipped under the "4105" name were in
# fact the 110 kV rule.  Keeping it as the default means existing nets, whose
# var_q spans 0..2, keep resolving to a three-variant code.
DEFAULT_GRID_CODE = VDE_AR_N_4120


# --- Backwards-compatible module-level constants -------------------------
# These mirror the default grid code so that existing imports keep working.
# New code should read the envelopes off a GridCode instead.
VQU_V_POINTS = DEFAULT_GRID_CODE.vqu_v_points
VQU_Q_MAX = DEFAULT_GRID_CODE.vqu_q_max
QP_P_HIGH = DEFAULT_GRID_CODE.qp_p_high
QP_P_LOW = DEFAULT_GRID_CODE.qp_p_low
VPU_V_CURTAIL = DEFAULT_GRID_CODE.vpu_v_curtail
VPU_V_MAX = DEFAULT_GRID_CODE.vpu_v_max
CPP_P_THRESHOLD_PU = DEFAULT_GRID_CODE.cpp_p_threshold_pu


def resolve_grid_code(grid_code=None) -> GridCode:
    """Return the :class:`GridCode` for ``grid_code``.

    Args:
        grid_code: ``None`` for :data:`DEFAULT_GRID_CODE`, a short name such
            as ``"4105"`` / ``"4110"`` / ``"4120"`` (``"VDE-AR-N 4110"`` is
            also accepted), or a :class:`GridCode` instance, which is
            returned unchanged.

    Returns:
        The resolved :class:`GridCode`.

    Raises:
        ValueError: If ``grid_code`` names an unknown grid code.

    Warns:
        ProvisionalGridCodeWarning: If the resolved grid code carries
            placeholder rather than normative parameters.
    """
    if grid_code is None:
        resolved = DEFAULT_GRID_CODE
    elif isinstance(grid_code, GridCode):
        resolved = grid_code
    else:
        key = str(grid_code).strip()
        # Accept "VDE-AR-N 4110" and "vde-ar-n-4110" as well as "4110".
        normalised = key.upper().replace(" ", "").replace("-", "")
        normalised = normalised.replace("VDEARN", "")
        resolved = GRID_CODES.get(key) or GRID_CODES.get(normalised)
        if resolved is None:
            raise ValueError(
                f"Unknown grid code {grid_code!r}. "
                f"Available: {sorted(GRID_CODES)}"
            )

    if resolved.provisional:
        warnings.warn(
            resolved.provisional_note
            or f"{resolved.title} carries placeholder parameters.",
            ProvisionalGridCodeWarning,
            stacklevel=3,
        )
    return resolved


def check_var_q(values, grid_code=None, context="") -> None:
    """Raise if any ``var_q`` value is not a variant the grid code defines.

    ``var_q`` selects one column of the capability table, so it has to be a
    whole number in range.  Integral floats are fine — pandas stores the
    column as ``float64`` whenever it contains NaN — but a fractional value
    is rejected rather than truncated, which would silently pick a
    neighbouring variant and change dispatch with no error anywhere.

    Args:
        values: Iterable of ``net.sgen.var_q`` entries; NaN is ignored.
        grid_code: Grid code to check against, as accepted by
            :func:`resolve_grid_code`.
        context: Optional prefix for the error message.

    Raises:
        ValueError: If any value is fractional or out of range.
    """
    code = resolve_grid_code(grid_code)
    arr = pd.Series(list(values)).dropna()
    if arr.empty:
        return
    fractional, out_of_range = [], []
    for raw in arr:
        try:
            val = float(raw)
        except (TypeError, ValueError):
            fractional.append(raw)
            continue
        if not val.is_integer():
            fractional.append(raw)
        elif not 0 <= int(val) < code.n_variants:
            out_of_range.append(int(val))
    prefix = f"{context}: " if context else ""
    if fractional:
        raise ValueError(
            f"{prefix}var_q values {sorted(map(str, set(fractional)))} are "
            f"not whole numbers. var_q selects a column of the "
            f"{code.title} capability table, so it must be an integer "
            f"0..{code.n_variants - 1}, not a value to round."
        )
    if out_of_range:
        raise ValueError(
            f"{prefix}var_q values {sorted(set(out_of_range))} are not "
            f"valid for {code.title}, which defines {code.n_variants} "
            f"variant(s) (var_q must be 0..{code.n_variants - 1}). Either "
            f"set var_q to a valid variant or select a grid code with more "
            f"variants."
        )


def bus_voltage_range(net, default=DEFAULT_V_RANGE_PU):
    """Widest voltage band any bus in ``net`` may operate in.

    Passed as ``x_range`` when building the Q(U) pieces.  A grid whose
    limits sit inside the grid code's own voltage span — the usual case —
    gets exactly the same pieces either way; a grid that allows wider
    excursions gets the hull instead of an extrapolated end segment, which
    is what stops the model from demanding reactive power the standard does
    not ask for.  For VDE-AR-N 4105 at 0.85 p.u. the extrapolation requires
    Q >= +0.329 Pn where the standard requires nothing at all.

    Args:
        net: pandapower network.
        default: Band to assume when the columns are absent or all NaN.

    Returns:
        ``(v_min, v_max)`` in per unit.
    """
    lo, hi = default
    if "min_vm_pu" in net.bus:
        col = net.bus["min_vm_pu"].astype(float)
        if col.notna().any():
            lo = min(lo, float(col.min()))
    if "max_vm_pu" in net.bus:
        col = net.bus["max_vm_pu"].astype(float)
        if col.notna().any():
            hi = max(hi, float(col.max()))
    return lo, hi


def resolve_qu_curve(spec, grid_code=None):
    """Turn a ``qu_deadband`` argument into a :class:`QVCurve`, or ``None``.

    Args:
        spec: ``None``/``False`` for no dead band (the capability *area* is
            used instead); ``True`` for the grid code's own QV plateau as
            the dead band; a ``(v_low, v_high)`` pair to set it explicitly;
            or a ready-made :class:`QVCurve`.
        grid_code: Code supplying the reactive limits and breakpoints, as
            accepted by :func:`resolve_grid_code`.

    Returns:
        The curve to pin Q to, or ``None`` to keep the area formulation.
    """
    if spec is None or spec is False:
        return None
    if isinstance(spec, QVCurve):
        return spec
    code = resolve_grid_code(grid_code)
    if spec is True:
        return code.deadband_curve()
    return code.deadband_curve(deadband=spec)


def attach_deadband_qu(
    model,
    name,
    keys,
    *,
    q_of,
    v_of,
    pn_of,
    variant_of,
    curve,
    v_bounds=None,
    pw_repn="INC",
):
    """Attach a dead-band Q(U) characteristic as a piecewise equality.

    Where an :class:`Envelope` *bounds* Q and leaves the optimiser free
    inside the band, this *pins* Q to the curve, reproducing what a Q(U)
    droop controller with a dead band actually does.  The graph of the curve
    is not a convex set, so this is a genuine integer formulation: it needs
    a MIP-capable solver (MindtPy, CBC, GLPK, Gurobi), not IPOPT alone.

    Creates, prefixed with ``name``: ``_IDX`` (index set), ``_v`` and
    ``_q_pu`` (auxiliary variables the piecewise block operates on),
    ``_v_link`` and ``_q_link`` (tying them to the model's own voltage and
    reactive power) and ``_pw`` (the equality itself).  The auxiliaries
    exist because Pyomo's ``Piecewise`` needs its input and output indexed
    alike, whereas voltage is indexed by bus and Q by generator.

    Args:
        model: Pyomo model to attach to.
        name: Prefix for every created component.
        keys: Index tuples identifying each Q-controlled element.
        q_of: ``key -> `` the reactive-power variable to pin (per unit).
        v_of: ``key -> `` the bus-voltage variable driving the curve.
        pn_of: ``key -> float``, installed active power (per unit).
        variant_of: ``key -> int``, the element's ``var_q``.
        curve: The :class:`QVCurve` to follow.
        v_bounds: ``(v_min, v_max)`` the model can reach.  The curve is
            padded to cover it, so a curve narrower than the bus limits
            saturates instead of making the model infeasible.
        pw_repn: Pyomo piecewise representation.  The default ``"INC"``
            (incremental) is pure-binary and works with any MIP solver;
            ``"SOS2"`` is tighter but needs solver SOS support, which CBC
            and Gurobi have and GLPK does not.

    Returns:
        The number of elements constrained.
    """
    import pyomo.environ as pyo

    keys = list(keys)
    if not keys:
        return 0
    if v_bounds is not None:
        curve = curve.padded(*v_bounds)
    v_pts = [float(x) for x in curve.v_points_pu]

    def _key(args):
        """Cache key for a resolved capability curve.

        Returns:
            A Pyomo expression.
        """
        return args[0] if len(args) == 1 else tuple(args)

    tupled = isinstance(keys[0], tuple)
    idx = pyo.Set(initialize=keys, dimen=len(keys[0]) if tupled else 1)
    model.add_component(f"{name}_IDX", idx)

    # Start at nominal voltage, on the curve.  Beyond helping the solver,
    # this guarantees the auxiliaries carry values even when a decomposition
    # solver returns without writing every subproblem variable back.
    v_start = min(max(1.0, v_pts[0]), v_pts[-1])
    v_aux = pyo.Var(idx, bounds=(v_pts[0], v_pts[-1]), initialize=v_start)
    q_aux = pyo.Var(
        idx,
        initialize=lambda _m, *args: float(
            curve.step(v_start, variant_of(_key(args)))
        ),
    )
    model.add_component(f"{name}_v", v_aux)
    model.add_component(f"{name}_q_pu", q_aux)

    def _link_v(_m, *args):
        """Tie the voltage slacks to the voltage deviation.

        Returns:
            A Pyomo expression.
        """
        k = _key(args)
        return v_aux[k] == v_of(k)

    def _link_q(_m, *args):
        """Tie the reactive slacks to the reactive dispatch.

        Returns:
            A Pyomo expression.
        """
        k = _key(args)
        return q_of(k) == q_aux[k] * pn_of(k)

    # One of the few places `rule=` has to stay: the component names are
    # built from `name` at run time, and a decorator can only name a
    # component after the function it decorates. See
    # docs/contributing-pyomo.md.
    model.add_component(f"{name}_v_link", pyo.Constraint(idx, rule=_link_v))
    model.add_component(f"{name}_q_link", pyo.Constraint(idx, rule=_link_q))

    def _f_rule(_m, *args):
        """One piece of the piecewise characteristic.

        Returns:
            A Pyomo expression.
        """
        *key_parts, x = args
        return float(curve.step(x, variant_of(_key(tuple(key_parts)))))

    model.add_component(
        f"{name}_pw",
        pyo.Piecewise(
            keys,
            q_aux,
            v_aux,
            pw_pts={k: list(v_pts) for k in keys},
            pw_constr_type="EQ",
            f_rule=_f_rule,
            pw_repn=pw_repn,
        ),
    )
    return len(keys)


def warn_if_curve_leaves_pq_area(
    curve, pq_area, variant=0, v_range=None, p_range=None, context=""
) -> bool:
    """Warn when a Q(U) characteristic cannot coexist with the Q(P) area.

    The two express different things and are imposed together.  The Q(P)
    area *bounds* Q from active power; a Q(U) characteristic *assigns* Q
    from voltage.  Where the assigned value falls outside the bound there is
    no feasible reactive power at all, and the solver reports a plain
    infeasibility with nothing pointing at the cause.

    The conflict is real but conditional — it bites only when voltage sits
    far from the dead band while active power is low.  For VDE-AR-N 4110
    with its default dead band, the curve demands +0.484 Pn at 0.90 p.u.,
    which the Q(P) area permits only at rated output.

    Args:
        curve: The :class:`QVCurve` being imposed.
        pq_area: The :class:`Envelope` bounding Q against ``P/Pn``.
        variant: Variant index to check.
        v_range: ``(v_min, v_max)`` the buses may reach.  Defaults to the
            curve's own span.
        p_range: ``(p_min, p_max)`` in ``P/Pn``.  Defaults to
            :data:`DEFAULT_P_RANGE_PU`.
        context: Prefix identifying the caller in the warning.

    Returns:
        ``True`` if a warning was emitted.
    """
    p_range = DEFAULT_P_RANGE_PU if p_range is None else p_range
    if v_range is None:
        v_range = (
            float(curve.v_points_pu[0]),
            float(curve.v_points_pu[-1]),
        )
    vs = np.linspace(float(v_range[0]), float(v_range[1]), 201)
    ps = np.linspace(float(p_range[0]), float(p_range[1]), 201)
    hi = np.array(
        [
            min(m * p + b for m, b in pq_area.upper_pieces(variant, p_range))
            for p in ps
        ]
    )
    lo = np.array(
        [
            max(m * p + b for m, b in pq_area.lower_pieces(variant, p_range))
            for p in ps
        ]
    )

    worst_v, worst_p = None, 0.0
    for v in vs:
        q = float(curve.step(v, variant))
        ok = (lo <= q + 1e-12) & (q - 1e-12 <= hi)
        if not ok.any():
            worst_v, worst_p = float(v), float("inf")
            break
        need = float(ps[ok][0])
        if need > worst_p:
            worst_v, worst_p = float(v), need
    if worst_v is None or worst_p <= float(p_range[0]) + 1e-12:
        return False

    prefix = f"{context}: " if context else ""
    where = (
        "no active power at all satisfies it"
        if worst_p == float("inf")
        else f"only at P >= {worst_p:.3f} Pn"
    )
    warnings.warn(
        f"{prefix}the Q(U) characteristic and the Q(P) capability area are "
        f"both imposed, and they disagree over part of the operating range. "
        f"At v = {worst_v:.4f} p.u. the characteristic assigns "
        f"Q = {float(curve.step(worst_v, variant)):+.4f} Pn, which the Q(P) "
        f"area permits {where}. If a bus reaches that voltage at lower "
        f"active power the model is infeasible, with nothing in the solver "
        f"output pointing here. Widen the dead band, raise the sgens' "
        f"minimum active power, or drop qu_deadband and use the Q(U) area.",
        QuCurveOutsidePqAreaWarning,
        stacklevel=3,
    )
    return True


def warn_if_outside_exact_range(
    envelope: Envelope, x_min, x_max, context=""
) -> bool:
    """Warn when the model can drive an envelope past its exact range.

    The linear-inequality form the Pyomo constraints use reproduces the
    envelope only between its outermost breakpoints; beyond them the bound
    is relaxed to its hull, so the model permits a little more reactive
    power than the standard requires.  Bus voltage limits inside the span
    — the usual 0.9/1.1 against a QV area reaching 0.87/1.15 — never trigger
    this.

    Args:
        envelope: The capability envelope about to be constrained.
        x_min: Smallest value the driving quantity can take.
        x_max: Largest value it can take.
        context: Prefix identifying the caller in the warning.

    Returns:
        ``True`` if a warning was emitted.
    """
    lo, hi = envelope.exact_range()
    below = x_min < lo - 1e-9
    above = x_max > hi + 1e-9
    if not (below or above):
        return False
    prefix = f"{context}: " if context else ""
    warnings.warn(
        f"{prefix}the operating range [{x_min:.4f}, {x_max:.4f}] leaves the "
        f"capability envelope's exact range [{lo:.4f}, {hi:.4f}]. Outside it "
        f"the bound is relaxed to its hull, so the model permits somewhat "
        f"more reactive power than the standard requires. Tighten the "
        f"operating bounds, or extend the envelope with an explicit "
        f"breakpoint at the operating limit, to keep it exact.",
        EnvelopeRangeWarning,
        stacklevel=3,
    )
    return True


def compute_q_curves(grid_code=None) -> pd.DataFrame:
    """Return the sloped-segment parameters of the Q(P) and Q(U) bounds.

    .. deprecated::
       Superseded by :attr:`GridCode.pq_area` and :attr:`GridCode.qv_area`,
       which carry the saturation shelves this representation cannot express.
       Kept for callers written against the pre-0.4.1 registry.

    Args:
        grid_code: Grid code to evaluate, as accepted by
            :func:`resolve_grid_code`.

    Returns a :class:`~pandas.DataFrame` indexed by variant with columns:

    * ``m_qv``      — Q(U) ramp slope (ΔQ/Pn per p.u. voltage)
    * ``b_qv_min``  — intercept of the lower Q(U) ramp
    * ``b_qv_max``  — intercept of the upper Q(U) ramp
    * ``m_qp_max``  — Q(P) upper ramp slope
    * ``b_qp_max``  — Q(P) upper ramp intercept
    * ``m_qp_min``  — Q(P) lower ramp slope
    * ``b_qp_min``  — Q(P) lower ramp intercept

    All Q values are normalised to the installed active power Pn.

    Warning:
        These describe the **sloped segment only**.  Evaluating them outside
        that segment extrapolates without limit — at nominal voltage the
        Q(U) pair spans roughly three times the grid-code range, and the
        Q(P) pair reaches several times Pn at full output.  This is what
        potpourri did before 0.4.1.  Use
        :meth:`GridCode.q_flexibility` for the actual capability.
    """
    code = resolve_grid_code(grid_code)

    def _steepest(pieces):
        """The ramp piece — the one the old two-point form encoded."""
        return min(pieces, key=lambda mb: mb[0])

    variants = range(code.n_variants)
    qv_lower = [_steepest(code.qv_area.lower_pieces(v)) for v in variants]
    qv_upper = [_steepest(code.qv_area.upper_pieces(v)) for v in variants]
    pq_upper = [
        max(code.pq_area.upper_pieces(v), key=lambda mb: mb[0])
        for v in variants
    ]
    pq_lower = [
        min(code.pq_area.lower_pieces(v), key=lambda mb: mb[0])
        for v in variants
    ]

    return pd.DataFrame(
        {
            "m_qv": [m for m, _ in qv_lower],
            "b_qv_min": [b for _, b in qv_lower],
            "b_qv_max": [b for _, b in qv_upper],
            "m_qp_max": [m for m, _ in pq_upper],
            "m_qp_min": [m for m, _ in pq_lower],
            "b_qp_max": [b for _, b in pq_upper],
            "b_qp_min": [b for _, b in pq_lower],
        }
    )
