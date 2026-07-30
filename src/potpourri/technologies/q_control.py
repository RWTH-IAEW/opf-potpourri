"""Grid-code reactive-power control parameters and Q-curve computation.

The German technical connection rules (TAR) are represented as selectable
:class:`GridCode` parameter sets rather than hard-coded constants, so a study
can target the rule that applies to its voltage level:

* :data:`VDE_AR_N_4105` — low voltage.
* :data:`VDE_AR_N_4110` — medium voltage.  **Provisional**, see below.

Select one by short name (``"4105"``, ``"4110"``) or by passing the
:class:`GridCode` itself.  The default is :data:`VDE_AR_N_4105`, so existing
callers are unaffected.

Consumed by :class:`~potpourri.technologies.pv.PV_multi_period`,
:class:`~potpourri.technologies.sgens.Sgens_multi_period` and the
single-period :class:`~potpourri.models.ACOPF_base.ACOPF`.

.. warning::
   :data:`VDE_AR_N_4110` currently carries the **VDE-AR-N 4105 values as a
   placeholder**.  Its normative medium-voltage parameters have not been
   entered yet, so results obtained with ``grid_code="4110"`` are *not*
   4110-compliant.  Selecting it emits a
   :class:`ProvisionalGridCodeWarning`.
"""

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd


class ProvisionalGridCodeWarning(UserWarning):
    """Raised when a grid code whose parameters are placeholders is used."""


@dataclass(frozen=True, eq=False)
class GridCode:
    """Reactive-power capability parameters for one grid code.

    Attributes:
        name: Short identifier, e.g. ``"4105"``.
        title: Full designation, e.g. ``"VDE-AR-N 4105"``.
        voltage_level: Voltage level the rule applies to.
        vqu_v_points: Q(U) voltage breakpoints [p.u.], shape (2, 2) as
            ``[[V1, V2], [V3, V4]]``.
        vqu_q_max: Q/Pn capability table, shape (2, n_variants).  Row 0 is
            the capacitive limit at low voltage, row 1 the inductive limit
            at high voltage; columns are the variants selected by
            ``net.sgen.var_q``.
        qp_p_high: Lower P/Pn breakpoint of the Q(P) characteristic.
        qp_p_low: Upper P/Pn breakpoint of the Q(P) characteristic.
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
    vqu_v_points: np.ndarray
    vqu_q_max: np.ndarray
    qp_p_high: float
    qp_p_low: float
    vpu_v_curtail: float
    vpu_v_max: float
    cpp_p_threshold_pu: float
    provisional: bool = False
    provisional_note: str = ""

    @property
    def n_variants(self) -> int:
        """Number of ``var_q`` variants this grid code defines."""
        return int(self.vqu_q_max.shape[1])


# --- VDE-AR-N 4105 (low voltage) -----------------------------------------
# Normative values.
VDE_AR_N_4105 = GridCode(
    name="4105",
    title="VDE-AR-N 4105",
    voltage_level="low voltage",
    # Normalised voltage breakpoints [p.u.]: [[V1, V2], [V3, V4]]
    vqu_v_points=np.array([[96, 103], [120, 127]]) / 110.0,
    # Rows: [capacitive limit at low voltage, inductive limit at high
    # voltage]; columns: variant [0, 1, 2]
    vqu_q_max=np.array([[0.48, 0.41, 0.33], [-0.23, -0.33, -0.41]]),
    qp_p_high=0.1,
    qp_p_low=0.2,
    vpu_v_curtail=1.06,  # VDE-AR-N 4105 section 8.5
    vpu_v_max=1.10,
    cpp_p_threshold_pu=0.2,
)


# --- VDE-AR-N 4110 (medium voltage) --------------------------------------
# PLACEHOLDER: these are the VDE-AR-N 4105 values, not 4110 figures.  The
# medium-voltage parameters still have to be taken from the standard (and,
# where the standard leaves them to the network operator, from the operator's
# parameterisation).  Until then, selecting this grid code warns.
_PROVISIONAL_4110_NOTE = (
    "VDE-AR-N 4110 currently reuses the VDE-AR-N 4105 (low-voltage) "
    "parameters as a placeholder. Its normative medium-voltage values have "
    "not been entered yet, so these results are NOT 4110-compliant. Replace "
    "the vqu_v_points / vqu_q_max / qp_* fields of VDE_AR_N_4110 in "
    "potpourri.technologies.q_control before relying on them."
)

VDE_AR_N_4110 = GridCode(
    name="4110",
    title="VDE-AR-N 4110",
    voltage_level="medium voltage",
    vqu_v_points=VDE_AR_N_4105.vqu_v_points,
    vqu_q_max=VDE_AR_N_4105.vqu_q_max,
    qp_p_high=VDE_AR_N_4105.qp_p_high,
    qp_p_low=VDE_AR_N_4105.qp_p_low,
    vpu_v_curtail=VDE_AR_N_4105.vpu_v_curtail,
    vpu_v_max=VDE_AR_N_4105.vpu_v_max,
    cpp_p_threshold_pu=VDE_AR_N_4105.cpp_p_threshold_pu,
    provisional=True,
    provisional_note=_PROVISIONAL_4110_NOTE,
)


GRID_CODES = {
    VDE_AR_N_4105.name: VDE_AR_N_4105,
    VDE_AR_N_4110.name: VDE_AR_N_4110,
}

DEFAULT_GRID_CODE = VDE_AR_N_4105


# --- Backwards-compatible module-level constants -------------------------
# These mirror the default (4105) grid code so that existing imports keep
# working.  New code should read the fields off a GridCode instead.
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
            as ``"4105"`` / ``"4110"`` (``"VDE-AR-N 4110"`` is also
            accepted), or a :class:`GridCode` instance, which is returned
            unchanged.

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


def compute_q_curves(grid_code=None) -> pd.DataFrame:
    """Return Q(P) and Q(U) slope/intercept parameters for a grid code.

    Computes a linearised piecewise-linear approximation of the grid-code
    reactive-power capability envelope.

    Args:
        grid_code: Grid code to evaluate, as accepted by
            :func:`resolve_grid_code`.  Defaults to :data:`VDE_AR_N_4105`.

    Returns a :class:`~pandas.DataFrame` indexed by variant (0, 1, 2) with
    columns:

    * ``m_qv``      — Q(U) slope (ΔQ/Pn per p.u. voltage)
    * ``b_qv_min``  — Q(U) lower intercept / Pn (capacitive, at V1)
    * ``b_qv_max``  — Q(U) upper intercept / Pn (inductive, at V3)
    * ``m_qp_max``  — Q(P) upper-bound slope
    * ``b_qp_max``  — Q(P) upper-bound intercept
    * ``m_qp_min``  — Q(P) lower-bound slope
    * ``b_qp_min``  — Q(P) lower-bound intercept

    All Q values are normalised to the installed active power Pn.

    Note:
        The Q(P) bounds are a single linear segment and are not clipped
        above ``qp_p_low``, so the envelope keeps widening beyond the
        reference point.  The Q(P) constraint is therefore only binding at
        low active power; the inverter S² circle and the cos(phi) cone
        provide the limit at high output.
    """
    code = resolve_grid_code(grid_code)

    x = code.vqu_v_points
    y = code.vqu_q_max
    m = (y[1] - y[0]) / (x[0, 1] - x[0, 0])
    b = np.array([y[0] - m * x[i, 0] for i in range(len(x))]).T

    p_range = code.qp_p_high - code.qp_p_low
    m_qp_max = (code.qp_p_high - y[0]) / p_range
    m_qp_min = (-code.qp_p_high - y[1]) / p_range
    b_qp_max = code.qp_p_high - m_qp_max * code.qp_p_high
    b_qp_min = -code.qp_p_high - m_qp_min * code.qp_p_high

    return pd.DataFrame(
        {
            "m_qv": m,
            "b_qv_min": b[:, 0],
            "b_qv_max": b[:, 1],
            "m_qp_max": m_qp_max,
            "m_qp_min": m_qp_min,
            "b_qp_max": b_qp_max,
            "b_qp_min": b_qp_min,
        }
    )
