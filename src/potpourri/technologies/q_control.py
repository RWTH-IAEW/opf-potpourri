"""VDE-AR-N 4105 / BDEW grid-code reactive-power control constants and
Q-curve computation.

Shared by :class:`~potpourri.technologies.pv.PV_multi_period`,
:class:`~potpourri.technologies.sgens.Sgens_multi_period`, and the
single-period :class:`~potpourri.models.ACOPF_base.ACOPF`.
"""

import numpy as np
import pandas as pd

# Normalised voltage breakpoints [p.u.]: [[V1, V2], [V3, V4]]
VQU_V_POINTS = np.array([[96, 103], [120, 127]]) / 110.0

# Q/P bounds table:
#   rows = [capacitive limit at low voltage, inductive limit at high voltage]
#   cols = variant [0, 1, 2]
VQU_Q_MAX = np.array([[0.48, 0.41, 0.33], [-0.23, -0.33, -0.41]])

# P/Pn breakpoints for the Q(P) piecewise characteristic
QP_P_HIGH = 0.1
QP_P_LOW = 0.2

# P(U) curtailment thresholds (VDE-AR-N 4105 §8.5)
VPU_V_CURTAIL = 1.06  # voltage above which curtailment begins (p.u.)
VPU_V_MAX = 1.10  # voltage at which active output reaches zero (p.u.)

# cos(φ)(P) profile: P/Pn below which Q = 0 (VDE-AR-N 4105)
CPP_P_THRESHOLD_PU = 0.2


def compute_q_curves() -> pd.DataFrame:
    """Return VDE-AR-N 4105 Q(P) and Q(U) slope/intercept parameters.

    Computes a linearised piecewise-linear approximation of the grid-code
    reactive-power capability envelope.

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
    """
    x = VQU_V_POINTS
    y = VQU_Q_MAX
    m = (y[1] - y[0]) / (x[0, 1] - x[0, 0])
    b = np.array([y[0] - m * x[i, 0] for i in range(len(x))]).T

    p_range = QP_P_HIGH - QP_P_LOW
    m_qp_max = (QP_P_HIGH - y[0]) / p_range
    m_qp_min = (-QP_P_HIGH - y[1]) / p_range
    b_qp_max = QP_P_HIGH - m_qp_max * QP_P_HIGH
    b_qp_min = -QP_P_HIGH - m_qp_min * QP_P_HIGH

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
