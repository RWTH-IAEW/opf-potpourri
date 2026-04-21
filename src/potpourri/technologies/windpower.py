"""Wind power mix-in: attaches wind generator sets, parameters, and
Q-control constraints to a multi-period model."""

import numpy as np
import pandas as pd
import pyomo.environ as pyo
from loguru import logger
from potpourri.technologies.sgens import Sgens_multi_period


# ---------------------------------------------------------------------------
# VDE-AR-N 4105 / BDEW grid-code Q-curve constants (medium-voltage, 110 kV)
#
# Voltage operating points [p.u. on 110 kV base]:
#   V1 = 96/110,  V2 = 103/110,  V3 = 120/110,  V4 = 127/110
#
# Q/P limits at V1 and V3 for three operating variants:
#   variant 0: Q/P_max = 0.48 @ V1,  Q/P_max = 0.33 @ V2,
#             Q/P_min = -0.41 @ V1
#   variant 1: Q/P_max = 0.41 @ V1,  etc.
#   variant 2: Q/P_max = 0.33 @ V1,  etc.
# ---------------------------------------------------------------------------

# Normalised voltage points [p.u.]: [[V1, V2], [V3, V4]]
_VQU_V_POINTS = np.array([[96, 103], [120, 127]]) / 110.0

# Q/P table: rows = voltage level (low, high), cols = variant (0, 1, 2)
_VQU_Q_MAX = np.array([[0.48, 0.41, 0.33], [-0.23, -0.33, -0.41]])

# Q/P boundaries at the P-axis breakpoints (0.1 Pn and 0.2 Pn)
_QP_P_BREAK_HIGH = 0.1  # upper P/Pn at which Q/P characteristic changes
_QP_P_BREAK_LOW = 0.2  # lower P/Pn at which Q/P characteristic changes (abs)

# Default Q/P ratio bounds for the simplified HC grid-code check
_QP_HC_MAX = 0.48  # maximum Q/P (capacitive) – variant 0 at low voltage
_QP_HC_MIN = -0.41  # minimum Q/P (inductive) – variant 0 at low voltage


class Windpower_multi_period(Sgens_multi_period):
    """Multi-period wind generator device module, extending sgen with
    Q-control and hosting-capacity (HC) support.

    The Q-control constraints implement the VDE-AR-N 4105 / BDEW grid-code
    voltage-reactive power characteristic for medium-voltage wind generators.

    Args:
        net: pandapower network.  If ``net.bus`` contains a ``windpot_p_mw``
            column, it is used as the active-power upper bound for HC
            generators.
        T: Number of time steps.
        scenario: Unused for wind (no penetration scenario); reserved for
            interface compatibility.
        sw_max_mva: Default apparent-power upper limit per HC wind generator
            (MVA).  Can be overridden per network via
            ``_calc_wind_opf_parameters``.
        sw_min_mva: Minimum apparent power for an *active* HC wind generator
            (MVA).  A generator with ``y=1`` must carry at least this much
            apparent power.
        qp_max: Maximum Q/P ratio (capacitive, positive) for the simplified
            HC grid-code Q-P constraint.  Default: 0.48 (VDE-AR-N 4105,
            variant 0).
        qp_min: Minimum Q/P ratio (inductive, negative) for the simplified
            HC grid-code Q-P constraint.  Default: -0.41 (VDE-AR-N 4105,
            variant 0).
    """

    def __init__(
        self,
        net,
        T=None,
        scenario=None,
        *,
        sw_max_mva: float = 10_000.0,
        sw_min_mva: float = 0.0,
        qp_max: float = _QP_HC_MAX,
        qp_min: float = _QP_HC_MIN,
    ):
        super().__init__(net, T, scenario)
        self._sw_max_mva = sw_max_mva
        self._sw_min_mva = sw_min_mva
        self.qp_max = qp_max
        self.qp_min = qp_min

        if "windpot_p_mw" in net.bus:
            self.static_generation_data["windpot"] = net.bus.windpot_p_mw[
                net.sgen.bus.values
            ].values
            self.static_generation_data["type"] = net.sgen.type.values

    def get_all(self, model):
        """No-op: wind generators are initialised via get_all_opf."""

    def get_all_opf(self, model):
        """Attach OPF sets and parameters for controllable wind generators."""
        self.get_opf_sets(model)
        self.get_opf_parameters(model)

    def get_opf_sets(self, model):
        """Define WIND_HC, WIND, and WINDc sets."""
        model.WIND_HC = pyo.Set(
            within=model.sG,
            initialize=self.static_generation_data.index[
                self.static_generation_data["wind_hc"]
                & self.static_generation_data.in_service
            ],
        )
        model.WIND = model.WIND_HC | pyo.Set(
            within=model.sG,
            initialize=self.static_generation_data.index[
                (self.static_generation_data["type"] == "Wind")
                & self.static_generation_data.in_service
            ],
        )
        model.WINDc = (
            model.WIND
            & model.sGc
            & pyo.Set(
                initialize=self.static_generation_data.index[
                    self.static_generation_data["var_q"].values != None  # noqa: E711
                ],
            )
        )
        return True

    def _calc_wind_opf_parameters(
        self,
        model,
        sw_max_mva: float | None = None,
        sw_min_mva: float | None = None,
    ):
        """Compute SWmax/SWmin apparent-power limits and Q(U) slope
        parameters for HC wind generators.

        Args:
            model: Pyomo model (requires ``model.WIND_HC`` and ``model.T``).
            sw_max_mva: Override the instance-level ``sw_max_mva`` for this
                call.  Defaults to the value passed at construction.
            sw_min_mva: Override the instance-level ``sw_min_mva`` for this
                call.  Defaults to the value passed at construction.
        """
        if sw_max_mva is None:
            sw_max_mva = self._sw_max_mva
        if sw_min_mva is None:
            sw_min_mva = self._sw_min_mva

        if "windpot_p_mw" in self.net.bus:
            self.static_generation_data["windpot"] = self.net.bus.windpot_p_mw[
                self.net.sgen.bus.values
            ].values

        wind_hc_set = np.arange(len(self.net.sgen))[
            self.net.sgen.wind_hc & self.net.sgen.in_service
        ]
        self.SWmax_data = pd.Series(sw_max_mva / self.baseMVA, wind_hc_set)
        self.SWmax_data_dict, self.SWmax_tuple = self.make_to_dict(
            model.WIND_HC, model.T, self.SWmax_data
        )
        self.SWmin_data = pd.Series(sw_min_mva / self.baseMVA, wind_hc_set)
        self.SWmin_data_dict, self.SWmin_tuple = self.make_to_dict(
            model.WIND_HC, model.T, self.SWmin_data
        )

        # Q(U) slopes from VDE-AR-N 4105 grid-code characteristic (variant 1)
        # Slope from low-voltage to high-voltage point; intercepts at V3 and V1
        self.m_qu_max = (
            (_QP_HC_MAX + abs(_VQU_Q_MAX[1, 0]))
            / (_VQU_V_POINTS[0, 0] - _VQU_V_POINTS[1, 0])
            * 1.0  # already normalised
        )
        self.qu_max = -self.m_qu_max * _VQU_V_POINTS[1, 0] + _QP_HC_MAX
        self.m_qu_min = (abs(_VQU_Q_MAX[0, 2]) + abs(_VQU_Q_MAX[1, 2])) / (
            _VQU_V_POINTS[0, 0] - _VQU_V_POINTS[1, 0]
        )
        self.qu_min = -self.m_qu_min * _VQU_V_POINTS[0, 0] + _VQU_Q_MAX[0, 2]
        return True

    def get_hc_acopf_parameters(self, model, net):
        """Attach SWmax, SWmin, and optional pWmax parameters for HC-ACOPF."""
        model.SWmax = pyo.Param(
            self.SWmax_tuple, initialize=self.SWmax_data_dict, mutable=True
        )
        model.SWmin = pyo.Param(
            self.SWmin_data_dict, initialize=self.SWmin_data_dict, mutable=True
        )

        if "windpot_p_mw" in self.net.bus:
            self.Windpot_data_dict, self.Windpot_tuple = self.make_to_dict(
                model.WIND_HC, model.T, self.static_generation_data["windpot"]
            )
            model.pWmax = pyo.Param(
                self.Windpot_tuple,
                initialize=self.Windpot_data_dict,
                mutable=True,
            )
        return True

    def get_hc_acopf_variables(self, model):
        """Attach binary HC placement variable y for each wind generator."""
        model.y = pyo.Var(
            self.Windpot_tuple, within=pyo.Binary, initialize=1.0
        )
        return True

    def get_opf_parameters(self, model):
        """Attach var_q and PsG_inst parameters for controllable wind
        generators."""
        model.var_q = pyo.Param(
            model.WINDc,
            model.T,
            initialize=self.static_generation_data["var_q"][model.WINDc],
        )
        model.PsG_inst = pyo.Param(
            model.WINDc,
            model.T,
            initialize=self.static_generation_data["p_inst"][model.WINDc],
        )
        return True

    def static_generation_wind_var_q(self, net):
        """Compute Q(P) and Q(U) characteristic slopes from VDE-AR-N 4105
        and populate ``static_generation_data`` Q limits.

        The characteristic is parameterised by three operating variants
        (``var_q`` column in ``net.sgen``).  Each variant selects a different
        Q/P intercept from the ``_VQU_Q_MAX`` table.
        """
        x = _VQU_V_POINTS  # shape (2, 2): [[V1, V2], [V3, V4]]
        y = _VQU_Q_MAX  # shape (2, 3): rows = voltage level, cols = variant

        m = (y[1] - y[0]) / (x[0, 1] - x[0, 0])
        b = np.array([y[0] - m * x[i, 0] for i in range(len(x))]).T

        p_range = _QP_P_BREAK_HIGH - _QP_P_BREAK_LOW
        m_qp_max = (_QP_P_BREAK_HIGH - y[0]) / p_range
        m_qp_min = (-_QP_P_BREAK_HIGH - y[1]) / p_range
        b_qp_max = _QP_P_BREAK_HIGH - m_qp_max * _QP_P_BREAK_HIGH
        b_qp_min = -_QP_P_BREAK_HIGH - m_qp_min * _QP_P_BREAK_HIGH

        self.q_limit_parameter = pd.DataFrame(
            {
                "m_qv": m,
                "b_qv_min": b[:, 0],
                "b_qv_max": b[:, 1],
                "m_qp_max": m_qp_max,
                "m_qp_min": m_qp_min,
                "b_qp_max": b_qp_max,
                "b_qp_min": b_qp_min,
            },
        )

        if "var_q" in self.net.sgen:
            self.static_generation_data["var_q"] = self.net.sgen.var_q.values
            sgens_var_q = self.static_generation_data.index[
                self.static_generation_data.var_q.notna()
            ]

            try:
                p_inst = self.net.sgen.p_inst_mw.values / self.baseMVA
            except AttributeError:
                logger.warning(
                    "No p_inst_mw attribute found in net.sgen. "
                    "Using p_mw as p_inst for wind generator power limits."
                )
                p_inst = self.static_generation_data["p"]

            self.static_generation_data["p_inst"] = p_inst

            self.static_generation_data["max_q"][sgens_var_q] = [
                y[0, int(self.static_generation_data.var_q[g])]
                * self.static_generation_data["p_inst"][g]
                for g in sgens_var_q
            ]
            self.static_generation_data["min_q"][sgens_var_q] = [
                y[1, int(self.static_generation_data.var_q[g])]
                * self.static_generation_data["p_inst"][g]
                for g in sgens_var_q
            ]

            self.static_generation_data["max_p"][sgens_var_q] = p_inst[
                sgens_var_q
            ]
            self.static_generation_data["min_p"][sgens_var_q] = (
                p_inst[sgens_var_q] * _QP_P_BREAK_HIGH
            )

        else:
            self.static_generation_data["var_q"] = None
            self.static_generation_data["p_inst"] = None

        if "wind_hc" in self.net.sgen:
            self.static_generation_data["wind_hc"] = (
                self.net.sgen.wind_hc.values
            )
        else:
            self.static_generation_data["wind_hc"] = False

    def get_objective(self, model):
        """Add a wind-maximisation objective that subtracts line losses."""

        def obj_wind_loss_rule(model):
            return (
                sum(model.psG[w] for w in model.WIND_HC)
                - sum(model.pLfrom[l] + model.pLto[l] for l in model.L)
                - sum(model.pThv[t] + model.pTlv[t] for t in model.TRANSF)
            )

        model.obj = pyo.Objective(rule=obj_wind_loss_rule, sense=pyo.maximize)

    def get_constraints(self, model, net):
        """Add Q(P) and Q(U) constraints for controllable wind and HC
        generators."""

        @model.Constraint(model.WINDc)
        def QW_pos(model, w):
            return (
                model.qsG[w]
                <= self.q_limit_parameter.b_qp_max[model.var_q[w]]
                * model.PsG_inst[w]
                + self.q_limit_parameter.m_qp_max[model.var_q[w]]
                * model.psG[w]
            )

        @model.Constraint(model.WINDc)
        def QW_neg(model, w):
            return (
                model.qsG[w]
                >= self.q_limit_parameter.b_qp_min[model.var_q[w]]
                * model.PsG_inst[w]
                + self.q_limit_parameter.m_qp_min[model.var_q[w]]
                * model.psG[w]
            )

        @model.Constraint(model.WINDc)
        def QV_min(model, w):
            for g, b in model.sGbs:
                if g == w:
                    return (
                        model.qsG[w]
                        >= (
                            self.q_limit_parameter.m_qv[model.var_q[w]]
                            * model.v[b]
                            + self.q_limit_parameter.b_qv_min[model.var_q[w]]
                        )
                        * model.PsG_inst[w]
                    )

        @model.Constraint(model.WINDc)
        def QV_max(model, w):
            for g, b in model.sGbs:
                if g == w:
                    return (
                        model.qsG[w]
                        <= (
                            self.q_limit_parameter.m_qv[model.var_q[w]]
                            * model.v[b]
                            + self.q_limit_parameter.b_qv_max[model.var_q[w]]
                        )
                        * model.PsG_inst[w]
                    )

        @model.Constraint(model.WIND_HC)
        def SW_max(model, w):
            return (
                model.psG[w] ** 2 + model.qsG[w] ** 2
                <= model.SWmax[w] ** 2 * model.y[w]
            )

        @model.Constraint(model.WIND_HC)
        def SW_min(model, w):
            return (
                model.psG[w] ** 2 + model.qsG[w] ** 2
                >= model.SWmin[w] ** 2 * model.y[w]
            )

        # Simplified HC Q-P bounds (VDE-AR-N 4105, variant 0)
        @model.Constraint(model.WIND_HC)
        def QW_min(model, w):
            return model.qsG[w] >= self.qp_min * model.psG[w]

        @model.Constraint(model.WIND_HC)
        def QW_max(model, w):
            return model.qsG[w] <= self.qp_max * model.psG[w]

        @model.Constraint(model.WIND_HC)
        def QU_min_hc(model, w):
            for g, b in model.sGbs:
                if g == w:
                    return (
                        model.qsG[w]
                        >= (self.m_qu_min * model.v[b] + self.qu_min)
                        * model.psG[w]
                    )

        @model.Constraint(model.WIND_HC)
        def QU_max_hc(model, w):
            for g, b in model.sGbs:
                if g == w:
                    return (
                        model.qsG[w]
                        <= (self.m_qu_max * model.v[b] + self.qu_max)
                        * model.psG[w]
                    )

        if "windpot_p_mw" in net.bus:

            @model.Constraint(model.WIND_HC)
            def PW_max(model, w):
                return model.psG[w] <= model.pWmax[w]

    def unfix_variables(self, model):
        """Unfix real and reactive power for all HC wind generators."""
        for w in model.WIND_HC:
            model.psG[w].unfix()
            model.qsG[w].unfix()

    def get_all_acopf(self, model):
        """No additional ACOPF components needed for wind (called via
        get_all_opf)."""
