"""Wind power mix-in: attaches wind generator sets, parameters, and
Q-control constraints to a multi-period model."""

import numpy as np
import pandas as pd
import pyomo.environ as pyo
from loguru import logger
from potpourri.technologies.q_control import (
    DEFAULT_GRID_CODE,
    DEFAULT_P_RANGE_PU,
    DEFAULT_WIND_SGEN_TYPES,
    bus_voltage_range,
    check_var_q,
    compute_q_curves,
    resolve_grid_code,
)
from potpourri.technologies.sgens import Sgens_multi_period


# ---------------------------------------------------------------------------
# Grid-code Q-curve parameters come from the registry in
# potpourri.technologies.q_control, not from a private copy.  This module used
# to carry its own capability table and its own reimplementation of the
# Q-curve maths, which meant `grid_code` never reached the wind or
# hosting-capacity paths and the numbers could drift from the registry
# silently.
#
# The simplified hosting-capacity check uses the widest envelope the selected
# code offers: the largest capacitive entry and the most negative inductive
# one.  For the default VDE-AR-N 4120 that is +0.484322 / -0.410775.
# ---------------------------------------------------------------------------


def _hc_q_bounds(code):
    """Return (max, min) Q/P for the simplified HC check of ``code``."""
    q = code.vqu_q_max
    return float(q[0].max()), float(q[1].min())


# Defaults for the HC keyword arguments below.  Derived from the default
# grid code so they track the registry: +0.484322 / -0.410775 for VDE-AR-N
# 4120, the widest capacitive and widest inductive variant it offers.
_DEFAULT_HC_Q_MAX, _DEFAULT_HC_Q_MIN = _hc_q_bounds(DEFAULT_GRID_CODE)


class Windpower_multi_period(Sgens_multi_period):
    """Multi-period wind generator device module, extending sgen with
    Q-control and hosting-capacity (HC) support.

    The Q-control constraints implement the capability area of the selected
    grid code (see :mod:`potpourri.technologies.q_control`), defaulting to
    VDE-AR-N 4120.

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
            HC grid-code Q-P constraint.  Default: 0.484322, the widest
            capacitive variant of VDE-AR-N 4120.
        qp_min: Minimum Q/P ratio (inductive, negative) for the simplified
            HC grid-code Q-P constraint.  Default: -0.410775, the widest
            inductive variant of the same code.
    """

    def __init__(
        self,
        net,
        T=None,
        scenario=None,
        *,
        sw_max_mva: float = 10_000.0,
        sw_min_mva: float = 0.0,
        qp_max: float = _DEFAULT_HC_Q_MAX,
        qp_min: float = _DEFAULT_HC_Q_MIN,
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
        # Match every SimBench wind spelling, not just the HV one; see
        # DEFAULT_WIND_SGEN_TYPES. Shared with the single-period path so the
        # two cannot drift apart.
        wind_types = getattr(self, "wind_sgen_types", DEFAULT_WIND_SGEN_TYPES)
        model.WIND = model.WIND_HC | pyo.Set(
            within=model.sG,
            initialize=self.static_generation_data.index[
                self.static_generation_data["type"].isin(wind_types)
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

        # Q(U) slopes from the selected grid code's characteristic.
        # Slope from low-voltage to high-voltage point; intercepts at V3 and V1
        code = resolve_grid_code(getattr(self, "grid_code", None))
        x = code.vqu_v_points
        y = code.vqu_q_max
        hc_max, _hc_min = _hc_q_bounds(code)
        # Narrowest capacitive / widest inductive variant is the last column.
        last = y.shape[1] - 1
        self.m_qu_max = (hc_max + abs(y[1, 0])) / (x[0, 0] - x[1, 0])
        self.qu_max = -self.m_qu_max * x[1, 0] + hc_max
        self.m_qu_min = (abs(y[0, last]) + abs(y[1, last])) / (
            x[0, 0] - x[1, 0]
        )
        self.qu_min = -self.m_qu_min * x[0, 0] + y[0, last]
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

    def static_generation_wind_var_q(self, net, grid_code=None):
        """Populate ``static_generation_data`` Q limits from the grid code.

        The capability is parameterised by operating variants (the ``var_q``
        column in ``net.sgen``); each selects one column of the selected grid
        code's capability table.  How many variants exist depends on the
        code — VDE-AR-N 4105 defines two, 4110 one and 4120 three — so
        ``var_q`` is validated against the resolved code.
        """
        code = resolve_grid_code(grid_code)
        self.grid_code = code
        self.q_limit_parameter = compute_q_curves(code)
        # Q/Pn capability table: row 0 capacitive, row 1 inductive;
        # columns are the var_q variants.
        y = code.vqu_q_max

        if "var_q" in self.net.sgen:
            self.static_generation_data["var_q"] = self.net.sgen.var_q.values
            sgens_var_q = self.static_generation_data.index[
                self.static_generation_data.var_q.notna()
            ]
            check_var_q(
                self.static_generation_data.var_q[sgens_var_q],
                code,
                context="net.sgen.var_q",
            )

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
                p_inst[sgens_var_q] * code.qp_p_high
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

        # Each grid-code bound is a piecewise-linear envelope, so it becomes
        # one inequality per affine piece: the upper bound is the pointwise
        # minimum of its pieces, the lower bound the pointwise maximum.  A
        # single line cannot express the saturation shelf.
        pq_area = self.grid_code.pq_area
        qv_area = self.grid_code.qv_area
        v_span = bus_voltage_range(self.net)
        sGbs_lookup = {g: b for (g, b) in model.sGbs}
        model.W_QP_PIECE = pyo.RangeSet(
            0, pq_area.max_pieces(DEFAULT_P_RANGE_PU) - 1
        )
        model.W_QU_PIECE = pyo.RangeSet(0, qv_area.max_pieces(v_span) - 1)

        @model.Constraint(model.WINDc, model.W_QP_PIECE)
        def QW_pos(model, w, k):
            pieces = pq_area.upper_pieces(
                int(pyo.value(model.var_q[w])), DEFAULT_P_RANGE_PU
            )
            if k >= len(pieces):
                return pyo.Constraint.Skip
            m, b = pieces[k]
            return model.qsG[w] <= m * model.psG[w] + b * model.PsG_inst[w]

        @model.Constraint(model.WINDc, model.W_QP_PIECE)
        def QW_neg(model, w, k):
            pieces = pq_area.lower_pieces(
                int(pyo.value(model.var_q[w])), DEFAULT_P_RANGE_PU
            )
            if k >= len(pieces):
                return pyo.Constraint.Skip
            m, b = pieces[k]
            return model.qsG[w] >= m * model.psG[w] + b * model.PsG_inst[w]

        @model.Constraint(model.WINDc, model.W_QU_PIECE)
        def QV_min(model, w, k):
            if w not in sGbs_lookup:
                return pyo.Constraint.Skip
            pieces = qv_area.lower_pieces(
                int(pyo.value(model.var_q[w])), v_span
            )
            if k >= len(pieces):
                return pyo.Constraint.Skip
            b_bus = sGbs_lookup[w]
            m, b = pieces[k]
            return model.qsG[w] >= (m * model.v[b_bus] + b) * model.PsG_inst[w]

        @model.Constraint(model.WINDc, model.W_QU_PIECE)
        def QV_max(model, w, k):
            if w not in sGbs_lookup:
                return pyo.Constraint.Skip
            pieces = qv_area.upper_pieces(
                int(pyo.value(model.var_q[w])), v_span
            )
            if k >= len(pieces):
                return pyo.Constraint.Skip
            b_bus = sGbs_lookup[w]
            m, b = pieces[k]
            return model.qsG[w] <= (m * model.v[b_bus] + b) * model.PsG_inst[w]

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

        # Simplified HC Q-P bounds: the widest band the grid code offers
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
