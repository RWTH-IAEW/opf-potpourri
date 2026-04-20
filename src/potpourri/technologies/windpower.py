"""Wind power mix-in: attaches wind generator sets, parameters, and
Q-control constraints to a multi-period model."""

import numpy as np
import pandas as pd
import pyomo.environ as pyo
from loguru import logger
from potpourri.technologies.sgens import Sgens_multi_period


# TODO: check if properly made multiperiod
class Windpower_multi_period(Sgens_multi_period):
    """Multi-period wind generator device module, extending sgen with
    Q-control and HC support."""

    def __init__(self, net, T=None, scenario=None):
        super().__init__(net, T, scenario)
        # taken from HC_ACOPF, formerly used in _calc_opf_parameters
        if "windpot_p_mw" in net.bus:
            self.static_generation_data["windpot"] = net.bus.windpot_p_mw[
                net.sgen.bus.values
            ].values
            # taken from ACOPF_base, formerly used in
            # get static_generation_reactive_power_limits()
            self.static_generation_data["type"] = net.sgen.type.values

    def get_all(self, model):
        """No-op placeholder; wind generators are initialised via
        get_all_opf."""

    def get_all_opf(self, model):
        """Attach OPF sets and parameters for controllable wind generators."""
        self.get_opf_sets(model)
        self.get_opf_parameters(model)

    def get_opf_sets(self, model):
        """Define WIND_HC, WIND, and WINDc sets for HC and Q-control
        formulations."""
        # generators for hc calculation
        model.WIND_HC = pyo.Set(
            within=model.sG,
            initialize=self.static_generation_data.index[
                self.static_generation_data["wind_hc"]
                & self.static_generation_data.in_service
            ],
        )
        # all wind generators
        model.WIND = model.WIND_HC | pyo.Set(
            within=model.sG,
            initialize=self.static_generation_data.index[
                (self.static_generation_data["type"] == "Wind")
                & self.static_generation_data.in_service
            ],
        )
        # controllable wind generators, not for hc calculation
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

    def _calc_wind_opf_parameters(self, model, SWmax=10000, SWmin=0):
        """Compute SWmax/SWmin limits and Q(U) slope parameters for HC wind
        generators."""
        if "windpot_p_mw" in self.net.bus:
            self.static_generation_data["windpot"] = self.net.bus.windpot_p_mw[
                self.net.sgen.bus.values
            ].values

        wind_hc_set = np.arange(len(self.net.sgen))[
            self.net.sgen.wind_hc & self.net.sgen.in_service
        ]
        self.SWmax_data = pd.Series(SWmax / self.baseMVA, wind_hc_set)
        self.SWmax_data_dict, self.SWmax_tuple = self.make_to_dict(
            model.WIND_HC, model.T, self.SWmax_data
        )
        self.SWmin_data = pd.Series(SWmin / self.baseMVA, wind_hc_set)
        self.SWmin_data_dict, self.SWmin_tuple = self.make_to_dict(
            model.WIND_HC, model.T, self.SWmin_data
        )

        self.m_qu_max = (0.48 + 0.23) / (96 - 103) * 110  # Variante 1
        self.qu_max = -self.m_qu_max * 120 / 110 + 0.48
        self.m_qu_min = (0.33 + 0.41) / (96 - 103) * 110  # Variante 3
        self.qu_min = -self.m_qu_min * 96 / 110 + 0.33
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
        """Attach binary HC variable y for each wind generator over time."""
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
        """Compute Q(P) and Q(U) characteristic slopes and populate
        static_generation_data limits."""
        x = np.array([[96, 103], [120, 127]]) / 110
        y = np.array([[0.48, 0.41, 0.33], [-0.23, -0.33, -0.41]])
        m = (y[1] - y[0]) / (x[0, 1] - x[0, 0])
        b = np.array([y[0] - m * x[i, 0] for i in range(len(x))]).T

        m_qp_max = (0.1 - y[0]) / (0.1 - 0.2)
        m_qp_min = (-0.1 - y[1]) / (0.1 - 0.2)
        b_qp_max = 0.1 - m_qp_max * 0.1
        b_qp_min = -0.1 - m_qp_min * 0.1

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
                p_inst[sgens_var_q] * 0.1
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
        """Add a wind maximisation objective that subtracts line losses."""

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

        # --- QU Variante 1 ---
        @model.Constraint(model.WIND_HC)
        def QW_min(model, w):
            return model.qsG[w] >= -0.41 * model.psG[w]

        @model.Constraint(model.WIND_HC)
        def QW_max(model, w):
            return model.qsG[w] <= 0.48 * model.psG[w]

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
        """Unfix real and reactive power variables for all HC wind
        generators."""
        for w in model.WIND_HC:
            model.psG[w].unfix()
            model.qsG[w].unfix()

    def get_all_acopf(self, model):
        """No-op placeholder; required by the multi-period mix-in interface."""
        pass
