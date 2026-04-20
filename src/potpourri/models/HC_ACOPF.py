"""Hosting Capacity AC OPF model for sizing wind generation in distribution
grids."""

import copy

import numpy as np
import pandas as pd
import pyomo.environ as pyo
from src.potpourri.models.ACOPF_base import ACOPF
import pandapower as pp


class HC_ACOPF(ACOPF):
    """Hosting Capacity AC OPF for wind generation integration studies.

    Extends ACOPF with binary variables y[w] ∈ {0, 1} indicating whether each
    wind generator is active. The default objective maximises total wind
    generation minus network losses.

    If net.sgen has no 'wind_hc' column, candidate wind generators are
    automatically placed at every bus that is not an external grid bus.

    Args:
        net: A pandapower network. Wind candidates identified via sgen.wind_hc.
    """

    def __init__(self, net):
        net = copy.deepcopy(net)
        if "wind_hc" not in net.sgen:
            buses_excl_extGrids = net.bus.loc[
                ~net.bus.index.isin(net.ext_grid.bus)
            ].index
            pp.create_sgens(net, buses_excl_extGrids, p_mw=0, wind_hc=True)
        net.sgen.wind_hc.fillna(False, inplace=True)

        super().__init__(net)

    def _calc_opf_parameters(self, SWmax=10000, SWmin=0):
        """Compute HC-OPF limit data from the network.

        Extends ACOPF._calc_opf_parameters() with apparent power bounds
        SWmax and SWmin for WIND_HC generators, and Q-U characteristic slopes
        for grid code variant 1.

        Args:
            SWmax: Maximum apparent power per wind generator (MVA).
            SWmin: Minimum apparent power per active wind generator (MVA).
        """
        super()._calc_opf_parameters()

        if "windpot_p_mw" in self.net.bus:
            self.static_generation_data["windpot"] = self.net.bus.windpot_p_mw[
                self.net.sgen.bus.values
            ].values

        wind_hc_set = np.arange(len(self.net.sgen))[
            self.net.sgen.wind_hc & self.net.sgen.in_service
        ]
        self.SWmax_data = pd.Series(SWmax / self.baseMVA, wind_hc_set)
        self.SWmin_data = pd.Series(SWmin / self.baseMVA, wind_hc_set)

        self.m_qu_max = (0.48 + 0.23) / (96 - 103) * 110  # pyo.Variante 1
        self.qu_max = -self.m_qu_max * 120 / 110 + 0.48
        self.m_qu_min = (0.33 + 0.41) / (96 - 103) * 110  # pyo.Variante 3
        self.qu_min = -self.m_qu_min * 96 / 110 + 0.33

    def add_OPF(self, **kwargs):
        """Attach HC-OPF variables and constraints to self.model.

        Extends ACOPF.add_OPF() with binary variable y[w], apparent power
        envelope (SWmin·y ≤ S² ≤ SWmax·y), Q-P bounds (variant 1), Q-U bounds
        from the grid-code voltage characteristic, an optional real power cap
        from net.bus.windpot_p_mw, and a default objective maximising wind
        generation minus network losses.

        Args:
            **kwargs: Forwarded to ACOPF.add_OPF().
        """
        super().add_OPF(**kwargs)

        self.model.name = "HC_ACOPF"

        self.model.SWmax = pyo.Param(
            self.model.WIND_HC,
            initialize=self.SWmax_data[self.model.WIND_HC],
            mutable=True,
        )
        self.model.SWmin = pyo.Param(
            self.model.WIND_HC,
            initialize=self.SWmin_data[self.model.WIND_HC],
            mutable=True,
        )

        if "windpot_p_mw" in self.net.bus:
            self.model.pWmax = pyo.Param(
                self.model.WIND_HC,
                initialize=self.static_generation_data["windpot"][
                    self.model.WIND_HC
                ],
                mutable=True,
            )

        self.model.y = pyo.Var(
            self.model.WIND_HC, within=pyo.Binary, initialize=1.0
        )

        def obj_wind_loss_rule(model):
            return (
                sum(model.psG[w] for w in model.WIND_HC)
                - sum(model.pLfrom[l] + model.pLto[l] for l in model.L)
                - sum(model.pThv[t] + model.pTlv[t] for t in model.TRANSF)
            )

        self.model.obj = pyo.Objective(
            rule=obj_wind_loss_rule, sense=pyo.maximize
        )

        def SW_max(model, w):
            return (
                model.psG[w] ** 2 + model.qsG[w] ** 2
                <= model.SWmax[w] ** 2 * model.y[w]
            )

        def SW_min(model, w):
            return (
                model.psG[w] ** 2 + model.qsG[w] ** 2
                >= model.SWmin[w] ** 2 * model.y[w]
            )

        self.model.SW_max_constraint = pyo.Constraint(
            self.model.WIND_HC, rule=SW_max
        )
        self.model.SW_min_constraint = pyo.Constraint(
            self.model.WIND_HC, rule=SW_min
        )

        # --- generator power ---
        for w in self.model.WIND_HC:
            self.model.psG[w].unfix()
            self.model.qsG[w].unfix()

        # --- QU pyo.Variante 1 ---
        def QW_min(model, w):
            return model.qsG[w] >= -0.41 * model.psG[w]

        def QW_max(model, w):
            return model.qsG[w] <= 0.48 * model.psG[w]

        self.model.QW_min_constraint = pyo.Constraint(
            self.model.WIND_HC, rule=QW_min
        )
        self.model.QW_max_constraint = pyo.Constraint(
            self.model.WIND_HC, rule=QW_max
        )

        hc_sGbs_lookup = {g: b for (g, b) in self.model.sGbs}

        def QU_min_hc(model, w):
            if w not in hc_sGbs_lookup:
                return pyo.Constraint.Skip
            b = hc_sGbs_lookup[w]
            return (
                model.qsG[w]
                >= (self.m_qu_min * model.v[b] + self.qu_min) * model.psG[w]
            )

        self.model.QU_min_hc_constraint = pyo.Constraint(
            self.model.WIND_HC, rule=QU_min_hc
        )

        def QU_max_hc(model, w):
            if w not in hc_sGbs_lookup:
                return pyo.Constraint.Skip
            b = hc_sGbs_lookup[w]
            return (
                model.qsG[w]
                <= (self.m_qu_max * model.v[b] + self.qu_max) * model.psG[w]
            )

        self.model.QU_max_hc_constraint = pyo.Constraint(
            self.model.WIND_HC, rule=QU_max_hc
        )

        if "windpot_p_mw" in self.net.bus:

            def PW_max(model, w):
                return model.psG[w] <= model.pWmax[w]

            self.model.PW_max_constraint = pyo.Constraint(
                self.model.WIND_HC, rule=PW_max
            )

    def add_loss_obj(self):
        """Replace default objective with a weighted wind-vs-loss objective.

        Maximises: ε · Σ psG[w]  +  (1 - ε) · (- Σ losses)

        The mutable parameter eps (default 1.0) can be updated for sensitivity
        analysis without rebuilding the model: model.eps.set_value(0.5).
        """
        self.model.eps = pyo.Param(
            domain=pyo.Reals, initialize=1.0, mutable=True
        )

        def objective_pwind_loss(model):
            return model.eps * sum(model.psG[w] for w in model.WIND_HC) + (
                1 - model.eps
            ) * (-sum(model.pLfrom[l] + model.pLto[l] for l in model.L))

        self.model.obj.deactivate()
        self.model.OBJ_with_loss = pyo.Objective(
            rule=objective_pwind_loss, sense=pyo.maximize
        )
