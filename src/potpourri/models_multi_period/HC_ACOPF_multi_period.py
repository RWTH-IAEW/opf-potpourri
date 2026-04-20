"""Multi-period Hosting Capacity AC OPF for wind generation integration studies."""

import copy

import numpy as np
import pandas as pd
from pyomo.environ import *
from src.potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period
from src.potpourri.technologies.windpower import Windpower_multi_period
import pandapower as pp

# TODO make multiperiod


class HC_ACOPF_multi_period(ACOPF_multi_period):
    """Multi-period hosting capacity AC OPF; delegates HC constraints to Windpower_multi_period."""

    def __init__(self, net, toT, fromT=None, pf=1):
        if 'wind_hc' not in net.sgen:
            net = copy.deepcopy(net)
            buses_excl_extGrids = net.bus.loc[~net.bus.index.isin(net.ext_grid.bus)].index

            pp.create_sgens(net, buses_excl_extGrids, p_mw=0, wind_hc=True)
        net.sgen.wind_hc.fillna(False, inplace=True)

        super().__init__(net, toT, fromT, pf)

        # noinspection PyProtectedMember

    def _calc_opf_parameters(self, SWmax=10000, SWmin=0):
        """Extend AC-OPF parameters with wind HC apparent power bounds via Windpower_multi_period."""
        super()._calc_opf_parameters()

        if 'windpot_p_mw' in self.net.bus:
            windpower_object = next(
                (obj for obj in self.flexibilities if isinstance(obj, Windpower_multi_period)), None
            )
            windpower_object._calc_wind_opf_parameters(self.net, SWmax=SWmax, SWmin=SWmin)

    def add_OPF(self, **kwargs):
        """Extend ACOPF.add_OPF() with HC wind constraints and objective via Windpower_multi_period."""
        super().add_OPF(**kwargs)

        self.model.name = "HC_ACOPF"

        if 'windpot_p_mw' in self.net.bus:
            windpower_object = next(
                (obj for obj in self.flexibilities if isinstance(obj, Windpower_multi_period)), None
            )
            windpower_object.get_hc_acopf_parameters(self.model, self.net)
            windpower_object.get_hc_acopf_variables(self.model)
            windpower_object.get_objective(self.model)
            windpower_object.get_constraints(self.model, self.net)
            windpower_object.unfix_variables(self.model)

    def add_loss_obj(self):
        """Replace default objective with weighted wind-vs-loss objective using mutable eps parameter."""
        self.model.eps = Param(domain=Reals, initialize=1., mutable=True)

        def objective_pwind_loss(model):
            return model.eps * sum(model.psG[w] for w in model.WIND_HC) + (1 - model.eps) * (
                - sum(model.pLfrom[l] + model.pLto[l] for l in model.L))

        self.model.obj_hc.deactivate()
        self.model.OBJ_with_loss = Objective(rule=objective_pwind_loss, sense=maximize)