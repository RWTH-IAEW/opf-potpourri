"""Multi-period DC OPF combining linearised DC power flow with operational limit constraints."""

from pyomo.environ import *
from src.potpourri.models_multi_period.DC_multi_period import DC_multi_period
from src.potpourri.models_multi_period.OPF_multi_period import OPF_multi_period

# TODO make multiperiod


class DCOPF_multi_period(DC_multi_period, OPF_multi_period):
    """Multi-period DC OPF model: linearised power flow with generator and thermal limit constraints."""

    def __init__(self, net, toT, fromT=None, pf=1):
        super().__init__(net, toT, fromT, pf)

        self.create_model()


    def create_model(self):
        """Build the multi-period DCOPF model and attach line/transformer thermal limit constraints."""
        super().create_model()
        self.model.name = "DCOPF"


        # --- line power limits ---
        @self.model.Constraint(self.model.L, self.model.T)
        def line_lim1_def(model, l, t):
            return model.pL[l, t] <= model.SLmax[l]

        @self.model.Constraint(self.model.L, self.model.T)
        def line_lim2_def(model, l, t):
            return model.pL[l,t] >= -model.SLmax[l]



        # --- power flow limits on transformer lines---
        @self.model.Constraint(self.model.TRANSF, self.model.T)
        def transf_lim1_def(model, l, t):
            return model.pLT[l,t] <= model.SLmaxT[l]
        @self.model.Constraint(self.model.TRANSF, self.model.T)
        def transf_lim2_def(model, l, t):
            return model.pLT[l, t] >= -(model.SLmaxT[l])



