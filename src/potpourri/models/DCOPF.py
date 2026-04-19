import pyomo.environ as pyo
from src.potpourri.models.DC import DC
from src.potpourri.models.OPF import OPF


class DCOPF(DC, OPF):
    def __init__(self, net):
        super().__init__(net)
        self.create_model()

    def create_model(self):
        super().create_model()
        self.model.name = "DCOPF"

    def add_OPF(self, **kwargs):
        """Attach DC-OPF sets, parameters, and constraints to self.model.

        Calls OPF.add_OPF() for generator/demand limits and line ratings, then
        adds DC thermal limit constraints using the correct pLfrom / pThv variables.
        """
        super().add_OPF(**kwargs)

        # --- line power limits (check sending end; DC is approximately lossless) ---
        def line_lim_upper(model, l):
            return model.pLfrom[l] <= model.SLmax[l]

        def line_lim_lower(model, l):
            return model.pLfrom[l] >= -model.SLmax[l]

        self.model.line_lim_from = pyo.Constraint(self.model.L, rule=line_lim_upper)
        self.model.line_lim_to = pyo.Constraint(self.model.L, rule=line_lim_lower)

        # --- transformer power limits ---
        def transf_lim_upper(model, l):
            return model.pThv[l] <= model.SLmaxT[l]

        def transf_lim_lower(model, l):
            return model.pThv[l] >= -model.SLmaxT[l]

        self.model.transf_lim1 = pyo.Constraint(self.model.TRANSF, rule=transf_lim_upper)
        self.model.transf_lim2 = pyo.Constraint(self.model.TRANSF, rule=transf_lim_lower)
