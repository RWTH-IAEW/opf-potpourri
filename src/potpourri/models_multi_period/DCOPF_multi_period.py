"""Multi-period DC OPF combining linearised DC power flow with operational
limit constraints."""

from pyomo.environ import *
from potpourri.models_multi_period.DC_multi_period import DC_multi_period
from potpourri.models_multi_period.OPF_multi_period import OPF_multi_period


class DCOPF_multi_period(DC_multi_period, OPF_multi_period):
    """Multi-period DC OPF: linearised power flow with generator and thermal
    limit constraints.

    Construction follows the single-period DCOPF pattern: ``__init__`` only
    builds the underlying DC power-flow model (KCL + KVL + angle relations
    over the time horizon). Thermal/operational constraints are attached
    explicitly with :meth:`add_OPF` (which delegates to
    :meth:`OPF_multi_period.add_OPF` for ``SLmax`` / ``SLmaxT`` / generator
    limits and then adds DC line / transformer flow bounds on top).
    """

    def __init__(self, net, toT, fromT=None, pf=1):
        super().__init__(net, toT, fromT, pf)
        self.model.name = "DCOPF"

    def add_OPF(self, **kwargs):
        """Attach OPF constraints — line / transformer ratings, generator
        and demand limits, plus DC apparent-power thermal limits."""
        OPF_multi_period.add_OPF(self, **kwargs)

        # --- line power limits (DC: lossless, check sending end only) ---
        @self.model.Constraint(self.model.L, self.model.T)
        def line_lim_upper(model, l, t):
            return model.pLfrom[l, t] <= model.SLmax[l]

        @self.model.Constraint(self.model.L, self.model.T)
        def line_lim_lower(model, l, t):
            return model.pLfrom[l, t] >= -model.SLmax[l]

        # --- transformer power limits ---
        @self.model.Constraint(self.model.TRANSF, self.model.T)
        def transf_lim_upper(model, l, t):
            return model.pThv[l, t] <= model.SLmaxT[l]

        @self.model.Constraint(self.model.TRANSF, self.model.T)
        def transf_lim_lower(model, l, t):
            return model.pThv[l, t] >= -model.SLmaxT[l]
