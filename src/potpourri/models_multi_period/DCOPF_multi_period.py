# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Multi-period DC OPF: linearised power flow plus limits."""

from pyomo.environ import *
from potpourri.models_multi_period.DC_multi_period import DC_multi_period
from potpourri.models_multi_period.OPF_multi_period import OPF_multi_period


class DCOPF_multi_period(DC_multi_period, OPF_multi_period):
    """Multi-period DC OPF: linearised flow plus limits.

    Multi-period DC OPF: linearised flow plus generator and thermal limit
    constraints.

    Construction follows the single-period DCOPF pattern: ``__init__`` only
    builds the underlying DC power-flow model (KCL + KVL + angle relations over
    the time horizon). Thermal/operational constraints are attached explicitly
    with :meth:`add_OPF` (which delegates to :meth:`OPF_multi_period.add_OPF`
    for ``SLmax`` / ``SLmaxT`` / generator limits and then adds DC line /
    transformer flow bounds on top).
    """

    def __init__(self, net, toT, fromT=None, pf=1):
        super().__init__(net, toT, fromT, pf)
        self.model.name = "DCOPF"

    def add_OPF(self, angle_limits: bool = False, **kwargs):
        """Add branch ratings, generator and demand limits.

        Attach OPF constraints — line / transformer ratings, generator
        and demand limits, plus DC apparent-power thermal limits.

        Args:
            angle_limits: When ``True``, enforce branch
                phase-angle-difference constraints
                ``angmin ≤ δ_from − δ_to ≤ angmax`` at every time step, read
                from ``net.line.angmin_degree`` / ``net.line.angmax_degree``
                and the transformer equivalent. Defaults to ``False``, as on
                the single-period :class:`~potpourri.models.DCOPF.DCOPF`.
            **kwargs: Forwarded to :meth:`OPF_multi_period.add_OPF`, which
                rejects unsupported names. In particular ``thermal_limit``
                is AC-only — the DC model is lossless and carries no
                reactive power, so its limit is a real-power bound with no
                current-versus-MVA distinction to make.
        """
        OPF_multi_period.add_OPF(self, **kwargs)

        # --- line power limits (DC: lossless, check sending end only) ---
        @self.model.Constraint(self.model.L, self.model.T)
        def line_lim_upper(model, l, t):
            """Upper branch-flow limit on line `l`.

            Args:
                model: The Pyomo model being built.
                l: Branch index.
                t: Time index.

            Returns:
                A Pyomo expression.
            """
            return model.pLfrom[l, t] <= model.SLmax[l]

        @self.model.Constraint(self.model.L, self.model.T)
        def line_lim_lower(model, l, t):
            """Lower branch-flow limit on line `l`.

            Args:
                model: The Pyomo model being built.
                l: Branch index.
                t: Time index.

            Returns:
                A Pyomo expression.
            """
            return model.pLfrom[l, t] >= -model.SLmax[l]

        # --- transformer power limits ---
        @self.model.Constraint(self.model.TRANSF, self.model.T)
        def transf_lim_upper(model, l, t):
            """Upper branch-flow limit on transformer `l`.

            Args:
                model: The Pyomo model being built.
                l: Branch index.
                t: Time index.

            Returns:
                A Pyomo expression.
            """
            return model.pThv[l, t] <= model.SLmaxT[l]

        @self.model.Constraint(self.model.TRANSF, self.model.T)
        def transf_lim_lower(model, l, t):
            """Lower branch-flow limit on transformer `l`.

            Args:
                model: The Pyomo model being built.
                l: Branch index.
                t: Time index.

            Returns:
                A Pyomo expression.
            """
            return model.pThv[l, t] >= -model.SLmaxT[l]

        # --- optional branch angle-difference limits ---
        if angle_limits:
            self._add_branch_angle_limits()
