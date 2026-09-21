# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Multi-period DC power flow mixin.

Adds linearised DC equations indexed over time steps.
"""

import numpy as np
import pandas as pd
from pyomo.environ import *
from potpourri.models_multi_period.basemodel_multi_period import (
    Basemodel_multi_period,
)

# TODO make multiperiod


class DC_multi_period(Basemodel_multi_period):
    """Multi-period linearised DC power flow model."""

    def __init__(self, net, toT, fromT=None, pf=1):
        super().__init__(net, toT, fromT, pf)

        x = self.net._ppc["branch"][:, 3].real
        BL = -1 / x
        trafo_start = len(self.net.line)
        trafo_end = trafo_start + len(self.net.trafo)
        imp_table = self.net.get("impedance")
        n_imp = (
            len(imp_table)
            if imp_table is not None and not imp_table.empty
            else 0
        )

        self.trafo_data = self.trafo_data.assign(
            **{"BLT_data": BL[trafo_start:trafo_end]}
        )

        # native lines + impedance rows (in that order) populate line_data
        if n_imp:
            line_idx_ppc = np.r_[
                np.arange(0, trafo_start),
                np.arange(trafo_end, trafo_end + n_imp),
            ]
        else:
            line_idx_ppc = np.arange(0, trafo_start)
        self.line_data["BL_data"] = BL[line_idx_ppc]

        ZN = self.net.bus.vn_kv**2 / self.baseMVA
        y_s_line = -1 / (
            self.net.line.x_ohm_per_km * self.net.line.length_km
        )  # according to matpower manual dc modeling
        bl_line = y_s_line * ZN[self.net.line.from_bus].values
        # Use a pd.Series so ``self.BL_data[self.model.L]`` indexing works
        # the same way it does for native lines + impedance synthetic IDs.
        line_index = list(self.net.line.index)
        if n_imp:
            bl_imp = BL[trafo_end : trafo_end + n_imp]
            self.BL_data = pd.Series(
                np.concatenate([bl_line, bl_imp]),
                index=line_index + [trafo_start + i for i in range(n_imp)],
            )
        else:
            self.BL_data = pd.Series(bl_line, index=line_index)

        self.create_model()

    def create_model(self):
        """Build the multi-period DC model, in place.

        Build the multi-period DC Pyomo model with susceptance parameters
        and time-variant KCL/KVL constraints.

        Conventions
        -----------
        * `BL`, `BLT`, `shift`, `GB` — single-period parameters indexed by
          branch / shunt (no time index).
        * `delta[bus, t]`, `pLfrom[l, t]`, `pLto[l, t]`, `pThv/pTlv[l, t]`,
          `psG/pG[g, t]`, `pD[d, t]` — time-indexed variables.
        * `deltaL[l, t]`, `deltaLT[l, t]` — time-indexed angle differences.

        The previous version had a dead ``self.T is None`` branch (``self.T``
        is always an int set by ``Basemodel_multi_period.__init__``), used the
        bare int ``self.T`` rather than the Pyomo Set ``self.model.T`` for
        constraint indexing, and accessed several single-period parameters
        with a spurious time index.
        """
        super().create_model()

        self.model.name = "DC"

        # --- single-period line / transformer parameters ---
        self.model.BL = Param(
            self.model.L, within=Reals, initialize=self.BL_data[self.model.L]
        )  # line + impedance series susceptance
        self.model.BLT = Param(
            self.model.TRANSF,
            within=Reals,
            initialize=self.trafo_data.BLT_data[self.model.TRANSF],
        )  # transformer series susceptance

        # --- time-indexed angle differences ---
        self.model.deltaL = Var(
            self.model.L, self.model.T, domain=Reals
        )  # angle difference across lines + impedance branches
        self.model.deltaLT = Var(
            self.model.TRANSF, self.model.T, domain=Reals
        )  # angle difference across transformers

        # --- KCL at each bus, per time step ---
        self.build_kcl()

        # --- KVL on lines + impedance branches ---
        @self.model.Constraint(self.model.L, self.model.T)
        def KVL_real_fromend(model, l, t):
            """Active power entering line `l` at its from bus.

            Args:
                model: The Pyomo model being built.
                l: Branch index.
                t: Time index.

            Returns:
                A Pyomo expression.
            """
            return model.pLfrom[l, t] == (-model.BL[l]) * model.deltaL[l, t]

        @self.model.Constraint(self.model.L, self.model.T)
        def KVL_real_toend(model, l, t):
            """Active power entering line `l` at its to bus.

            Args:
                model: The Pyomo model being built.
                l: Branch index.
                t: Time index.

            Returns:
                A Pyomo expression.
            """
            return model.pLto[l, t] == (model.BL[l]) * model.deltaL[l, t]

        # --- KVL on transformers ---
        @self.model.Constraint(self.model.TRANSF, self.model.T)
        def KVL_trans_fromend(model, l, t):
            """Active power entering transformer `l` at its HV bus.

            Args:
                model: The Pyomo model being built.
                l: Branch index.
                t: Time index.

            Returns:
                A Pyomo expression.
            """
            return model.pThv[l, t] == (-model.BLT[l]) * model.deltaLT[l, t]

        @self.model.Constraint(self.model.TRANSF, self.model.T)
        def KVL_trans_toend(model, l, t):
            """Active power entering transformer `l` at its LV bus.

            Args:
                model: The Pyomo model being built.
                l: Branch index.
                t: Time index.

            Returns:
                A Pyomo expression.
            """
            return model.pTlv[l, t] == (model.BLT[l]) * model.deltaLT[l, t]

        # --- angle-difference identities ---
        @self.model.Constraint(self.model.L, self.model.T)
        def phase_angle_diff1(model, l, t):
            """Upper bound on a branch's angle difference.

            Args:
                model: The Pyomo model being built.
                l: Branch index.
                t: Time index.

            Returns:
                A Pyomo expression.
            """
            return (
                model.deltaL[l, t]
                == model.delta[model.A[l, 1], t]
                - model.delta[model.A[l, 2], t]
            )

        @self.model.Constraint(self.model.TRANSF, self.model.T)
        def phase_angle_diff2(model, l, t):
            """Lower bound on a branch's angle difference.

            Args:
                model: The Pyomo model being built.
                l: Branch index.
                t: Time index.

            Returns:
                A Pyomo expression.
            """
            return (
                model.deltaLT[l, t]
                == model.delta[model.AT[l, 1], t]
                - model.delta[model.AT[l, 2], t]
                - model.shift[l]
            )

    #: The DC formulation has a single real-power balance and no reactive one.
    KCL_CONSTRAINTS = ("KCL_def",)

    def _kcl_def_rule(self, model, b, t):
        """Active-power balance at bus `b`.

        Args:
            model: The Pyomo model being built.
            b: Bus index (a ppc bus number).
            t: Time index.

        Returns:
            A Pyomo expression.
        """
        kcl = sum(
            model.psG[g, t] for g in model.sG if (g, b) in model.sGbs
        ) + sum(model.pG[g, t] for g in model.G if (g, b) in model.Gbs) == sum(
            model.pD[d, t] for d in model.D if (b, d) in model.Dbs
        ) + sum(
            model.pLfrom[l, t] for l in model.L if model.A[l, 1] == b
        ) + sum(model.pLto[l, t] for l in model.L if model.A[l, 2] == b) + sum(
            model.pThv[l, t] for l in model.TRANSF if model.AT[l, 1] == b
        ) + sum(
            model.pTlv[l, t] for l in model.TRANSF if model.AT[l, 2] == b
        ) + sum(
            model.GB[s, t] for s in model.SHUNT if (b, s) in model.SHUNTbs
        ) + self.KCL_flexibility(model, b, t)
        if isinstance(kcl, bool):
            return Constraint.Skip
        return kcl
