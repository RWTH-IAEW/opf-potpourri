# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Multi-period AC power flow mixin.

Adds full AC equations with voltage magnitudes over time.
"""

import numpy as np

from potpourri.models.basemodel import branch_charging_admittance
from pyomo.environ import *
from potpourri.models_multi_period.basemodel_multi_period import (
    Basemodel_multi_period,
)
from potpourri.technologies.demand import Demand_multi_period


class AC_multi_period(Basemodel_multi_period):
    """Multi-period AC power flow, indexed over time steps."""

    def __init__(self, net, toT, fromT=None, pf=1):
        super().__init__(net, toT, fromT, pf)

        self.BB_data = (
            -self.net.shunt.q_mvar * self.net.shunt.step / self.baseMVA
        )

        # line and transformer admittances
        r = self.net._ppc["branch"][:, 2].real
        x = self.net._ppc["branch"][:, 3].real
        y = branch_charging_admittance(self.net._ppc["branch"])
        gt_ik = r / (r**2 + x**2)
        bt_ik = -x / (r**2 + x**2)
        BiiT = bt_ik + y.imag / 2
        BikT = -bt_ik
        GiiT = gt_ik + y.real / 2
        GikT = -gt_ik
        trafo_start = len(self.net.line)
        trafo_end = trafo_start + len(self.net.trafo)
        imp_table = self.net.get("impedance")
        n_imp = (
            len(imp_table)
            if imp_table is not None and not imp_table.empty
            else 0
        )
        # Mirror of single-period AC: include impedance rows (slice
        # [trafo_end : trafo_end + n_imp]) alongside native lines in
        # line_data so model.L sees them.
        if n_imp:
            line_idx_ppc = np.r_[
                np.arange(0, trafo_start),
                np.arange(trafo_end, trafo_end + n_imp),
            ]
        else:
            line_idx_ppc = np.arange(0, trafo_start)

        self.line_data = self.line_data.assign(
            **{
                "Bii_data": BiiT[line_idx_ppc],
                "Bik_data": BikT[line_idx_ppc],
                "Gii_data": GiiT[line_idx_ppc],
                "Gik_data": GikT[line_idx_ppc],
            }
        )
        self.trafo_data = self.trafo_data.assign(
            **{
                "BiiT_data": BiiT[trafo_start:trafo_end],
                "BikT_data": BikT[trafo_start:trafo_end],
                "GiiT_data": GiiT[trafo_start:trafo_end],
                "GikT_data": GikT[trafo_start:trafo_end],
            }
        )

        self.create_model()

    def create_model(self):
        """Build the multi-period AC model, in place.

        Build the multi-period AC Pyomo model with admittance parameters
        and AC KCL/KVL constraints.
        """
        super().create_model()
        self.model.name = "AC"

        # shunt
        self.BB_data_dict, self.BB_tuple = self.make_to_dict(
            self.model.SHUNT,
            self.model.T,
            self.BB_data[self.model.SHUNT],
            False,
        )
        self.model.BB = Param(
            self.BB_tuple, within=Reals, initialize=self.BB_data_dict
        )  # shunt susceptance

        # derived line parameters
        self.Bii_data_dict, self.Bii_tuple = self.make_to_dict(
            self.model.L,
            self.model.T,
            self.line_data.Bii_data[self.model.L],
            False,
        )
        self.Bik_data_dict, self.Bik_tuple = self.make_to_dict(
            self.model.L,
            self.model.T,
            self.line_data.Bik_data[self.model.L],
            False,
        )
        self.Gii_data_dict, self.Gii_tuple = self.make_to_dict(
            self.model.L,
            self.model.T,
            self.line_data.Gii_data[self.model.L],
            False,
        )
        self.Gik_data_dict, self.Gik_tuple = self.make_to_dict(
            self.model.L,
            self.model.T,
            self.line_data.Gik_data[self.model.L],
            False,
        )

        self.model.Bii = Param(
            self.Bii_tuple, within=Reals, initialize=self.Bii_data_dict
        )
        self.model.Bik = Param(
            self.Bik_tuple, within=Reals, initialize=self.Bik_data_dict
        )
        self.model.Gii = Param(
            self.Gii_tuple, within=Reals, initialize=self.Gii_data_dict
        )
        self.model.Gik = Param(
            self.Gik_tuple, within=Reals, initialize=self.Gik_data_dict
        )

        self.model.BiiT = Param(
            self.model.TRANSF,
            within=Reals,
            initialize=self.trafo_data.BiiT_data[self.model.TRANSF],
        )
        self.model.BikT = Param(
            self.model.TRANSF,
            within=Reals,
            initialize=self.trafo_data.BikT_data[self.model.TRANSF],
        )
        self.model.GiiT = Param(
            self.model.TRANSF,
            within=Reals,
            initialize=self.trafo_data.GiiT_data[self.model.TRANSF],
        )
        self.model.GikT = Param(
            self.model.TRANSF,
            within=Reals,
            initialize=self.trafo_data.GikT_data[self.model.TRANSF],
        )
        # create instance of demand
        demand_object = next(
            (
                obj
                for obj in self.flexibilities
                if isinstance(obj, Demand_multi_period)
            ),
            None,
        )
        demand_object.get_all_ac(self.model)

        # external grid voltage
        self.model.v_b0 = Param(
            self.model.b0,
            within=NonNegativeReals,
            initialize=self.bus_data.v_m[self.model.b0],
        )

        # time dependent control variables

        # --- control variables --- stay multiperiod
        self.model.qLfrom = Var(
            self.model.L, self.model.T, domain=Reals
        )  # reactive power injected at b onto line
        self.model.qLto = Var(
            self.model.L, self.model.T, domain=Reals
        )  # reactive power injected at b' onto line
        self.model.qThv = Var(
            self.model.TRANSF, self.model.T, domain=Reals
        )  # reactive power injected at b onto transformer
        self.model.qTlv = Var(
            self.model.TRANSF, self.model.T, domain=Reals
        )  # reactive power injected at b' onto transformer
        self.model.v = Var(
            self.model.B, self.model.T, domain=NonNegativeReals, initialize=1.0
        )  # voltage magnitude at bus b, rad

        # correct?
        self.model.qG = Var(self.model.G, self.model.T, domain=Reals)

        # --- Kirchoff's current law at each bus b ---
        self.build_kcl()

        # --- Kirchoff's voltage law on each line ---
        def KVL_real_fromend(model, l, t):
            r"""Active power entering line `l` at its from bus, at time `t`.

            The time-indexed twin of `AC.create_model`'s `KVL_real_fromend`;
            see [`potpourri.models.AC`][potpourri.models.AC] for the
            $\pi$-model relation and the sign conventions, which are identical
            here.

            Args:
                model: The Pyomo model being built.
                l: Line index from `model.L`.
                t: Time index from `model.T`.

            Returns:
                A Pyomo equality expression defining `pLfrom[l, t]`.
            """
            return model.pLfrom[l, t] == model.Gii[(l, t)] * (
                model.v[model.A[l, 1], t] ** 2
            ) + model.v[model.A[l, 1], t] * model.v[model.A[l, 2], t] * (
                model.Bik[l, t]
                * sin(
                    model.delta[model.A[l, 1], t]
                    - model.delta[model.A[l, 2], t]
                )
                + model.Gik[l, t]
                * cos(
                    model.delta[model.A[l, 1], t]
                    - model.delta[model.A[l, 2], t]
                )
            )

        def KVL_real_toend(model, l, t):
            r"""Active power entering line `l` at its to bus, at time `t`.

            The time-indexed twin of `AC.create_model`'s `KVL_real_toend`; see
            [`potpourri.models.AC`][potpourri.models.AC] for the $\pi$-model
            relation and the sign conventions, which are identical here.

            Args:
                model: The Pyomo model being built.
                l: Line index from `model.L`.
                t: Time index from `model.T`.

            Returns:
                A Pyomo equality expression defining `pLto[l, t]`.
            """
            return model.pLto[l, t] == model.Gii[l, t] * (
                model.v[model.A[l, 2], t] ** 2
            ) + model.v[model.A[l, 1], t] * model.v[model.A[l, 2], t] * (
                model.Bik[l, t]
                * sin(
                    model.delta[model.A[l, 2], t]
                    - model.delta[model.A[l, 1], t]
                )
                + model.Gik[l, t]
                * cos(
                    model.delta[model.A[l, 2], t]
                    - model.delta[model.A[l, 1], t]
                )
            )

        def KVL_reactive_fromend(model, l, t):
            r"""Reactive power entering line `l` at its from bus, at time `t`.

            The time-indexed twin of `AC.create_model`'s
            `KVL_reactive_fromend`; see
            [`potpourri.models.AC`][potpourri.models.AC] for the $\pi$-model
            relation and the sign conventions, which are identical here.

            Args:
                model: The Pyomo model being built.
                l: Line index from `model.L`.
                t: Time index from `model.T`.

            Returns:
                A Pyomo equality expression defining `qLfrom[l, t]`.
            """
            return model.qLfrom[l, t] == -model.Bii[l, t] * (
                model.v[model.A[l, 1], t] ** 2
            ) + model.v[model.A[l, 1], t] * model.v[model.A[l, 2], t] * (
                model.Gik[l, t]
                * sin(
                    model.delta[model.A[l, 1], t]
                    - model.delta[model.A[l, 2], t]
                )
                - model.Bik[l, t]
                * cos(
                    model.delta[model.A[l, 1], t]
                    - model.delta[model.A[l, 2], t]
                )
            )

        def KVL_reactive_toend(model, l, t):
            r"""Reactive power entering line `l` at its to bus, at time `t`.

            The time-indexed twin of `AC.create_model`'s `KVL_reactive_toend`;
            see [`potpourri.models.AC`][potpourri.models.AC] for the
            $\pi$-model relation and the sign conventions, which are identical
            here.

            Args:
                model: The Pyomo model being built.
                l: Line index from `model.L`.
                t: Time index from `model.T`.

            Returns:
                A Pyomo equality expression defining `qLto[l, t]`.
            """
            return model.qLto[l, t] == -model.Bii[l, t] * (
                model.v[model.A[l, 2], t] ** 2
            ) + model.v[model.A[l, 1], t] * model.v[model.A[l, 2], t] * (
                model.Gik[l, t]
                * sin(
                    model.delta[model.A[l, 2], t]
                    - model.delta[model.A[l, 1], t]
                )
                - model.Bik[l, t]
                * cos(
                    model.delta[model.A[l, 2], t]
                    - model.delta[model.A[l, 1], t]
                )
            )

        self.model.KVL_real_from = Constraint(
            self.model.L, self.model.T, rule=KVL_real_fromend
        )
        self.model.KVL_real_to = Constraint(
            self.model.L, self.model.T, rule=KVL_real_toend
        )
        self.model.KVL_reactive_from = Constraint(
            self.model.L, self.model.T, rule=KVL_reactive_fromend
        )
        self.model.KVL_reactive_to = Constraint(
            self.model.L, self.model.T, rule=KVL_reactive_toend
        )

        # --- Kirchoff's voltage law on each transformer line ---
        def KVL_real_fromendTransf(model, l, t):
            """Active power entering transformer `l` at its HV bus.

            The time-indexed twin of `AC.create_model`'s
            `KVL_real_fromendTransf`. The tap ratio and phase shift are
            time-independent, so only the voltages and angles carry the extra
            index.

            Args:
                model: The Pyomo model being built.
                l: Transformer index from `model.TRANSF`.
                t: Time index from `model.T`.

            Returns:
                A Pyomo equality expression defining `pThv[l, t]`.
            """
            if model.shift[l]:
                return model.pThv[l, t] == model.GiiT[l] / model.Tap[
                    l, t
                ] ** 2 * (model.v[model.AT[l, 1], t] ** 2) + model.v[
                    model.AT[l, 1], t
                ] * model.v[model.AT[l, 2], t] / model.Tap[l, t] * (
                    model.GikT[l]
                    * cos(
                        model.delta[model.AT[l, 1], t]
                        - model.delta[model.AT[l, 2], t]
                        - model.shift[l]
                    )
                    + model.BikT[l]
                    * sin(
                        model.delta[model.AT[l, 1], t]
                        - model.delta[model.AT[l, 2], t]
                        - model.shift[l]
                    )
                )

            return model.pThv[l, t] == model.GiiT[l] / model.Tap[l, t] ** 2 * (
                model.v[model.AT[l, 1], t] ** 2
            ) + model.v[model.AT[l, 1], t] * model.v[
                model.AT[l, 2], t
            ] / model.Tap[l, t] * (
                model.GikT[l]
                * cos(
                    model.delta[model.AT[l, 1], t]
                    - model.delta[model.AT[l, 2], t]
                )
                + model.BikT[l]
                * sin(
                    model.delta[model.AT[l, 1], t]
                    - model.delta[model.AT[l, 2], t]
                )
            )

        def KVL_real_toendTransf(model, l, t):
            """Active power entering transformer `l` at its LV bus.

            The time-indexed twin of `AC.create_model`'s
            `KVL_real_toendTransf`. The tap ratio and phase shift are
            time-independent, so only the voltages and angles carry the extra
            index.

            Args:
                model: The Pyomo model being built.
                l: Transformer index from `model.TRANSF`.
                t: Time index from `model.T`.

            Returns:
                A Pyomo equality expression defining `pTlv[l, t]`.
            """
            if model.shift[l]:
                return model.pTlv[l, t] == model.GiiT[l] * (
                    model.v[model.AT[l, 2], t] ** 2
                ) + model.v[model.AT[l, 1], t] * model.v[
                    model.AT[l, 2], t
                ] / model.Tap[l, t] * (
                    model.BikT[l]
                    * sin(
                        model.delta[model.AT[l, 2], t]
                        - model.delta[model.AT[l, 1], t]
                        + model.shift[l]
                    )
                    + model.GikT[l]
                    * cos(
                        model.delta[model.AT[l, 2], t]
                        - model.delta[model.AT[l, 1], t]
                        + model.shift[l]
                    )
                )

            return model.pTlv[l, t] == model.GiiT[l] * (
                model.v[model.AT[l, 2], t] ** 2
            ) + model.v[model.AT[l, 1], t] * model.v[
                model.AT[l, 2], t
            ] / model.Tap[l, t] * (
                model.BikT[l]
                * sin(
                    model.delta[model.AT[l, 2], t]
                    - model.delta[model.AT[l, 1], t]
                )
                + model.GikT[l]
                * cos(
                    model.delta[model.AT[l, 2], t]
                    - model.delta[model.AT[l, 1], t]
                )
            )

        def KVL_reactive_fromendTransf(model, l, t):
            """Reactive power entering transformer `l` at its HV bus.

            The time-indexed twin of `AC.create_model`'s
            `KVL_reactive_fromendTransf`. The tap ratio and phase shift are
            time-independent, so only the voltages and angles carry the extra
            index.

            Args:
                model: The Pyomo model being built.
                l: Transformer index from `model.TRANSF`.
                t: Time index from `model.T`.

            Returns:
                A Pyomo equality expression defining `qThv[l, t]`.
            """
            if model.shift[l]:
                return model.qThv[l, t] == -model.BiiT[l] / model.Tap[
                    l, t
                ] ** 2 * (model.v[model.AT[l, 1], t] ** 2) + model.v[
                    model.AT[l, 1], t
                ] * model.v[model.AT[l, 2], t] / model.Tap[l, t] * (
                    -model.BikT[l]
                    * cos(
                        model.delta[model.AT[l, 1], t]
                        - model.delta[model.AT[l, 2], t]
                        - model.shift[l]
                    )
                    + model.GikT[l]
                    * sin(
                        model.delta[model.AT[l, 1], t]
                        - model.delta[model.AT[l, 2], t]
                        - model.shift[l]
                    )
                )

            return model.qThv[l, t] == -model.BiiT[l] / model.Tap[
                l, t
            ] ** 2 * (model.v[model.AT[l, 1], t] ** 2) + model.v[
                model.AT[l, 1], t
            ] * model.v[model.AT[l, 2], t] / model.Tap[l, t] * (
                -model.BikT[l]
                * cos(
                    model.delta[model.AT[l, 1], t]
                    - model.delta[model.AT[l, 2], t]
                )
                + model.GikT[l]
                * sin(
                    model.delta[model.AT[l, 1], t]
                    - model.delta[model.AT[l, 2], t]
                )
            )

        def KVL_reactive_toendTransf(model, l, t):
            """Reactive power entering transformer `l` at its LV bus.

            The time-indexed twin of `AC.create_model`'s
            `KVL_reactive_toendTransf`. The tap ratio and phase shift are
            time-independent, so only the voltages and angles carry the extra
            index.

            Args:
                model: The Pyomo model being built.
                l: Transformer index from `model.TRANSF`.
                t: Time index from `model.T`.

            Returns:
                A Pyomo equality expression defining `qTlv[l, t]`.
            """
            if model.shift[l]:
                return model.qTlv[l, t] == -model.BiiT[l] * (
                    model.v[model.AT[l, 2], t] ** 2
                ) + model.v[model.AT[l, 1], t] * model.v[
                    model.AT[l, 2], t
                ] / model.Tap[l, t] * (
                    -model.BikT[l]
                    * cos(
                        model.delta[model.AT[l, 2], t]
                        - model.delta[model.AT[l, 1], t]
                        + model.shift[l]
                    )
                    + model.GikT[l]
                    * sin(
                        model.delta[model.AT[l, 2], t]
                        - model.delta[model.AT[l, 1], t]
                        + model.shift[l]
                    )
                )

            return model.qTlv[l, t] == -model.BiiT[l] * (
                model.v[model.AT[l, 2], t] ** 2
            ) + model.v[model.AT[l, 1], t] * model.v[
                model.AT[l, 2], t
            ] / model.Tap[l, t] * (
                -model.BikT[l]
                * cos(
                    model.delta[model.AT[l, 2], t]
                    - model.delta[model.AT[l, 1], t]
                )
                + model.GikT[l]
                * sin(
                    model.delta[model.AT[l, 2], t]
                    - model.delta[model.AT[l, 1], t]
                )
            )

        # Constraint Definitions for Pyomo, needed for the model
        self.model.KVL_real_fromTransf = Constraint(
            self.model.TRANSF, self.model.T, rule=KVL_real_fromendTransf
        )
        self.model.KVL_real_toTransf = Constraint(
            self.model.TRANSF, self.model.T, rule=KVL_real_toendTransf
        )
        self.model.KVL_reactive_fromTransf = Constraint(
            self.model.TRANSF, self.model.T, rule=KVL_reactive_fromendTransf
        )
        self.model.KVL_reactive_toTransf = Constraint(
            self.model.TRANSF, self.model.T, rule=KVL_reactive_toendTransf
        )

        # --- reactive generator power limits ---
        for g in self.model.sG:
            for t in self.model.T:
                self.model.qsG[(g, t)].fix(
                    self.model.QsG[(g, t)]
                )  # reactive power of static generators fixed

        # --- reactive demand limits ---
        for d in self.model.D:
            for t in self.model.T:
                self.model.qD[(d, t)].fix(self.model.QD[(d, t)])

        # --- reference bus voltage constraint ---
        for b0 in self.model.b0:
            for t in self.model.T:
                self.model.v[b0, t].fix(self.model.v_b0[b0])

    def _kcl_real_rule(self, model, b, t):
        """Active-power balance at bus `b` and time `t`.

        Generation minus storage equals demand, the branch injections at that
        bus and the shunt term, per time step. Device contributions registered
        through `register_kcl_real` are added here, which is what lets a
        flexible device attach after the power-flow equations were built.

        Args:
            model: The Pyomo model being built.
            b: Bus index from `model.B` (a ppc bus number).
            t: Time index from `model.T`.

        Returns:
            A Pyomo equality expression, or `Constraint.Skip` where every term
            at that bus is constant.
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
            model.GB[s, t] * model.v[b, t] ** 2
            for s in model.SHUNT
            if (b, s) in model.SHUNTbs and model.GB[s, t] != 0
        ) + self.KCL_flexibility(model, b, t)
        if isinstance(kcl, bool):
            return Constraint.Skip
        return kcl

    def _kcl_reactive_rule(self, model, b, t):
        """Reactive-power balance at bus `b` and time `t`.

        The reactive counterpart of `_kcl_real_rule`, including any terms
        registered via `register_kcl_reactive`.

        Args:
            model: The Pyomo model being built.
            b: Bus index from `model.B` (a ppc bus number).
            t: Time index from `model.T`.

        Returns:
            A Pyomo equality expression, or `Constraint.Skip` where every term
            at that bus is constant.
        """
        kcl = sum(
            model.qsG[g, t] for g in model.sG if (g, b) in model.sGbs
        ) + sum(model.qG[g, t] for g in model.G if (g, b) in model.Gbs) == sum(
            model.qD[d, t] for d in model.D if (b, d) in model.Dbs
        ) + sum(
            model.qLfrom[l, t] for l in model.L if model.A[l, 1] == b
        ) + sum(model.qLto[l, t] for l in model.L if model.A[l, 2] == b) + sum(
            model.qThv[l, t] for l in model.TRANSF if model.AT[l, 1] == b
        ) + sum(
            model.qTlv[l, t] for l in model.TRANSF if model.AT[l, 2] == b
        ) - sum(
            model.BB[s, t] * model.v[b, t] ** 2
            for s in model.SHUNT
            if (b, s) in model.SHUNTbs and model.BB[s, t] != 0
        ) + self.KCL_flexibility(model, b, t, reactive=True)
        if isinstance(kcl, bool):
            return Constraint.Skip
        return kcl

    # build_kcl / rebuild_kcl / KCL_flexibility come from
    # Basemodel_multi_period, shared with the LPAC and DC power flows.
