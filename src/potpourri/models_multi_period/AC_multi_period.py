"""Multi-period AC power flow mixin: adds full AC equations with voltage
magnitudes over time."""

from pyomo.environ import *
from potpourri.models_multi_period.basemodel_multi_period import (
    Basemodel_multi_period,
)
from potpourri.technologies.demand import Demand_multi_period
from potpourri.technologies.evs import EV_multi_period


class AC_multi_period(Basemodel_multi_period):
    """Multi-period AC power flow model with full AC equations indexed over
    time steps."""

    def __init__(self, net, toT, fromT=None, pf=1, num_vehicles=None):
        super().__init__(net, toT, fromT, pf, num_vehicles)

        self.BB_data = (
            -self.net.shunt.q_mvar * self.net.shunt.step / self.baseMVA
        )

        # line and transformer admittances
        r = self.net._ppc["branch"][:, 2].real
        x = self.net._ppc["branch"][:, 3].real
        y = self.net._ppc["branch"][:, 4] * 1j
        gt_ik = r / (r**2 + x**2)
        bt_ik = -x / (r**2 + x**2)
        BiiT = bt_ik + y.imag / 2
        BikT = -bt_ik
        GiiT = gt_ik + y.real / 2
        GikT = -gt_ik
        trafo_start = len(self.net.line)
        trafo_end = trafo_start + len(self.net.trafo)

        self.line_data = self.line_data.assign(
            **{
                "Bii_data": BiiT[:trafo_start],
                "Bik_data": BikT[:trafo_start],
                "Gii_data": GiiT[:trafo_start],
                "Gik_data": GikT[:trafo_start],
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
        """Build the multi-period AC Pyomo model with admittance parameters
        and AC KCL/KVL constraints."""
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
        # extended by EV charging/discharging
        def KCL_real_def(model, b, t):
            kcl = sum(
                model.psG[g, t] for g in model.sG if (g, b) in model.sGbs
            ) + sum(
                model.pG[g, t] for g in model.G if (g, b) in model.Gbs
            ) == sum(
                model.pD[d, t] for d in model.D if (b, d) in model.Dbs
            ) + sum(
                model.pLfrom[l, t] for l in model.L if model.A[l, 1] == b
            ) + sum(
                model.pLto[l, t] for l in model.L if model.A[l, 2] == b
            ) + sum(
                model.pThv[l, t] for l in model.TRANSF if model.AT[l, 1] == b
            ) + sum(
                model.pTlv[l, t] for l in model.TRANSF if model.AT[l, 2] == b
            ) + sum(
                model.GB[s] * model.v[b, t] ** 2
                for s in model.SHUNT
                if (b, s) in model.SHUNTbs and model.GB[s] != 0
            ) + self.KCL_flexibility(model, b, t)
            if isinstance(kcl, bool):
                return Constraint.Skip
            return kcl

        def KCL_reactive_def(model, b, t):
            kcl = sum(
                model.qsG[g, t] for g in model.sG if (g, b) in model.sGbs
            ) + sum(
                model.qG[g, t] for g in model.G if (g, b) in model.Gbs
            ) == sum(
                model.qD[d, t] for d in model.D if (b, d) in model.Dbs
            ) + sum(
                model.qLfrom[l, t] for l in model.L if model.A[l, 1] == b
            ) + sum(
                model.qLto[l, t] for l in model.L if model.A[l, 2] == b
            ) + sum(
                model.qThv[l, t] for l in model.TRANSF if model.AT[l, 1] == b
            ) + sum(
                model.qTlv[l, t] for l in model.TRANSF if model.AT[l, 2] == b
            ) - sum(
                model.BB[s] * model.v[b, t] ** 2
                for s in model.SHUNT
                if (b, s) in model.SHUNTbs and model.BB[s] != 0
            )
            if isinstance(kcl, bool):
                return Constraint.Skip
            return kcl

        self.model.KCL_real = Constraint(
            self.model.B, self.model.T, rule=KCL_real_def
        )
        self.model.KCL_reactive = Constraint(
            self.model.B, self.model.T, rule=KCL_reactive_def
        )

        # --- Kirchoff's voltage law on each line ---
        def KVL_real_fromend(model, l, t):
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

    def KCL_flexibility(self, model, b, t):
        """Return the net power injection from flexible assets (EVs) at bus b
        at time t."""
        power_sum = 0
        self.EV_object = next(
            (
                obj
                for obj in self.flexibilities
                if isinstance(obj, EV_multi_period)
            ),
            None,
        )

        if self.EV_object is not None:
            power_sum += sum(
                model.p_opf[v, t]
                for v in model.veh
                if (b, model.veh_cps[v, t]) in model.Bcp
            )
            # += sum(
            #     model.p_discharging[v, t] for v in model.veh
            #     if (b, model.veh_cps[v, t]) in model.Bcp
            # )

        return power_sum
