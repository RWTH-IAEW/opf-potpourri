"""AC power flow mixin: adds full AC equations (voltage magnitudes, reactive
power) to Basemodel."""

import numpy as np
import pyomo.environ as pyo
from potpourri.models.basemodel import Basemodel


class AC(Basemodel):
    """AC power flow model for distribution network analysis.

    Extends Basemodel with full AC power flow equations including voltage
    magnitudes and reactive power. Adds real and reactive KCL at each bus,
    and real and reactive KVL on each line and transformer.

    Args:
        net: A pandapower network compatible with pp.runpp().
    """

    def __init__(self, net):
        super().__init__(net)

        self.BB_data = (
            -self.net.shunt.q_mvar * self.net.shunt.step / self.baseMVA
        )

        # --- line and transformer admittances (symmetric π model) ---
        # The branch series admittance is y_s = 1/(r + jx) = g + jb. We split
        # the branch into a from/to shunt of half the total line-charging
        # susceptance b_c and a mutual series term y_s.
        #
        # ASSUMPTION (D11 in the formulation audit): we read y from
        # ``_ppc["branch"][:, 4]`` and treat it as purely imaginary
        # (``y = j·b_c``), so no branch-shunt CONDUCTANCE is modelled. This
        # matches the standard MATPOWER convention (BR_B is susceptance
        # only; MATPOWER has no BR_G column) and is fine for PGLib cases
        # and pandapower's standard ``_ppc`` build. Magnetising-loss
        # transformers with a non-zero g_m would be silently approximated
        # by g_m = 0 here.
        r = self.net._ppc["branch"][:, 2].real
        x = self.net._ppc["branch"][:, 3].real
        y = self.net._ppc["branch"][:, 4] * 1j  # j·b_c (no branch shunt G)
        gt_ik = r / (r**2 + x**2)  # series conductance g
        bt_ik = -x / (r**2 + x**2)  # series susceptance b (b<0 inductive)
        BiiT = bt_ik + y.imag / 2  # self susceptance: b + b_c/2
        BikT = -bt_ik  # mutual susceptance: -b
        GiiT = gt_ik + y.real / 2  # self conductance: g + g_c/2 (g_c==0)
        GikT = -gt_ik  # mutual conductance: -g
        trafo_start = len(self.net.line)
        trafo_end = trafo_start + len(self.net.trafo)
        imp_table = self.net.get("impedance")
        n_imp = (
            len(imp_table)
            if imp_table is not None and not imp_table.empty
            else 0
        )
        # _ppc['branch'] layout (pandapower convention):
        #   [0 : trafo_start)              → lines
        #   [trafo_start : trafo_end)      → trafos
        #   [trafo_end : trafo_end + n_imp)→ impedance branches
        # We treat impedance rows as additional lines in the model — see
        # the matching block in Basemodel.__init__.
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

        # generator and external grids voltage set points
        self.generation_data["v"] = self.net._ppc["gen"][:, 5]

        qsg = (
            self.net.sgen.q_mvar.fillna(0)
            * self.net.sgen.scaling
            / self.baseMVA
        ).values
        self.static_generation_data["q"] = qsg

        self.QD_data = (
            self.net.load.q_mvar * self.net.load.scaling / self.baseMVA
        )

        self.create_model()

    def create_model(self):
        """Build the Pyomo ConcreteModel with AC power flow constraints.

        Adds line and transformer admittance parameters (Bii, Bik, Gii, Gik),
        voltage magnitude variables, reactive power variables (qsG, qD, qLfrom,
        qLto, qThv, qTlv, qG), and KCL/KVL constraints for real and reactive
        power on every bus, line, and transformer.
        """
        super().create_model()

        self.model.name = "AC"

        # shunt
        self.model.BB = pyo.Param(
            self.model.SHUNT,
            within=pyo.Reals,
            initialize=self.BB_data[self.model.SHUNT],
        )  # shunt susceptance

        # derived line pyo.Parameters
        self.model.Bii = pyo.Param(
            self.model.L,
            within=pyo.Reals,
            initialize=self.line_data.Bii_data[self.model.L],
        )
        self.model.Bik = pyo.Param(
            self.model.L,
            within=pyo.Reals,
            initialize=self.line_data.Bik_data[self.model.L],
        )
        self.model.Gii = pyo.Param(
            self.model.L,
            within=pyo.Reals,
            initialize=self.line_data.Gii_data[self.model.L],
        )
        self.model.Gik = pyo.Param(
            self.model.L,
            within=pyo.Reals,
            initialize=self.line_data.Gik_data[self.model.L],
        )

        ## derived transformer pyo.Parameters
        self.model.BiiT = pyo.Param(
            self.model.TRANSF,
            within=pyo.Reals,
            initialize=self.trafo_data.BiiT_data[self.model.TRANSF],
        )
        self.model.BikT = pyo.Param(
            self.model.TRANSF,
            within=pyo.Reals,
            initialize=self.trafo_data.BikT_data[self.model.TRANSF],
        )
        self.model.GiiT = pyo.Param(
            self.model.TRANSF,
            within=pyo.Reals,
            initialize=self.trafo_data.GiiT_data[self.model.TRANSF],
        )
        self.model.GikT = pyo.Param(
            self.model.TRANSF,
            within=pyo.Reals,
            initialize=self.trafo_data.GikT_data[self.model.TRANSF],
        )

        # reactive generation
        self.model.QsG = pyo.Param(
            self.model.sG,
            initialize=self.static_generation_data["q"][self.model.sG],
        )

        self.model.v_bPV = pyo.Param(
            self.model.bPV,
            within=pyo.NonNegativeReals,
            initialize=self.bus_data.v_m[self.model.bPV],
        )

        # reactive demand
        self.model.QD = pyo.Param(
            self.model.D, initialize=self.QD_data[self.model.D]
        )

        # external grid voltage
        self.model.v_b0 = pyo.Param(
            self.model.b0,
            within=pyo.NonNegativeReals,
            initialize=self.bus_data.v_m[self.model.b0],
        )

        # --- control pyo.Variables ---
        self.model.qsG = pyo.Var(
            self.model.sG, domain=pyo.Reals
        )  # reactive power of static generators

        self.model.qD = pyo.Var(
            self.model.D, domain=pyo.Reals
        )  # reactive power absorbed by demand

        self.model.qLfrom = pyo.Var(
            self.model.L, domain=pyo.Reals
        )  # reactive power injected at b onto line
        self.model.qLto = pyo.Var(
            self.model.L, domain=pyo.Reals
        )  # reactive power injected at b' onto line
        self.model.qThv = pyo.Var(
            self.model.TRANSF, domain=pyo.Reals
        )  # reactive power injected at b onto transformer
        self.model.qTlv = pyo.Var(
            self.model.TRANSF, domain=pyo.Reals
        )  # reactive power injected at b' onto transformer

        v_init = self.bus_data.v_m.fillna(1.0).to_dict()
        self.model.v = pyo.Var(
            self.model.B,
            domain=pyo.NonNegativeReals,
            bounds=(0.0, 2.0),
            initialize=v_init,
        )  # voltage magnitude at bus b (p.u.)

        self.model.qG = pyo.Var(self.model.G, domain=pyo.Reals)

        # --- Kirchoff's current law at each bus b ---
        def KCL_real_def(model, b):
            kcl = sum(
                model.psG[g] for g in model.sG if (g, b) in model.sGbs
            ) + sum(model.pG[g] for g in model.G if (g, b) in model.Gbs) - sum(
                model.pSTOR[s] for s in model.STOR if model.STOR_bus[s] == b
            ) == sum(
                model.pD[d] for d in model.D if (b, d) in model.Dbs
            ) + sum(
                model.pLfrom[l] for l in model.L if model.A[l, 1] == b
            ) + sum(
                model.pLto[l] for l in model.L if model.A[l, 2] == b
            ) + sum(
                model.pThv[l] for l in model.TRANSF if model.AT[l, 1] == b
            ) + sum(
                model.pTlv[l] for l in model.TRANSF if model.AT[l, 2] == b
            ) + sum(
                model.GB[s] * model.v[b] ** 2
                for s in model.SHUNT
                if (b, s) in model.SHUNTbs and model.GB[s] != 0
            )
            if isinstance(kcl, bool):
                return pyo.Constraint.Skip
            return kcl

        def KCL_reactive_def(model, b):
            kcl = sum(
                model.qsG[g] for g in model.sG if (g, b) in model.sGbs
            ) + sum(model.qG[g] for g in model.G if (g, b) in model.Gbs) + sum(
                model.qSTOR[s] for s in model.STOR if model.STOR_bus[s] == b
            ) == sum(
                model.qD[d] for d in model.D if (b, d) in model.Dbs
            ) + sum(
                model.qLfrom[l] for l in model.L if model.A[l, 1] == b
            ) + sum(
                model.qLto[l] for l in model.L if model.A[l, 2] == b
            ) + sum(
                model.qThv[l] for l in model.TRANSF if model.AT[l, 1] == b
            ) + sum(
                model.qTlv[l] for l in model.TRANSF if model.AT[l, 2] == b
            ) - sum(
                model.BB[s] * model.v[b] ** 2
                for s in model.SHUNT
                if (b, s) in model.SHUNTbs and model.BB[s] != 0
            )
            if isinstance(kcl, bool):
                return pyo.Constraint.Skip
            return kcl

        self.model.KCL_real = pyo.Constraint(self.model.B, rule=KCL_real_def)
        self.model.KCL_reactive = pyo.Constraint(
            self.model.B, rule=KCL_reactive_def
        )

        # --- Kirchoff's voltage law on each line ---
        def KVL_real_fromend(model, l):
            return model.pLfrom[l] == model.Gii[l] * (
                model.v[model.A[l, 1]] ** 2
            ) + model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                model.Bik[l]
                * pyo.sin(
                    model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]
                )
                + model.Gik[l]
                * pyo.cos(
                    model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]
                )
            )

        def KVL_real_toend(model, l):
            return model.pLto[l] == model.Gii[l] * (
                model.v[model.A[l, 2]] ** 2
            ) + model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                model.Bik[l]
                * pyo.sin(
                    model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]
                )
                + model.Gik[l]
                * pyo.cos(
                    model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]
                )
            )

        def KVL_reactive_fromend(model, l):
            return model.qLfrom[l] == -model.Bii[l] * (
                model.v[model.A[l, 1]] ** 2
            ) + model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                model.Gik[l]
                * pyo.sin(
                    model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]
                )
                - model.Bik[l]
                * pyo.cos(
                    model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]
                )
            )

        def KVL_reactive_toend(model, l):
            return model.qLto[l] == -model.Bii[l] * (
                model.v[model.A[l, 2]] ** 2
            ) + model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                model.Gik[l]
                * pyo.sin(
                    model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]
                )
                - model.Bik[l]
                * pyo.cos(
                    model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]
                )
            )

        self.model.KVL_real_from = pyo.Constraint(
            self.model.L, rule=KVL_real_fromend
        )
        self.model.KVL_real_to = pyo.Constraint(
            self.model.L, rule=KVL_real_toend
        )
        self.model.KVL_reactive_from = pyo.Constraint(
            self.model.L, rule=KVL_reactive_fromend
        )
        self.model.KVL_reactive_to = pyo.Constraint(
            self.model.L, rule=KVL_reactive_toend
        )

        # --- Kirchoff's voltage law on each transformer line ---
        def KVL_real_fromendTransf(model, l):
            if model.shift[l]:
                return model.pThv[l] == model.GiiT[l] / model.Tap[l] ** 2 * (
                    model.v[model.AT[l, 1]] ** 2
                ) + model.v[model.AT[l, 1]] * model.v[
                    model.AT[l, 2]
                ] / model.Tap[l] * (
                    model.GikT[l]
                    * pyo.cos(
                        model.delta[model.AT[l, 1]]
                        - model.delta[model.AT[l, 2]]
                        - model.shift[l]
                    )
                    + model.BikT[l]
                    * pyo.sin(
                        model.delta[model.AT[l, 1]]
                        - model.delta[model.AT[l, 2]]
                        - model.shift[l]
                    )
                )

            return model.pThv[l] == model.GiiT[l] / model.Tap[l] ** 2 * (
                model.v[model.AT[l, 1]] ** 2
            ) + model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[
                l
            ] * (
                model.GikT[l]
                * pyo.cos(
                    model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]
                )
                + model.BikT[l]
                * pyo.sin(
                    model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]
                )
            )

        def KVL_real_toendTransf(model, l):
            if model.shift[l]:
                return model.pTlv[l] == model.GiiT[l] * (
                    model.v[model.AT[l, 2]] ** 2
                ) + model.v[model.AT[l, 1]] * model.v[
                    model.AT[l, 2]
                ] / model.Tap[l] * (
                    model.BikT[l]
                    * pyo.sin(
                        model.delta[model.AT[l, 2]]
                        - model.delta[model.AT[l, 1]]
                        + model.shift[l]
                    )
                    + model.GikT[l]
                    * pyo.cos(
                        model.delta[model.AT[l, 2]]
                        - model.delta[model.AT[l, 1]]
                        + model.shift[l]
                    )
                )

            return model.pTlv[l] == model.GiiT[l] * (
                model.v[model.AT[l, 2]] ** 2
            ) + model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[
                l
            ] * (
                model.BikT[l]
                * pyo.sin(
                    model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]
                )
                + model.GikT[l]
                * pyo.cos(
                    model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]
                )
            )

        def KVL_reactive_fromendTransf(model, l):
            if model.shift[l]:
                return model.qThv[l] == -model.BiiT[l] / model.Tap[l] ** 2 * (
                    model.v[model.AT[l, 1]] ** 2
                ) + model.v[model.AT[l, 1]] * model.v[
                    model.AT[l, 2]
                ] / model.Tap[l] * (
                    -model.BikT[l]
                    * pyo.cos(
                        model.delta[model.AT[l, 1]]
                        - model.delta[model.AT[l, 2]]
                        - model.shift[l]
                    )
                    + model.GikT[l]
                    * pyo.sin(
                        model.delta[model.AT[l, 1]]
                        - model.delta[model.AT[l, 2]]
                        - model.shift[l]
                    )
                )

            return model.qThv[l] == -model.BiiT[l] / model.Tap[l] ** 2 * (
                model.v[model.AT[l, 1]] ** 2
            ) + model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[
                l
            ] * (
                -model.BikT[l]
                * pyo.cos(
                    model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]
                )
                + model.GikT[l]
                * pyo.sin(
                    model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]
                )
            )

        def KVL_reactive_toendTransf(model, l):
            if model.shift[l]:
                return model.qTlv[l] == -model.BiiT[l] * (
                    model.v[model.AT[l, 2]] ** 2
                ) + model.v[model.AT[l, 1]] * model.v[
                    model.AT[l, 2]
                ] / model.Tap[l] * (
                    -model.BikT[l]
                    * pyo.cos(
                        model.delta[model.AT[l, 2]]
                        - model.delta[model.AT[l, 1]]
                        + model.shift[l]
                    )
                    + model.GikT[l]
                    * pyo.sin(
                        model.delta[model.AT[l, 2]]
                        - model.delta[model.AT[l, 1]]
                        + model.shift[l]
                    )
                )

            return model.qTlv[l] == -model.BiiT[l] * (
                model.v[model.AT[l, 2]] ** 2
            ) + model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] / model.Tap[
                l
            ] * (
                -model.BikT[l]
                * pyo.cos(
                    model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]
                )
                + model.GikT[l]
                * pyo.sin(
                    model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]
                )
            )

        self.model.KVL_real_fromTransf = pyo.Constraint(
            self.model.TRANSF, rule=KVL_real_fromendTransf
        )
        self.model.KVL_real_toTransf = pyo.Constraint(
            self.model.TRANSF, rule=KVL_real_toendTransf
        )
        self.model.KVL_reactive_fromTransf = pyo.Constraint(
            self.model.TRANSF, rule=KVL_reactive_fromendTransf
        )
        self.model.KVL_reactive_toTransf = pyo.Constraint(
            self.model.TRANSF, rule=KVL_reactive_toendTransf
        )

        # --- reactive demand limits ---
        for d in self.model.D:
            self.model.qD[d].fix(self.model.QD[d])

        # --- generator voltage operating point ---
        def v_bPV_setpoint_rule(model, b):
            return model.v[b] == model.v_bPV[b]

        self.model.v_bPV_setpoint = pyo.Constraint(
            self.model.bPV, rule=v_bPV_setpoint_rule
        )

        # --- reference bus voltage pyo.Constraint ---
        for b0 in self.model.b0:
            self.model.v[b0].fix(self.model.v_b0[b0])
