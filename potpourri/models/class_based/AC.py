from pyomo.environ import *
from potpourri.models.class_based.basemodel import Basemodel
import numpy as np


class AC(Basemodel):
    def __init__(self, net):
        super().__init__(net)

        self.BB_data = - self.net.shunt.q_mvar * net.shunt.step / self.baseMVA

        ZN = net.bus.vn_kv ** 2 / self.baseMVA

        b_ik_data = -net.line.x_ohm_per_km / (
                (net.line.x_ohm_per_km ** 2 + net.line.r_ohm_per_km ** 2) * net.line.length_km) * ZN[
                        net.line.from_bus].values
        b_ik_C_data = 2 * np.pi * net.f_hz * net.line.c_nf_per_km * 1e-9 * net.line.length_km * ZN[
            net.line.from_bus].values
        g_ik_data = net.line.r_ohm_per_km / (
                (net.line.x_ohm_per_km ** 2 + net.line.r_ohm_per_km ** 2) * net.line.length_km) * ZN[
                        net.line.from_bus].values
        self.Bii_data = b_ik_data + 0.5 * b_ik_C_data
        self.Bik_data = -b_ik_data
        self.Gii_data = g_ik_data
        self.Gik_data = -g_ik_data

        r_k = net.trafo.vkr_percent / 100 * (net.sn_mva / net.trafo.sn_mva)
        z_k = net.trafo.vk_percent / 100 * (net.sn_mva / net.trafo.sn_mva)
        x_k = np.sqrt(z_k ** 2 - r_k ** 2)
        gt_ik = r_k / z_k ** 2
        bt_ik = -x_k / z_k ** 2
        yt_m = net.trafo.i0_percent / 100  # magnetising addmittance
        gt_m = net.trafo.pfe_kw / (net.trafo.sn_mva * 1000) * net.sn_mva / net.trafo.sn_mva

        # TODO: calculation of magnetising addmittance trafo
        # bt_m = np.sqrt(yt_m ** 2 - gt_m ** 2)
        bt_m = 0
        Zt_ref = net.trafo.vn_lv_kv ** 2 * net.sn_mva / net.trafo.sn_mva

        self.BiiT_data = (bt_ik + 0.5 * bt_m) * ZN[net.trafo.lv_bus].values / Zt_ref
        self.BikT_data = -bt_ik * ZN[net.trafo.lv_bus].values / Zt_ref
        self.GiiT_data = gt_ik * ZN[net.trafo.lv_bus].values / Zt_ref
        self.GikT_data = -gt_ik * ZN[net.trafo.lv_bus].values / Zt_ref

        self.get_generator_reactive_data()
        self.get_demand_reactive_data()

    def get_generator_reactive_data(self):
        # generation set points
        self.QG_data = self.net.sgen.q_mvar / self.baseMVA

        # use active power for reactive power limits, if no reactive power given for any sgen
        if self.net.sgen.q_mvar.sum() == 0:
            lim_q = self.net.sgen.p_mw
        else:
            lim_q = self.net.sgen.q_mvar

        # add rows with reactive generation limits if not existing
        if 'max_q_mvar' not in self.net.sgen:
            self.net.sgen['max_q_mvar'] = lim_q

        if 'min_q_mvar' not in self.net.sgen:
            self.net.sgen['min_q_mvar'] = -lim_q

        self.QGmax_data = self.net.sgen.max_q_mvar.fillna(lim_q) / self.baseMVA
        self.QGmin_data = self.net.sgen.min_q_mvar.fillna(-lim_q) / self.baseMVA

    def get_demand_reactive_data(self):
        # reactive power demand
        self.QD_data = self.net.load.q_mvar / self.baseMVA

        # use active power for reactive power limits, if no reactive power given for any sgen
        if self.net.load.q_mvar.sum() == 0:
            lim_q = self.net.load.p_mw
        else:
            lim_q = self.net.load.q_mvar

        # add rows with reactive generation limits if not existing
        if 'max_q_mvar' not in self.net.load:
            self.net.load['max_q_mvar'] = lim_q

        if 'min_q_mvar' not in self.net.load:
            self.net.load['min_q_mvar'] = -lim_q

        # demand limits for loads
        self.QDmax_data = self.net.load.max_q_mvar.fillna(self.net.load.q_mvar) / self.baseMVA
        self.QDmin_data = self.net.load.min_q_mvar.fillna(0) / self.baseMVA

    def create_model(self):
        super().create_model()

        # shunt
        self.model.BB = Param(self.model.SHUNT, within=Reals, initialize=self.BB_data)  # shunt susceptance

        # derived line parameters
        self.model.Bii = Param(self.model.L, within=Reals, initialize=self.Bii_data)
        self.model.Bik = Param(self.model.L, within=Reals, initialize=self.Bik_data)
        self.model.Gii = Param(self.model.L, within=Reals, initialize=self.Gii_data)
        self.model.Gik = Param(self.model.L, within=Reals, initialize=self.Gik_data)

        ## derived transformer parameters
        self.model.BiiT = Param(self.model.TRANSF, within=Reals, initialize=self.BiiT_data)
        self.model.BikT = Param(self.model.TRANSF, within=Reals, initialize=self.BikT_data)
        self.model.GiiT = Param(self.model.TRANSF, within=Reals, initialize=self.GiiT_data)
        self.model.GikT = Param(self.model.TRANSF, within=Reals, initialize=self.GikT_data)

        # reactive generation
        self.model.QG = Param(self.model.G, initialize=self.QG_data)
        self.model.QGmax = Param(self.model.G, initialize=self.QGmax_data)
        self.model.QGmin = Param(self.model.G, initialize=self.QGmin_data)

        # reactive demand
        self.model.QD = Param(self.model.D, initialize=self.QD_data)
        self.model.QDmax = Param(self.model.D, initialize=self.QDmax_data)
        self.model.QDmin = Param(self.model.D, initialize=self.QDmin_data)

        # --- control variables ---
        self.model.qG = Var(self.model.G, domain=Reals) # reactive generator power

        self.model.qD = Var(self.model.D, domain=Reals)  # reactive power absorbed by demand

        self.model.qLfrom = Var(self.model.L, domain=Reals)  # reactive power injected at b onto line
        self.model.qLto = Var(self.model.L, domain=Reals)  # reactive power injected at b' onto line
        self.model.qLfromT = Var(self.model.TRANSF, domain=Reals)  # reactive power injected at b onto transformer
        self.model.qLtoT = Var(self.model.TRANSF, domain=Reals)  # reactive power injected at b' onto transformer

        self.model.v = Var(self.model.B, domain=NonNegativeReals, initialize=1.0)  # voltage magnitude at bus b, rad

        self.model.qeG = Var(self.model.eG, domain=Reals)

        # --- Kirchoff's current law at each bus b ---
        def KCL_real_def(model, b):
            return sum(model.pG[g] for g in model.G if (b, g) in model.Gbs) + \
                sum(model.peG[g] for g in model.eG if (b, g) in model.eGbs) == \
                sum(model.pD[d] for d in model.D if (b, d) in model.Dbs) + \
                sum(model.pLfrom[l] for l in model.L if model.A[l, 1] == b) + \
                sum(model.pLto[l] for l in model.L if model.A[l, 2] == b) + \
                sum(model.pLfromT[l] for l in model.TRANSF if model.AT[l, 1] == b) + \
                sum(model.pLtoT[l] for l in model.TRANSF if model.AT[l, 2] == b) + \
                sum(model.GB[s] * model.v[b] ** 2 for s in model.SHUNT if (b, s) in model.SHUNTbs)

        def KCL_reactive_def(model, b):
            return sum(model.qG[g] for g in model.G if (b, g) in model.Gbs) + \
                sum(model.qeG[g] for g in model.eG if (b, g) in model.eGbs) == \
                sum(model.qD[d] for d in model.D if (b, d) in model.Dbs) + \
                sum(model.qLfrom[l] for l in model.L if model.A[l, 1] == b) + \
                sum(model.qLto[l] for l in model.L if model.A[l, 2] == b) + \
                sum(model.qLfromT[l] for l in model.TRANSF if model.AT[l, 1] == b) + \
                sum(model.qLtoT[l] for l in model.TRANSF if model.AT[l, 2] == b) - \
                sum(model.BB[s] * model.v[b] ** 2 for s in model.SHUNT if (b, s) in model.SHUNTbs)

        self.model.KCL_real = Constraint(self.model.B, rule=KCL_real_def)
        self.model.KCL_reactive = Constraint(self.model.B, rule=KCL_reactive_def)

        # --- Kirchoff's voltage law on each line ---
        def KVL_real_fromend(model, l):
            return model.pLfrom[l] == model.Gii[l] * (model.v[model.A[l, 1]] ** 2) + \
                model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                        model.Bik[l] * sin(model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]) +
                        model.Gik[l] * cos(model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]))

        def KVL_real_toend(model, l):
            return model.pLto[l] == model.Gii[l] * (model.v[model.A[l, 2]] ** 2) + \
                model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                        model.Bik[l] * sin(model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]) +
                        model.Gik[l] * cos(model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]))

        def KVL_reactive_fromend(model, l):
            return model.qLfrom[l] == -model.Bii[l] * (model.v[model.A[l, 1]] ** 2) + \
                model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                        model.Gik[l] * sin(model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]) -
                        model.Bik[l] * cos(model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]))

        def KVL_reactive_toend(model, l):
            return model.qLto[l] == -model.Bii[l] * (model.v[model.A[l, 2]] ** 2) + \
                model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                        model.Gik[l] * sin(model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]) -
                        model.Bik[l] * cos(model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]))

        self.model.KVL_real_from = Constraint(self.model.L, rule=KVL_real_fromend)
        self.model.KVL_real_to = Constraint(self.model.L, rule=KVL_real_toend)
        self.model.KVL_reactive_from = Constraint(self.model.L, rule=KVL_reactive_fromend)
        self.model.KVL_reactive_to = Constraint(self.model.L, rule=KVL_reactive_toend)

        # --- reactive generator power limits ---
        def reactive_power_bounds(model, g):
            if g in model.Gc:
                return model.QGmin[g], model.qG[g], model.QGmax[g]
            else:
                model.qG[g].fix(model.QG[g])
                return Constraint.Skip

        self.model.QG_Constraint = Constraint(self.model.G, rule=reactive_power_bounds)

        # --- reactive demand limits ---
        def reactive_demand_bounds(model, d):
            if d in model.Dc:
                return model.QDmin[d], model.qD[d], model.QDmax[d]
            else:
                model.qD[d].fix(model.QD[d])
                return Constraint.Skip
        self.model.QD_Constraint = Constraint(self.model.D, rule=reactive_demand_bounds)

        # --- reference bus voltage constraint ---
        def ref_bus_def(model, b, g):
            return model.v[b] == 1

        self.model.refbus_v = Constraint(self.model.eGbs, rule=ref_bus_def)

