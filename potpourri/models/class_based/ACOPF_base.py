from pyomo.environ import *
from potpourri.models.class_based.AC import AC
from potpourri.models.class_based.OPF import OPF
import pandas as pd


class ACOPF(AC, OPF):
    def __init__(self, net):
        super().__init__(net)

        if 'max_vm_pu' not in self.net.bus:
            self.net.bus['max_vm_pu'] = 1.1
        if 'min_vm_pu' not in net.bus:
            self.net.bus['min_vm_pu'] = 0.9
        self.vmax_data = self.net.bus.max_vm_pu.values
        self.vmin_data = self.net.bus.min_vm_pu.values

        self.get_generator_reactive_data()
        self.get_demand_reactive_data()
        self.get_generator_v_limits()

    def get_generator_reactive_data(self):
        # use active power for reactive power limits, if no reactive power given for any sgen
        if self.generators.q_mvar.sum() == 0:
            lim_q = abs(self.generators.p_mw)
        else:
            lim_q = abs(self.generators.q_mvar)

        # add rows with reactive generation limits if not existing
        if 'max_q_mvar' not in self.generators:
            self.generators['max_q_mvar'] = lim_q

        if 'min_q_mvar' not in self.generators:
            self.generators['min_q_mvar'] = -lim_q

        self.QGmax_data = self.generators.max_q_mvar.fillna(lim_q) / self.baseMVA
        self.QGmin_data = self.generators.min_q_mvar.fillna(-lim_q) / self.baseMVA

    def get_generator_v_limits(self):
        if self.gen_controllable_set.any():
            # gG_controllable_ind = self.gen_controllable_set.intersection(self.gen_set)
            # if gG_controllable_ind.any():
            if 'max_vm_pu' not in self.generators:
                self.v_max_data = pd.Series(
                    self.net.bus.max_vm_pu[self.generators.bus[self.gen_controllable_set]].values,
                    self.gen_controllable_set)
            else:
                self.v_max_data = self.net.gen.max_vm_pu.fillna(
                    pd.Series(self.net.bus.max_vm_pu[self.net.gen.bus].values, self.net.gen.index))

            if 'min_vm_pu' not in self.generators:
                self.v_min_data = pd.Series(self.net.bus.min_vm_pu[self.generators.bus[self.gen_controllable_set]].values, self.gen_controllable_set)
            else:
                self.v_min_data = self.net.gen.min_vm_pu.fillna(
                        pd.Series(self.net.bus.min_vm_pu[self.net.gen.bus[self.gen_controllable_set]].values, self.net.gen.index))
        else:
            self.v_max_data = pd.Series()
            self.v_min_data = pd.Series()

    def get_demand_reactive_data(self):
        # reactive power demand
        # use active power for reactive power limits, if no reactive power given for any sgen
        if self.net.load.q_mvar.sum() == 0:
            lim_q = abs(self.net.load.p_mw)
        else:
            lim_q = abs(self.net.load.q_mvar)

        # add rows with reactive generation limits if not existing
        if 'max_q_mvar' not in self.net.load:
            self.net.load['max_q_mvar'] = lim_q

        if 'min_q_mvar' not in self.net.load:
            self.net.load['min_q_mvar'] = -lim_q

        # demand limits for loads
        self.QDmax_data = self.net.load.max_q_mvar.fillna(self.net.load.q_mvar) / self.baseMVA
        self.QDmin_data = self.net.load.min_q_mvar.fillna(0) / self.baseMVA

    def add_OPF(self):
        super().add_OPF()

        self.model.name = "ACOPF"

        self.model.Vmax = Param(self.model.B, within=NonNegativeReals,
                                initialize=self.vmax_data[self.model.B])  # max voltage (p.u.)
        self.model.Vmin = Param(self.model.B, within=NonNegativeReals,
                                initialize=self.vmin_data[self.model.B])  # min voltage (p.u.)

        # reactive generation
        self.model.QGmax = Param(self.model.G, initialize=self.QGmax_data[self.model.G])
        self.model.QGmin = Param(self.model.G, initialize=self.QGmin_data[self.model.G])

        # reactive demand
        self.model.QDmax = Param(self.model.D, initialize=self.QDmax_data[self.model.D])
        self.model.QDmin = Param(self.model.D, initialize=self.QDmin_data[self.model.D])

        # TODO: combine with model.Vmax and model.Vmin
        self.model.v_gG_max = Param(self.model.Gc, within=NonNegativeReals, initialize=self.v_max_data[self.model.Gc])
        self.model.v_gG_min = Param(self.model.Gc, within=NonNegativeReals, initialize=self.v_min_data[self.model.Gc])

        # --- cost function ---
        def objective(model):
            obj = sum(
                model.c2[g] * (model.baseMVA * model.pG[g]) ** 2 + model.c1[g] * model.baseMVA * model.pG[g] +
                model.c0[g] for g in model.G) + \
                  sum(model.VOLL[d] * (model.PD[d] - model.pD[d]) * model.baseMVA for d in model.D)
            return obj

        self.model.OBJ = Objective(rule=objective, sense=minimize)

        # --- line power limits ---
        def line_lim_from_def(model, l):
            return model.pLfrom[l] ** 2 + model.qLfrom[l] ** 2 <= model.SLmax[l] ** 2 * model.v[model.A[l, 1]] ** 2

        def line_lim_to_def(model, l):
            return model.pLto[l] ** 2 + model.qLto[l] ** 2 <= model.SLmax[l] ** 2 * model.v[model.A[l, 2]] ** 2

        self.model.line_lim_from = Constraint(self.model.L, rule=line_lim_from_def)
        self.model.line_lim_to = Constraint(self.model.L, rule=line_lim_to_def)

        # --- power flow limits on transformer lines---
        def transf_lim1_def(model, l):
            return model.pThv[l] ** 2 + model.qThv[l] ** 2 <= model.SLmaxT[l] ** 2 * model.v[model.AT[l, 1]] ** 2

        def transf_lim2_def(model, l):
            return model.pTlv[l] ** 2 + model.qTlv[l] ** 2 <= model.SLmaxT[l] ** 2 * model.v[model.AT[l, 2]] ** 2

        self.model.transf_lim1 = Constraint(self.model.TRANSF, rule=transf_lim1_def)
        self.model.transf_lim2 = Constraint(self.model.TRANSF, rule=transf_lim2_def)

        # --- reactive generator power limits ---
        def reactive_power_bounds(model, g):
            model.qG[g].unfix()
            return model.QGmin[g], model.qG[g], model.QGmax[g]

        self.model.QGc_Constraint = Constraint(self.model.Gc, rule=reactive_power_bounds)

        # --- reactive demand limits ---
        def reactive_demand_bounds(model, d):
            model.qD[d].unfix()
            return model.QDmin[d], model.qD[d], model.QDmax[d]

        self.model.QD_Constraint = Constraint(self.model.Dc, rule=reactive_demand_bounds)

        # --- voltage constraints ---
        def v_bounds(model, b):
            # voltage limits from controllable generators
            for g in model.gG.intersection(model.Gc):
                if (b, g) in model.Gbs:
                    model.v[b].unfix()
                    return model.v_gG_min[g], model.v[b], model.v_gG_max[g]

            for g in model.gG:
                if (b, g) in model.Gbs:
                    return Constraint.Skip

            return model.Vmin[b], model.v[b], model.Vmax[b]

        self.model.v_constraint = Constraint(self.model.B, rule=v_bounds)
