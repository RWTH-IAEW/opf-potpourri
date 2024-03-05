from pyomo.environ import *
from potpourri.models.class_based.AC import AC
from potpourri.models.class_based.OPF import OPF
import numpy as np


class ACOPF(AC, OPF):
    def __init__(self, net):
        super().__init__(net)

    def _calc_opf_parameters(self):
        super()._calc_opf_parameters()

        max_vm_pu, min_vm_pu = self.get_v_limits()
        self.v_limits = (max_vm_pu, min_vm_pu)

        self.static_generation_reactive_power_limits()
        self.generation_reactive_power_limits()
        self.get_demand_reactive_data()

    def static_generation_reactive_power_limits(self):
        if 'controllable' in self.net.sgen:
            self.static_generation_data['controllable'] = self.net.sgen.controllable.values
        else:
            self.static_generation_data['controllable'] = False

        lim_q = abs(self.net.sgen.q_mvar) / self.baseMVA
        if 'max_q_mvar' in self.net.sgen:
            self.static_generation_data['max_q'] = self.net.sgen.max_q_mvar.fillna(lim_q).values / self.baseMVA
        else:
            self.static_generation_data['max_q'] = lim_q / self.baseMVA

        if 'min_q_mvar' in self.net.sgen:
            self.static_generation_data['min_q'] = self.net.sgen.min_q_mvar.fillna(-lim_q).values / self.baseMVA
        else:
            self.static_generation_data['min_q'] = -lim_q / self.baseMVA

        if 'wind_hc' in self.net.sgen:
            self.static_generation_data['wind_hc'] = self.net.sgen.wind_hc.values
            self.static_generation_data['max_q'][self.static_generation_data['wind_hc']] = 0.48 * self.static_generation_data['p'][self.static_generation_data['wind_hc']]
            self.static_generation_data['min_q'][self.static_generation_data['wind_hc']] = -0.23 * self.static_generation_data['p'][self.static_generation_data['wind_hc']]
        else:
            self.static_generation_data['wind_hc'] = False

        self.static_generation_data['type'] = self.net.sgen.type.values
        if any(self.static_generation_data['type'] == 'Wind'):
            self.static_generation_data['max_q'][self.static_generation_data['type'] == 'Wind'] = 0.48 * self.static_generation_data['p'][self.static_generation_data['type'] == 'Wind']
            self.static_generation_data['min_q'][self.static_generation_data['type'] == 'Wind'] = -0.23 * self.static_generation_data['p'][self.static_generation_data['type'] == 'Wind']



    def generation_reactive_power_limits(self):
        max_q = np.full(len(self.generation_data), 1e9) / self.baseMVA
        min_q = np.full(len(self.generation_data), -1e9) / self.baseMVA

        for element, (f, t) in self.net._gen_order.items():
            if 'max_q_mvar' in self.net[element]:
                max_q[f:t] = self.net[element].max_q_mvar.fillna(1e9) / self.baseMVA
            if 'min_q_mvar' in self.net[element]:
                min_q[f:t] = self.net[element].min_q_mvar.fillna(-1e9) / self.baseMVA

        self.generation_data['max_q'] = max_q
        self.generation_data['min_q'] = min_q

    def get_v_limits(self):
        if 'max_vm_pu' in self.net.bus:
            max_vm_pu = self.net.bus.max_vm_pu.values
        else:
            max_vm_pu = np.full(len(self.net.bus.index), 1.1)

        if 'min_vm_pu' in self.net.bus:
            min_vm_pu = self.net.bus.min_vm_pu.values
        else:
            min_vm_pu = np.full(len(self.net.bus.index), 0.9)

        if any(self.net.gen.index):
            self.add_generator_v_limits(max_vm_pu, min_vm_pu)

        return max_vm_pu, min_vm_pu

    def add_generator_v_limits(self, max_vm_pu, min_vm_pu):
        # check max_vm_pu / min_vm_pu bus limit violation by gens
        gen_buses = self.bus_lookup[self.net.gen.bus.values]
        if "max_vm_pu" in self.net["gen"].columns:
            v_max_bound = max_vm_pu[gen_buses] < self.net["gen"]["max_vm_pu"].values
            if np.any(v_max_bound):
                bound_gens = self.net["gen"].index.values[v_max_bound]
                print("gen max_vm_pu > bus max_vm_pu for gens {}. "
                      "Setting bus limit for these gens.".format(bound_gens))
                # set only vm of gens which do not violate the limits
                max_vm_pu[gen_buses[~v_max_bound]] = self.net["gen"]["max_vm_pu"].values[~v_max_bound]
            else:
                # set vm of all gens
                max_vm_pu[gen_buses] = self.net["gen"]["max_vm_pu"].values

        if "min_vm_pu" in self.net["gen"].columns:
            v_min_bound = self.net["gen"]["min_vm_pu"].values < min_vm_pu[gen_buses]
            if np.any(v_min_bound):
                bound_gens = self.net["gen"].index.values[v_min_bound]
                print("gen min_vm_pu < bus min_vm_pu for gens {}. "
                      "Setting bus limit for these gens.".format(bound_gens))
                # set only vm of gens which do not violate the limits
                min_vm_pu[gen_buses[~v_min_bound]] = self.net["gen"]["min_vm_pu"].values[~v_min_bound]
            else:
                # set vm of all gens
                min_vm_pu[gen_buses] = self.net["gen"]["min_vm_pu"].values

        if 'controllable' in self.net.gen:
            controllable = self.net["gen"]["controllable"].values
            not_controllable = ~controllable.astype(bool)

            # get voltage setpoints for not controllable generators
            if np.any(not_controllable):
                bus = self.net["gen"]["bus"].values[not_controllable]
                vm_pu = self.net["gen"]["vm_pu"].values[not_controllable]

                not_controllable_buses = self.bus_lookup[bus]
                max_vm_pu[not_controllable_buses] = vm_pu
                min_vm_pu[not_controllable_buses] = vm_pu

        return max_vm_pu, min_vm_pu

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

    def hc_wind_generation(self):
            if 'var_q' in self.net.sgen:
                self.static_generation_data['var_q'] = self.net.sgen.q_var
            else:
                var = 0
                self.static_generation_data['var_q'] = var

            x = np.array([[96, 103], [120, 127]]) / 110
            y = np.array([[0.48, 0.41, 0.33], [-0.23, -0.33, -0.41]])
            m = (y[1] - y[0]) / (x[0, 1] - x[0, 0])
            b = np.array([y[0] - m * x[i, 0] for i in range(len(x))]).T
            self.m = m[0]
            self.b_min = b[0, 0]
            self.b_max = b[0, 1]
            self.static_generation_data['wind_slope'][self.wind_hc_set] = m[self.static_generation_data['var_q']][
                self.wind_hc_set]
            self.static_generation_data['wind_intercept_min'] = b[self.static_generation_data['var_q'], 0]
            self.static_generation_data['wind_intercept_max'] = b[self.static_generation_data['var_q'], 1]


    def add_OPF(self, **kwargs):
        super().add_OPF(**kwargs)

        self.model.name = "ACOPF"

        # --- sets ---
        self.model.WIND_HC = Set(within=self.model.sG, initialize=self.static_generation_data.index[
            self.static_generation_data['wind_hc'] & self.static_generation_data.in_service])
        self.model.WIND = self.model.WIND_HC | Set(within=self.model.sG, initialize=self.static_generation_data.index[(self.static_generation_data['type'] == 'Wind') & self.static_generation_data.in_service])
        self.model.WINDc = self.model.WIND & self.model.sGc

        # voltage limits
        self.model.Vmax = Param(self.model.B, within=NonNegativeReals,
                                initialize=self.v_limits[0][self.model.B])  # max voltage (p.u.)
        self.model.Vmin = Param(self.model.B, within=NonNegativeReals,
                                initialize=self.v_limits[1][self.model.B])  # min voltage (p.u.)

        # generation reactive power limits
        self.model.QGmax = Param(self.model.G, initialize=self.generation_data['max_q'][self.model.G])
        self.model.QGmin = Param(self.model.G, initialize=self.generation_data['min_q'][self.model.G])

        # static generation reactive power limits
        self.model.QsGmax = Param(self.model.sGc, within=Reals, initialize=self.static_generation_data['max_q'][self.model.sGc])
        self.model.QsGmin = Param(self.model.sGc, within=Reals, initialize=self.static_generation_data['min_q'][self.model.sGc])

        # reactive demand
        self.model.QDmax = Param(self.model.D, initialize=self.QDmax_data[self.model.D])
        self.model.QDmin = Param(self.model.D, initialize=self.QDmin_data[self.model.D])

        # --- cost function ---
        # def objective(model):
        #     obj = sum(
        #         model.c2[g] * (model.baseMVA * model.pG[g]) ** 2 + model.c1[g] * model.baseMVA * model.pG[g] +
        #         model.c0[g] for g in model.G) + \
        #           sum(model.VOLL[d] * (model.PD[d] - model.pD[d]) * model.baseMVA for d in model.D)
        #     return obj

        # self.model.OBJ = Objective(rule=objective, sense=minimize)

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

        # --- static generation reactive power limits ---
        def static_generation_reactive_power_bounds(model, g):
            model.qsG[g].unfix()
            return model.QsGmin[g], model.qsG[g], model.QsGmax[g]

        self.model.QsG_Constraint = Constraint(self.model.sGc, rule=static_generation_reactive_power_bounds)

        # --- reactive generator power limits ---
        def reactive_power_bounds(model, g):
            model.qG[g].unfix()
            return model.QGmin[g], model.qG[g], model.QGmax[g]

        self.model.QG_Constraint = Constraint(self.model.G, rule=reactive_power_bounds)

        # --- reactive demand limits ---
        def reactive_demand_bounds(model, d):
            model.qD[d].unfix()
            return model.QDmin[d], model.qD[d], model.QDmax[d]

        self.model.QD_Constraint = Constraint(self.model.Dc, rule=reactive_demand_bounds)

        # --- voltage constraints ---
        self.model.v_bPV_setpoint.deactivate()

        def v_bounds(model, b):
            return model.Vmin[b], model.v[b], model.Vmax[b]

        self.model.v_constraint = Constraint(self.model.B, rule=v_bounds)

        # --- wind generation q requirements variant 3---
        def QW_pos(model, w):
            return model.qsG[w] <= -0.28 * model.PsG[w] + 3.8 * model.psG[w]

        def QW_neg(model, w):
            return model.qsG[w] >= 0.03 * model.PsG[w] - 1.3 * model.psG[w]

        self.model.QW_pos_constraint = Constraint(self.model.WINDc, rule=QW_pos)
        self.model.QW_neg_constraint = Constraint(self.model.WINDc, rule=QW_neg)

        #
        x = np.array([[96, 103], [120, 127]]) / 110
        y = np.array([[0.48, 0.41, 0.33], [-0.23, -0.33, -0.41]])
        m = (y[1] - y[0]) / (x[0, 1] - x[0, 0])
        b = np.array([y[0] - m * x[i, 0] for i in range(len(x))]).T
        self.m = m[0]
        self.b_min = b[0, 0]
        self.b_max = b[0, 1]
        def QU_min(model, w):
            for (g, b) in model.sGbs:
                if g == w:
                    return model.qsG[w] >= (self.m * model.v[b] + self.b_min) * model.PsG[w]
        self.model.QU_min_constraint = Constraint(self.model.WINDc, rule=QU_min)

        def QU_max(model, w):
            for (g, b) in model.sGbs:
                if g == w:
                    return model.qsG[w] <= (self.m * model.v[b] + self.b_max) * model.PsG[w]
        self.model.QU_max_constraint = Constraint(self.model.WINDc, rule=QU_max)
