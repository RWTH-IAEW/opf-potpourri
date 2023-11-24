# ==================================================================
# DCOPF.py.mod
# PYOMO model file of "DC" optimal power flow problem (DCOPF.py)
# This formulation uses the standard "DC" linearization of the AC power flow equations
# ---Author---
# W. Bukhsh,
# wbukhsh@gmail.com
# OATS
# Copyright (c) 2017 by W Bukhsh, Glasgow, Scotland
# OATS is distributed under the GNU GENERAL PUBLIC LICENSE v3. (see LICENSE file for details).
# ==================================================================

# ==========Import==========
from __future__ import division
from pyomo.environ import *
import pandas as pd
import numpy as np


# ==========================

class HCOPF:
    def __init__(self, net):
        self.net = net

        # --- Set indices ---
        self.bus_set = net.bus.index
        self.sgen_set = net.sgen.index
        self.demand_set = net.load.index
        self.demand_neg_set = net.load[net.load.p_mw < 0].index
        self.line_set = net.line.index
        self.shunt_set = net.shunt.index
        self.trafo_set = net.trafo.index

        self.bus_gen_set = list(zip(net.sgen.bus, self.sgen_set))
        self.bus_demand_set = list(zip(net.load.bus, self.demand_set))
        self.bus_shunt_set = list(zip(net.shunt.bus, self.shunt_set))

        self.bus_line_dict = dict(
            zip(list(zip(self.line_set, [1] * len(self.line_set))) + list(zip(self.line_set, [2] * len(self.line_set))),
                pd.concat([net.line.from_bus, net.line.to_bus])))
        self.bus_trafo_dict = dict(zip(list(zip(self.trafo_set, [1] * len(self.trafo_set))) + list(
            zip(self.trafo_set, [2] * len(self.trafo_set))), pd.concat([net.trafo.hv_bus, net.trafo.lv_bus])))

        # --- Param Data ---
        self.baseMVA = net.sn_mva

        self.demand_data = net.load.p_mw / self.baseMVA

        self.PGmax_data = net.sgen.p_mw / self.baseMVA
        self.PGmin_data = 0

        max_load = net.line.max_loading_percent.values / 100. if "max_loading_percent" in net.line else 1.
        vr = net.bus.loc[net.line["from_bus"].values, "vn_kv"].values * np.sqrt(3.)
        max_i_ka = net.line.max_i_ka.values
        df = net.line.df.values
        self.SLmax_data = max_load * max_i_ka * df * net.line.parallel.values * vr / self.baseMVA

        max_load_T = net.trafo.max_loading_percent.values / 100. if "max_loading_percent" in net.trafo else 1.
        sn_mva = net.trafo.sn_mva.values
        df_T = net.trafo.df.values
        self.SLmaxT_data = max_load_T * sn_mva * df_T * net.trafo.parallel.values / self.baseMVA

        self.GB_data = net.shunt.p_mw * net.shunt.step / self.baseMVA

        self.c0_data = np.random.randint(100, size=len(net.sgen.index))
        self.c1_data = np.random.randint(100, size=len(net.sgen.index))
        self.c2_data = np.random.randint(100, size=len(net.sgen.index))

        # ---AC---
        self.QD_data = net.load.p_mw.values / self.baseMVA

        if net.sgen.q_mvar.sum() != 0:
            self.QGmax_data = net.sgen.q_mvar / self.baseMVA
        else:
            self.QGmax_data = net.sgen.p_mw / self.baseMVA
        self.QGmin_data = -self.QGmax_data

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
        bt_m = np.sqrt(yt_m ** 2 - gt_m ** 2)
        Zt_ref = net.trafo.vn_lv_kv ** 2 * net.sn_mva / net.trafo.sn_mva

        self.BiiT_data = (bt_ik + 0.5 * bt_m) * ZN[net.trafo.lv_bus].values / Zt_ref
        self.BikT_data = -bt_ik * ZN[net.trafo.lv_bus].values / Zt_ref
        self.GiiT_data = gt_ik * ZN[net.trafo.lv_bus].values / Zt_ref
        self.GikT_data = -gt_ik * ZN[net.trafo.lv_bus].values / Zt_ref

        self.BB_data = net.shunt.q_mvar * net.shunt.step / self.baseMVA

        self.Vmax_data = net.bus.max_vm_pu
        self.Vmin_data = net.bus.min_vm_pu

        self.create_model()

    def get_generator_data(self):
        if 'controllable' in self.net.sgen:
            self.sgen_controllable_set = self.net.sgen.index[self.net.sgen.controllable]
            self.sgen_static_set = self.net.sgen.index[self.net.sgen.controllable == False]
        else:
            self.sgen_controllable_set = None
            self.sgen_static_set = self.net.sgen.index

        if 'max_p_mw' in self.net.sgen:
            self.sgen_controllable_PGmax = self.net.sgen.max_p_mw[self.sgen_controllable_set]
        else:
            print("Maximum active power injection for controllable sgens not defined, using active power instead")
        if 'min_p_mw' in self.net:
            self.sgen_controllable_PGmin = self.net.sgen.min_p_mw[self.sgen_controllable_set]
        if 'max_q_mvar' in self.net.sgen:
            self.sgen_controllable_QGmax = self.net.sgen.max_q_mvar[self.sgen_controllable_set]
        if 'min_q_mvar' in self.net.sgen:
            self.sgen_controllable_QGmin = self.net.sgen.min_q_mvar[self.sgen_controllable_set]



    def create_model(self, value_of_lost_load: int = 100000):
        self.model = ConcreteModel()

        # --- SETS ---
        self.model.B = Set(initialize=self.bus_set)
        self.model.G = Set(initialize=self.sgen_set)
        self.model.WIND = Set()
        self.model.D = Set(initialize=self.demand_set)
        self.model.DNeg = Set(initialize=self.demand_neg_set)
        self.model.L = Set(initialize=self.line_set)
        self.model.SHUNT = Set(initialize=self.shunt_set)
        self.model.b0 = Set(within=self.model.B)
        self.model.LE = Set(initialize=[1, 2])
        self.model.TRANSF = Set(initialize=self.trafo_set)

        # generators, buses, loads linked to each bus b
        self.model.Gbs = Set(within=self.model.B * self.model.G,
                             initialize=self.bus_gen_set)  # set of generator-bus mapping
        self.model.Wbs = Set(within=self.model.B * self.model.WIND)  # set of wind-bus mapping
        self.model.Dbs = Set(within=self.model.B * self.model.D,
                             initialize=self.bus_demand_set)  # set of demand-bus mapping
        self.model.SHUNTbs = Set(within=self.model.B * self.model.SHUNT,
                                 initialize=self.bus_shunt_set)  # set of shunt-bus mapping

        # --- parameters ---
        # line matrix
        self.model.A = Param(self.model.L * self.model.LE, initialize=self.bus_line_dict)  # bus-line matrix
        self.model.AT = Param(self.model.TRANSF * self.model.LE,
                              initialize=self.bus_trafo_dict)  # bus-transformer matrix

        # demands
        self.model.PD = Param(self.model.D, within=Reals, initialize=self.demand_data)  # real power demand
        self.model.VOLL = Param(self.model.D, within=Reals, initialize=value_of_lost_load)  # value of lost load

        self.model.QD = Param(self.model.D, within=Reals, initialize=self.QD_data)  # reactive power demand

        # generators
        self.model.PGmax = Param(self.model.G, within=NonNegativeReals,
                                 initialize=self.PGmax_data)  # max real power of generator, p.u.
        self.model.PGmin = Param(self.model.G, within=Reals,
                                 initialize=self.PGmin_data)  # min real power of generator, p.u.
        self.model.WGmax = Param(self.model.WIND, within=NonNegativeReals)  # max real power of wind generator, p.u.
        self.model.WGmin = Param(self.model.WIND, within=NonNegativeReals)  # min real power of wind generator, p.u.

        self.model.QGmax = Param(self.model.G, within=NonNegativeReals,
                                 initialize=self.QGmax_data)  # max reactive power of generator
        self.model.QGmin = Param(self.model.G, within=Reals,
                                 initialize=self.QGmin_data)  # min reactive power of generator
        self.model.WGQmax = Param(self.model.WIND, within=NonNegativeReals)  # max reactive power of wind generator
        self.model.WGQmin = Param(self.model.WIND, within=Reals)  # min reactive power of wind generator

        # lines and transformer chracteristics and ratings
        self.model.SLmax = Param(self.model.L, within=NonNegativeReals,
                                 initialize=self.SLmax_data)  # real power line limit
        self.model.SLmaxT = Param(self.model.TRANSF, within=NonNegativeReals,
                                  initialize=self.SLmaxT_data)  # real power transformer limit

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

        # buses
        self.model.Vmax = Param(self.model.B, within=NonNegativeReals, initialize=self.Vmax_data)  # max voltage angle
        self.model.Vmin = Param(self.model.B, within=NonNegativeReals, initialize=self.Vmin_data)  # min voltage angle

        # shunt
        self.model.GB = Param(self.model.SHUNT, within=Reals, initialize=self.GB_data)  # shunt conductance
        self.model.BB = Param(self.model.SHUNT, within=Reals, initialize=self.BB_data)  # shunt susceptance

        # cost data
        self.model.c2 = Param(self.model.G, within=NonNegativeReals,
                              initialize=self.c2_data)  # generator cost coefficient c2 (*pG^2)
        self.model.c1 = Param(self.model.G, within=NonNegativeReals,
                              initialize=self.c1_data)  # generator cost coefficient c1 (*pG)
        self.model.c0 = Param(self.model.G, within=NonNegativeReals,
                              initialize=self.c0_data)  # generator cost coefficient c0

        self.model.baseMVA = Param(within=NonNegativeReals, initialize=self.baseMVA)  # base MVA

        # --- control variables ---
        def pG_bounds_rule(model, g):
            return (model.PGmin[g], model.PGmax[g])

        self.model.pG = Var(self.model.G, domain=Reals, bounds=pG_bounds_rule)  # real power generation

        def pW_bounds_rule(model, w):
            return (model.WGmin[w], model.WGmax[w])

        self.model.pW = Var(self.model.WIND, domain=Reals, bounds=pW_bounds_rule)  # real power generation from wind
        self.model.pD = Var(self.model.D, domain=Reals)  # real power demand delivered
        self.model.pLfrom = Var(self.model.L, domain=Reals)  # real power injected at b onto line
        self.model.pLto = Var(self.model.L, domain=Reals)  # real power injected at b' onto line
        self.model.pLfromT = Var(self.model.TRANSF, domain=Reals)  # real power injected at b onto transformer
        self.model.pLtoT = Var(self.model.TRANSF, domain=Reals)  # real power injected at b' onto transformer

        def qG_bounds_rule(model, g):
            return (model.QGmin[g], model.QGmax[g])

        self.model.qG = Var(self.model.G, domain=Reals, bounds=qG_bounds_rule)  # reactive power output of generator

        def qW_bounds_rule(model, g):
            return (model.WGQmin[g], model.WGQmax[g])

        self.model.qW = Var(self.model.WIND, domain=Reals, bounds=qW_bounds_rule)  # reactive power generation from wind
        self.model.qD = Var(self.model.D, domain=Reals)  # reactive power absorbed by demand
        self.model.qLfrom = Var(self.model.L, domain=Reals)  # reactive power injected at b onto line
        self.model.qLto = Var(self.model.L, domain=Reals)  # reactive power injected at b' onto line
        self.model.qLfromT = Var(self.model.TRANSF, domain=Reals)  # reactive power injected at b onto transformer
        self.model.qLtoT = Var(self.model.TRANSF, domain=Reals)  # reactive power injected at b' onto transformer

        self.model.delta = Var(self.model.B, domain=Reals, initialize=0.0)  # voltage phase angle at bus b, rad

        def v_bounds_rule(model, b):
            return (model.Vmin[b], model.Vmax[b])

        self.model.v = Var(self.model.B, domain=NonNegativeReals, initialize=1.0,
                           bounds=v_bounds_rule)  # voltage magnitude at bus b, rad
        self.model.alpha = Var(self.model.D, initialize=1.0, domain=NonNegativeReals,
                               bounds=(0, 1))  # proportion to supply of load d

        # --- cost function ---
        def objective(model):
            obj = sum(
                model.c2[g] * (model.baseMVA * model.pG[g]) ** 2 + model.c1[g] * model.baseMVA * model.pG[g] + model.c0[
                    g] for g in model.G) + \
                  sum(model.VOLL[d] * (model.PD[d] - model.pD[d]) * model.baseMVA for d in model.D)
            return obj

        self.model.OBJ = Objective(rule=objective, sense=minimize)

        # --- Kirchoff's current law at each bus b ---
        def KCL_real_def(model, b):
            return sum(model.pG[g] for g in model.G if (b, g) in model.Gbs) + \
                sum(model.pW[w] for w in model.WIND if (b, w) in model.Wbs) == \
                sum(model.pD[d] for d in model.D if (b, d) in model.Dbs) + \
                sum(model.pLfrom[l] for l in model.L if model.A[l, 1] == b) + \
                sum(model.pLto[l] for l in model.L if model.A[l, 2] == b) + \
                sum(model.pLfromT[l] for l in model.TRANSF if model.AT[l, 1] == b) + \
                sum(model.pLtoT[l] for l in model.TRANSF if model.AT[l, 2] == b) + \
                sum(model.GB[s] * model.v[b] ** 2 for s in model.SHUNT if (b, s) in model.SHUNTbs)

        def KCL_reactive_def(model, b):
            return sum(model.qG[g] for g in model.G if (b, g) in model.Gbs) + \
                sum(model.qW[w] for w in model.WIND if (b, w) in model.Wbs) == \
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
                model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (model.Bik[l] * sin(model.delta[model.A[l, 1]] - \
                                                                                      model.delta[model.A[l, 2]]) +
                                                                   model.Gik[l] * cos(
                            model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]))

        def KVL_real_toend(model, l):
            return model.pLto[l] == model.Gii[l] * (model.v[model.A[l, 2]] ** 2) + \
                model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (model.Bik[l] * sin(model.delta[model.A[l, 2]] - \
                                                                                      model.delta[model.A[l, 1]]) +
                                                                   model.Gik[l] * cos(
                            model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]))

        def KVL_reactive_fromend(model, l):
            return model.qLfrom[l] == -model.Bii[l] * (model.v[model.A[l, 1]] ** 2) + \
                model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (model.Gik[l] * sin(model.delta[model.A[l, 1]] - \
                                                                                      model.delta[model.A[l, 2]]) -
                                                                   model.Bik[l] * cos(
                            model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]))

        def KVL_reactive_toend(model, l):
            return model.qLto[l] == (-model.Bii[l] * (model.v[model.A[l, 2]] ** 2) + \
                                     model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                                             model.Gik[l] * sin(model.delta[model.A[l, 2]] - \
                                                                model.delta[model.A[l, 1]]) - model.Bik[l] * cos(
                                         model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]])))

        self.model.KVL_real_from = Constraint(self.model.L, rule=KVL_real_fromend)
        self.model.KVL_real_to = Constraint(self.model.L, rule=KVL_real_toend)
        self.model.KVL_reactive_from = Constraint(self.model.L, rule=KVL_reactive_fromend)
        self.model.KVL_reactive_to = Constraint(self.model.L, rule=KVL_reactive_toend)

        # --- Kirchoff's voltage law on each transformer line ---
        def KVL_real_fromendTransf(model, l):
            return model.pLfromT[l] == model.GiiT[l] * (model.v[model.AT[l, 1]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (model.BikT[l] * sin(model.delta[model.AT[l, 1]] - \
                                                                                         model.delta[model.AT[l, 2]]) +
                                                                     model.GikT[l] * cos(
                            model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]))

        def KVL_real_toendTransf(model, l):
            return model.pLtoT[l] == model.GiiT[l] * (model.v[model.AT[l, 2]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (model.BikT[l] * sin(model.delta[model.AT[l, 2]] - \
                                                                                         model.delta[model.AT[l, 1]]) +
                                                                     model.GikT[l] * cos(
                            model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]))

        def KVL_reactive_fromendTransf(model, l):
            return model.qLfromT[l] == -model.BiiT[l] * (model.v[model.AT[l, 1]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (model.GikT[l] * sin(model.delta[model.AT[l, 1]] - \
                                                                                         model.delta[model.AT[l, 2]]) -
                                                                     model.BikT[l] * cos(
                            model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]))

        def KVL_reactive_toendTransf(model, l):
            return model.qLtoT[l] == -model.BiiT[l] * (model.v[model.AT[l, 2]] ** 2) + \
                model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (model.GikT[l] * sin(model.delta[model.AT[l, 2]] - \
                                                                                         model.delta[model.AT[l, 1]]) -
                                                                     model.BikT[l] * cos(
                            model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]))

        self.model.KVL_real_fromTransf = Constraint(self.model.TRANSF, rule=KVL_real_fromendTransf)
        self.model.KVL_real_toTransf = Constraint(self.model.TRANSF, rule=KVL_real_toendTransf)
        self.model.KVL_reactive_fromTransf = Constraint(self.model.TRANSF, rule=KVL_reactive_fromendTransf)
        self.model.KVL_reactive_toTransf = Constraint(self.model.TRANSF, rule=KVL_reactive_toendTransf)

        # --- demand and load shedding ---
        def Load_Shed_real(model, d):
            return model.pD[d] == model.alpha[d] * model.PD[d]

        def Load_Shed_reactive(model, d):
            return model.qD[d] == model.alpha[d] * model.QD[d]

        def alpha_FixNegDemands(model, d):
            return model.alpha[d] == 1

        self.model.LoadShed_real = Constraint(self.model.D, rule=Load_Shed_real)
        self.model.LoadShed_reactive = Constraint(self.model.D, rule=Load_Shed_reactive)
        self.model.alphaFix = Constraint(self.model.DNeg, rule=alpha_FixNegDemands)

        # --- line power limits ---
        def line_lim1_def(model, l):
            return model.pLfrom[l] ** 2 + model.qLfrom[l] ** 2 <= model.SLmax[l] ** 2

        def line_lim2_def(model, l):
            return model.pLto[l] ** 2 + model.qLto[l] ** 2 <= model.SLmax[l] ** 2

        self.model.line_lim1 = Constraint(self.model.L, rule=line_lim1_def)
        self.model.line_lim2 = Constraint(self.model.L, rule=line_lim2_def)

        # --- power flow limits on transformer lines---
        def transf_lim1_def(model, l):
            return model.pLfromT[l] ** 2 + model.qLfromT[l] ** 2 <= model.SLmaxT[l] ** 2

        def transf_lim2_def(model, l):
            return model.pLtoT[l] ** 2 + model.qLtoT[l] ** 2 <= model.SLmaxT[l] ** 2

        self.model.transf_lim1 = Constraint(self.model.TRANSF, rule=transf_lim1_def)
        self.model.transf_lim2 = Constraint(self.model.TRANSF, rule=transf_lim2_def)

        # --- reference bus constraint ---
        def ref_bus_def(model, b):
            return model.delta[b] == 0

        self.model.refbus = Constraint(self.model.b0, rule=ref_bus_def)

    def solve(self):
        solver = SolverFactory('ipopt')
        self.results = solver.solve(self.model, load_solutions=True)
