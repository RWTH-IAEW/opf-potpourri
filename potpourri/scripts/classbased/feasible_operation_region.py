#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Calculation of the feasible operation region

example of usage:


Institut für Elektrische Anlagen und Netze, Digitalisierung und Energiewirtschaft (IAEW)
(c) 2023, Steffen Kortmann
"""

# ==========Import==========
from __future__ import division

import copy

from pyomo.environ import *
import pandapower as pp
import pickle

from potpourri.models.class_based.HC_ACOPF import HC_ACOPF


# ==========================

def feasible_operation_region(net: pp.pandapowerNet) -> ConcreteModel:
    model = ConcreteModel()

    def _create_sets(model: ConcreteModel, net: pp.pandapowerNet):
        # --- sets ---
        # buses, generators, loads, lines, sections
        model.B = Set(initialize=net.bus.index)  # set of buses
        model.G = Set(initialize=net.gen.index)  # set of generators
        model.WIND = Set(initialize=net.sgen.index)  # set of wind generators
        model.D = Set(initialize=net.load.index)  # set of demands
        model.DNeg = Set(initialize=net.load.index)  # set of demands
        model.L = Set(initialize=net.line.index)  # set of lines
        model.SHUNT = Set(initialize=net.shunt.index)  # set of shunts
        model.LE = Set()  # line-to and from ends set (1,2)
        model.TRANSF = Set(initialize=net.trafo.index)  # set of transformers
        model.b0 = Set(initialize=net.ext_grid.bus, within=model.B)  # set of reference buses

        # generators, buses, loads linked to each bus b
        model.Gbs = Set(within=model.B * model.G)  # generator-bus mapping
        model.Dbs = Set(within=model.B * model.D)  # demand-bus mapping
        model.Wbs = Set(within=model.B * model.WIND)  # wind-bus mapping
        model.SHUNTbs = Set(within=model.B * model.SHUNT)  # shunt-bus mapping

        return model

    model = _create_sets(model, net)

    def _create_parameters(model: ConcreteModel, net: pp.pandapowerNet):
        # --- parameters ---
        # line matrix
        model.A = Param(model.L * model.LE, within=pyomo.core.Any)  # bus-line
        model.AT = Param(model.TRANSF * model.LE, within=pyomo.core.Any)  # bus-transformer

        # demands
        model.PD = Param(model.D,
                         initialize={load_idx: net.load.p_mw[load_idx] for load_idx in net.load.index},
                         within=pyomo.core.Reals)  # real power demand
        model.QD = Param(model.D,
                         initialize={load_idx: net.load.q_mvar[load_idx] for load_idx in net.load.index},
                         within=pyomo.core.Reals)  # reactive power demand
        model.VOLL = Param(model.D,
                           initialize={load_idx: 1000 for load_idx in net.load.index},
                           within=pyomo.core.Reals)  # value of lost load

        # generators
        model.PGmax = Param(model.G,
                            initialize={gen_idx: net.gen.max_p_mw[gen_idx] for gen_idx in net.gen.index},
                            within=pyomo.core.NonNegativeReals)  # max real power of generator
        model.PGmin = Param(model.G,
                            initialize={gen_idx: 0 for gen_idx in net.gen.index},
                            within=pyomo.core.Reals)  # min real power of generator
        model.QGmax = Param(model.G,
                            initialize={gen_idx: net.gen.max_q_mvar[gen_idx] for gen_idx in net.gen.index},
                            within=pyomo.core.NonNegativeReals)  # max reactive power of generator
        model.QGmin = Param(model.G,
                            initialize={gen_idx: net.gen.min_q_mvar[gen_idx] for gen_idx in net.gen.index},
                            within=pyomo.core.Reals)  # min reactive power of generator

        # wind generators
        model.WGmax = Param(model.WIND,
                            initialize={sgen_idx: net.sgen.max_p_mw[sgen_idx] for sgen_idx in net.sgen.index},
                            within=pyomo.core.NonNegativeReals)  # max real power of wind generator
        model.WGmin = Param(model.WIND,
                            initialize={sgen_idx: 0 for sgen_idx in net.sgen.index},
                            within=pyomo.core.NonNegativeReals)  # min real power of wind generator
        model.WGQmax = Param(model.WIND,
                             initialize={sgen_idx: net.sgen.max_q_mvar[sgen_idx] for sgen_idx in net.sgen.index},
                             within=pyomo.core.NonNegativeReals)  # max reactive power of wind generator
        model.WGQmin = Param(model.WIND,
                             initialize={sgen_idx: net.sgen.min_q_mvar[sgen_idx] for sgen_idx in net.sgen.index},
                             within=pyomo.core.Reals)  # min reactive power of wind generator

        # lines
        model.SLmax = Param(model.L,
                            initialize={line_idx: net.line.max_i_ka[line_idx] for line_idx in net.line.index},
                            within=pyomo.core.NonNegativeReals)  # max real power limit on flow in a line
        model.GL = Param(model.L, within=pyomo.core.Reals)
        model.BL = Param(model.L, within=pyomo.core.Reals)
        model.BC = Param(model.L, within=pyomo.core.Reals)

        # emergency ratings
        model.SLmax_E = Param(model.L, within=pyomo.core.NonNegativeReals)  # max emergency real power flow limit
        model.SLmaxT_E = Param(model.TRANSF, within=pyomo.core.NonNegativeReals)  # max emergency real power flow limit

        # transformers
        model.Tap = Param(model.TRANSF, within=pyomo.core.NonNegativeReals)  # turns ratio of a transformer
        model.TapLB = Param(model.TRANSF, within=pyomo.core.NonNegativeReals)  # lower bound on turns ratio
        model.TapUB = Param(model.TRANSF, within=pyomo.core.NonNegativeReals)  # upper bound on turns ratio
        model.Deltashift = Param(model.TRANSF)  # phase shift of transformer, rad
        model.DeltashiftLB = Param(model.TRANSF)  # lower bound on phase shift of transformer, rad
        model.DeltashiftUB = Param(model.TRANSF)  # upper bound on phase shift of transformer, rad

        model.SLmaxT = Param(model.TRANSF, within=pyomo.core.NonNegativeReals)  # max real power flow limit
        model.GLT = Param(model.TRANSF, within=pyomo.core.Reals)
        model.BLT = Param(model.TRANSF, within=pyomo.core.Reals)

        # derived line parameters
        model.G11 = Param(model.L, within=pyomo.core.Reals)
        model.G12 = Param(model.L, within=pyomo.core.Reals)
        model.G21 = Param(model.L, within=pyomo.core.Reals)
        model.G22 = Param(model.L, within=pyomo.core.Reals)
        model.B11 = Param(model.L, within=pyomo.core.Reals)
        model.B12 = Param(model.L, within=pyomo.core.Reals)
        model.B21 = Param(model.L, within=pyomo.core.Reals)
        model.B22 = Param(model.L, within=pyomo.core.Reals)
        ## derived transformer parameters
        model.G11T = Param(model.TRANSF, within=pyomo.core.Reals)
        model.G12T = Param(model.TRANSF, within=pyomo.core.Reals)
        model.G21T = Param(model.TRANSF, within=pyomo.core.Reals)
        model.G22T = Param(model.TRANSF, within=pyomo.core.Reals)
        model.B11T = Param(model.TRANSF, within=pyomo.core.Reals)
        model.B12T = Param(model.TRANSF, within=pyomo.core.Reals)
        model.B21T = Param(model.TRANSF, within=pyomo.core.Reals)
        model.B22T = Param(model.TRANSF, within=pyomo.core.Reals)

        # buses
        model.Vmax = Param(model.B, within=pyomo.core.NonNegativeReals)  # max voltage angle
        model.Vmin = Param(model.B, within=pyomo.core.NonNegativeReals)  # min voltage angle

        # shunt
        model.GB = Param(model.SHUNT, within=pyomo.core.Reals)  # shunt conductance
        model.BB = Param(model.SHUNT, within=pyomo.core.Reals)  # shunt susceptance

        # cost
        model.c2 = Param(model.G, within=pyomo.core.NonNegativeReals)  # generator cost coefficient c2 (*pG^2)
        model.c1 = Param(model.G, within=pyomo.core.NonNegativeReals)  # generator cost coefficient c1 (*pG)
        model.c0 = Param(model.G, within=pyomo.core.NonNegativeReals)  # generator cost coefficient c0

        model.baseMVA = Param(within=pyomo.core.NonNegativeReals)  # base MVA

        # constants
        model.eps = Param(within=pyomo.core.NonNegativeReals)

        return model

    model = _create_parameters(model, net)

    def _create_variables(model: ConcreteModel, net: pp.pandapowerNet):
        # --- variables ---
        model.pG = Var(model.G, domain=pyomo.core.NonNegativeReals)  # real power output of generator
        model.qG = Var(model.G, domain=pyomo.core.Reals)  # reactive power output of generator
        model.pW = Var(model.WIND, domain=pyomo.core.Reals)  # real power generation from wind
        model.qW = Var(model.WIND, domain=pyomo.core.Reals)  # reactive power generation from wind
        model.pD = Var(model.D, domain=pyomo.core.Reals)  # real power absorbed by demand
        model.qD = Var(model.D, domain=pyomo.core.Reals)  # reactive power absorbed by demand
        model.pLfrom = Var(model.L, domain=pyomo.core.Reals)  # real power injected at b onto line
        model.pLto = Var(model.L, domain=pyomo.core.Reals)  # real power injected at b' onto line
        model.qLfrom = Var(model.L, domain=pyomo.core.Reals)  # reactive power injected at b onto line
        model.qLto = Var(model.L, domain=pyomo.core.Reals)  # reactive power injected at b' onto line
        model.pLfromT = Var(model.TRANSF, domain=pyomo.core.Reals)  # real power injected at b onto transformer
        model.pLtoT = Var(model.TRANSF, domain=pyomo.core.Reals)  # real power injected at b' onto transformer
        model.qLfromT = Var(model.TRANSF, domain=pyomo.core.Reals)  # reactive power injected at b onto transformer
        model.qLtoT = Var(model.TRANSF, domain=pyomo.core.Reals)  # reactive power injected at b' onto transformer

        # model.deltaL = Var(model.L, domain= Reals) # angle difference across lines
        model.delta = Var(model.B, domain=pyomo.core.Reals, initialize=0.0)  # voltage phase angle at bus b, rad
        model.v = Var(model.B, domain=pyomo.core.NonNegativeReals, initialize=1.0)  # voltage magnitude at bus b, rad
        model.alpha = Var(model.D, initialize=1.0, domain=pyomo.core.NonNegativeReals)  # proportion to supply of load d

    # --- cost function ---
    def objective(model):
        obj = sum(
            model.c2[g] * (model.baseMVA * model.pG[g]) ** 2 + model.c1[g] * model.baseMVA * model.pG[g] + model.c0[g]
            for g in model.G) + \
              sum(model.VOLL[d] * (1 - model.alpha[d]) * model.baseMVA * model.PD[d] for d in model.D)
        return obj

    model.OBJ = Objective(rule=objective, sense=minimize)

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

    model.KCL_real = Constraint(model.B, rule=KCL_real_def)
    model.KCL_reactive = Constraint(model.B, rule=KCL_reactive_def)

    # --- Kirchoff's voltage law on each line ---
    def KVL_real_fromend(model, l):
        return model.pLfrom[l] == model.G11[l] * (model.v[model.A[l, 1]] ** 2) + \
            model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (model.B12[l] * sin(model.delta[model.A[l, 1]] - \
                                                                                  model.delta[model.A[l, 2]]) +
                                                               model.G12[l] * cos(
                        model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]))

    def KVL_real_toend(model, l):
        return model.pLto[l] == model.G22[l] * (model.v[model.A[l, 2]] ** 2) + \
            model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (model.B21[l] * sin(model.delta[model.A[l, 2]] - \
                                                                                  model.delta[model.A[l, 1]]) +
                                                               model.G21[l] * cos(
                        model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]]))

    def KVL_reactive_fromend(model, l):
        return model.qLfrom[l] == -model.B11[l] * (model.v[model.A[l, 1]] ** 2) + \
            model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (model.G12[l] * sin(model.delta[model.A[l, 1]] - \
                                                                                  model.delta[model.A[l, 2]]) -
                                                               model.B12[l] * cos(
                        model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]))

    def KVL_reactive_toend(model, l):
        return model.qLto[l] == (-model.B22[l] * (model.v[model.A[l, 2]] ** 2) + \
                                 model.v[model.A[l, 1]] * model.v[model.A[l, 2]] * (
                                         model.G21[l] * sin(model.delta[model.A[l, 2]] - \
                                                            model.delta[model.A[l, 1]]) - model.B21[l] * cos(
                                     model.delta[model.A[l, 2]] - model.delta[model.A[l, 1]])))

    model.KVL_real_from = Constraint(model.L, rule=KVL_real_fromend)
    model.KVL_real_to = Constraint(model.L, rule=KVL_real_toend)
    model.KVL_reactive_from = Constraint(model.L, rule=KVL_reactive_fromend)
    model.KVL_reactive_to = Constraint(model.L, rule=KVL_reactive_toend)

    # --- Kirchoff's voltage law on each transformer line ---
    def KVL_real_fromendTransf(model, l):
        return model.pLfromT[l] == model.G11T[l] * (model.v[model.AT[l, 1]] ** 2) + \
            model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (model.B12T[l] * sin(model.delta[model.AT[l, 1]] - \
                                                                                     model.delta[model.AT[l, 2]]) +
                                                                 model.G12T[l] * cos(
                        model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]))

    def KVL_real_toendTransf(model, l):
        return model.pLtoT[l] == model.G22T[l] * (model.v[model.AT[l, 2]] ** 2) + \
            model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (model.B21T[l] * sin(model.delta[model.AT[l, 2]] - \
                                                                                     model.delta[model.AT[l, 1]]) +
                                                                 model.G21T[l] * cos(
                        model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]))

    def KVL_reactive_fromendTransf(model, l):
        return model.qLfromT[l] == -model.B11T[l] * (model.v[model.AT[l, 1]] ** 2) + \
            model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (model.G12T[l] * sin(model.delta[model.AT[l, 1]] - \
                                                                                     model.delta[model.AT[l, 2]]) -
                                                                 model.B12T[l] * cos(
                        model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]))

    def KVL_reactive_toendTransf(model, l):
        return model.qLtoT[l] == -model.B22T[l] * (model.v[model.AT[l, 2]] ** 2) + \
            model.v[model.AT[l, 1]] * model.v[model.AT[l, 2]] * (model.G21T[l] * sin(model.delta[model.AT[l, 2]] - \
                                                                                     model.delta[model.AT[l, 1]]) -
                                                                 model.B21T[l] * cos(
                        model.delta[model.AT[l, 2]] - model.delta[model.AT[l, 1]]))

    model.KVL_real_fromTransf = Constraint(model.TRANSF, rule=KVL_real_fromendTransf)
    model.KVL_real_toTransf = Constraint(model.TRANSF, rule=KVL_real_toendTransf)
    model.KVL_reactive_fromTransf = Constraint(model.TRANSF, rule=KVL_reactive_fromendTransf)
    model.KVL_reactive_toTransf = Constraint(model.TRANSF, rule=KVL_reactive_toendTransf)

    # --- generator power limits ---
    def Real_Power_Max(model, g):
        return model.pG[g] <= model.PGmax[g]

    def Real_Power_Min(model, g):
        return model.pG[g] >= model.PGmin[g]

    def Reactive_Power_Max(model, g):
        return model.qG[g] <= model.QGmax[g]

    def Reactive_Power_Min(model, g):
        return model.qG[g] >= model.QGmin[g]

    model.PGmaxC = Constraint(model.G, rule=Real_Power_Max)
    model.PGminC = Constraint(model.G, rule=Real_Power_Min)
    model.QGmaxC = Constraint(model.G, rule=Reactive_Power_Max)
    model.QGminC = Constraint(model.G, rule=Reactive_Power_Min)

    # ---wind generator power limits ---
    def Wind_Real_Power_Max(model, w):
        return model.pW[w] <= model.WGmax[w]

    def Wind_Real_Power_Min(model, w):
        return model.pW[w] >= model.WGmin[w]

    def Wind_Reactive_Power_Max(model, w):
        return model.qW[w] <= model.WGQmax[w]

    def Wind_Reactive_Power_Min(model, w):
        return model.qW[w] >= model.WGQmin[w]

    model.WGmaxC = Constraint(model.WIND, rule=Wind_Real_Power_Max)
    model.WGminC = Constraint(model.WIND, rule=Wind_Real_Power_Min)
    model.WGQmaxC = Constraint(model.WIND, rule=Wind_Reactive_Power_Max)
    model.WGQminC = Constraint(model.WIND, rule=Wind_Reactive_Power_Min)

    # --- demand and load shedding ---
    def Load_Shed_real(model, d):
        return model.pD[d] == model.alpha[d] * model.PD[d]

    def Load_Shed_reactive(model, d):
        return model.qD[d] == model.alpha[d] * model.QD[d]

    def alpha_FixNegDemands(model, d):
        return model.alpha[d] == 1

    model.LoadShed_real = Constraint(model.D, rule=Load_Shed_real)
    model.LoadShed_reactive = Constraint(model.D, rule=Load_Shed_reactive)
    model.alphaFix = Constraint(model.DNeg, rule=alpha_FixNegDemands)

    def alpha_BoundUB(model, d):
        return model.alpha[d] <= 1

    def alpha_BoundLB(model, d):
        return model.alpha[d] >= 0

    model.alphaBoundUBC = Constraint(model.D, rule=alpha_BoundUB)
    model.alphaBoundLBC = Constraint(model.D, rule=alpha_BoundLB)

    # --- line power limits ---
    def line_lim1_def(model, l):
        return model.pLfrom[l] ** 2 + model.qLfrom[l] ** 2 <= model.SLmax[l] ** 2

    def line_lim2_def(model, l):
        return model.pLto[l] ** 2 + model.qLto[l] ** 2 <= model.SLmax[l] ** 2

    model.line_lim1 = Constraint(model.L, rule=line_lim1_def)
    model.line_lim2 = Constraint(model.L, rule=line_lim2_def)

    # --- power flow limits on transformer lines---
    def transf_lim1_def(model, l):
        return model.pLfromT[l] ** 2 + model.qLfromT[l] ** 2 <= model.SLmaxT[l] ** 2

    def transf_lim2_def(model, l):
        return model.pLtoT[l] ** 2 + model.qLtoT[l] ** 2 <= model.SLmaxT[l] ** 2

    model.transf_lim1 = Constraint(model.TRANSF, rule=transf_lim1_def)
    model.transf_lim2 = Constraint(model.TRANSF, rule=transf_lim2_def)

    # --- voltage constraints ---
    def bus_max_voltage(model, b):
        return model.v[b] <= model.Vmax[b]

    def bus_min_voltage(model, b):
        return model.v[b] >= model.Vmin[b]

    model.Vmaxc = Constraint(model.B, rule=bus_max_voltage)
    model.Vminc = Constraint(model.B, rule=bus_min_voltage)

    # --- reference bus constraint ---
    def ref_bus_def(model, b):
        return model.delta[b] == 0

    model.refbus = Constraint(model.b0, rule=ref_bus_def)

    return model


def run_feasible_operation_region(self):
    print("Run feasible operation region")
    import numpy as np

    theta_values = np.linspace(0, 2 * np.pi, 8)  # One degree resolution
    boundary_P_values = [[], [], []]
    boundary_Q_values = [[], [], []]
    boundary_U_values = []

    # Set objective to maximize absolute power for the given angle
    self.model.obj = Objective(
        expr=sum(self.model.peG[b0]**2 + self.model.qeG[b0]**2 for b0 in self.model.eG),
        sense=maximize)

    # Constraint for the given power factor angle
    self.model.tan_theta = Param(mutable=True, initialize=0.0)

    def constrain_pf(model, b0):
        return model.qeG[b0] == self.model.tan_theta * model.peG[b0]

    self.model.power_factor_constraint = Constraint(self.model.eG, rule=constrain_pf)

    for theta in theta_values:
        # Constraint for the given power factor angle
        tan_theta = np.tan(theta)
        self.model.tan_theta = tan_theta

        # Solve
        self.solve()

        for i in self.model.eG:
            print(value(self.model.peG[i]))

        for b0 in self.model.eG:
            boundary_P_values[b0].append(value(self.model.peG[b0]))
            boundary_Q_values[b0].append(value(self.model.qeG[b0]))

    return boundary_P_values, boundary_Q_values

def for_setpoint_based(opf):
    alpha_beta = [(1, 0), (1, 1), (0, 1), (-1, 0), (-1, -1), (0, -1), (1, -1), (-1, 1)]
    opf.model.alpha_beta = Set(initialize=alpha_beta)
    opf.model.alpha = Param(initialize=0, mutable=True)
    opf.model.beta = Param(initialize=0, mutable=True)

    def setpoint_based(model):
        return sum(model.peG[g]*model.alpha + model.qeG[g]*model.beta for g in model.eG)
    opf.model.obj = Objective(expr=setpoint_based, sense=maximize)

    boundary_P_values = [[], [], []]
    boundary_Q_values = [[], [], []]

    for (alpha, beta) in alpha_beta:
        opf.model.alpha = alpha
        opf.model.beta = beta
        opf.solve()
        print(value(opf.model.obj))

        for i in opf.model.eG:
            print(value(opf.model.peG[i]))

        for b0 in opf.model.eG:
            boundary_P_values[b0].append(value(opf.model.peG[b0]))
            boundary_Q_values[b0].append(value(opf.model.qeG[b0]))

    return boundary_P_values, boundary_Q_values


if __name__ == "__main__":
    net = pp.networks.create_cigre_network_mv()

    with open('C:/Users/f.lohse/PycharmProjects/potpourri/potpourri/data/simbench_hv_grid_with_potential_pkl.pkl',
              'rb') as f:
        net = pickle.load(f)

    # net = sb.get_simbench_net("1-HV-mixed--0-no_sw")

    case = 'lW'

    factors = net.loadcases.loc[case]
    net.load.p_mw *= factors['pload']
    net.load.q_mvar *= factors['qload']
    net.sgen.scaling[net.sgen.type == 'Wind'] = factors['Wind_p']
    net.sgen.scaling[net.sgen.type == 'PV'] = factors['PV_p']
    net.sgen.scaling[(net.sgen.type != 'Wind') & (net.sgen.type != 'Solar')] = factors['RES_p']
    net.ext_grid.vm_pu = factors['Slack_vm']

    hc = HC_ACOPF(net, SWmin=10)
    hc.solve()
    hc.add_OPF()
    hc.solve(solver='mindtpy', mip_solver='gurobi')

    net_wind = copy.deepcopy(hc.net)
    net_wind.sgen['controllable'] = True
    net_wind.load['controllable'] = True

    hc_for = HC_ACOPF(net_wind)
    hc_for.solve()
    hc_for.add_OPF()
    hc_for.fix_vars('y')
    hc_for.model.OBJ.deactivate()

    # p, q = for_setpoint_based(hc_for)
    p, q = run_feasible_operation_region(hc_for)
    # model = feasible_operation_region(net)
    # fig, ax = plt.subplots()
    # for g in hc_for.model.eG:
    #     ax.plot(p[g], q[g], '.-', label = 'Slack #' + str(g))
    # ax.legend()
    # ax.set_xlabel('P [MW]')
    # ax.set_ylabel('Q [MVAr]')