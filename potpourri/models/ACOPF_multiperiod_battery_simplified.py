#==================================================================
# ACOPF.mod
# PYOMO model file of "AC" optimal power flow problem (ACOPF)
# ---Author---
# W. Bukhsh,
# wbukhsh@gmail.com
# OATS
# Copyright (c) 2015 by W Bukhsh, Glasgow, Scotland
# OATS is distributed under the GNU GENERAL PUBLIC LICENSE license. (see LICENSE file for details).
#==========Import==========
from __future__ import division
from pyomo.environ import *
#==========================

model = AbstractModel()

# --- sets ---
# buses, generators, loads, lines, sections
model.B      = Set()  # set of buses
model.G      = Set()  # set of generators
model.WIND   = Set()  # set of wind generators
model.D      = Set()  # set of demands
model.DNeg   = Set()  # set of demands
model.L      = Set()  # set of lines
model.SHUNT  = Set()  # set of shunts
model.LE     = Set()  # line-to and from ends set (1,2)
model.TRANSF = Set()  # set of transformers
model.b0     = Set(within=model.B)  # set of reference buses
model.T      = Set() # set of time
model.Tred   = Set()
model.BATTERY = Set() # set of batteries
model.Tend    = Param(within=NonNegativeReals) #number of periods
model.S       = Set()  # Set of solar power plants 
#Tend = 0
#model.Tend = (Tend +1 for i in model.T)
# generators, buses, loads linked to each bus b
model.Gbs = Set(within=model.B * model.G)    # generator-bus mapping
model.Dbs = Set(within=model.B * model.D)    # demand-bus mapping
model.Wbs = Set(within=model.B * model.WIND) # wind-bus mapping
model.SHUNTbs = Set(within=model.B * model.SHUNT)# shunt-bus mapping
model.Bbs = Set(within=model.B * model.BATTERY) # battery-bus mapping
model.Sbs = Set(within=model.B * model.S) #solar-bus mapping

# --- parameters ---
# line matrix
model.A = Param(model.L*model.LE,within=Any)       # bus-line
model.AT = Param(model.TRANSF*model.LE,within=Any) # bus-transformer

# demands
model.PD = Param(model.D, model.T, within=Reals)  # real power demand
model.QD = Param(model.D, within=Reals)  # reactive power demand
model.VOLL    = Param(model.D, within=Reals) #value of lost load

# generators
model.PGmax    = Param(model.G, within=NonNegativeReals) # max real power of generator
model.PGmin    = Param(model.G, within=Reals)            # min real power of generator
model.QGmax    = Param(model.G, within=NonNegativeReals) # max reactive power of generator
model.QGmin    = Param(model.G, within=Reals)            # min reactive power of generator

#wind generators
model.WGmax    = Param(model.WIND,  within=NonNegativeReals) # max real power of wind generator
model.WGmin    = Param(model.WIND,  within=NonNegativeReals) # min real power of wind generator
model.WGQmax   = Param(model.WIND,  within=NonNegativeReals) # max reactive power of wind generator
model.WGQmin   = Param(model.WIND,  within=Reals)            # min reactive power of wind generator

# batteries
model.PCmax = Param(model.BATTERY, within=NonNegativeReals) # max real power for charging of batteries
model.PCmin = Param(model.BATTERY, within=NonNegativeReals) # min real power for charging of batteries
model.PDmax = Param(model.BATTERY, within=NonNegativeReals) # max real power for discharging of batteries
model.PDmin = Param(model.BATTERY, within=NonNegativeReals) # min real power for discharging of batteries
model.EMAX = Param(model.BATTERY, within=NonNegativeReals) # max SoC level
model.EMIN = Param(model.BATTERY, within=NonNegativeReals) # min SoC level
model.E0 = Param(model.BATTERY, within=NonNegativeReals) # backlog energy
model.nchar = Param(model.BATTERY, within=NonNegativeReals) # charging efficiency 
model.ndis = Param(model.BATTERY, within=PositiveReals) #discharging efficiency
#model.z = Var(model.BATTERY, model.T, within=Binary)
#model.y = Var(model.BATTERY, model.T, within=Binary)

#solar power plant
model.PS = Param(model.S, model.T, within=NonNegativeReals)
model.QS = Param(model.S, model.T, within=Reals)

# lines
model.SLmax = Param(model.L, within=NonNegativeReals) # max real power limit on flow in a line
model.GL = Param(model.L, within=Reals)
model.BL = Param(model.L, within=Reals)
model.BC = Param(model.L, within=Reals)

#emergency ratings
model.SLmax_E = Param(model.L, within=NonNegativeReals)       # max emergency real power flow limit
model.SLmaxT_E = Param(model.TRANSF, within=NonNegativeReals) # max emergency real power flow limit

#transformers
model.Tap          = Param(model.TRANSF, within=NonNegativeReals)  # turns ratio of a transformer
model.TapLB        = Param(model.TRANSF, within=NonNegativeReals)  # lower bound on turns ratio
model.TapUB        = Param(model.TRANSF, within=NonNegativeReals)  # upper bound on turns ratio
model.Deltashift   = Param(model.TRANSF) #  phase shift of transformer, rad
model.DeltashiftLB = Param(model.TRANSF) #  lower bound on phase shift of transformer, rad
model.DeltashiftUB = Param(model.TRANSF) #  upper bound on phase shift of transformer, rad

model.SLmaxT = Param(model.TRANSF, within=NonNegativeReals) # max real power flow limit
model.GLT    = Param(model.TRANSF, within=Reals)
model.BLT    = Param(model.TRANSF, within=Reals)

# derived line parameters
model.G11 = Param(model.L, within=Reals)
model.G12 = Param(model.L, within=Reals)
model.G21 = Param(model.L, within=Reals)
model.G22 = Param(model.L, within=Reals)
model.B11 = Param(model.L, within=Reals)
model.B12 = Param(model.L, within=Reals)
model.B21 = Param(model.L, within=Reals)
model.B22 = Param(model.L, within=Reals)
## derived transformer parameters
model.G11T = Param(model.TRANSF, within=Reals)
model.G12T = Param(model.TRANSF, within=Reals)
model.G21T = Param(model.TRANSF, within=Reals)
model.G22T = Param(model.TRANSF, within=Reals)
model.B11T = Param(model.TRANSF, within=Reals)
model.B12T = Param(model.TRANSF, within=Reals)
model.B21T = Param(model.TRANSF, within=Reals)
model.B22T = Param(model.TRANSF, within=Reals)

# buses
model.Vmax = Param(model.B, within=NonNegativeReals) #  max voltage angle
model.Vmin = Param(model.B, within=NonNegativeReals) #  min voltage angle

#shunt
model.GB = Param(model.SHUNT, within=Reals) #  shunt conductance
model.BB = Param(model.SHUNT, within=Reals) #  shunt susceptance

# cost
model.c2 = Param(model.G, within=NonNegativeReals)# generator cost coefficient c2 (*pG^2)
model.c1 = Param(model.G, within=NonNegativeReals)# generator cost coefficient c1 (*pG)
model.c0 = Param(model.G, within=NonNegativeReals)# generator cost coefficient c0
model.c3 = Param(model.BATTERY, within=NonNegativeReals) #battery cost coefficient 
model.c4 = Param(model.S, within=NonNegativeReals)

model.baseMVA = Param(within=NonNegativeReals)# base MVA

#constants
model.eps = Param(within=NonNegativeReals)

# --- variables ---
# TODO: Bengisu, make time-variant
model.pG       = Var(model.G, model.T, domain= NonNegativeReals)# real power output of generator
model.qG       = Var(model.G, model.T, domain= Reals)# reactive power output of generator
model.pW       = Var(model.WIND, model.T, domain= Reals) #real power generation from wind
model.qW       = Var(model.WIND, model.T, domain= Reals) #reactive power generation from wind
model.pD       = Var(model.D, model.T, domain= Reals)# real power absorbed by demand
model.qD       = Var(model.D, model.T, domain= Reals)# reactive power absorbed by demand
model.pLfrom   = Var(model.L, model.T, domain= Reals) # real power injected at b onto line
model.pLto     = Var(model.L, model.T, domain= Reals) # real power injected at b' onto line
model.qLfrom   = Var(model.L, model.T, domain= Reals) # reactive power injected at b onto line
model.qLto     = Var(model.L, model.T, domain= Reals) # reactive power injected at b' onto line
model.pLfromT  = Var(model.TRANSF, model.T, domain= Reals) # real power injected at b onto transformer
model.pLtoT    = Var(model.TRANSF, model.T, domain= Reals) # real power injected at b' onto transformer
model.qLfromT  = Var(model.TRANSF, model.T, domain= Reals) # reactive power injected at b onto transformer
model.qLtoT    = Var(model.TRANSF, model.T, domain= Reals) # reactive power injected at b' onto transformer
model.pChar    = Var(model.BATTERY, model.T, domain = NonNegativeReals) # real charging power of battery
model.pDis     = Var(model.BATTERY, model.T, domain=NonNegativeReals) # real discharging power of battery
model.qChar    = Var(model.BATTERY, model.T, domain= NonNegativeReals) # reactive charging power of battery
model.qDis     = Var(model.BATTERY, model.T, domain=NonNegativeReals) # reactive discharging power of battery
model.e        = Var(model.BATTERY, model.T, domain=NonNegativeReals) # state of charge in battery

#model.deltaL = Var(model.L, domain= Reals) # angle difference across lines
# TODO: Bengisu, make time-variant
model.delta  = Var(model.B, model.T, domain= Reals, initialize=0.0) # voltage phase angle at bus b, rad
model.v      = Var(model.B, model.T, domain= NonNegativeReals, initialize=1.0) # voltage magnitude at bus b, rad
# TODO: Bengisu, make not time-variant
model.alpha  = Var(model.D, initialize=1.0, domain= NonNegativeReals)# proportion to supply of load d

# --- cost function ---
'''
def objective(model):
    # obj = sum((model.baseMVA*model.pG[g,t])**2+model.baseMVA*model.pG[g,t] for g in model.G for t in model.T)
    obj = sum(model.baseMVA*model.pG[g,t] for g in model.G for t in model.T)
    return obj
model.OBJ = Objective(rule=objective, sense=minimize)
'''
def objective(model):
    #TODO: sum(model.pChar[a,t] for a in model.BATTERY for t in model.T) +\

    obj = sum(model.c2[g]*(model.baseMVA*model.pG[g,t])**2+model.c1[g]*model.baseMVA*model.pG[g,t]+ model.c0[g] for g in model.G for t in model.T)+\
    sum(model.VOLL[d]*(1-model.alpha[d])*model.baseMVA*model.PD[d,t] for d in model.D for t in model.T)+\
    sum(model.c3[a]*model.pDis[a,t]*model.baseMVA for a in model.BATTERY for t in model.T) +\
    sum(model.c4[s]*model.baseMVA*model.PS[s,t] for s in model.S for t in model.T)
    return obj
model.OBJ = Objective(rule=objective, sense=minimize)

# --- Kirchoff's current law at each bus b ---
# TODO: Bengisu, make time-variant


def KCL_real_def(model, b, t):
    return sum(model.pG[g, t] for g in model.G if (b,g) in model.Gbs) +\
    sum(model.pW[w,t] for w in model.WIND if (b,w) in model.Wbs)+\
    sum(model.PS[s,t] for s in model.S if (b,s) in model.Sbs)+\
    sum(model.pDis[a,t] for a in model.BATTERY if (b,a) in model.Bbs)==\
    sum(model.pChar[a,t] for a in model.BATTERY if (b,a) in model.Bbs)+\
    sum(model.pD[d,t] for d in model.D if (b,d) in model.Dbs)+\
    sum(model.pLfrom[l,t] for l in model.L if model.A[l,1]==b)+ \
    sum(model.pLto[l,t] for l in model.L if model.A[l,2]==b)+\
    sum(model.pLfromT[l,t] for l in model.TRANSF if model.AT[l,1]==b)+ \
    sum(model.pLtoT[l,t] for l in model.TRANSF if model.AT[l,2]==b)+\
    sum(model.GB[s]*model.v[b,t]**2 for s in model.SHUNT if (b,s) in model.SHUNTbs) 
def KCL_reactive_def(model, b, t):
    return sum(model.qG[g,t] for g in model.G if (b,g) in model.Gbs) + \
    sum(model.QS[s,t] for s in model.S if (b,s) in model.Sbs)+\
    sum(model.qW[w,t] for w in model.WIND if (b,w) in model.Wbs) == \
    sum(model.qD[d,t] for d in model.D if (b,d) in model.Dbs)+\
    sum(model.qLfrom[l,t] for l in model.L if model.A[l,1]==b)+ \
    sum(model.qLto[l,t] for l in model.L if model.A[l,2]==b)+\
    sum(model.qLfromT[l,t] for l in model.TRANSF if model.AT[l,1]==b)+ \
    sum(model.qLtoT[l,t] for l in model.TRANSF if model.AT[l,2]==b)-\
    sum(model.BB[s]*model.v[b,t]**2 for s in model.SHUNT if (b,s) in model.SHUNTbs)
model.KCL_real     = Constraint(model.B, model.T, rule=KCL_real_def)
model.KCL_reactive = Constraint(model.B, model.T, rule=KCL_reactive_def)

'''
def KCL_real_def(model, b, t):
    return sum(model.pG[g, t] for g in model.G if (b,g) in model.Gbs) +\
    sum(model.pW[w,t] for w in model.WIND if (b,w) in model.Wbs)==\
    sum(model.pD[d,t] for d in model.D if (b,d) in model.Dbs)+\
    sum(model.pLfrom[l,t] for l in model.L if model.A[l,1]==b)+ \
    sum(model.pLto[l,t] for l in model.L if model.A[l,2]==b)+\
    sum(model.pLfromT[l,t] for l in model.TRANSF if model.AT[l,1]==b)+ \
    sum(model.pLtoT[l,t] for l in model.TRANSF if model.AT[l,2]==b)+\
    sum(model.GB[s]*model.v[b,t]**2 for s in model.SHUNT if (b,s) in model.SHUNTbs)
def KCL_reactive_def(model, b, t):
    return sum(model.qG[g,t] for g in model.G if (b,g) in model.Gbs) +\
    sum(model.qW[w,t] for w in model.WIND if (b,w) in model.Wbs)== \
    sum(model.qD[d,t] for d in model.D if (b,d) in model.Dbs)+\
    sum(model.qLfrom[l,t] for l in model.L if model.A[l,1]==b)+ \
    sum(model.qLto[l,t] for l in model.L if model.A[l,2]==b)+\
    sum(model.qLfromT[l,t] for l in model.TRANSF if model.AT[l,1]==b)+ \
    sum(model.qLtoT[l,t] for l in model.TRANSF if model.AT[l,2]==b)-\
    sum(model.BB[s]*model.v[b,t]**2 for s in model.SHUNT if (b,s) in model.SHUNTbs)
model.KCL_real     = Constraint(model.B, model.T, rule=KCL_real_def)
model.KCL_reactive = Constraint(model.B, model.T, rule=KCL_reactive_def)
'''
# --- Kirchoff's voltage law on each line ---
# TODO: Bengisu, make time-variant
# battery systems in voltage law relevance?

def KVL_real_fromend(model,l,t):
    return model.pLfrom[l,t] == model.G11[l]*(model.v[model.A[l,1],t]**2)+\
    model.v[model.A[l,1],t]*model.v[model.A[l,2],t]*(model.B12[l]*sin(model.delta[model.A[l,1],t]-\
    model.delta[model.A[l,2],t])+model.G12[l]*cos(model.delta[model.A[l,1],t]-model.delta[model.A[l,2],t]))
def KVL_real_toend(model,l,t):
    return model.pLto[l,t] == model.G22[l]*(model.v[model.A[l,2],t]**2)+\
    model.v[model.A[l,1],t]*model.v[model.A[l,2],t]*(model.B21[l]*sin(model.delta[model.A[l,2],t]-\
    model.delta[model.A[l,1],t])+model.G21[l]*cos(model.delta[model.A[l,2],t]-model.delta[model.A[l,1],t]))
def KVL_reactive_fromend(model,l,t):
    return model.qLfrom[l,t] == -model.B11[l]*(model.v[model.A[l,1],t]**2)+\
    model.v[model.A[l,1],t]*model.v[model.A[l,2],t]*(model.G12[l]*sin(model.delta[model.A[l,1],t]-\
    model.delta[model.A[l,2],t])-model.B12[l]*cos(model.delta[model.A[l,1],t]-model.delta[model.A[l,2],t]))
def KVL_reactive_toend(model,l,t):
    return model.qLto[l,t] == (-model.B22[l]*(model.v[model.A[l,2],t]**2)+\
    model.v[model.A[l,1],t]*model.v[model.A[l,2],t]*(model.G21[l]*sin(model.delta[model.A[l,2],t]-\
    model.delta[model.A[l,1],t])-model.B21[l]*cos(model.delta[model.A[l,2],t]-model.delta[model.A[l,1],t])))
model.KVL_real_from     = Constraint(model.L, model.T, rule=KVL_real_fromend)
model.KVL_real_to       = Constraint(model.L, model.T, rule=KVL_real_toend)
model.KVL_reactive_from = Constraint(model.L, model.T, rule=KVL_reactive_fromend)
model.KVL_reactive_to   = Constraint(model.L, model.T, rule=KVL_reactive_toend)

# --- Kirchoff's voltage law on each transformer line ---
# TODO: Bengisu, make time-variant

def KVL_real_fromendTransf(model,l,t):
    return model.pLfromT[l,t] == model.G11T[l]*(model.v[model.AT[l,1],t]**2)+\
    model.v[model.AT[l,1],t]*model.v[model.AT[l,2],t]*(model.B12T[l]*sin(model.delta[model.AT[l,1],t]-\
    model.delta[model.AT[l,2],t])+model.G12T[l]*cos(model.delta[model.AT[l,1],t]-model.delta[model.AT[l,2],t]))
def KVL_real_toendTransf(model,l,t):
    return model.pLtoT[l,t] == model.G22T[l]*(model.v[model.AT[l,2],t]**2)+\
    model.v[model.AT[l,1],t]*model.v[model.AT[l,2],t]*(model.B21T[l]*sin(model.delta[model.AT[l,2],t]-\
    model.delta[model.AT[l,1],t])+model.G21T[l]*cos(model.delta[model.AT[l,2],t]-model.delta[model.AT[l,1],t]))
def KVL_reactive_fromendTransf(model,l,t):
    return model.qLfromT[l,t] == -model.B11T[l]*(model.v[model.AT[l,1],t]**2)+\
    model.v[model.AT[l,1],t]*model.v[model.AT[l,2],t]*(model.G12T[l]*sin(model.delta[model.AT[l,1],t]-\
    model.delta[model.AT[l,2],t])-model.B12T[l]*cos(model.delta[model.AT[l,1],t]-model.delta[model.AT[l,2],t]))
def KVL_reactive_toendTransf(model,l,t):
    return model.qLtoT[l,t] == -model.B22T[l]*(model.v[model.AT[l,2],t]**2)+\
    model.v[model.AT[l,1],t]*model.v[model.AT[l,2],t]*(model.G21T[l]*sin(model.delta[model.AT[l,2],t]-\
    model.delta[model.AT[l,1],t])-model.B21T[l]*cos(model.delta[model.AT[l,2],t]-model.delta[model.AT[l,1],t]))
model.KVL_real_fromTransf     = Constraint(model.TRANSF, model.T, rule=KVL_real_fromendTransf)
model.KVL_real_toTransf       = Constraint(model.TRANSF, model.T, rule=KVL_real_toendTransf)
model.KVL_reactive_fromTransf = Constraint(model.TRANSF, model.T, rule=KVL_reactive_fromendTransf)
model.KVL_reactive_toTransf   = Constraint(model.TRANSF, model.T, rule=KVL_reactive_toendTransf)

# --- generator power limits ---
# TODO: Bengisu, make time-variant
def Real_Power_Max(model,g, t):
    return model.pG[g, t] <= model.PGmax[g]
def Real_Power_Min(model,g, t):
    return model.pG[g, t] >= model.PGmin[g]
def Reactive_Power_Max(model,g, t):
    return model.qG[g, t] <= model.QGmax[g]
def Reactive_Power_Min(model,g, t):
    return model.qG[g, t] >= model.QGmin[g]
model.PGmaxC = Constraint(model.G, model.T, rule=Real_Power_Max)
model.PGminC = Constraint(model.G, model.T, rule=Real_Power_Min)
model.QGmaxC = Constraint(model.G, model.T, rule=Reactive_Power_Max)
model.QGminC = Constraint(model.G, model.T, rule=Reactive_Power_Min)

# ---wind generator power limits ---
# TODO: Bengisu, make time-variant
def Wind_Real_Power_Max(model,w, t):
    return model.pW[w,t] <= model.WGmax[w]
def Wind_Real_Power_Min(model,w,t):
    return model.pW[w,t] >= model.WGmin[w]
def Wind_Reactive_Power_Max(model,w,t):
    return model.qW[w,t] <= model.WGQmax[w]
def Wind_Reactive_Power_Min(model,w,t):
    return model.qW[w,t] >= model.WGQmin[w]
model.WGmaxC  = Constraint(model.WIND, model.T, rule=Wind_Real_Power_Max)
model.WGminC  = Constraint(model.WIND, model.T, rule=Wind_Real_Power_Min)
model.WGQmaxC = Constraint(model.WIND, model.T, rule=Wind_Reactive_Power_Max)
model.WGQminC = Constraint(model.WIND, model.T, rule=Wind_Reactive_Power_Min)


# --- demand and load shedding ---
# TODO: Bengisu, make time-variant
def Load_Shed_real(model,d,t):
    return model.pD[d,t] == model.alpha[d]*model.PD[d,t]
def Load_Shed_reactive(model,d,t):
    return model.qD[d,t] == model.alpha[d]*model.QD[d]
def alpha_FixNegDemands(model,d):
    return model.alpha[d] == 1

model.LoadShed_real     = Constraint(model.D, model.T, rule=Load_Shed_real)
model.LoadShed_reactive = Constraint(model.D, model.T, rule=Load_Shed_reactive)
model.alphaFix          = Constraint(model.DNeg, rule=alpha_FixNegDemands)


def alpha_BoundUB(model,d):
    return model.alpha[d] <= 1
def alpha_BoundLB(model,d):
    return model.alpha[d] >= 0
model.alphaBoundUBC = Constraint(model.D, rule=alpha_BoundUB)
model.alphaBoundLBC = Constraint(model.D, rule=alpha_BoundLB)

# --- line power limits ---
# TODO: Bengisu, make time-variant
def line_lim1_def(model,l,t):
    return model.pLfrom[l,t]**2+model.qLfrom[l,t]**2 <= model.SLmax[l]**2
def line_lim2_def(model,l,t):
    return model.pLto[l,t]**2+model.qLto[l,t]**2 <= model.SLmax[l]**2
model.line_lim1 = Constraint(model.L, model.T, rule=line_lim1_def)
model.line_lim2 = Constraint(model.L, model.T, rule=line_lim2_def)

# --- power flow limits on transformer lines---
# TODO: Bengisu, make time-variant
def transf_lim1_def(model,l,t):
    return model.pLfromT[l,t]**2+model.qLfromT[l,t]**2 <= model.SLmaxT[l]**2
def transf_lim2_def(model,l,t):
    return model.pLtoT[l,t]**2+model.qLtoT[l,t]**2 <= model.SLmaxT[l]**2
model.transf_lim1 = Constraint(model.TRANSF, model.T, rule=transf_lim1_def)
model.transf_lim2 = Constraint(model.TRANSF, model.T, rule=transf_lim2_def)

# --- voltage constraints ---
# TODO: Bengisu, make time-variant
def bus_max_voltage(model,b,t):
    return model.v[b,t] <= model.Vmax[b]
def bus_min_voltage(model,b,t):
    return model.v[b,t] >= model.Vmin[b]
model.Vmaxc = Constraint(model.B, model.T, rule=bus_max_voltage)
model.Vminc = Constraint(model.B, model.T, rule=bus_min_voltage)

# --- reference bus constraint ---
# TODO: Bengisu, make time-variant
def ref_bus_def(model,b,t):
    return model.delta[b,t]==0
model.refbus = Constraint(model.b0, model.T, rule=ref_bus_def)


# ---battery power and energy limits
'''
#--- EXACT FORMULATION ---
def Battery_Charge_Real_Power_Max(model, a, t):
    return model.pChar[a,t] <= model.PCmax[a]*model.z[a,t]
def Battery_Charge_Real_Power_Min(model, a, t):
    return model.pChar[a,t] >= model.PCmin[a]*model.z[a,t]
def Battery_Discharge_Real_Power_Max(model, a, t):
    return model.pDis[a,t] <= model.PDmax[a]*model.y[a,t]
def Battery_Discharge_Real_Power_Min(model, a ,t):
    return model.pDis[a,t] >= model.PDmin[a]*model.y[a,t]
def complementary_rule(model, a, t):
    return model.y[a,t] + model.z[a,t] == 1
def Battery_Energy_Max(model, a ,t):
    return model.e[a,t] <= model.EMAX[a]
def Battery_Energy_Min(model, a, t):
    return model.e[a,t] >= model.EMIN[a]
'''

 # --- SIMPLIFIED FORMULATION---
def Battery_Charge_Real_Power_Max(model, a, t):
    return model.pChar[a,t] <= model.PCmax[a]
def Battery_Charge_Real_Power_Min(model, a, t):
    return model.pChar[a,t] >= model.PCmin[a]
def Battery_Discharge_Real_Power_Max(model, a, t):
    return model.pDis[a,t] <= model.PDmax[a]
def Battery_Discharge_Real_Power_Min(model, a ,t):
    return model.pDis[a,t] >= model.PDmin[a]
def Battery_Energy_Max(model, a ,t):
    return model.e[a,t] <= model.EMAX[a]
def Battery_Energy_Min(model, a, t):
    return model.e[a,t] >= model.EMIN[a]

'''
#---EXTENDED MODEL---
def Battery_Charge_Real_Power_Max(model, a, t):
    return model.pChar[a,t] <= model.PCmax[a]
def Battery_Charge_Real_Power_Min(model, a, t):
    return model.pChar[a,t] >= model.PCmin[a]
def Battery_Discharge_Real_Power_Max(model, a, t):
    return model.pDis[a,t] <= model.PDmax[a]
def Battery_Discharge_Real_Power_Min(model, a ,t):
    return model.pDis[a,t] >= model.PDmin[a]
def Battery_Energy_Max(model, a ,t):
    return model.e[a,t] <= model.EMAX[a]
def Battery_Energy_Min(model, a, t):
    return model.e[a,t] >= model.EMIN[a]
def Charge_Constraint(model, a, t):
    return model.pChar[a,t] <= ((model.EMAX[a]-model.E0[a]) / model.nchar[a])
def Discharge_Constraint(model, a ,t):
    return model.pDis[a,t] <= ((model.E0[a]-model.EMIN[a])*model.ndis[a])
def Discharge_Charge_Constraint(model, a, t):
    return model.pDis[a,t] <= (model.PDmax[a] - ((model.PDmax[a]/model.PCmax[a])*model.pChar[a,t]))
'''
model.Battery_Charge_Max = Constraint(model.BATTERY, model.T, rule=Battery_Charge_Real_Power_Max)
model.Battery_Charge_Min = Constraint(model.BATTERY, model.T, rule=Battery_Charge_Real_Power_Min)
model.Battery_Discharge_Max = Constraint(model.BATTERY, model.T, rule=Battery_Discharge_Real_Power_Max)
model.Battery_Discharge_Min = Constraint(model.BATTERY, model.T, rule=Battery_Discharge_Real_Power_Min)
model.Battery_EMax = Constraint(model.BATTERY, model.T, rule=Battery_Energy_Max)
model.Battery_Emin = Constraint(model.BATTERY, model.T, rule=Battery_Energy_Min)
#model.complementary = Constraint(model.BATTERY, model.T, rule=complementary_rule)
#model.Charge_Con = Constraint(model.BATTERY, model.T, rule=Charge_Constraint)
#model.Discharge_Con = Constraint(model.BATTERY, model.T, rule=Discharge_Constraint)
#model.Discharge_Charge = Constraint(model.BATTERY, model.T, rule=Discharge_Charge_Constraint)


# --- battery energy level ---
def energy_level_battery(model, a, t):
    if t==0:
        return model.e[a,0] == model.E0[a]
    else:
        return model.e[a,t] == model.e[a,t-1] + model.nchar[a]*model.pChar[a,t-1] - model.pDis[a,t-1]/model.ndis[a]
model.SoC = Constraint(model.BATTERY, model.T, rule=energy_level_battery)
# burada t mi yoksa t-1 mi yapmam gerekiyor emin olamadim

'''
#battery starting
def  battery_start(model, a):
  return model.e[a,0] == model.E0[a]
model.bt_start = Constraint(model.BATTERY, rule=battery_start)
'''
#start-end battery loop
def battery_loop(model, a):
   return model.e[a, model.Tend.value] == model.e[a, 0]
model.loop = Constraint(model.BATTERY, rule=battery_loop)


