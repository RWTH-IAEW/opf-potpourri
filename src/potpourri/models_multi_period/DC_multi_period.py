from pyomo.environ import *
from potpourri.models_multi_period.basemodel_multi_period import Basemodel_multi_period
#TODO make multiperiod

class DC_multi_period(Basemodel_multi_period):
    def __init__(self, net, toT, fromT=None):
        super().__init__(net, toT, fromT)

        x = self.net._ppc["branch"][:, 3].real
        BL = -1 / x
        trafo_start = len(self.net.line)
        trafo_end = trafo_start + len(self.net.trafo)

        self.trafo_data = self.trafo_data.assign(**{"BLT_data": BL[trafo_start:trafo_end]})
        # self.line_data = pd.DataFrame(BL[:trafo_start], columns=["BL_data"])
        self.line_data['BL_data'] = BL[:trafo_start]

        ZN = self.net.bus.vn_kv ** 2 / self.baseMVA
        y_s = - 1 / (self.net.line.x_ohm_per_km * self.net.line.length_km)  # according to matpower manual dc modeling

        self.BL_data = y_s * ZN[self.net.line.from_bus].values

        # self.BLT_data = pd.Series(-1 / self.trafo_parameters["x"] / self.trafo_parameters["ratio"][0], self.trafo_set)

        self.create_model()

    def create_model(self):
        super().create_model()

        self.model.name = "DC"

        # lines and transformer chracteristics
        self.model.BL = Param(self.model.L, within=Reals,
                              initialize=self.BL_data[self.model.L])  # susceptance of a line
        self.model.BLT = Param(self.model.TRANSF, within=Reals,
                               initialize=self.trafo_data.BLT_data[self.model.TRANSF])  # susceptance of a transformer

        # --- Variables ---
        self.model.deltaL = Var(self.model.L, self.model.T, domain=Reals)  # angle difference across lines
        self.model.deltaLT = Var(self.model.TRANSF, self.model.T, domain=Reals)  # angle difference across transformers

        if self.T == None:
            @self.model.Constraint(self.model.B)
            def KCL_def(model, b):
                kcl = (sum(self.model.psG[g] for g in self.model.sG if (g, b) in self.model.sGbs) +
                       sum(self.model.pG[g] for g in self.model.G if (g, b) in self.model.Gbs) ==
                       sum(self.model.pD[d] for d in self.model.D if (b, d) in self.model.Dbs) +
                       sum(self.model.pLfrom[l] for l in self.model.L if self.model.A[l, 1] == b) +
                       sum(self.model.pLto[l] for l in self.model.L if self.model.A[l, 2] == b) +
                       sum(self.model.pThv[l] for l in self.model.TRANSF if self.model.AT[l, 1] == b) +
                       sum(self.model.pTlv[l] for l in self.model.TRANSF if self.model.AT[l, 2] == b) +
                       sum(self.model.GB[s] for s in self.model.SHUNT if (b, s) in self.model.SHUNTbs))
                if isinstance(kcl, bool):
                    return Constraint.Skip
                return kcl == 0

            # --- Kirchoff's voltage law at each line and transformer---
            @self.model.Constraint(self.model.L)
            def KVL_real_fromend(model, l):
                return model.pLfrom[l] == (-model.BL[l]) * model.deltaL[l]

            @self.model.Constraint(self.model.L)
            def KVL_real_toend(model, l):
                return model.pLto[l] == (model.BL[l]) * model.deltaL[l]

            @self.model.Constraint(self.model.TRANSF, self.model.T)
            def KVL_trans_fromend(model, l):
                return model.pThv[l] == (-model.BLT[l]) * (model.deltaLT[l])

            @self.model.Constraint(self.model.TRANSF)
            def KVL_trans_toend(model, l):
                return model.pTlv[l] == (model.BLT[l]) * (model.deltaLT[l])

            # --- phase angle constraints ---
            @self.model.Constraint(self.model.L)
            def phase_angle_diff1(model, l):
                return model.deltaL[l] == model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]

            # --- phase angle constraints ---
            @self.model.Constraint(self.model.TRANSF)
            def phase_angle_diff2(model, l):
                return model.deltaLT[l] == model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]] - model.shift[l]

        else:
            # --- Kirchoff's current law at each bus b time variant---
            @self.model.Constraint(self.model.B, self.T)
            def KCL_def(model, b, t):
                kcl = (sum(self.model.psG[g, t] for g in self.model.sG if (g, b) in self.model.sGbs) +
                       sum(self.model.pG[g, t] for g in self.model.G if (g, b) in self.model.Gbs) ==
                       sum(self.model.pD[d, t] for d in self.model.D if (b, d) in self.model.Dbs) +
                       sum(self.model.pLfrom[l, t] for l in self.model.L if self.model.A[l, 1] == b) +
                       sum(self.model.pLto[l, t] for l in self.model.L if self.model.A[l, 2] == b) +
                       sum(self.model.pThv[l, t] for l in self.model.TRANSF if self.model.AT[l, 1] == b) +
                       sum(self.model.pTlv[l, t] for l in self.model.TRANSF if self.model.AT[l, 2] == b) +
                       sum(self.model.GB[s, t] for s in self.model.SHUNT if (b, s) in self.model.SHUNTbs))
                if isinstance(kcl, bool):
                    return Constraint.Skip
                return kcl == 0

            # --- Kirchoff's voltage law at each line and transformer---
            @self.model.Constraint(self.model.L, self.model.T)
            def KVL_real_fromend(model, l, t):
                return model.pLfrom[l, t] == (-model.BL[l, t]) * model.deltaL[l]

            @self.model.Constraint(self.model.L, self.model.T)
            def KVL_real_toend(model, l, t):
                return model.pLto[l, t] == (model.BL[l, t]) * model.deltaL[l, t]

            @self.model.Constraint(self.model.TRANSF, self.model.T)
            def KVL_trans_fromend(model, l, t):
                return model.pThv[l, t] == (-model.BLT[l, t]) * (model.deltaLT[l, t])

            @self.model.Constraint(self.model.TRANSF, self.model.T)
            def KVL_trans_toend(model, l, t):
                return model.pTlv[l, t] == (model.BLT[l, t]) * (model.deltaLT[l, t])

            # --- phase angle constraints ---
            @self.model.Constraint(self.model.L, self.model.T)
            def phase_angle_diff1(model, l, t):
                return model.deltaL[l, t] == model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]

            # --- phase angle constraints ---
            @self.model.Constraint(self.model.TRANSF, self.model.T)
            def phase_angle_diff2(model, l, t):
                return model.deltaLT[l, t] == model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]] - \
                    model.shift[l, t]
