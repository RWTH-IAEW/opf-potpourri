import pandas as pd
import pandapower as pp
import copy
from pyomo.environ import *
from potpourri.models.class_based.basemodel import Basemodel


class DC(Basemodel):
    def __init__(self, net):
        super().__init__(net)

        x = self.net._ppc["branch"][:, 3].real
        BL = -1/x
        trafo_start = len(self.net.line)
        trafo_end = trafo_start + len(self.net.trafo)

        self.trafo_data = self.trafo_data.assign(**{"BLT_data": BL[trafo_start:trafo_end]})
        # self.line_data = pd.DataFrame(BL[:trafo_start], columns=["BL_data"])
        self.line_data['BL_data'] = BL[:trafo_start]

        ZN = self.net.bus.vn_kv ** 2 / self.baseMVA
        y_s = - 1 / (self.net.line.x_ohm_per_km * self.net.line.length_km)      # according to matpower manual dc modeling

        self.BL_data = y_s * ZN[self.net.line.from_bus].values

        # self.BLT_data = pd.Series(-1 / self.trafo_parameters["x"] / self.trafo_parameters["ratio"][0], self.trafo_set)

        self.create_model()

    def create_model(self):
        super().create_model()

        self.model.name = "DC"

        # lines and transformer chracteristics
        self.model.BL = Param(self.model.L, within=Reals, initialize=self.BL_data[self.model.L])  # susceptance of a line
        self.model.BLT = Param(self.model.TRANSF, within=Reals,
                               initialize=self.trafo_data.BLT_data[self.model.TRANSF])  # susceptance of a transformer

        # --- Variables ---
        self.model.deltaL = Var(self.model.L, domain=Reals)  # angle difference across lines
        self.model.deltaLT = Var(self.model.TRANSF, domain=Reals)  # angle difference across transformers

        # --- Kirchoff's current law at each bus b ---
        def KCL_def(model, b):
            kcl = (sum(model.pG[g] for g in model.G if (b, g) in model.Gbs) +
                    sum(model.peG[g] for g in model.eG if (b, g) in model.eGbs) ==
                    sum(model.pD[d] for d in model.D if (b, d) in model.Dbs) +
                    sum(model.pLfrom[l] for l in model.L if model.A[l, 1] == b) +
                    sum(model.pLto[l] for l in model.L if model.A[l, 2] == b) +
                    sum(model.pThv[l] for l in model.TRANSF if model.AT[l, 1] == b) +
                    sum(model.pTlv[l] for l in model.TRANSF if model.AT[l, 2] == b) +
                    sum(model.GB[s] for s in model.SHUNT if (b, s) in model.SHUNTbs))
            if isinstance(kcl, bool):
                return Constraint.Skip
            return kcl

        self.model.KCL_const = Constraint(self.model.B, rule=KCL_def)

        # --- Kirchoff's voltage law at each line and transformer---
        def KVL_real_fromend(model, l):
            return model.pLfrom[l] == (-model.BL[l]) * model.deltaL[l]

        def KVL_real_toend(model, l):
            return model.pLto[l] == (model.BL[l]) * model.deltaL[l]

        self.model.KVL_real_from = Constraint(self.model.L, rule=KVL_real_fromend)
        self.model.KVL_real_to = Constraint(self.model.L, rule=KVL_real_toend)

        def KVL_trans_fromend(model, l):
            return model.pThv[l] == (-model.BLT[l]) * (model.deltaLT[l])
        def KVL_trans_toend(model, l):
            return model.pTlv[l] == (model.BLT[l]) * (model.deltaLT[l])

        self.model.KVL_trans_from = Constraint(self.model.TRANSF, rule=KVL_trans_fromend)
        self.model.KVL_trans_to = Constraint(self.model.TRANSF, rule=KVL_trans_toend)

        # --- phase angle constraints ---
        def phase_angle_diff1(model, l):
            return model.deltaL[l] == model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]

        self.model.phase_diff1 = Constraint(self.model.L, rule=phase_angle_diff1)

        # --- phase angle constraints ---
        def phase_angle_diff2(model, l):
            return model.deltaLT[l] == model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]] - model.shift[l]

        self.model.phase_diff2 = Constraint(self.model.TRANSF, rule=phase_angle_diff2)
