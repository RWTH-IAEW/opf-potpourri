from pyomo.environ import *
from potpourri.models.class_based.basemodel import Basemodel
import numpy as np


class DC(Basemodel):
    def __init__(self, net):
        super().__init__(net)

        ZN = self.net.bus.vn_kv ** 2 / self.baseMVA
        bik = - self.net.line.x_ohm_per_km / (
                    (self.net.line.x_ohm_per_km ** 2 + self.net.line.r_ohm_per_km ** 2) * self.net.line.length_km)

        self.BL_data = bik * ZN[self.net.line.from_bus].values

        r_k = self.net.trafo.vkr_percent / 100 * (self.net.sn_mva / self.net.trafo.sn_mva)
        z_k = self.net.trafo.vk_percent / 100 * (self.net.sn_mva / self.net.trafo.sn_mva)
        x_k = np.sqrt(z_k ** 2 - r_k ** 2)

        self.BLT_data = - x_k / z_k ** 2

    def create_model(self):
        super().create_model()
        # lines and transformer chracteristics
        self.model.BL = Param(self.model.L, within=Reals, initialize=self.BL_data)  # susceptance of a line
        self.model.BLT = Param(self.model.TRANSF, within=Reals,
                               initialize=self.BLT_data)  # susceptance of a transformer

        # --- Variables ---
        self.model.deltaL = Var(self.model.L, domain=Reals)  # angle difference across lines
        self.model.deltaLT = Var(self.model.TRANSF, domain=Reals)  # angle difference across transformers

        # --- Kirchoff's current law at each bus b ---
        def KCL_def(model, b):
            return (sum(model.pG[g] for g in model.G if (b, g) in model.Gbs) +
                    sum(model.peG[g] for g in model.eG if (b, g) in model.eGbs) ==
                    sum(model.pD[d] for d in model.D if (b, d) in model.Dbs) +
                    sum(model.pLfrom[l] for l in model.L if model.A[l, 1] == b) +
                    sum(model.pLto[l] for l in model.L if model.A[l, 2] == b) +
                    sum(model.pLfromT[l] for l in model.TRANSF if model.AT[l, 1] == b) +
                    sum(model.pLtoT[l] for l in model.TRANSF if model.AT[l, 2] == b) +
                    sum(model.GB[s] for s in model.SHUNT if (b, s) in model.SHUNTbs))

        self.model.KCL_const = Constraint(self.model.B, rule=KCL_def)

        def KVL_real_fromend(model, l):
            return model.pLfrom[l] == (-model.BL[l]) * model.deltaL[l]

        def KVL_real_toend(model, l):
            return model.pLto[l] == (model.BL[l]) * model.deltaL[l]

        self.model.KVL_real_from = Constraint(self.model.L, rule=KVL_real_fromend)
        self.model.KVL_real_to = Constraint(self.model.L, rule=KVL_real_toend)

        # --- phase angle constraints ---
        def phase_angle_diff1(model, l):
            return model.deltaL[l] == model.delta[model.A[l, 1]] - \
                model.delta[model.A[l, 2]]

        self.model.phase_diff1 = Constraint(self.model.L, rule=phase_angle_diff1)

        # --- phase angle constraints ---
        def phase_angle_diff2(model, l):
            return model.deltaLT[l] == model.delta[model.AT[l, 1]] - \
                model.delta[model.AT[l, 2]]

        self.model.phase_diff2 = Constraint(self.model.TRANSF, rule=phase_angle_diff2)
