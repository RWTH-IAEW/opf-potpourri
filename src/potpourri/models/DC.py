"""DC power flow mixin: adds linearised DC equations (voltage angles only)
to Basemodel."""

import numpy as np
import pandas as pd
import pyomo.environ as pyo
from potpourri.models.basemodel import Basemodel


class DC(Basemodel):
    """Linearised DC power flow model for distribution network analysis.

    Extends Basemodel with lossless DC power flow equations. Voltage
    magnitudes are assumed 1.0 p.u. and reactive power is ignored, reducing
    the power flow to a linear system in voltage angles.

    Args:
        net: A pandapower network compatible with pp.runpp().
        dc_convention: Which linearisation the DC flow follows.
            ``"matpower"`` (default): branch susceptance ``-1/x`` and the
            transformer phase shift in the angle difference, the textbook DC
            power flow and MATPOWER's ``makeBdc``. ``"powermodels"``: the
            series susceptance of the full impedance at the nominal tap
            position, ``-x/(r² + x²)``, and no phase shift, i.e.
            ``p = -b (θ_from − θ_to)`` as in PowerModels.jl's
            ``DCPPowerModel`` that produced the PGLib-OPF DC reference values
            (neither tap nor shift enter its DC flow; pandapower's referral
            of a transformer's impedance to the tapped LV voltage is undone
            for the same reason). The two
            agree where r ≪ x and no phase shifters exist, and differ by
            several percent on high-r/x networks such as the RTE cases.
    """

    DC_CONVENTIONS = ("matpower", "powermodels")

    def __init__(self, net, dc_convention: str = "matpower"):
        if dc_convention not in self.DC_CONVENTIONS:
            raise ValueError(
                f"dc_convention must be one of {self.DC_CONVENTIONS}, "
                f"got {dc_convention!r}"
            )
        self.dc_convention = dc_convention
        super().__init__(net)

        r = self.net._ppc["branch"][:, 2].real.copy()
        x = self.net._ppc["branch"][:, 3].real.copy()
        trafo_start = len(self.net.line)
        trafo_end = trafo_start + len(self.net.trafo)
        if dc_convention == "powermodels":
            # pandapower refers a transformer's series impedance to the
            # tapped LV voltage when the tap sits on the LV side, so r and x
            # carry a factor (vn_trafo_lv / vn_lv_kv)². PowerModels reads the
            # MATPOWER reactance, which no tap enters; undo the factor.
            scale = self._lv_tap_impedance_scale()
            r[trafo_start:trafo_end] /= scale
            x[trafo_start:trafo_end] /= scale
            BL = -x / (r**2 + x**2)
        else:
            BL = -1 / x
        imp_table = self.net.get("impedance")
        n_imp = (
            len(imp_table)
            if imp_table is not None and not imp_table.empty
            else 0
        )

        self.trafo_data = self.trafo_data.assign(
            **{"BLT_data": BL[trafo_start:trafo_end]}
        )

        # DC susceptance for lines: −1/x (matpower DC convention).
        # For native lines we recompute from net.line for traceability;
        # for impedance branches (which carry per-unit r/x referenced to
        # impedance.sn_mva, in actual ohms per the from-bus), we read the
        # per-unit value pandapower already wrote into _ppc['branch'].
        if dc_convention == "powermodels":
            bl_line = BL[:trafo_start]
        else:
            ZN = self.net.bus.vn_kv**2 / self.baseMVA
            y_s_line = -1 / (
                self.net.line.x_ohm_per_km * self.net.line.length_km
            )  # according to matpower manual dc modeling
            bl_line = y_s_line * ZN[self.net.line.from_bus].values
        # Use a pd.Series keyed by the model.L indices (native line indices +
        # synthetic impedance indices) so out-of-service line rows — whose
        # indices are skipped in model.L — don't shift the lookup.
        line_index = list(self.net.line.index)
        if n_imp:
            bl_imp = BL[trafo_end : trafo_end + n_imp]
            self.BL_data = pd.Series(
                np.concatenate([bl_line, bl_imp]),
                index=line_index + [trafo_start + i for i in range(n_imp)],
            )
            line_idx_ppc = np.r_[
                np.arange(0, trafo_start),
                np.arange(trafo_end, trafo_end + n_imp),
            ]
        else:
            self.BL_data = pd.Series(bl_line, index=line_index)
            line_idx_ppc = np.arange(0, trafo_start)
        self.line_data["BL_data"] = BL[line_idx_ppc]

        self.create_model()

    def _lv_tap_impedance_scale(self):
        """Factor pandapower applied to each transformer's series impedance
        for a tap on the LV side: ``(vn_trafo_lv / vn_lv_kv)²`` per
        ``net.trafo`` row, 1 where the tap sits on the HV side or is absent.
        """
        trafo = self.net.trafo
        if trafo.empty:
            return np.ones(0)
        rated = trafo["vn_lv_kv"].to_numpy(dtype=float)
        try:
            from pandapower.build_branch import _calc_tap_from_dataframe

            _, tapped, _ = _calc_tap_from_dataframe(self.net, trafo)
        except Exception:  # noqa: BLE001 — no _options yet, or an old API
            tapped = rated.copy()
            if "tap_pos" in trafo:
                pos = trafo["tap_pos"].to_numpy(dtype=float)
                neutral = trafo.get(
                    "tap_neutral", pd.Series(0.0, index=trafo.index)
                ).to_numpy(dtype=float)
                step = trafo.get(
                    "tap_step_percent", pd.Series(0.0, index=trafo.index)
                ).to_numpy(dtype=float)
                on_lv = (trafo["tap_side"] == "lv").to_numpy()
                steps = np.nan_to_num(pos) - np.nan_to_num(neutral)
                ratio = 1.0 + steps * np.nan_to_num(step) / 100.0
                tapped = np.where(on_lv, rated * ratio, rated)
        return (np.asarray(tapped, dtype=float) / rated) ** 2

    def create_model(self):
        """Build the Pyomo ConcreteModel with DC power flow constraints.

        Adds line susceptance parameters (BL, BLT), angle difference variables
        (deltaL, deltaLT), real-power KCL at each bus, and KVL constraints
        linking power flows to angle differences. All quantities in per-unit.
        """
        super().create_model()

        self.model.name = "DC"

        # lines and transformer chracteristics
        self.model.BL = pyo.Param(
            self.model.L,
            within=pyo.Reals,
            initialize=self.BL_data[self.model.L],
        )  # susceptance of a line
        self.model.BLT = pyo.Param(
            self.model.TRANSF,
            within=pyo.Reals,
            initialize=self.trafo_data.BLT_data[self.model.TRANSF],
        )  # susceptance of a transformer

        # --- pyo.Variables ---
        self.model.deltaL = pyo.Var(
            self.model.L, domain=pyo.Reals
        )  # angle difference across lines
        self.model.deltaLT = pyo.Var(
            self.model.TRANSF, domain=pyo.Reals
        )  # angle difference across transformers

        # --- Kirchoff's current law at each bus b ---
        def KCL_def(model, b):
            kcl = sum(
                model.psG[g] for g in model.sG if (g, b) in model.sGbs
            ) + sum(
                model.pG[g] for g in model.G if (g, b) in model.Gbs
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
                model.GB[s] for s in model.SHUNT if (b, s) in model.SHUNTbs
            )
            if isinstance(kcl, (bool, np.bool_)):
                return pyo.Constraint.Skip
            return kcl

        self.model.KCL_const = pyo.Constraint(self.model.B, rule=KCL_def)

        # --- Kirchoff's voltage law at each line and transformer---
        def KVL_real_fromend(model, l):
            return model.pLfrom[l] == (-model.BL[l]) * model.deltaL[l]

        def KVL_real_toend(model, l):
            return model.pLto[l] == (model.BL[l]) * model.deltaL[l]

        self.model.KVL_real_from = pyo.Constraint(
            self.model.L, rule=KVL_real_fromend
        )
        self.model.KVL_real_to = pyo.Constraint(
            self.model.L, rule=KVL_real_toend
        )

        def KVL_trans_fromend(model, l):
            return model.pThv[l] == (-model.BLT[l]) * (model.deltaLT[l])

        def KVL_trans_toend(model, l):
            return model.pTlv[l] == (model.BLT[l]) * (model.deltaLT[l])

        self.model.KVL_trans_from = pyo.Constraint(
            self.model.TRANSF, rule=KVL_trans_fromend
        )
        self.model.KVL_trans_to = pyo.Constraint(
            self.model.TRANSF, rule=KVL_trans_toend
        )

        # --- phase angle pyo.Constraints ---
        def phase_angle_diff1(model, l):
            return (
                model.deltaL[l]
                == model.delta[model.A[l, 1]] - model.delta[model.A[l, 2]]
            )

        self.model.phase_diff1 = pyo.Constraint(
            self.model.L, rule=phase_angle_diff1
        )

        # --- phase angle pyo.Constraints ---
        # PowerModels' DC flow is p = -b (θ_from − θ_to): the transformer
        # phase shift does not enter it (nor does the tap).
        use_shift = self.dc_convention != "powermodels"

        def phase_angle_diff2(model, l):
            diff = model.delta[model.AT[l, 1]] - model.delta[model.AT[l, 2]]
            if use_shift:
                diff = diff - model.shift[l]
            return model.deltaLT[l] == diff

        self.model.phase_diff2 = pyo.Constraint(
            self.model.TRANSF, rule=phase_angle_diff2
        )
