"""PV mix-in: attaches PV generation sets, parameters, variables, and
power-bound constraints to a multi-period model."""

import numpy as np
import pyomo.environ as pyo
from potpourri.technologies.flexibility import Flexibility_multi_period


class PV_multi_period(Flexibility_multi_period):
    """Multi-period PV device module with scenario-based penetration levels.

    PV units are placed randomly at a fraction of non-slack buses.  Their
    available generation follows a time-varying upper bound taken from the
    SimBench renewables profile; curtailment is modelled implicitly (the
    optimiser may dispatch below the available potential).

    Penetration percentages are based on:
    *Anforderungen an aktuelle Verteilnetze und deren zukuenftige
    Versorgungsaufgabe*.

    Args:
        net: pandapower network with simbench profiles.
        T: Number of time steps (must match the model's time horizon).
        scenario: Predefined penetration scenario 0–3.  Ignored when
            *penetration* is supplied.

            =========  =================
            scenario   % of non-slack buses
            =========  =================
            0          13.4 %
            1          22.4 %
            2          24.4 %
            3          25.4 %
            =========  =================

        penetration: Percentage of non-slack buses that receive a PV unit
            (0–100).  Overrides *scenario* when given.
        profile_column: Column name from ``net.pv_load_profiles`` (the
            SimBench renewables table) to use as the PV generation profile.
            Defaults to ``"PV5"``.
        pv_pmin: Minimum PV output in per-unit (curtailment lower bound).
            Defaults to ``0.0`` (full curtailment allowed).

    Example::

        pv = PV_multi_period(net, T=96, scenario=0, profile_column="PV3")
        pv.get_all(model)
    """

    SCENARIO_PENETRATION: dict[int, float] = {
        0: 13.4,
        1: 22.4,
        2: 24.4,
        3: 25.4,
    }

    def __init__(
        self,
        net,
        T=None,
        scenario=None,
        *,
        penetration: float | None = None,
        profile_column: str = "PV5",
        pv_pmin: float = 0.0,
    ):
        super().__init__(net, T, scenario)
        self.net = net

        if penetration is not None:
            self.pv_percentage = float(penetration)
        elif scenario is not None:
            self.pv_percentage = self.SCENARIO_PENETRATION[scenario]
        else:
            raise ValueError(
                "Provide either scenario (0–3) or an explicit penetration "
                "percentage via the penetration= argument."
            )

        # PV profiles are stored as positive generation values in SimBench;
        # negate so that injection into the network is positive in the model.
        pv_profiles = -self.net.pv_load_profiles
        if profile_column not in pv_profiles.columns:
            available = list(pv_profiles.columns)
            raise ValueError(
                f"profile_column '{profile_column}' not found in "
                f"net.pv_load_profiles. Available columns: {available}"
            )
        self.pv_load_profile = pv_profiles[profile_column]

        self.pv_pmax = self.pv_load_profile
        self.pv_pmin = pv_pmin

        num_indexes = round(
            len(self.buses_excl_extGrids) * self.pv_percentage / 100
        )
        self.random_indexes = np.random.choice(
            self.buses_excl_extGrids, num_indexes, replace=False
        )

    def get_all(self, model):
        """Attach PV sets, parameters, variables, constraints, and unfix
        variables."""
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.get_all_constraints(model)
        self.unfix_variables(model)

    def unfix_variables(self, model):
        """Unfix pPV variables for all PV units over the full time horizon."""
        for t in model.T:
            for pv in model.PV:
                model.pPV[pv, t].unfix()

    def get_sets(self, model):
        """Define PV and PV_bus sets from randomly placed PV units."""
        super().get_sets(model)
        model.PV = pyo.Set(initialize=list(range(len(self.random_indexes))))
        model.PV_bus = pyo.Set(initialize=list(enumerate(self.random_indexes)))
        return True

    def get_parameters(self, model):
        """Attach time-indexed PV_Pmax and PV_Pmin parameters."""
        self.PV_Pmax_dict, self.PV_Pmax_tuple = self.make_to_dict(
            model.PV, model.T, self.pv_pmax
        )
        self.PV_Pmin_dict, self.PV_Pmin_tuple = self.make_to_dict(
            model.PV, model.T, self.pv_pmin, False
        )

        model.PV_Pmax = pyo.Param(
            self.PV_Pmax_tuple, within=pyo.Reals, initialize=self.PV_Pmax_dict
        )
        model.PV_Pmin = pyo.Param(
            self.PV_Pmin_tuple, within=pyo.Reals, initialize=self.PV_Pmin_dict
        )

    def get_variables(self, model):
        """Create pPV variable initialised from the load profile."""
        self.pPV_data_dict, self.pPV_tuple = self.make_to_dict(
            model.PV, model.T, self.pv_load_profile
        )
        model.pPV = pyo.Var(
            self.pPV_tuple, within=pyo.Reals, initialize=self.pPV_data_dict
        )

    def get_all_constraints(self, model):
        """Add real-power bound constraints for all PV units over all time
        steps."""

        @model.Constraint(model.PV, model.T)
        def PV_real_power_bounds(model, pv, t):
            return model.PV_Pmax[pv, t], model.pPV[pv, t], model.PV_Pmin[pv, t]

    def get_all_acopf(self, model):
        """No additional ACOPF components needed for PV."""

    def get_all_ac(self, model):
        """No additional AC components needed for PV."""

    def get_all_opf(self, model):
        """No additional OPF components needed for PV."""
