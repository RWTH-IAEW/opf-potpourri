"""Battery mix-in: attaches battery storage sets, parameters, variables,
and constraints to a multi-period model."""

import numpy as np
import pyomo.environ as pyo
from potpourri.technologies.flexibility import Flexibility_multi_period


class Battery_multi_period(Flexibility_multi_period):
    """Multi-period battery storage device module.

    Randomly places batteries at a fraction of non-slack buses and attaches
    the corresponding Pyomo Sets, Parameters, Variables, and Constraints to
    an existing multi-period model.

    Args:
        net: pandapower network with simbench profiles (already passed to
            ``ACOPF_multi_period``).
        T: Number of time steps (must match the model's time horizon).
        scenario: Predefined penetration scenario 0–3.  Ignored when
            *penetration* is supplied explicitly.

            =========  =================
            scenario   % of non-slack buses
            =========  =================
            0          1.0 %
            1          7.9 %
            2          9.9 %
            3          10.6 %
            =========  =================

        penetration: Percentage of non-slack buses that receive a battery
            (0–100).  Overrides *scenario* when given.
        power_pu: Symmetric charge/discharge power limit in per-unit on the
            system base ``net.sn_mva``.
        soc_max: Maximum state of charge (p.u. of capacity, 0–1).
        soc_min: Minimum state of charge (p.u. of capacity, 0–1).
        capacity_pu_h: Battery energy capacity in per-unit power · hours
            (i.e. ``capacity_MWh / net.sn_mva``).
        efficiency: One-way charge/discharge efficiency (0–1).
        initial_soc_fraction: Initial SOC expressed as a fraction of
            *soc_max* (0–1).  Default 0.5 = 50 % of maximum capacity.

    Example::

        battery = Battery_multi_period(
            net, T=96, scenario=1
        )
        # or with explicit parameters:
        battery = Battery_multi_period(
            net, T=96,
            penetration=15.0,   # 15 % of buses
            power_pu=0.01,
            capacity_pu_h=0.025,
            efficiency=0.95,
        )
        battery.get_all(model)
    """

    # Default penetration levels per scenario (percentage of non-slack buses)
    SCENARIO_PENETRATION: dict[int, float] = {0: 1.0, 1: 7.9, 2: 9.9, 3: 10.6}

    def __init__(
        self,
        net,
        T=None,
        scenario=None,
        *,
        penetration: float | None = None,
        power_pu: float = 0.006,
        soc_max: float = 1.0,
        soc_min: float = 0.2,
        capacity_pu_h: float = 0.015,
        efficiency: float = 0.9,
        initial_soc_fraction: float = 0.5,
    ):
        super().__init__(net, T, scenario)

        if penetration is not None:
            self.bat_percentage = float(penetration)
        elif scenario is not None:
            self.bat_percentage = self.SCENARIO_PENETRATION[scenario]
        else:
            raise ValueError(
                "Provide either scenario (0–3) or an explicit penetration "
                "percentage via the penetration= argument."
            )

        # randomly select buses for battery placement
        num_indexes = round(
            len(self.buses_excl_extGrids) * self.bat_percentage / 100
        )
        self.random_indexes = np.random.choice(
            self.buses_excl_extGrids, num_indexes, replace=False
        )

        self.bat_power = power_pu
        self.bat_soc_max = soc_max
        self.bat_soc_min = soc_min
        self.bat_cap = capacity_pu_h
        self.bat_efficiency = efficiency
        self.bat_initial_soc_fraction = initial_soc_fraction

    def get_all(self, model):
        """Attach battery sets, parameters, variables, and constraints to the
        model."""
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.get_all_constraints(model)

    def get_sets(self, model):
        """Define BAT and BAT_bus sets from randomly placed batteries."""
        super().get_sets(model)
        model.BAT = pyo.Set(initialize=list(range(len(self.random_indexes))))
        model.BAT_bus = pyo.Set(
            initialize=list(enumerate(self.random_indexes))
        )
        return True

    def get_parameters(self, model):
        """Attach battery power, SOC bounds, capacity, and efficiency
        parameters."""
        model.BAT_Pmax = pyo.Param(
            model.BAT, within=pyo.Reals, initialize=self.bat_power
        )
        model.BAT_Pmin = pyo.Param(
            model.BAT, within=pyo.Reals, initialize=-self.bat_power
        )
        model.BAT_SOCmax = pyo.Param(
            model.BAT, within=pyo.Reals, initialize=self.bat_soc_max
        )
        model.BAT_SOCmin = pyo.Param(
            model.BAT, within=pyo.Reals, initialize=self.bat_soc_min
        )
        model.BAT_Cap = pyo.Param(
            model.BAT, within=pyo.Reals, initialize=self.bat_cap
        )
        model.BAT_Eff = pyo.Param(
            model.BAT, within=pyo.Reals, initialize=self.bat_efficiency
        )
        model.BAT_SOC_init = pyo.Param(
            model.BAT,
            within=pyo.Reals,
            initialize={
                b: self.bat_initial_soc_fraction * self.bat_soc_max
                for b in range(len(self.random_indexes))
            },
        )
        return True

    def get_variables(self, model):
        """Create BAT_P (power) and BAT_SOC (state of charge) variables."""
        model.BAT_P = pyo.Var(model.BAT, model.T, within=pyo.Reals)
        model.BAT_SOC = pyo.Var(model.BAT, model.T, within=pyo.Reals)
        return True

    def get_all_constraints(self, model):
        """Add power-bound, SOC-bound, and SOC-update constraints."""

        def bat_power_rule(model, b, t):
            return model.BAT_Pmin[b], model.BAT_P[b, t], model.BAT_Pmax[b]

        model.bat_power_con = pyo.Constraint(
            model.BAT, model.T, rule=bat_power_rule
        )

        def bat_soc_rule(model, b, t):
            if t == model.T.at(1):
                return model.BAT_SOC[b, t] == model.BAT_SOC_init[b]
            return (
                model.BAT_SOCmin[b],
                model.BAT_SOC[b, t],
                model.BAT_SOCmax[b],
            )

        model.bat_soc_con = pyo.Constraint(
            model.BAT, model.T, rule=bat_soc_rule
        )

        def bat_soc_update_rule(model, b, t):
            if t == model.T.at(1):
                return pyo.Constraint.Skip
            return (
                model.BAT_SOC[b, t]
                == model.BAT_SOC[b, t - 1]
                + model.deltaT
                * (model.BAT_P[b, t] * model.BAT_Eff[b])
                / model.BAT_Cap[b]
            )

        model.bat_soc_update_con = pyo.Constraint(
            model.BAT, model.T, rule=bat_soc_update_rule
        )
        return True

    def get_all_acopf(self, model):
        """No additional ACOPF components needed for batteries."""

    def get_all_ac(self, model):
        """No additional AC components needed for batteries."""

    def get_all_opf(self, model):
        """No additional OPF components needed for batteries."""
