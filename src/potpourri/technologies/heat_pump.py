# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Heat pump mix-in.

Attaches heat pump sets, parameters, variables, and constraints to a
multi-period model.
"""

import numpy as np
import pyomo.environ as pyo
from potpourri.technologies.flexibility import Flexibility_multi_period


class Heatpump_multi_period(Flexibility_multi_period):
    """Multi-period heat pump device module with thermal building model.

    Models the electro-thermal dynamics of a building heated by a heat pump:
    each heat pump consumes electrical power ``hp_p`` and raises the indoor
    temperature according to a first-order thermal model.  Heat loss is
    proportional to the load profile.

    Args:
        net: pandapower network with simbench profiles.
        T: Number of time steps (must match the model's time horizon).
        scenario: Predefined penetration scenario 0–3.  Ignored when
            *penetration* is supplied explicitly.

            =========  =================
            scenario   % of non-slack buses
            =========  =================
            0          6.3 %
            1          26.4 %
            2          44.6 %
            3          59.9 %
            =========  =================

        penetration: Percentage of non-slack buses that receive a heat pump
            (0–100).  Overrides *scenario* when given.
        power_max_pu: Maximum electrical input power in per-unit on the
            system base ``net.sn_mva``.
        cop: Coefficient of performance (heat output / electrical input).
            Becomes a Pyomo parameter ``HP_CoP`` indexed over the HP set,
            so it can be updated after model construction.
        temp_max_c: Upper indoor temperature bound in °C (or K offset).
        temp_min_c: Lower indoor temperature bound in °C (or K offset).
        avg_house_size_m3: Average house volume in m³ for the thermal capacity
            calculation.
        wall_thickness_m: Wall thickness in m used for heat-loss scaling.
        qloss_max: Scaling target for heat loss — the maximum heat-loss value
            (in p.u. power) that is mapped to the peak electrical load.
        seed: Seed for the placement draw.  Defaults to
            :data:`~potpourri.technologies.flexibility.DEFAULT_PLACEMENT_SEED`,
            so repeated runs equip the same buses.  Vary it to sample
            placements.
        rng: An existing :class:`numpy.random.Generator`, taking precedence
            over *seed*.  Thread one through a Monte Carlo sweep for
            independent scenarios from a run that replays exactly.

    The thermal capacity of the building is derived from the supplied house
    geometry and material constants.  To override the calculated value, set
    ``self.heat_cap`` after construction but before calling ``get_all``.

    Example::

        hp = Heatpump_multi_period(
            net, T=96, scenario=1,
            cop=3.5,
            temp_max_c=22,
            temp_min_c=18,
        )
        hp.get_all(model)
    """

    SCENARIO_PENETRATION: dict[int, float] = {
        0: 6.3,
        1: 26.4,
        2: 44.6,
        3: 59.9,
    }

    # Physical constants for the default building thermal model
    _AIR_DENSITY_KG_M3: float = 1.2
    _AIR_SPECIFIC_HEAT_J_KGK: float = 1005.0
    _CONCRETE_DENSITY_KG_M3: float = 2400.0
    _CONCRETE_SPECIFIC_HEAT_KJ_KGK: float = 0.84

    def __init__(
        self,
        net,
        T=None,
        scenario=None,
        *,
        penetration: float | None = None,
        power_max_pu: float = 0.005,
        cop: float = 4.0,
        temp_max_c: float = 20.0,
        temp_min_c: float = 15.0,
        avg_house_size_m3: float = 500.0,
        wall_thickness_m: float = 0.2,
        qloss_max: float = 0.01,
        seed=None,
        rng=None,
    ):
        super().__init__(net, T, scenario, seed=seed, rng=rng)

        if penetration is not None:
            self.hp_percentage = float(penetration)
        elif scenario is not None:
            self.hp_percentage = self.SCENARIO_PENETRATION[scenario]
        else:
            raise ValueError(
                "Provide either scenario (0–3) or an explicit penetration "
                "percentage via the penetration= argument."
            )

        # Drawn from this device's own generator, so the placement is
        # reproducible (see Flexibility_multi_period.rng).
        self.random_indexes = self.draw_placement(self.hp_percentage)

        self.temp_max = temp_max_c
        self.temp_min = temp_min_c
        self.hp_power_max = power_max_pu
        self.hp_cop = cop

        # --- Thermal building model ---
        # Wall area ≈ 1/3 of volume (cube approximation)
        wall_area_m2 = avg_house_size_m3 / 3.0
        air_mass_kg = avg_house_size_m3 * self._AIR_DENSITY_KG_M3
        concrete_mass_kg = (
            4.0
            * wall_area_m2
            * wall_thickness_m
            * self._CONCRETE_DENSITY_KG_M3
        )
        # Thermal capacity in MWh/K (convert from J/K via / 3_600_000)
        self.heat_cap = (
            air_mass_kg * self._AIR_SPECIFIC_HEAT_J_KGK
            + concrete_mass_kg * self._CONCRETE_SPECIFIC_HEAT_KJ_KGK * 1000.0
        ) / 3_600_000.0

        # --- Heat-loss profile proportional to load ---
        # Both the reference peak and the profile it scales are per-unit, so
        # heat_load comes out in p.u. power like hp_p and HP_Pmax. Scaling a
        # per-unit factor onto an MW profile made Qloss wrong by baseMVA, which
        # fed straight into the building temperature update.
        self.Qloss_max = qloss_max
        load_profile_pu = net.profiles[("load", "p_mw")][0] / self.baseMVA
        self.max_load = float(np.max(load_profile_pu))
        self.heat_scaling_fac = self.Qloss_max / self.max_load
        self.heat_load = self.heat_scaling_fac * load_profile_pu

    def get_all(self, model):
        """Attach the heat pump and couple it to the balance.

        Attach heat pump sets, parameters, variables, and constraints, and
        couple the electrical power into the nodal balance.
        """
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.get_all_constraints(model)
        self.couple_to_power_balance(model)

    def couple_to_power_balance(self, model):
        """Register the heat pump's electrical demand with the nodal balance.

        ``hp_p`` is the electrical input power, bounded in ``[0, HP_Pmax]``, so
        it is a pure load and enters the balance with a ``+``. Heat pumps have
        no reactive-power variable, so nothing goes on the reactive balance.
        """
        self.register_kcl_real(
            model, self.bus_term(model, model.HP_bus, "hp_p", sign=1.0)
        )

    def get_sets(self, model):
        """Define HP and HP_bus sets from randomly placed heat pumps."""
        super().get_sets(model)
        model.HP = pyo.Set(initialize=list(range(len(self.random_indexes))))
        model.HP_bus = pyo.Set(initialize=list(enumerate(self.random_indexes)))
        return True

    def get_parameters(self, model):
        """Attach the power, temperature, CoP and heat-loss data.

        Attach power limits, temperature bounds, CoP, thermal capacity,
        and heat-loss profile.
        """
        model.HP_Pmax = pyo.Param(
            model.HP, within=pyo.Reals, initialize=self.hp_power_max
        )
        model.HP_Pmin = pyo.Param(model.HP, within=pyo.Reals, initialize=0.0)
        model.TempMax = pyo.Param(
            model.HP, within=pyo.Reals, initialize=self.temp_max
        )
        model.TempMin = pyo.Param(
            model.HP, within=pyo.Reals, initialize=self.temp_min
        )
        # CoP as a mutable per-unit Pyomo parameter so it can be updated
        # after model construction (e.g. for seasonal COP sensitivity).
        model.HP_CoP = pyo.Param(
            model.HP,
            within=pyo.NonNegativeReals,
            initialize=self.hp_cop,
            mutable=True,
        )
        # Thermal capacity of the building (MWh/K) as a Pyomo parameter so
        # it can be adjusted per heat-pump instance if needed.
        model.HP_ThermCap = pyo.Param(
            model.HP,
            within=pyo.NonNegativeReals,
            initialize=self.heat_cap,
            mutable=True,
        )

        self.Qloss_data_dict, self.Qloss_tuple = self.make_to_dict(
            model.HP, model.T, self.heat_load, True
        )
        model.Qloss = pyo.Param(
            self.Qloss_tuple, initialize=self.Qloss_data_dict
        )
        return True

    def get_variables(self, model):
        """Attach the electrical-power and indoor-temperature variables.

        Create hp_p (electrical power) and temp (indoor temperature)
        variables.

        ``temp`` starts mid-band rather than at Pyomo's default of 0, which
        sits outside ``[temp_min_c, temp_max_c]`` and so begins the solve on
        an infeasible point. ``hp_p`` starts at 0, inside its own bounds.
        """
        model.hp_p = pyo.Var(
            model.HP, model.T, domain=pyo.Reals, initialize=0.0
        )
        model.temp = pyo.Var(
            model.HP,
            model.T,
            domain=pyo.Reals,
            initialize=0.5 * (self.temp_min + self.temp_max),
        )
        return True

    def get_all_constraints(self, model):
        """Add the power, temperature and thermal-update constraints."""

        @model.Constraint(model.HP, model.T)
        def hp_power_con(model, h, t):
            """Bound the heat pump's electrical power.

            Args:
                h: Heat-pump index.
                model: The Pyomo model being built.
                t: Time index.

            Returns:
                A Pyomo expression.
            """
            return model.HP_Pmin[h], model.hp_p[h, t], model.HP_Pmax[h]

        @model.Constraint(model.HP, model.T)
        def hp_temp_con(model, h, t):
            """Keep the indoor temperature inside its band.

            Args:
                h: Heat-pump index.
                model: The Pyomo model being built.
                t: Time index.

            Returns:
                A Pyomo expression.
            """
            return model.TempMin[h], model.temp[h, t], model.TempMax[h]

        @model.Constraint(model.HP, model.T)
        def hp_temp_update_con(model, h, t):
            """Carry the indoor temperature from one step to the next.

            The thermal analogue of a storage balance: the building's heat
            capacity couples the time steps, so the heat pump can shift
            consumption the way a battery shifts energy.

            Args:
                h: Heat-pump index.
                model: The Pyomo model being built.
                t: Time index.

            Returns:
                A Pyomo expression.
            """
            if t == model.T.at(1):
                return pyo.Constraint.Skip
            return (
                model.temp[h, t]
                == model.temp[h, t - 1]
                + model.deltaT
                * (model.hp_p[h, t] * model.HP_CoP[h] - model.Qloss[h, t])
                / model.HP_ThermCap[h]
            )

        return True

    def get_all_opf(self, model):
        """No additional OPF components needed for heat pumps."""

    def get_all_ac(self, model):
        """No additional AC components needed for heat pumps."""

    def get_all_acopf(self, model):
        """No additional ACOPF components needed for heat pumps."""
