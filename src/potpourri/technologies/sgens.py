"""Static generator (sgen) mix-in: attaches sgen profiles and OPF limits to
a multi-period model."""

import numpy as np
import pyomo.environ as pyo
from potpourri.technologies.flexibility import Flexibility_multi_period


class Sgens_multi_period(Flexibility_multi_period):
    """Multi-period static generator device module; reads sgen profiles from
    net.profiles."""

    def __init__(self, net, T=None, scenario=None):
        super().__init__(net, T, scenario)

        # TODO: check whether q values are provided in simbench networks
        sgen_bus = self.bus_lookup[self.net.sgen.bus.values]
        self.static_generation_data = {
            "p": self.net.profiles[("sgen", "p_mw")],
            "q": self.net.profiles[("sgen", "q_mvar")],
            "in_service": self.net.sgen.in_service.values,
            "bus": sgen_bus,
        }
        self.static_generation_data["gen_bus"] = list(
            enumerate(self.static_generation_data["bus"])
        )

    def get_all(self, model):
        """Attach sets, parameters, and variables; fix all sgen variables to
        profile values."""
        self.get_sets(model)
        self.get_parameters(model)
        self.get_variables(model)
        self.fix_variables(model)

    def get_all_opf(self, model):
        """Attach OPF sets, parameters, and real-power bound constraints for
        controllable sgens."""
        self.get_opf_sets(model)
        self.get_opf_parameters(model)
        self.get_all_Constraints_opf(model)

    def get_all_acopf(self, model):
        self.get_all_Constraints_acopf(model)

    def get_sets(self, model):
        super().get_sets(model)
        # list, no set, because list is ordered data source, set is not
        self.sgens_in_service_list = np.where(
            self.static_generation_data["in_service"]
        )[0].tolist()
        model.sG = pyo.Set(
            initialize=self.sgens_in_service_list
        )  # static generators
        model.sGbs = pyo.Set(
            within=model.sG * model.B,
            initialize=self.static_generation_data["gen_bus"],
        )  # set of static generator-bus mapping
        return True

    def get_opf_sets(self, model):
        # list, no set, because list is ordered data source, set is not
        self.sgens_controllable_list = np.where(
            self.static_generation_data["controllable"]
        )[0].tolist()

        model.sGc = pyo.Set(
            within=model.sG,
            initialize=[
                g
                for g in self.sgens_controllable_list
                if g in self.sgens_in_service_list
            ],
        )  # static generators that are controllable and in service

    def get_parameters(self, model):
        self.PsG_data_dict, self.PsG_tuple = self.make_to_dict(
            model.sG, model.T, self.static_generation_data["p"]
        )
        self.QsG_data_dict, self.QsG_tuple = self.make_to_dict(
            model.sG, model.T, self.static_generation_data["q"]
        )
        # --- Parameters ---
        model.PsG = pyo.Param(self.PsG_tuple, initialize=self.PsG_data_dict)
        # reactive generation
        model.QsG = pyo.Param(self.QsG_tuple, initialize=self.QsG_data_dict)
        return True

    def get_opf_parameters(self, model):
        # static generation real power limits
        model.sPGmax = pyo.Param(
            self.PsGmax_tuple, initialize=self.PsGmax_data_dict
        )
        model.sPGmin = pyo.Param(
            self.PsGmin_tuple, initialize=self.PsGmin_data_dict
        )

        # Only for AC OPF, maybe need to be moved into other function?
        # static generation reactive power limits
        # static generation reactive power limits
        model.QsGmax = pyo.Param(
            self.QsGmax_tuple,
            within=pyo.Reals,
            initialize=self.QsGmax_data_dict,
            mutable=True,
        )
        model.QsGmin = pyo.Param(
            self.QsGmin_tuple,
            within=pyo.Reals,
            initialize=self.QsGmin_data_dict,
            mutable=True,
        )

    def get_all_Constraints_opf(self, model):
        # psG Constraint
        @model.Constraint(model.sGc, model.T)
        def static_generation_real_power_bounds(model, g, t):
            model.psG[(g, t)].unfix()
            return (
                model.sPGmin[(g, t)],
                model.psG[(g, t)],
                model.sPGmax[(g, t)],
            )

    def get_variables(self, model):
        # --- Variables ---
        model.psG = pyo.Var(
            self.PsG_tuple, domain=pyo.NonNegativeReals
        )  # real static generator power
        model.qsG = pyo.Var(
            self.QsG_tuple, domain=pyo.Reals
        )  # reactive power of static generators
        return True

    def unfix_variables(self, model):
        # unfix the static generation values
        for g in model.sG:
            for t in model.T:
                model.psG[(g, t)].unfix()
        return True

    def fix_variables(self, model):
        # fix the static generation values
        for g in model.sG:
            for t in model.T:
                model.psG[(g, t)].fix(model.PsG[(g, t)])

        for g in model.sG:
            for t in model.T:
                model.qsG[(g, t)].fix(model.QsG[(g, t)])

    #     return True

    def static_generation_real_power_limits(self, model):
        if "controllable" in self.net.sgen:
            self.static_generation_data["controllable"] = (
                self.net.sgen.controllable.values
            )
        else:
            self.static_generation_data["controllable"] = np.full(
                len(self.net.sgen), False
            )

        self.PsGmax_data_dict, self.PsGmax_tuple = self.make_to_dict(
            model.sG, model.T, self.static_generation_data["p"]
        )
        self.PsGmin_data_dict, self.PsGmin_tuple = self.make_to_dict(
            model.sG, model.T, 0, False
        )

    def static_generation_reactive_power_limits(self, model):
        if "controllable" in self.net.sgen:
            self.static_generation_data["controllable"] = (
                self.net.sgen.controllable.values
            )
        else:
            self.static_generation_data["controllable"] = np.full(
                len(self.net.sgen), False
            )

        self.QsGmax_data_dict, self.QsGmax_tuple = self.make_to_dict(
            model.sG, model.T, abs(self.static_generation_data["q"])
        )
        self.QsGmin_data_dict, self.QsGmin_tuple = self.make_to_dict(
            model.sG, model.T, -abs(self.static_generation_data["q"])
        )
        # self.static_generation_wind_var_q( self.net)
        self.static_generation_data["type"] = self.net.sgen.type.values

    def get_all_Constraints_acopf(self, model):

        # QsG_Constraint
        @model.Constraint(model.sGc, model.T)
        def static_generation_reactive_power_bounds(model, g, t):
            model.qsG[(g, t)].unfix()
            return (
                model.QsGmin[(g, t)],
                model.qsG[(g, t)],
                model.QsGmax[(g, t)],
            )
