"""Shunt mix-in: attaches shunt conductance sets and parameters to a multi-period model."""

from pyomo.environ import *
from src.potpourri.technologies.flexibility import Flexibility_multi_period


class Shunts_multi_period(Flexibility_multi_period):
    """Multi-period shunt device module; reads shunt data from the pandapower network."""

    def __init__(self, net, T=None, scenario=None):
        super().__init__(net, T, scenario)

        self.shunt_set = self.net.shunt.index[self.net.shunt.in_service]
        self.bus_shunt_set = list(zip(self.bus_lookup[self.net.shunt.bus[self.shunt_set].values], self.shunt_set))

    def get_all(self, model):
        """Attach shunt sets and parameters to the model."""
        self.get_sets(model)
        self.get_parameters(model)

    def get_all_opf(self, model):
        """No-op placeholder for the OPF mix-in interface."""

    def get_all_acopf(self, model):
        pass

    def get_sets(self, model):
        """Define SHUNT and SHUNTbs sets from in-service shunts."""
        super().get_sets(model)
        model.SHUNT = Set(initialize=self.shunt_set)  # set of shunts
        model.SHUNTbs = Set(within=model.B * model.SHUNT,
                            initialize=self.bus_shunt_set)  # set of shunt-bus mapping
        return True

    def get_parameters(self, model):
        """Attach shunt conductance parameter GB."""
        model.GB = Param(model.SHUNT, within=Reals,
                         initialize=self.GB_data[model.SHUNT])  # shunt conductance
        return True
