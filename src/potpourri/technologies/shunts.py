# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Shunt mix-in: shunt conductance over a horizon."""

import pyomo.environ as pyo
from potpourri.technologies.flexibility import Flexibility_multi_period


class Shunts_multi_period(Flexibility_multi_period):
    """Multi-period shunts, read from the network.

    Multi-period shunt device module; reads shunt data from the pandapower
    network.
    """

    def __init__(self, net, T=None, scenario=None):
        super().__init__(net, T, scenario)

        self.shunt_set = self.net.shunt.index[self.net.shunt.in_service]
        self.bus_shunt_set = list(
            zip(
                self.bus_lookup[self.net.shunt.bus[self.shunt_set].values],
                self.shunt_set,
            )
        )

    def get_all(self, model):
        """Attach shunt sets and parameters to the model."""
        self.get_sets(model)
        self.get_parameters(model)

    def get_all_opf(self, model):
        """No-op placeholder for the OPF mix-in interface."""

    def get_all_acopf(self, model):
        """Attach the AC OPF layer for this device.

        Adds the reactive-power bounds on top of the active-power ones.

        Returns:
                None. The components are added to `model` in place.
        """
        pass

    def get_sets(self, model):
        """Define SHUNT and SHUNTbs sets from in-service shunts."""
        super().get_sets(model)
        model.SHUNT = pyo.Set(initialize=self.shunt_set)  # set of shunts
        model.SHUNTbs = pyo.Set(
            within=model.B * model.SHUNT, initialize=self.bus_shunt_set
        )  # set of shunt-bus mapping
        return True

    def get_parameters(self, model):
        """Attach the shunt conductance parameter GB, indexed by (shunt, step).

        Time-indexed to match ``BB``, the susceptance, which
        ``AC_multi_period`` has always declared over ``SHUNT × T``. The two
        disagreed: ``GB`` was declared over ``SHUNT`` alone while the result
        mapper read ``GB[s, t]``, and the power-balance rules read ``BB[s]``.
        Each consumer therefore used the wrong form for one of the two, and any
        network carrying a shunt raised ``KeyError`` on model construction.
        Every SimBench network has no shunts, so ``for s in model.SHUNT`` never
        iterated and neither path was ever reached.

        The conductance is constant over the horizon today; the index is there
        so a switched shunt bank can vary with time without another mismatch.
        """
        self.GB_data_dict, self.GB_tuple = self.make_to_dict(
            model.SHUNT, model.T, self.GB_data[model.SHUNT], False
        )
        model.GB = pyo.Param(
            self.GB_tuple, within=pyo.Reals, initialize=self.GB_data_dict
        )  # shunt conductance, per (shunt, time step)
        return True
