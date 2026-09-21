# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Base mix-in class for all multi-period flexibility device modules
(batteries, heat pumps, PV, demand, etc.)."""

import pandas as pd
import pyomo.environ as pyo
from math import pi
import numpy as np

#: Seed used for device placement when the caller supplies neither *seed* nor
#: *rng*. Fixed rather than random so that the shipped examples, the test suite
#: and any study reproduce out of the box; matches the default in
#: ``research/lin_opf/scenarios.py``.
DEFAULT_PLACEMENT_SEED = 42


class Flexibility_multi_period:
    """Base class for technology mix-in objects that attach Pyomo components
    to a multi-period model.

    Reads network topology and profile data from *net* in ``__init__``.
    Subclasses implement ``get_all()``, ``get_sets()``, ``get_parameters()``,
    ``get_variables()``, and constraint methods.

    Args:
        net: pandapower network the device reads its data from.
        T: Number of time steps, where the subclass needs it.
        scenario: Predefined penetration scenario, where the subclass uses one.
        seed: Seed for this device's placement draw. Defaults to
            :data:`DEFAULT_PLACEMENT_SEED`.
        rng: An existing :class:`numpy.random.Generator` to draw from, which
            takes precedence over *seed*. Pass one generator through a whole
            Monte Carlo sweep to get independent scenarios from a run that
            replays exactly.

    Placement draws from ``self.rng``, never from the global NumPy state, so a
    device's buses depend only on its own seed and not on whatever else in the
    process happened to draw a random number first.
    """

    def __init__(self, net, T=None, scenario=None, *, seed=None, rng=None):
        self.rng = (
            rng
            if rng is not None
            else np.random.default_rng(
                DEFAULT_PLACEMENT_SEED if seed is None else seed
            )
        )
        self.net = net

        # buses that are not ext_grids (eligible for technology placement)
        self.buses_excl_extGrids = self.net.bus.loc[
            ~self.net.bus.index.isin(net.ext_grid.bus)
        ].index

        bus_set = self.net._ppc["bus"][:, [0, 1, 7, 8]]
        bus_set[:, -1] *= pi / 180
        self.bus_data = pd.DataFrame(
            bus_set[:, 1:],
            index=bus_set[:, 0].astype(int),
            columns=["type", "v_m", "v_a_rad"],
        )
        self.bus_lookup = self.net._pd2ppc_lookups["bus"]
        # ppc bus number carrying each pandapower bus, and the subset of ppc
        # buses that a pandapower bus maps onto.  Auxiliary ppc buses inserted
        # for node-node switches are absent from the latter: they have no
        # pandapower row, so no per-bus user data (voltage limits) exists for
        # them.  Consumed by get_sets() to build the Bpd set.
        self.pd_bus_to_ppc = self.bus_lookup[self.net.bus.index.values]
        self.ppc_buses_with_pd = pd.Index(
            sorted({int(b) for b in self.pd_bus_to_ppc})
        )

        self.baseMVA = self.net.sn_mva

        self.PD_data = self.net.profiles[("load", "p_mw")] / self.baseMVA
        self.QD_data = self.net.profiles[("load", "q_mvar")] / self.baseMVA

        self.GB_data = self.net.shunt.p_mw * self.net.shunt.step / self.baseMVA

    def draw_placement(self, percentage):
        """Draw the buses that receive a device, from this device's generator.

        Args:
            percentage: Share of eligible (non-slack) buses to equip, 0–100.

        Returns:
            ndarray of pandapower bus indices, drawn without replacement.

        Drawing from ``self.rng`` rather than ``np.random`` is what makes a
        study reproducible: the global generator is shared with every other
        library in the process, so an unseeded draw gave a different
        placement — and different results — on every run.
        """
        count = round(len(self.buses_excl_extGrids) * percentage / 100)
        return self.rng.choice(self.buses_excl_extGrids, count, replace=False)

    # --- power-balance coupling ------------------------------------------
    #
    # A device receives only the Pyomo model, never the model wrapper, so the
    # registry of coupling terms lives on the model. The power-flow classes
    # read it back through KCL_flexibility(). Terms are plain callables
    # ``(model, b, t) -> expression`` so they are re-evaluated whenever the
    # balance is rebuilt, which is what lets a device attach after the
    # power-flow equations were first constructed.
    _KCL_REAL_ATTR = "_kcl_real_terms"
    _KCL_REACTIVE_ATTR = "_kcl_reactive_terms"

    @classmethod
    def kcl_terms(cls, model, reactive=False):
        """Return the model's list of registered coupling terms."""
        attr = cls._KCL_REACTIVE_ATTR if reactive else cls._KCL_REAL_ATTR
        if not hasattr(model, attr):
            setattr(model, attr, [])
        return getattr(model, attr)

    def _claim_registration(self, model, kind):
        """Refuse a second registration of this device on the same model.

        Registering twice would add the device's injection to the balance
        twice — a silent doubling of its power, which no constraint would
        catch.
        """
        attr = "_kcl_registered"
        if not hasattr(model, attr):
            setattr(model, attr, set())
        # Key on the device object, not id(self): a recycled id could otherwise
        # make a fresh device look already-registered.
        key = (self, kind)
        registered = getattr(model, attr)
        if key in registered:
            raise RuntimeError(
                f"{type(self).__name__} is already coupled to this model's "
                f"{kind} power balance. Registering again would double-count "
                f"its power. Call get_all(model) once per device."
            )
        registered.add(key)

    def register_kcl_real(self, model, term):
        """Register a real-power contribution to the nodal balance.

        Args:
            term: ``(model, b, t) -> expression`` giving this device's net
                real power at ppc bus *b* and time *t*, in the **load sign
                convention**: positive is consumption, negative is injection.
                The nodal balance places it on the consumption side, so this
                matches ``pD`` rather than ``psG``.

        Raises:
            RuntimeError: If this device already registered a real-power term
                on this model.
        """
        self._claim_registration(model, "real")
        self.kcl_terms(model).append(term)

    def register_kcl_reactive(self, model, term):
        """Register a reactive-power contribution to the nodal balance.

        Same signature and sign convention as :meth:`register_kcl_real`:
        positive is consumption (inductive), negative is injection.

        Raises:
            RuntimeError: If this device already registered a reactive-power
                term on this model.
        """
        self._claim_registration(model, "reactive")
        self.kcl_terms(model, reactive=True).append(term)

    def bus_term(self, model, buses, var, sign=1.0):
        """Build a coupling term summing *var* over the devices at each bus.

        Args:
            model: The Pyomo model, used to resolve the device-bus set.
            buses: Mapping of device index to **pandapower** bus index, as the
                device's ``*_bus`` set holds it.
            var: Name of the model component indexed by ``(device, t)``.
            sign: ``+1`` if the variable is already in the load convention,
                ``-1`` if it is a generator-convention injection.

        The device placement sets hold pandapower bus indices, while the
        balance is indexed over **ppc** bus numbers — the two differ on grids
        where pandapower inserts auxiliary nodes for node-node switches. The
        mapping goes through ``bus_lookup``, the same way the sgen and load
        bus sets are built.
        """
        # Iterate the members, not dict(buses): dict() on a scalar Pyomo Set
        # yields {None: <the set>} rather than its (device, bus) pairs.
        by_ppc_bus = {}
        for device, pd_bus in list(buses):
            ppc_bus = int(self.bus_lookup[int(pd_bus)])
            by_ppc_bus.setdefault(ppc_bus, []).append(device)

        def term(model, b, t):
            devices = by_ppc_bus.get(b)
            if not devices:
                return 0
            component = getattr(model, var)
            return sign * sum(component[d, t] for d in devices)

        return term

    def get_sets(self, model):
        """Initialise (or re-initialise) the bus set B on the Pyomo model
        from network topology.

        Also defines ``Bpd``, the buses that a pandapower bus maps onto —
        everything in ``B`` except the auxiliary nodes pandapower inserts for
        node-node switches.  Per-bus user data (voltage limits) exists only
        for those, so constraints derived from ``net.bus`` are indexed over
        ``Bpd``.  On grids without auxiliary nodes ``Bpd == B``.
        """
        if hasattr(model, "B"):
            model.del_component(model.B)
        model.B = pyo.Set(initialize=self.bus_data.index)
        if hasattr(model, "Bpd"):
            model.del_component(model.Bpd)
        model.Bpd = pyo.Set(
            within=model.B,
            initialize=getattr(self, "ppc_buses_with_pd", self.bus_data.index),
        )
        return True

    def make_to_dict(self, model_obj, model_time, data, time_dependent=True):
        """Convert data into a ``{(object_index, time_index): value}`` dict
        and a matching list of index tuples for Pyomo Param initialisation.

        Args:
            model_obj: Iterable of object indices (e.g. a Pyomo Set).
            model_time: Iterable of time indices (e.g. ``model.T``).
            data: Data values.  Accepted types:

                * ``np.ndarray`` — indexed by object index.
                * ``pd.Series`` — indexed by time (``time_dependent=True``)
                  or by object index (``time_dependent=False``).
                * ``dict`` — ``{obj: value}`` (time-independent) or
                  ``{obj: {t: value}}`` (time-dependent).
                * ``list`` — same layout as ndarray.
                * ``float`` / ``int`` — scalar, broadcast to all (obj, t).
                * ``pd.DataFrame`` — calls ``.to_dict()`` internally.

            time_dependent: When ``True`` (default) the same value is
                repeated for every time step for each object.  When
                ``False`` the data varies only over objects (not time).

        Returns:
            tuple: ``(data_dict, tuple_list)`` where *data_dict* maps
            ``(obj, t)`` → value and *tuple_list* is the ordered list of
            index pairs.
        """
        tuple_list = [(o, t) for o in model_obj for t in model_time]

        if time_dependent:
            if isinstance(data, np.ndarray):
                data_dict = {(o, t): data[o] for o, t in tuple_list}
            elif isinstance(data, dict):
                # Support both flat {obj: val} and nested {obj: {t: val}}
                first_val = next(iter(data.values())) if data else None
                if isinstance(first_val, dict):
                    data_dict = {(o, t): data[o][t] for o, t in tuple_list}
                else:
                    data_dict = {(o, t): data[o] for o, t in tuple_list}
            elif isinstance(data, pd.Series):
                data_dict = {(o, t): data[t] for o, t in tuple_list}
            elif isinstance(data, list):
                data_dict = {(o, t): data[t] for o, t in tuple_list}
            elif isinstance(data, (int, float)):
                data_dict = {(o, t): data for o, t in tuple_list}
            else:
                raw = data.to_dict()
                data_dict = {(o, t): raw[o][t] for o, t in tuple_list}
        else:
            if isinstance(data, (int, float)) and data == 0:
                data_dict = {(o, t): 0 for o, t in tuple_list}
            elif isinstance(data, np.ndarray):
                data_dict = {(o, t): data[o] for o, t in tuple_list}
            elif isinstance(data, pd.Series):
                data_dict = {(o, t): data[o] for o, t in tuple_list}
            elif isinstance(data, dict):
                data_dict = {(o, t): data[o] for o, t in tuple_list}
            elif isinstance(data, list):
                data_dict = {(o, t): data[o] for o, t in tuple_list}
            elif isinstance(data, (int, float)):
                data_dict = {(o, t): data for o, t in tuple_list}
            else:
                raw = data.to_dict()
                data_dict = {(o, t): raw[o] for o, t in tuple_list}

        return data_dict, tuple_list
