# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

r"""Bus-to-element incidence lists for the nodal balance constructors.

The nodal balance of every power-flow formulation sums, at each bus, the
generators, loads, branch ends, storage units and shunts attached to it.
Written as a scan of every element for every bus that is
$O(|B| \cdot |E|)$ Pyomo set and parameter look-ups; the lists built here
are computed **once per balance construction** and make it
$O(|B| + |E|)$.

The lists are derived from the model's own sets and parameters (`Gbs`,
`sGbs`, `Dbs`, `SHUNTbs`, `A`, `AT`, `STOR_bus`) at the moment
[`bus_incidence`][potpourri.models.incidence.bus_incidence] is called, so a
rebuild after a component was added still sees the current state. Within a
bus every list follows the **iteration order of the element's own set**
(`model.G`, `model.L`, ...), which is the order the former scans visited
the elements in; the expressions built from them are therefore term for
term the ones the scans produced.
"""

from contextlib import contextmanager
from dataclasses import dataclass, field


@dataclass
class BusIncidence:
    """Elements attached to each bus, keyed by ppc bus number.

    Every attribute maps a bus to a list of element indices in the iteration
    order of the element's set. A bus with no such element is absent, so
    read the lists with `.get(b, ())`.

    Attributes:
        sG: Static generators of `model.sG` at the bus (via `sGbs`).
        G: Generators and external grids of `model.G` (via `Gbs`).
        STOR: Storage units of `model.STOR` (via `STOR_bus`).
        D: Loads of `model.D` (via `Dbs`).
        L_from: Lines of `model.L` whose from end (`A[l, 1]`) is the bus.
        L_to: Lines of `model.L` whose to end (`A[l, 2]`) is the bus.
        T_hv: Transformers of `model.TRANSF` with `AT[l, 1]` at the bus.
        T_lv: Transformers of `model.TRANSF` with `AT[l, 2]` at the bus.
        SHUNT: Shunts of `model.SHUNT` at the bus (via `SHUNTbs`).
    """

    sG: dict = field(default_factory=dict)
    G: dict = field(default_factory=dict)
    STOR: dict = field(default_factory=dict)
    D: dict = field(default_factory=dict)
    L_from: dict = field(default_factory=dict)
    L_to: dict = field(default_factory=dict)
    T_hv: dict = field(default_factory=dict)
    T_lv: dict = field(default_factory=dict)
    SHUNT: dict = field(default_factory=dict)


def _from_pairs(pairs, index_set, element_first):
    """Group a (bus, element) relation by bus, ordered like `index_set`.

    Args:
        pairs: The relation, an iterable of 2-tuples.
        index_set: The element's Pyomo set; fixes the order within a bus
            and drops pairs whose element is not a member.
        element_first: `True` if a pair reads `(element, bus)`, `False` if
            it reads `(bus, element)`.

    Returns:
        A dict mapping each bus to its list of elements.
    """
    position = {e: i for i, e in enumerate(index_set)}
    grouped = {}
    for first, second in pairs:
        element, bus = (first, second) if element_first else (second, first)
        if element in position:
            grouped.setdefault(bus, []).append(element)
    for elements in grouped.values():
        elements.sort(key=position.__getitem__)
    return grouped


def _from_ends(index_set, param, end):
    """Group elements by the bus stored in `param[element, end]`.

    Args:
        index_set: The element's Pyomo set, iterated in order.
        param: The bus-element matrix (`model.A` or `model.AT`).
        end: 1 for the from / HV end, 2 for the to / LV end.

    Returns:
        A dict mapping each bus to its list of elements.
    """
    grouped = {}
    for e in index_set:
        grouped.setdefault(param[e, end], []).append(e)
    return grouped


def bus_incidence(model, storage=False):
    """Compute the bus -> element lists of a power-flow model.

    Reads `model.Gbs`, `sGbs`, `Dbs`, `SHUNTbs`, `A`, `AT` and, with
    `storage`, `STOR_bus` as they are now. The model is not modified.

    Args:
        model: A Pyomo model carrying the single- or multi-period base sets
            (`G`, `sG`, `D`, `L`, `TRANSF`, `SHUNT` and the relations
            above).
        storage: Also fill `STOR` from `model.STOR` / `model.STOR_bus`.
            Only the single-period balances have a storage term, so the
            multi-period models leave this off.

    Returns:
        A [`BusIncidence`][potpourri.models.incidence.BusIncidence].
    """
    stor = {}
    if storage:
        stor = _from_ends_scalar(model.STOR, model.STOR_bus)
    return BusIncidence(
        sG=_from_pairs(model.sGbs, model.sG, element_first=True),
        G=_from_pairs(model.Gbs, model.G, element_first=True),
        STOR=stor,
        D=_from_pairs(model.Dbs, model.D, element_first=False),
        L_from=_from_ends(model.L, model.A, 1),
        L_to=_from_ends(model.L, model.A, 2),
        T_hv=_from_ends(model.TRANSF, model.AT, 1),
        T_lv=_from_ends(model.TRANSF, model.AT, 2),
        SHUNT=_from_pairs(model.SHUNTbs, model.SHUNT, element_first=False),
    )


def _from_ends_scalar(index_set, param):
    """Group elements by the bus stored in `param[element]`.

    Args:
        index_set: The element's Pyomo set, iterated in order.
        param: A one-index parameter holding the bus of each element.

    Returns:
        A dict mapping each bus to its list of elements.
    """
    grouped = {}
    for e in index_set:
        grouped.setdefault(param[e], []).append(e)
    return grouped


#: Attribute under which `build_kcl` parks the lists for its rules.
_CACHE_ATTR = "_kcl_incidence"


def cached_incidence(owner, model, storage=False):
    """Return the incidence lists a running `build_kcl` prepared.

    The balance rules are called once per bus (and time step). While
    `build_kcl` runs it stores one
    [`BusIncidence`][potpourri.models.incidence.BusIncidence] on the owner
    through [`incidence_scope`][potpourri.models.incidence.incidence_scope];
    this returns it. Called outside a build (a rule invoked directly), it
    computes the lists afresh, which is slow per call but always right.

    Args:
        owner: The model object (`self` of the rule) holding the cache.
        model: The Pyomo model the rule is building on.
        storage: Passed to `bus_incidence` when nothing is cached.

    Returns:
        A [`BusIncidence`][potpourri.models.incidence.BusIncidence] for
        `model`.
    """
    cached = getattr(owner, _CACHE_ATTR, None)
    if cached is not None and cached[0] is model:
        return cached[1]
    return bus_incidence(model, storage=storage)


@contextmanager
def incidence_scope(owner, model, storage=False):
    """Compute the incidence lists once for the duration of a build.

    Args:
        owner: The model object whose balance rules call
            [`cached_incidence`][potpourri.models.incidence.cached_incidence].
        model: The Pyomo model being built.
        storage: Include the storage lists (single-period balances).

    Yields:
        The [`BusIncidence`][potpourri.models.incidence.BusIncidence].
    """
    inc = bus_incidence(model, storage=storage)
    setattr(owner, _CACHE_ATTR, (model, inc))
    try:
        yield inc
    finally:
        setattr(owner, _CACHE_ATTR, None)


def python_floats(series):
    """Turn a float Series into a `{label: float}` dict of Python floats.

    A Pyomo `Param` initialised from a pandas Series keeps numpy scalars,
    and every `np.float64 * <Pyomo expression>` then goes through numpy's
    `__array_ufunc__` dispatch before Pyomo sees it, which dominates the
    time spent writing the branch equations of a large grid. Python floats
    skip that detour and hold the same IEEE value, so the expressions and
    everything computed from them are unchanged.

    Args:
        series: A numeric Series indexed by the Param's index set.

    Returns:
        A dict mapping each label to `float(value)`, in the Series' order.
    """
    return dict(zip(series.index, series.astype(float).tolist()))
