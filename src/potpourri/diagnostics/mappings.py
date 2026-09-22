# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Translate Pyomo indices back into pandapower objects.

Every diagnostic in this package ends with a sentence a power-system
engineer can act on, and that only works if `line_lim_from[23]` can be
turned into `net.line[7] "Cable 7"` reliably. Four index spaces are in play
and none of them can be assumed equal to another:

| space | what it indexes |
|---|---|
| caller's `net` | the network the user handed to the model |
| `model.net` | potpourri's preprocessed copy, with buses renumbered |
| ppc | pandapower's internal solver structure |
| Pyomo | the index sets on the Pyomo model |

**Buses move; other elements do not.** `Basemodel.__init__` runs
`preprocess_grid`, which merges buses joined by closed zero-impedance
bus-bus switches and then renumbers what is left to `0..n-1`. Line, load,
sgen and transformer indices survive that untouched, so only buses need a
map back to the caller. `preprocess_grid` records one under
`BUS_ORIGIN_KEY`; the multi-period models do no preprocessing at all, so
there the caller's indices are already the model's.

**`model.B` is in ppc numbering**, not pandapower numbering, and it can
contain auxiliary buses that pandapower invented for switch handling and
that have no row in `net.bus` at all. Those resolve to `None` rather than
to a wrong bus.

**`model.L` spans two tables.** potpourri models `net.impedance` rows as
lines, appending them after the real lines, so a Pyomo line index resolves
to either a `line` or an `impedance` row. The resolution is positional, not
by index value, because the synthetic indices the model gives impedance
rows can collide with real line indices on a network whose line indices are
non-contiguous.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from potpourri.models.basemodel import BUS_ORIGIN_KEY


@dataclass(frozen=True)
class ElementRef:
    """One pandapower object, named the way a user would refer to it.

    Attributes:
        table: pandapower table name, e.g. `"bus"` or `"line"`. `None` when
            an index could not be resolved to any object.
        index: Index into that table, in the *caller's* numbering.
        name: The object's `name` column, where it has one.
    """

    table: str | None
    index: int | None
    name: str | None = None

    def __str__(self) -> str:
        """Render as `net.line[7] "Cable 7"`, or a clear unknown."""
        if self.table is None or self.index is None:
            return "<unmapped>"
        label = f"net.{self.table}[{self.index}]"
        return f'{label} "{self.name}"' if self.name else label

    @property
    def is_resolved(self) -> bool:
        """Whether this points at a real pandapower object."""
        return self.table is not None and self.index is not None


@dataclass
class IndexMap:
    """Resolves Pyomo indices to `ElementRef`s for one model.

    Build it with `IndexMap.from_model`. Every lookup returns an
    `ElementRef`; unresolvable indices come back unresolved rather than
    raising, because a diagnostic run must survive a partially built model.

    Attributes:
        bus: Pyomo bus index (ppc numbering) to caller's bus index.
        line: Pyomo `model.L` index to `(table, caller index)`, where table
            is `"line"` or `"impedance"`.
        aux_buses: Pyomo bus indices that are auxiliary ppc nodes with no
            pandapower row.
    """

    bus: dict[int, int] = field(default_factory=dict)
    line: dict[int, tuple[str, int]] = field(default_factory=dict)
    aux_buses: set[int] = field(default_factory=set)
    merged_buses: dict[int, list[int]] = field(default_factory=dict)
    _bus_model_index: dict[int, int] = field(default_factory=dict)
    _net: Any = None

    # --- construction -------------------------------------------------

    @classmethod
    def from_model(cls, model_obj) -> IndexMap:
        """Build the map for a potpourri model object.

        Args:
            model_obj: A `Basemodel` or multi-period equivalent, after
                construction. `add_OPF` need not have run.

        Returns:
            An `IndexMap`. Spaces the model does not expose stay empty, so
            a lookup against them returns an unresolved `ElementRef`.
        """
        net = getattr(model_obj, "net", None)
        out = cls(_net=net)
        if net is None:
            return out
        out._build_bus_map(model_obj, net)
        out._build_line_map(model_obj, net)
        return out

    def _build_bus_map(self, model_obj, net) -> None:
        """Compose ppc -> model bus -> caller bus."""
        # model bus (pandapower numbering inside model.net) -> caller bus
        origin = net.get(BUS_ORIGIN_KEY) if hasattr(net, "get") else None
        if isinstance(origin, dict) and origin:
            # The map is caller -> model and several caller buses can share
            # one model bus, because a closed zero-impedance bus-bus switch
            # fuses them. Report the lowest of each group and remember the
            # rest, so a diagnostic can say which buses were merged instead
            # of naming an arbitrary one of them.
            groups: dict[int, list[int]] = {}
            for caller, model_bus in origin.items():
                groups.setdefault(int(model_bus), []).append(int(caller))
            model_to_caller = {}
            for model_bus, callers in groups.items():
                callers.sort()
                model_to_caller[model_bus] = callers[0]
                if len(callers) > 1:
                    self.merged_buses[callers[0]] = callers[1:]
        else:
            # Multi-period models do not preprocess, so the model's bus
            # numbering is already the caller's.
            model_to_caller = {int(b): int(b) for b in net.bus.index}

        lookup = getattr(model_obj, "bus_lookup", None)
        if lookup is None:
            self.bus = model_to_caller
            return

        for pp_bus in net.bus.index:
            try:
                ppc_bus = int(lookup[pp_bus])
            except (KeyError, IndexError, TypeError):
                continue
            caller = model_to_caller.get(int(pp_bus))
            if caller is not None:
                self.bus[ppc_bus] = caller
                self._bus_model_index[ppc_bus] = int(pp_bus)

        model = getattr(model_obj, "model", None)
        if model is not None and hasattr(model, "B"):
            self.aux_buses = {
                int(b) for b in model.B if int(b) not in self.bus
            }

    def _build_line_map(self, model_obj, net) -> None:
        """Resolve `model.L` positionally across line and impedance rows."""
        line_data = getattr(model_obj, "line_data", None)
        if line_data is None:
            return
        n_line = len(net.line.index)
        line_index = list(net.line.index)
        imp = net.get("impedance")
        imp_index = list(imp.index) if imp is not None and len(imp) else []

        for position, model_index in enumerate(line_data.index):
            if position < n_line:
                self.line[int(model_index)] = (
                    "line",
                    int(line_index[position]),
                )
            else:
                offset = position - n_line
                if offset < len(imp_index):
                    self.line[int(model_index)] = (
                        "impedance",
                        int(imp_index[offset]),
                    )

    # --- lookups ------------------------------------------------------

    def _named(self, table: str, index: int) -> ElementRef:
        """Attach the object's `name` column if the table carries one."""
        name = None
        net = self._net
        if net is not None and table in net:
            frame = net[table]
            if (
                "name" in getattr(frame, "columns", [])
                and index in frame.index
            ):
                raw = frame.at[index, "name"]
                name = None if raw is None or raw != raw else str(raw)
        return ElementRef(table, index, name)

    def bus_ref(self, pyomo_index) -> ElementRef:
        """The pandapower bus behind a Pyomo bus index.

        Args:
            pyomo_index: An index from `model.B`, in ppc numbering.

        Returns:
            The caller's bus, or an unresolved reference when the index is
            an auxiliary ppc node with no pandapower row.
        """
        caller = self.bus.get(int(pyomo_index))
        if caller is None:
            return ElementRef(None, None)
        # The name lives in the model's own bus table, under the model's
        # numbering; the index reported is the caller's. Looking the
        # caller's index up in the model's table would silently name a
        # different bus, or none.
        model_index = self._bus_model_index.get(int(pyomo_index))
        name = None
        net = self._net
        if net is not None and model_index is not None and "bus" in net:
            frame = net["bus"]
            if "name" in frame.columns and model_index in frame.index:
                raw = frame.at[model_index, "name"]
                name = None if raw is None or raw != raw else str(raw)
        return ElementRef("bus", caller, name)

    def line_ref(self, pyomo_index) -> ElementRef:
        """The pandapower line or impedance behind a Pyomo `model.L` index.

        Args:
            pyomo_index: An index from `model.L`.

        Returns:
            A reference into `net.line` or `net.impedance`.
        """
        entry = self.line.get(int(pyomo_index))
        if entry is None:
            return ElementRef(None, None)
        table, index = entry
        return self._named(table, index)

    def element_ref(self, table: str, pyomo_index) -> ElementRef:
        """A pandapower object whose index the model did not change.

        Transformers, loads, sgens, storage, generators and external grids
        keep the caller's indices through preprocessing, so this is a
        lookup rather than a translation. It still verifies the row exists,
        so a stale index is reported as unmapped instead of being printed
        as if it were real.

        Args:
            table: pandapower table name, e.g. `"trafo"`.
            pyomo_index: The Pyomo index, which equals the pandapower one.

        Returns:
            A reference to that object, unresolved if the row is absent.
        """
        net = self._net
        index = int(pyomo_index)
        if net is None or table not in net or index not in net[table].index:
            return ElementRef(None, None)
        return self._named(table, index)
