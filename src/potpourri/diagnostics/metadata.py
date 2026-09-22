# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""What each Pyomo component in a potpourri model physically means.

Diagnostics need to know that `line_lim_from` is a from-side thermal limit
on a line and that `v` is a bus voltage magnitude in per unit. Deriving that
from the component name at the point of use would mean string-matching
scattered through the package, which breaks the first time a name changes.
This module holds it in one place instead.

The registry describes the components as they are named today. Those names
are public API — tests and scripts index into `model.line_lim_from` — so
nothing here renames anything; it annotates. A component that is not
registered still works everywhere in this package, it simply gets no
physical interpretation, which is why every lookup returns `None` rather
than raising.

Units follow the model: voltages and angles in per unit and radians,
powers in per unit on `baseMVA`. Converting to MW/MVAr/kA for the user is
the reporting layer's job, not this one's.
"""

from __future__ import annotations

from dataclasses import dataclass

from potpourri.diagnostics.report import DiagnosticCategory


@dataclass(frozen=True)
class ConstraintMeta:
    """The physical meaning of one constraint family.

    Attributes:
        category: Which `DiagnosticCategory` a violation belongs to.
        element: pandapower table the index refers to, e.g. `"line"`.
            `None` where the constraint is not per element.
        description: One phrase naming what the constraint enforces.
        side: `"from"` or `"to"` for branch-end constraints.
        squared: True when the constraint is written on squared
            quantities, so a residual is in units of power squared and has
            to be converted before it is shown to anyone.
    """

    category: DiagnosticCategory
    element: str | None
    description: str
    side: str | None = None
    squared: bool = False


@dataclass(frozen=True)
class VariableMeta:
    """The physical meaning of one variable family.

    Attributes:
        quantity: What the variable represents, e.g.
            `"voltage_magnitude"`.
        unit: Unit as the model carries it, e.g. `"p.u."`.
        element: pandapower table the index refers to.
        per_unit: Whether the value is per unit on `baseMVA` and therefore
            needs scaling before being reported in MW or MVAr.
    """

    quantity: str
    unit: str
    element: str | None
    per_unit: bool = False


#: Constraint families, keyed by Pyomo component name. Covers the AC, DC and
#: LPAC single-period models and their multi-period twins, whose components
#: carry the same names with a time index appended.
CONSTRAINTS: dict[str, ConstraintMeta] = {
    # nodal balance
    "KCL_real": ConstraintMeta(
        DiagnosticCategory.POWER_BALANCE, "bus", "active power balance"
    ),
    "KCL_reactive": ConstraintMeta(
        DiagnosticCategory.POWER_BALANCE, "bus", "reactive power balance"
    ),
    "KCL_const": ConstraintMeta(
        DiagnosticCategory.POWER_BALANCE, "bus", "active power balance"
    ),
    "KCL_def": ConstraintMeta(
        DiagnosticCategory.POWER_BALANCE, "bus", "active power balance"
    ),
    # branch flow definitions (pi-model, despite the historical KVL names)
    "KVL_real_from": ConstraintMeta(
        DiagnosticCategory.POWER_BALANCE, "line", "active branch flow", "from"
    ),
    "KVL_real_to": ConstraintMeta(
        DiagnosticCategory.POWER_BALANCE, "line", "active branch flow", "to"
    ),
    "KVL_reactive_from": ConstraintMeta(
        DiagnosticCategory.POWER_BALANCE,
        "line",
        "reactive branch flow",
        "from",
    ),
    "KVL_reactive_to": ConstraintMeta(
        DiagnosticCategory.POWER_BALANCE, "line", "reactive branch flow", "to"
    ),
    "KVL_real_fromTransf": ConstraintMeta(
        DiagnosticCategory.POWER_BALANCE,
        "trafo",
        "active transformer flow",
        "from",
    ),
    "KVL_real_toTransf": ConstraintMeta(
        DiagnosticCategory.POWER_BALANCE,
        "trafo",
        "active transformer flow",
        "to",
    ),
    "KVL_reactive_fromTransf": ConstraintMeta(
        DiagnosticCategory.POWER_BALANCE,
        "trafo",
        "reactive transformer flow",
        "from",
    ),
    "KVL_reactive_toTransf": ConstraintMeta(
        DiagnosticCategory.POWER_BALANCE,
        "trafo",
        "reactive transformer flow",
        "to",
    ),
    # thermal limits
    "line_lim_from": ConstraintMeta(
        DiagnosticCategory.THERMAL,
        "line",
        "thermal limit",
        "from",
        squared=True,
    ),
    "line_lim_to": ConstraintMeta(
        DiagnosticCategory.THERMAL, "line", "thermal limit", "to", squared=True
    ),
    "transf_lim1": ConstraintMeta(
        DiagnosticCategory.THERMAL,
        "trafo",
        "thermal limit",
        "from",
        squared=True,
    ),
    "transf_lim2": ConstraintMeta(
        DiagnosticCategory.THERMAL,
        "trafo",
        "thermal limit",
        "to",
        squared=True,
    ),
    "thermal_facet_line_from": ConstraintMeta(
        DiagnosticCategory.THERMAL, "line", "polygonal thermal limit", "from"
    ),
    "thermal_facet_line_to": ConstraintMeta(
        DiagnosticCategory.THERMAL, "line", "polygonal thermal limit", "to"
    ),
    # voltage
    "v_pyo": ConstraintMeta(
        DiagnosticCategory.VOLTAGE, "bus", "voltage magnitude limits"
    ),
    "v_constraint": ConstraintMeta(
        DiagnosticCategory.VOLTAGE, "bus", "voltage magnitude limits"
    ),
    "v_fixed": ConstraintMeta(
        DiagnosticCategory.VOLTAGE, "bus", "fixed voltage setpoint"
    ),
    "v_bPV_setpoint": ConstraintMeta(
        DiagnosticCategory.VOLTAGE, "bus", "PV-bus voltage setpoint"
    ),
    # generator and sgen capability
    "PsG_Constraint": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY,
        "sgen",
        "active power limits",
    ),
    "QsG_pyo": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY,
        "sgen",
        "reactive power limits",
    ),
    "PG_Constraint": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY,
        "ext_grid",
        "active power limits",
    ),
    "QG_pyo": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY,
        "ext_grid",
        "reactive power limits",
    ),
    "PD_Constraint": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY, "load", "active demand limits"
    ),
    "QD_pyo": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY,
        "load",
        "reactive demand limits",
    ),
    "sgen_inverter_s2": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY,
        "sgen",
        "inverter apparent-power circle",
        squared=True,
    ),
    "sgen_cos_phi_upper": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY, "sgen", "cos(phi) cone"
    ),
    "sgen_cos_phi_lower": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY, "sgen", "cos(phi) cone"
    ),
    "sgen_fixed_cos_phi": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY, "sgen", "fixed cos(phi)"
    ),
    "sgen_pu_curtail": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY, "sgen", "P(U) curtailment"
    ),
    "sgen_cpp": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY,
        "sgen",
        "cos(phi)(P) characteristic",
    ),
    "QW_pos_pyo": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY, "sgen", "wind Q(P) capability"
    ),
    "QW_neg_pyo": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY, "sgen", "wind Q(P) capability"
    ),
    "QU_min_pyo": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY, "sgen", "wind Q(U) capability"
    ),
    "QU_max_pyo": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY, "sgen", "wind Q(U) capability"
    ),
    "PV_QP_pos": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY, "sgen", "PV Q(P) capability"
    ),
    "PV_QP_neg": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY, "sgen", "PV Q(P) capability"
    ),
    "PV_QU_min": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY, "sgen", "PV Q(U) capability"
    ),
    "PV_QU_max": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY, "sgen", "PV Q(U) capability"
    ),
    "inv_capability_facet": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY,
        "sgen",
        "inverter capability envelope",
    ),
    # angles
    "line_angle_diff": ConstraintMeta(
        DiagnosticCategory.NETWORK, "line", "phase-angle difference limit"
    ),
    "trafo_angle_diff": ConstraintMeta(
        DiagnosticCategory.NETWORK, "trafo", "phase-angle difference limit"
    ),
    "phase_diff1": ConstraintMeta(
        DiagnosticCategory.NETWORK, "line", "phase-angle difference limit"
    ),
    "phase_diff2": ConstraintMeta(
        DiagnosticCategory.NETWORK, "line", "phase-angle difference limit"
    ),
    # transformer taps
    "Tap_linear_constr": ConstraintMeta(
        DiagnosticCategory.BOUNDS, "trafo", "continuous tap-ratio limits"
    ),
    "Tap_pos_constr": ConstraintMeta(
        DiagnosticCategory.BOUNDS, "trafo", "tap-position limits"
    ),
    "Tap_discrete_constr": ConstraintMeta(
        DiagnosticCategory.BOUNDS, "trafo", "discrete tap positions"
    ),
    # storage
    "stor_soc_update": ConstraintMeta(
        DiagnosticCategory.STORAGE, "storage", "state-of-charge balance"
    ),
    "stor_soc_bounds": ConstraintMeta(
        DiagnosticCategory.STORAGE, "storage", "state-of-charge limits"
    ),
    "stor_chg_limit": ConstraintMeta(
        DiagnosticCategory.STORAGE, "storage", "charging power limit"
    ),
    "stor_dis_limit": ConstraintMeta(
        DiagnosticCategory.STORAGE, "storage", "discharging power limit"
    ),
    "stor_no_simul": ConstraintMeta(
        DiagnosticCategory.STORAGE,
        "storage",
        "no simultaneous charge/discharge",
    ),
    "stor_inverter_cap": ConstraintMeta(
        DiagnosticCategory.STORAGE,
        "storage",
        "inverter capability",
        squared=True,
    ),
    "bat_soc_con": ConstraintMeta(
        DiagnosticCategory.STORAGE, "storage", "battery state-of-charge limits"
    ),
    "bat_soc_update_con": ConstraintMeta(
        DiagnosticCategory.STORAGE, "storage", "battery energy balance"
    ),
    "bat_terminal_soc_con": ConstraintMeta(
        DiagnosticCategory.STORAGE, "storage", "terminal state of charge"
    ),
    "bat_power_con": ConstraintMeta(
        DiagnosticCategory.STORAGE, "storage", "battery power limits"
    ),
    # hosting capacity
    "SW_max": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY,
        "sgen",
        "hosting-capacity dispatch within installed rating",
        squared=True,
    ),
    "hc_size_upper": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY,
        "sgen",
        "hosting-capacity sizing",
    ),
    "hc_size_lower": ConstraintMeta(
        DiagnosticCategory.GENERATOR_CAPABILITY,
        "sgen",
        "hosting-capacity sizing",
    ),
}


#: Variable families, keyed by Pyomo component name.
VARIABLES: dict[str, VariableMeta] = {
    "v": VariableMeta("voltage_magnitude", "p.u.", "bus"),
    "delta": VariableMeta("voltage_angle", "rad", "bus"),
    "pLfrom": VariableMeta("active_flow", "p.u.", "line", per_unit=True),
    "pLto": VariableMeta("active_flow", "p.u.", "line", per_unit=True),
    "qLfrom": VariableMeta("reactive_flow", "p.u.", "line", per_unit=True),
    "qLto": VariableMeta("reactive_flow", "p.u.", "line", per_unit=True),
    "pThv": VariableMeta("active_flow", "p.u.", "trafo", per_unit=True),
    "pTlv": VariableMeta("active_flow", "p.u.", "trafo", per_unit=True),
    "qThv": VariableMeta("reactive_flow", "p.u.", "trafo", per_unit=True),
    "qTlv": VariableMeta("reactive_flow", "p.u.", "trafo", per_unit=True),
    "psG": VariableMeta("active_power", "p.u.", "sgen", per_unit=True),
    "qsG": VariableMeta("reactive_power", "p.u.", "sgen", per_unit=True),
    "pG": VariableMeta("active_power", "p.u.", "ext_grid", per_unit=True),
    "qG": VariableMeta("reactive_power", "p.u.", "ext_grid", per_unit=True),
    "pD": VariableMeta("active_power", "p.u.", "load", per_unit=True),
    "qD": VariableMeta("reactive_power", "p.u.", "load", per_unit=True),
    "Tap": VariableMeta("tap_ratio", "-", "trafo"),
    "SOC": VariableMeta("state_of_charge", "p.u.", "storage"),
    "pSTOR": VariableMeta("active_power", "p.u.", "storage", per_unit=True),
}


def constraint_meta(component_name: str) -> ConstraintMeta | None:
    """Look up a constraint family's meaning.

    Args:
        component_name: The Pyomo component's `local_name`, without any
            index.

    Returns:
        Its `ConstraintMeta`, or `None` when the family is not registered.
        Callers must treat `None` as "no physical interpretation available"
        rather than as an error.
    """
    return CONSTRAINTS.get(component_name)


def variable_meta(component_name: str) -> VariableMeta | None:
    """Look up a variable family's meaning.

    Args:
        component_name: The Pyomo component's `local_name`, without any
            index.

    Returns:
        Its `VariableMeta`, or `None` when the family is not registered.
    """
    return VARIABLES.get(component_name)


def element_table_for(component_name: str) -> str | None:
    """The pandapower table a component's index refers to.

    Args:
        component_name: A constraint or variable component name.

    Returns:
        The table name, or `None` when unknown. Checks constraints first,
        then variables, since the two namespaces do not overlap.
    """
    meta = CONSTRAINTS.get(component_name) or VARIABLES.get(component_name)
    return meta.element if meta else None
