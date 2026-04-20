# Flexible Devices API

All device classes live in `src/potpourri/technologies/` and inherit from `Flexibility_multi_period`.

---

## Flexibility_multi_period

```
potpourri.technologies.flexibility.Flexibility_multi_period
```

Abstract base class for all modular flexible device extensions. Provides common data extraction from the network and the standard interface that all device classes must implement.

### `__init__(net, T=None, scenario=None)`

**Parameters:**

- `net` — pandapower network
- `T` — list of time step indices; if `None`, defaults to `[0]`
- `scenario` — scenario identifier (device-specific meaning)

Extracts `baseMVA`, bus lookup table, and time-series profiles from `net.profiles`.

### Standard interface

All device subclasses implement:

| Method | Called by | Description |
|---|---|---|
| `get_sets(model)` | `get_all(model)` | Add Pyomo Sets |
| `get_parameters(model)` | `get_all(model)` | Add Pyomo Params |
| `get_variables(model)` | `get_all(model)` | Add Pyomo Vars |
| `get_all_constraints(model)` | `get_all(model)` | Add Pyomo Constraints |
| `get_all(model)` | `create_model()` | Calls all four above |
| `get_all_opf(model)` | `OPF_multi_period.add_OPF()` | OPF-specific additions |
| `get_all_acopf(model)` | `ACOPF_multi_period.add_OPF()` | AC-OPF specific additions |

### `make_to_dict(model_obj, model_time, data, time_dependent=True)`

Converts data arrays or dicts into Pyomo parameter initialisation format.

**Parameters:**

- `model_obj` — iterable of object indices
- `model_time` — iterable of time indices
- `data` — numpy array, pandas Series, or dict
- `time_dependent` — if `True`, output is keyed by `(obj, t)`; if `False`, keyed by `obj` only

**Returns:** `(data_dict, tuple_list)`

---

## Battery_multi_period

```
potpourri.technologies.battery.Battery_multi_period
```

Models grid-connected battery storage with SOC dynamics.

**Inherits:** `Flexibility_multi_period`

### `__init__(net, T=None, scenario=None)`

Randomly places batteries on network buses. The fraction of buses with batteries is controlled by `scenario` (0 = 2020 scenario, 1 = 2030, etc.).

**Hardcoded defaults:**

| Name | Value | Unit |
|---|---|---|
| `bat_power` | 0.006 | p.u. (MW/baseMVA) |
| `bat_cap` | 0.015 | p.u. (MWh/baseMVA) |
| `bat_soc_max` | 1.0 | — |
| `bat_soc_min` | 0.2 | — |
| `bat_eff` | 0.9 | — |

**Pyomo components added:**

- Sets: `BAT`, `BAT_bus`
- Params: `BAT_Pmax[b]`, `BAT_Pmin[b]`, `BAT_SOCmax[b]`, `BAT_SOCmin[b]`, `BAT_Cap[b]`, `BAT_Eff[b]`
- Vars: `BAT_P[b, t]`, `BAT_SOC[b, t]`
- Constraints: power bounds, SOC bounds, SOC dynamics, initial SOC

---

## Heatpump_multi_period

```
potpourri.technologies.heat_pump.Heatpump_multi_period
```

Models thermal loads with an electric heat pump and building thermal mass.

**Inherits:** `Flexibility_multi_period`

### `__init__(net, T=None, scenario=None)`

Uses a simplified building model with heat capacity, thermal mass, and heat loss parameters to compute heating demand. The heat pump operates with a fixed COP.

**Pyomo components added:**

- Sets: heat pump indices, bus-HP mapping
- Params: `hp_power`, COP, thermal mass, heat loss coefficients
- Vars: `hp_p[h, t]` (electrical power), `temp[h, t]` (indoor temperature)
- Constraints: thermal dynamics, temperature bounds, power limits

---

## PV_multi_period

```
potpourri.technologies.pv.PV_multi_period
```

Models controllable PV generators with generation profiles from `net.profiles`. Percentages of PV penetration are configurable.

**Inherits:** `Flexibility_multi_period`

---

## Windpower_multi_period

```
potpourri.technologies.windpower.Windpower_multi_period
```

Extends `Sgens_multi_period` with grid-code Q-curve constraints for wind generators. Reactive power control variants (0–2) correspond to different Q-P and Q-U slope/intercept parameters.

**Inherits:** `Sgens_multi_period`

### `static_generation_wind_var_q(net)`

Reads `sgen.var_q` to assign grid-code reactive power bounds for each wind generator. Same logic as `ACOPF.static_generation_wind_var_q()`.

### `get_objective(model)`

Adds a wind generation maximisation objective to the model.

---

## Demand_multi_period

```
potpourri.technologies.demand.Demand_multi_period
```

Handles time-indexed load power (P and Q) from SimBench profiles. Fixes load variables to profile values in base operation; adds controllable load bounds in OPF mode.

**Inherits:** `Flexibility_multi_period`

### `get_demand_real_power_data(model, max_p_mw=None, min_p_mw=None)`

Reads load real power bounds and attaches `PDmax`, `PDmin` parameters and `PD_Constraint` to the model.

### `get_demand_reactive_data(model, max_q_mvar=None, min_q_mvar=None)`

Reads load reactive power bounds and attaches `QDmax`, `QDmin` parameters and `QD_Constraint` to the model.

---

## Sgens_multi_period

```
potpourri.technologies.sgens.Sgens_multi_period
```

Handles time-indexed static generator (sgen) real and reactive power from profiles. In OPF mode, adds controllable sgen power bounds.

**Inherits:** `Flexibility_multi_period`

### `static_generation_real_power_limits(model)`

Reads `sgen.max_p_mw` / `min_p_mw` and adds `sPGmax`, `sPGmin` params and `PsG_Constraint` to the model.

### `static_generation_reactive_power_limits(model)`

Reads `sgen.max_q_mvar` / `min_q_mvar` and adds `QsGmax`, `QsGmin` params and `QsG_Constraint` to the model.

---

## Generator_multi_period

```
potpourri.technologies.generator.Generator_multi_period
```

Handles time-indexed external grid and synchronous generator real and reactive power. Provides OPF limits for real power and ACOPF limits for reactive power.

**Inherits:** `Flexibility_multi_period`

### `generation_real_power_limits_opf(model)`

Reads generator real power limits and attaches `PGmax`, `PGmin` params and `PG_Constraint` to the model.

### `generation_reactive_power_limits_acopf()`

Reads generator reactive power limits into `generation_data['max_q']` and `['min_q']`.
