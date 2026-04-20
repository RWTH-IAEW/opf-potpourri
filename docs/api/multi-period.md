# Multi-Period Models API

Multi-period OPF model classes live in `src/potpourri/models_multi_period/`. Device/technology mix-in classes (batteries, EVs, heat pumps, PV, wind, demand, generators) live in `src/potpourri/technologies/`.

---

## Basemodel_multi_period

```
potpourri.models_multi_period.basemodel_multi_period.Basemodel_multi_period
```

Extends `Basemodel` with a time dimension and a modular framework for flexible devices. Integrates SimBench time-series profiles for loads and generators.

**Inherits:** `Basemodel`

### `__init__(net, toT, fromT=None, pf=1, num_vehicles=None)`

**Parameters:**

| Parameter | Description |
|---|---|
| `net` | SimBench pandapower network with `net.profiles` |
| `toT` | End time step (exclusive, 0-indexed) |
| `fromT` | Start time step (inclusive); defaults to `0` |
| `pf` | Power factor for reactive power initialisation (default `1`) |
| `num_vehicles` | Number of EV profiles to load; activates `EV_multi_period` if set |

The constructor instantiates the following flexibility objects and appends them to `self.flexibilities`:

- `Demand_multi_period`
- `Shunts_multi_period`
- `Sgens_multi_period`
- `Generator_multi_period`
- `EV_multi_period` (if `num_vehicles` is set)
- `Windpower_multi_period` (if `net.bus` has a `windpot_p_mw` column)

### `create_model()`

Builds the Pyomo `ConcreteModel` with time-indexed sets and variables. Calls `flex.get_all(model)` for each flexibility object. Adds:

- `model.T = Set(range(fromT, toT))` — time index
- `model.deltaT = 0.25` — time step duration (hours, fixed at 15 min)
- Time-indexed variables: `delta[b, t]`, `v[b, t]`, `pLfrom[l, t]`, etc.

### `make_to_dict(model_obj, model_time, data, time_dependent=True)`

Converts load/generation data (numpy array, pandas Series, or dict) into a Pyomo-compatible dict keyed by `(object_index, time_index)` tuples.

**Returns:** `(data_dict, tuple_list)`

### `make_to_tuple(model_obj, model_time)`

Creates a list of `(object_index, time_index)` tuples for use as Pyomo set initialisers.

**Returns:** `tuple_list`

---

## ACOPF_multi_period

```
potpourri.models_multi_period.ACOPF_multi_period.ACOPF_multi_period
```

Full multi-period AC OPF model. Inherits from `AC_multi_period` and `OPF_multi_period`.

**Inherits:** `AC_multi_period`, `OPF_multi_period`

### `add_OPF(**kwargs)`

Calls `flex.get_all_acopf(model)` for each flexibility object, then adds time-indexed constraints:

- `Vmin[b] ≤ v[b, t] ≤ Vmax[b]` for all buses and time steps
- Apparent power limits on lines and transformers

### Objective methods

| Method | Objective |
|---|---|
| `add_voltage_deviation_objective()` | `min Σ_t Σ_b (v[b,t] - 1)²` |
| `add_minimize_power_objective()` | `min Σ_d,t pD[d,t]` |
| `add_generation_objective()` | `min Σ_g,t pG[g,t]²` |
| `add_weighted_generation_objective()` | `min 4·Σ pG + Σ p_discharging + Σ psG` |
| `add_charging_power_obj()` | Maximise total EV charging power |
| `add_discharging_power_obj()` | Maximise total EV discharging (V2G) |
| `add_arbitrage_objective()` | `min Σ_v,t p_opf[v,t] · p_da[t] · Δt/60` |

---

## OPF_multi_period

```
potpourri.models_multi_period.OPF_multi_period.OPF_multi_period
```

Adds OPF operational constraints (power limits, thermal limits, tap changers) to the multi-period base model.

**Inherits:** `Basemodel_multi_period`

### `add_OPF(**kwargs)`

Adds non-time-dependent thermal limit parameters (`SLmax[l]`, `SLmaxT[t]`) and calls `flex.get_all_opf(model)` for each flexibility object.

### `add_tap_changer_linear()`

Enables continuous transformer tap optimisation (same as single-period OPF).

### `add_tap_changer_discrete()`

Enables discrete tap position optimisation (requires MIP solver).

---

## AC_multi_period

```
potpourri.models_multi_period.AC_multi_period.AC_multi_period
```

Adds time-indexed AC power flow equations (KCL real/reactive, KVL real/reactive on lines and transformers) to the multi-period base model. All constraints are indexed over `model.T`.

**Inherits:** `Basemodel_multi_period`
