# Single-Period Models API

All single-period models live in `src/potpourri/models/`.

---

## Basemodel

```
potpourri.models.basemodel.Basemodel
```

A Pyomo-based optimisation model for power system analysis. Creates and solves optimisation models based on a given pandapower network. Supports creation of sets, parameters, and variables, as well as solving and post-processing.

**Attributes:**

- `net` — deep copy of the input pandapower network
- `model` — Pyomo `ConcreteModel` created during construction
- `results` — solver result object after calling `solve()`
- `baseMVA` — system apparent power base (MVA)

### `__init__(net)`

Initialises the Basemodel with a deep copy of the given network. Runs `pp.runpp(net)` to obtain an initial power flow solution, then extracts bus, line, transformer, generator, load, and shunt data into internal DataFrames. Calls `create_model()`.

**Parameters:**

- `net` — pandapower network

### `create_model()`

Creates a Pyomo model based on the input network. The model includes:

**Sets:**

| Set | Description |
|---|---|
| `B` | All buses |
| `b0` | Slack (reference) buses |
| `bPV` | PV buses (voltage-controlled generators) |
| `G` | Controllable generators (ext_grid) |
| `sG` | Static generators (sgen) |
| `D` | Loads |
| `L` | Lines |
| `TRANSF` | Transformers |
| `SHUNT` | Shunts |

**Key parameters:** `A[l, 1/2]` (line-bus incidence), `AT[t, 1/2]` (trafo-bus incidence), `PG[g]`, `PsG[g]`, `PD[d]`, `shift[t]`, `Tap[t]`, `baseMVA`

**Key variables:** `delta[b]` (voltage angle, rad), `pG[g]`, `psG[g]`, `pD[d]`, `pLfrom[l]`, `pLto[l]`, `pThv[t]`, `pTlv[t]`

### `solve(solver, to_net, print_solver_output, mip_solver, max_iter, time_limit, init_strategy, neos_opt)`

Solves the optimisation model using the specified solver.

**Parameters:**

| Parameter | Default | Description |
|---|---|---|
| `solver` | `'ipopt'` | Solver name: `'ipopt'`, `'mindtpy'`, `'neos'` |
| `to_net` | `True` | Write results back to `self.net.res_*` after solving |
| `print_solver_output` | `False` | Stream solver output to stdout |
| `mip_solver` | `'gurobi'` | MIP sub-solver for MindtPy |
| `max_iter` | `None` | Maximum solver iterations |
| `time_limit` | `600` | Solver time limit (seconds) |
| `init_strategy` | `'rNLP'` | MindtPy initialisation strategy |
| `neos_opt` | `'ipopt'` | Solver to request from NEOS |

### `change_vals(key, value)`

Changes the value of a Pyomo component in the model. `key` is the component name (string); `value` is a dict or scalar.

### `fix_vars(key, value=None)`

Fixes a Pyomo variable to `value`. If `value` is `None`, fixes each index to its current value.

### `unfix_vars(key, value=None)`

Unfixes a Pyomo variable, restoring it as a free decision variable. Optionally sets the initial value to `value`.

---

## AC

```
potpourri.models.AC.AC
```

Extends `Basemodel` with full AC power flow equations. Adds voltage magnitude variables and real/reactive KCL and KVL constraints at every bus, line, and transformer.

**Inherits:** `Basemodel`

**Additional variables:** `v[b]` (voltage magnitude, p.u.), `qsG[g]`, `qG[g]`, `qD[d]`, `qLfrom[l]`, `qLto[l]`, `qThv[t]`, `qTlv[t]`

**Additional parameters:** `Bii[l]`, `Bik[l]`, `Gii[l]`, `Gik[l]` (line admittances), `BiiT[t]`, `BikT[t]`, `GiiT[t]`, `GikT[t]` (transformer admittances), `BB[s]` (shunt susceptance), `QsG[g]`, `QD[d]`, `v_b0[b]`, `v_bPV[b]`

**Constraints added:**

- `KCL_real` — real power balance at each bus
- `KCL_reactive` — reactive power balance at each bus
- `KVL_real_from / KVL_real_to` — real power flow on lines (both ends)
- `KVL_reactive_from / KVL_reactive_to` — reactive power flow on lines
- `KVL_real_fromTransf / KVL_real_toTransf` — real power on transformers
- `KVL_reactive_fromTransf / KVL_reactive_toTransf` — reactive power on transformers

---

## DC

```
potpourri.models.DC.DC
```

Extends `Basemodel` with linearised DC power flow equations. Voltage magnitudes are fixed at 1.0 p.u. and reactive power is ignored.

**Inherits:** `Basemodel`

**Additional parameters:** `BL[l]` (line susceptance), `BLT[t]` (transformer susceptance)

**Additional variables:** `deltaL[l]` (angle difference on lines), `deltaLT[t]` (angle difference on transformers)

**Constraints added:**

- `KCL_const` — real power balance at each bus
- `KVL_real_from / KVL_real_to` — DC power flow on lines
- `KVL_trans_from / KVL_trans_to` — DC power flow on transformers
- `phase_diff1` — angle difference definition for lines
- `phase_diff2` — angle difference definition for transformers (includes phase shift)

---

## OPF

```
potpourri.models.OPF.OPF
```

Mixin that adds operational limit constraints to a power flow model. Intended for use via multiple inheritance alongside `AC` or `DC`.

**Inherits:** `Basemodel`

### `generation_real_power_limits()`

Reads generator real power limits from `net` into `generation_data`. Populates `generation_data['max_p']` and `['min_p']` (per-unit). Non-controllable generators are pinned to their current `p_mw` setpoint.

### `static_generation_real_power_limits()`

Reads static generator power limits from `net.sgen`. Populates `static_generation_data['max_p']`, `['min_p']`, and `['controllable']` (per-unit). Defaults: `max = p_mw`, `min = 0`.

### `get_demand_real_power_data()`

Reads load real power bounds from `net.load`. Populates `self.PDmax_data` and `self.PDmin_data` (per-unit). Falls back to `p_mw` / `0` if columns are absent.

### `_calc_opf_parameters(**kwargs)`

Computes all OPF limit data from the network. Calculates line apparent power limits (`SLmax`) and transformer limits (`SLmaxT`), then reads generator, static generator, and demand limits.

### `add_OPF(**kwargs)`

Attaches OPF sets, parameters, and constraints to `self.model`. Calls `_calc_opf_parameters()`, then adds:

- Sets: `sGc`, `Dc`
- Params: `PGmax/min`, `sPGmax/min`, `PDmax/min`, `SLmax`, `SLmaxT`
- Constraints: `PsG_Constraint`, `PG_Constraint`, `PD_Constraint`

### `add_tap_changer_linear()`

Enables continuous transformer tap ratio optimisation. Unfixes `Tap` variables and adds `[Tap_min, Tap_max]` bounds from `net.trafo`.

### `add_tap_changer_discrete()`

Enables discrete tap ratio optimisation via integer variable `Tap_pos`. Suitable for use with MIP/MINLP solvers (MindtPy, Gurobi).

---

## ACOPF

```
potpourri.models.ACOPF_base.ACOPF
```

Full AC Optimal Power Flow model. Combines AC power flow physics with operational limit constraints via multiple inheritance. Adds voltage bounds, reactive power bounds, apparent power thermal limits, and optional wind Q-curve constraints.

**Inherits:** `AC`, `OPF`

### `add_OPF(**kwargs)`

Extends `OPF.add_OPF()` with:

- Bus voltage bounds: `Vmin[b] ≤ v[b] ≤ Vmax[b]`
- Line apparent power limits: `pLfrom² + qLfrom² ≤ SLmax² · v²`
- Transformer apparent power limits
- Reactive power bounds for static generators, generators, and controllable loads
- Wind Q-curve constraints (Q-P and Q-U) for sgens with `var_q` set

### `add_voltage_deviation_objective()`

Sets objective to minimise sum of squared voltage deviations from 1 p.u.:

```
min  Σ_{b ∉ b0} (v[b] - 1)²  +  Σ_{b ∈ b0} (v[b] - v_b0[b])²
```

### `add_reactive_power_flow_objective()`

Sets objective to minimise total squared reactive generation:

```
min  Σ_g qsG[g]²
```

### `get_v_limits()`

Reads voltage bounds from `net.bus`. Returns `(max_vm_pu, min_vm_pu)` arrays. Generator-level limits override bus limits if stricter.

### `static_generation_wind_var_q()`

Computes Q-P and Q-U characteristic limits for wind generators based on grid code variants (0–2). Populates `q_limit_parameter` with slope/intercept values.

---

## HC_ACOPF

```
potpourri.models.HC_ACOPF.HC_ACOPF
```

Hosting Capacity AC OPF for wind generation integration studies. Extends `ACOPF` with binary variables `y[w] ∈ {0, 1}` indicating whether each wind generator is active. The default objective maximises total wind generation minus network losses.

**Inherits:** `ACOPF`

If `net.sgen` has no `wind_hc` column, candidate wind generators are automatically placed at every bus that is not an external grid bus.

### `_calc_opf_parameters(SWmax=10000, SWmin=0)`

Extends `ACOPF._calc_opf_parameters()` with apparent power bounds `SWmax` and `SWmin` for WIND_HC generators, and Q-U characteristic slopes for grid code variant 1.

**Parameters:**

- `SWmax` — maximum apparent power per wind generator (MVA)
- `SWmin` — minimum apparent power per active wind generator (MVA)

### `add_OPF(**kwargs)`

Extends `ACOPF.add_OPF()` with:

- Binary variable `y[w]` (active/inactive) for each WIND_HC generator
- Apparent power envelope: `SWmin · y ≤ S² ≤ SWmax · y`
- Q-P bounds (variant 1): `-0.41 · psG ≤ qsG ≤ 0.48 · psG`
- Q-U bounds from grid-code voltage characteristic
- Optional real power limit from `net.bus.windpot_p_mw` if present
- Default objective: maximise wind generation minus network losses

### `add_loss_obj()`

Replaces the default objective with a weighted wind-vs-loss objective:

```
max  ε · Σ psG[w]  +  (1 - ε) · (- Σ losses)
```

The mutable parameter `eps` (default 1.0) can be updated for sensitivity analysis without rebuilding the model:

```python
hc.model.eps.set_value(0.5)
```

---

## pyo_to_net

```
potpourri.models.pyo_to_net.pyo_sol_to_net_res(net, model)
```

Extracts the Pyomo solution from `model` and writes it to the pandapower result tables on `net`. Called automatically by `Basemodel.solve()` when `to_net=True`.

**Populated tables:**

- `net.res_bus` — `vm_pu`, `va_degree`, `p_mw`, `q_mvar`
- `net.res_line` — `p_from_mw`, `p_to_mw`, `q_from_mvar`, `q_to_mvar`, `i_ka`, `loading_percent`
- `net.res_trafo` — `p_hv_mw`, `p_lv_mw`, `q_hv_mvar`, `q_lv_mvar`, `loading_percent`, `tap`
- `net.res_sgen` — `p_mw`, `q_mvar`
- `net.res_gen` — `p_mw`, `q_mvar`
- `net.res_load` — `p_mw`, `q_mvar`
- `net.res_shunt` — `p_mw`, `q_mvar`

---

## init_pyo_from_pp_res

```
potpourri.models.init_pyo_from_pp_res.init_pyo_from_dcpp(net, model)
```

Warm-starts Pyomo model variables from a prior pandapower power flow result. Sets initial values for `delta[b]`, `v[b]` (AC only), `pLfrom[l]`, `pLto[l]`, `pThv[t]`, `pTlv[t]`, and generator power variables. Useful for improving IPOPT convergence on difficult AC OPF instances.
