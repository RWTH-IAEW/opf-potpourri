# Class Architecture

This page describes how `potpourri`'s classes relate to each other, what Pyomo
components each layer contributes, and how the model state changes as you
compose or extend classes.  All power quantities live in the **per-unit system**
on `net.sn_mva` (=`model.baseMVA`).

---

## Single-Period Models

### Inheritance diagram

```mermaid
classDiagram
    direction TB

    class Basemodel {
        +net : pandapowerNet
        +model : ConcreteModel
        +results : SolverResults
        +baseMVA : float
        +bus_data, line_data, trafo_data
        ── Sets ──
        +B  b0  bPV
        +G  sG  eG  gG  D  L  TRANSF  SHUNT
        +Dbs  Gbs  sGbs
        ── Parameters ──
        +A  AT  baseMVA  shift  delta_b0
        +PD  PG  PsG  GB
        ── Variables ──
        +delta[b]  pG[g]  psG[sg]  pD[d]
        +pLfrom[l]  pLto[l]  pThv[tr]  pTlv[tr]  Tap[tr]
        +solve(solver)
        +add_storage()
        +change_vals()  fix_vars()  unfix_vars()
    }

    class AC {
        <<power flow — AC>>
        ── Parameters ──
        +Gii[l]  Bii[l]  Gik[l]  Bik[l]
        +GiiT[tr]  BiiT[tr]  GikT[tr]  BikT[tr]
        +BB[shunt]  v_b0  QsG[sg]  QD[d]  v_bPV
        ── Variables ──
        +v[b]  qG[g]  qsG[sg]  qD[d]
        +qLfrom[l]  qLto[l]  qThv[tr]  qTlv[tr]
        ── Constraints ──
        +KCL_real[b]  KCL_reactive[b]
        +KVL_real_from[l]  KVL_real_to[l]
        +KVL_reactive_from[l]  KVL_reactive_to[l]
        +KVL_*Transf[tr]  v_bPV_setpoint[b]
    }

    class DC {
        <<power flow — DC linearised>>
        ── Parameters ──
        +BL[l]  BLT[tr]
        ── Variables ──
        +deltaL[l]  deltaLT[tr]
        ── Constraints ──
        +KCL_const[b]
        +KVL_real_from[l]  KVL_real_to[l]
        +KVL_trans_from[tr]  KVL_trans_to[tr]
        +phase_diff1[l]  phase_diff2[tr]
    }

    class OPF {
        <<mixin — operational limits>>
        ── Sets ──
        +sGc  Dc
        ── Parameters ──
        +sPGmax[sg]  sPGmin[sg]
        +PGmax[g]  PGmin[g]
        +PDmax[d]  PDmin[d]
        +SLmax[l]  SLmaxT[tr]
        ── Constraints ──
        +PsG_Constraint[sg]
        +PG_Constraint[g]
        +PD_Constraint[d]
        +add_OPF()
        +add_tap_changer_linear()
        +add_tap_changer_discrete()
    }

    class ACOPF {
        ── Sets ──
        +WIND  WINDc  sGnc  Bfix  WIND_HC
        ── Parameters ──
        +Vmax[b]  Vmin[b]
        +QGmax[g]  QGmin[g]
        +QsGmax[sg]  QsGmin[sg]  (mutable)
        +QDmax[d]  QDmin[d]
        +var_q[w]  PsG_inst[w]
        ── Constraints ──
        +line_lim_from[l]  line_lim_to[l]
        +transf_lim1[tr]  transf_lim2[tr]
        +v_pyo[b]
        +QsG_pyo[sg]  QG_pyo[g]  QD_pyo[d]
        +QW_pos_pyo[w]  QW_neg_pyo[w]
        +QU_min_pyo[w]  QU_max_pyo[w]
        ── Objectives (add separately) ──
        +obj_v_deviation
        +obj_reactive
        +add_OPF()
        +add_voltage_deviation_objective()
        +add_reactive_power_flow_objective()
        +add_active_change_objective()
    }

    class DCOPF {
        ── Constraints ──
        +line_lim_from[l]  line_lim_to[l]
        +transf_lim1[tr]  transf_lim2[tr]
        +add_OPF()
    }

    class HC_ACOPF {
        ── Sets ──
        +WIND_HC (wind HC generators)
        ── Parameters ──
        +SWmax[w]  SWmin[w]  (mutable)
        +eps  (mutable weight)
        ── Variables ──
        +y[w] : Binary
        ── Constraints ──
        +SW_max_constraint[w]
        +SW_min_constraint[w]
        +QW_min_constraint[w]  QW_max_constraint[w]
        +QU_min_hc_constraint[w]  QU_max_hc_constraint[w]
        ── Objective ──
        +obj (maximise hosted wind generation)
        +add_OPF()
        +add_loss_obj()
    }

    Basemodel <|-- AC : extends
    Basemodel <|-- DC : extends
    Basemodel <|-- OPF : extends
    AC        <|-- ACOPF : multiple inheritance
    OPF       <|-- ACOPF : multiple inheritance
    DC        <|-- DCOPF : multiple inheritance
    OPF       <|-- DCOPF : multiple inheritance
    ACOPF     <|-- HC_ACOPF : extends
```

> **Multiple inheritance** — `ACOPF` inherits from both `AC` and `OPF` using Python's
> MRO.  `Basemodel` is the common root so `model` and `net` are never duplicated.
> The same pattern applies to `DCOPF(DC, OPF)`.

---

### Pyomo components at each layer

The table below lists **every Pyomo component** (Set, Param, Var, Constraint,
Objective) that is introduced at each class level.  Lower classes inherit
everything listed above them.

| Layer | Sets | Parameters | Variables | Constraints | Objectives |
|---|---|---|---|---|---|
| **Basemodel** | `B` `b0` `bPV` `G` `sG` `eG` `gG` `D` `L` `TRANSF` `SHUNT` `STOR` | `A` `AT` `baseMVA` `shift` `delta_b0` `PD` `PG` `PsG` `GB` | `delta[b]` `pG[g]` `psG[sg]` `pD[d]` `pLfrom[l]` `pLto[l]` `pThv[tr]` `pTlv[tr]` `Tap[tr]` | — | — |
| **+AC** | — | `Gii` `Bii` `Gik` `Bik` `GiiT` `BiiT` `GikT` `BikT` `BB` `QD` `QsG` `v_b0` `v_bPV` | `v[b]` `qG[g]` `qsG[sg]` `qD[d]` `qLfrom[l]` `qLto[l]` `qThv[tr]` `qTlv[tr]` | `KCL_real[b]` `KCL_reactive[b]` `KVL_real_from/to[l]` `KVL_reactive_from/to[l]` `KVL_*Transf[tr]` `v_bPV_setpoint[b]` | — |
| **+DC** | — | `BL[l]` `BLT[tr]` | `deltaL[l]` `deltaLT[tr]` | `KCL_const[b]` `KVL_real_from/to[l]` `KVL_trans_from/to[tr]` `phase_diff1[l]` `phase_diff2[tr]` | — |
| **+OPF** | `sGc` `Dc` | `sPGmax` `sPGmin` `PGmax` `PGmin` `PDmax` `PDmin` `SLmax[l]` `SLmaxT[tr]` | — | `PsG_Constraint[sg]` `PG_Constraint[g]` `PD_Constraint[d]` | — |
| **ACOPF = AC+OPF** | `WIND` `WINDc` `sGnc` `Bfix` `WIND_HC` | `Vmax[b]` `Vmin[b]` `QGmax` `QGmin` `QsGmax` `QsGmin` `QDmax` `QDmin` `var_q` `PsG_inst` | — | `line_lim_from/to[l]` `transf_lim1/2[tr]` `v_pyo[b]` `QsG_pyo[sg]` `QG_pyo[g]` `QD_pyo[d]` `QW_*_pyo[w]` `QU_*_pyo[w]` | `obj_v_deviation` `obj_reactive` |
| **DCOPF = DC+OPF** | — | — | — | `line_lim_from/to[l]` `transf_lim1/2[tr]` | — |
| **HC_ACOPF ⊂ ACOPF** | `WIND_HC` | `SWmax[w]` `SWmin[w]` `eps` | `y[w]` (binary) | `SW_max/min_constraint[w]` `QW_min/max_constraint[w]` `QU_min/max_hc_constraint[w]` | `obj` (maximise HC) |

---

### Choosing between single-period models

| Goal | Class to use | Key components gained |
|---|---|---|
| Power flow only, full AC equations | `AC` | `v[b]`, `KCL_real/reactive[b]`, `KVL_*[l]` |
| Power flow only, linear DC | `DC` | `deltaL[l]`, `KCL_const[b]` — **no** `v`, **no** reactive |
| AC OPF (voltage & thermal limits) | `ACOPF` | everything in AC + `v_pyo[b]`, `line_lim_from/to[l]`, `QsG_pyo[sg]` + objectives |
| Fast screening, DC OPF | `DCOPF` | everything in DC + `line_lim_from/to[l]` — **no** voltage variables |
| Hosting-capacity study | `HC_ACOPF` | everything in ACOPF + binary `y[w]` + `SW_*_constraint[w]` |

> **Note on the DC model** — The DC approximation fixes all voltage magnitudes
> to 1.0 p.u. implicitly; `Basemodel` does *not* declare a `v` variable for
> `DC` or `DCOPF` models.  Checking `hasattr(model, 'v')` returns `False`.

---

## Multi-Period Models

### Inheritance diagram

```mermaid
classDiagram
    direction TB

    class Basemodel_multi_period {
        +net : pandapowerNet
        +model : ConcreteModel
        +results : SolverResults
        +fromT  toT  T : int
        +profiles : dict
        +flexibilities : list
        ── Sets ──
        +T  B  b0  bPV  L  TRANSF  LE
        ── Parameters ──
        +deltaT  A  AT  shift  delta_b0  baseMVA
        ── Variables ──
        +delta[b,t]
        +pLfrom[l,t]  pLto[l,t]
        +pThv[tr,t]   pTlv[tr,t]  Tap[tr,t]
        +solve(solver)
        +make_to_dict()  make_to_tuple()
        +change_vals()  fix_vars()  unfix_vars()
    }

    class AC_multi_period {
        <<power flow — AC, time-indexed>>
        ── Parameters ──
        +Gii[l]  Bii[l]  Gik[l]  Bik[l]
        +GiiT[tr]  BiiT[tr]  GikT[tr]  BikT[tr]
        +BB[shunt]  v_b0
        ── Variables ──
        +v[b,t]  qG[g,t]
        +qLfrom[l,t]  qLto[l,t]
        +qThv[tr,t]   qTlv[tr,t]
        ── Constraints ──
        +KCL_real[b,t]  KCL_reactive[b,t]
        +KVL_real_from[l,t]  KVL_real_to[l,t]
        +KVL_reactive_from[l,t]  KVL_reactive_to[l,t]
        +KVL_*Transf[tr,t]
    }

    class OPF_multi_period {
        <<mixin — operational limits>>
        ── Parameters ──
        +SLmax[l]  SLmaxT[tr]  (mutable, time-independent)
        +add_OPF()
        +add_tap_changer_linear()
        +add_tap_changer_discrete()
    }

    class ACOPF_multi_period {
        +v_limits
        ── Parameters ──
        +Vmax[b]  Vmin[b]  (mutable, time-independent)
        ── Constraints ──
        +line_lim_from[l,t]  line_lim_to[l,t]
        +transf_lim1[tr,t]   transf_lim2[tr,t]
        +v_constraint[b,t]
        ── Objectives (add separately) ──
        +obj_v_deviation
        +Objective (demand minimisation)
        +obj  (generation minimisation)
        +add_OPF()
        +add_voltage_deviation_objective()
        +add_minimize_power_objective()
        +add_generation_objective()
        +add_weighted_generation_objective()
    }

    Basemodel_multi_period <|-- AC_multi_period  : extends
    Basemodel_multi_period <|-- OPF_multi_period : extends
    AC_multi_period  <|-- ACOPF_multi_period : multiple inheritance
    OPF_multi_period <|-- ACOPF_multi_period : multiple inheritance
```

---

### Key differences from single-period

| Aspect | Single-period | Multi-period |
|---|---|---|
| Variable index | `v[b]` | `v[b, t]` — bus × time |
| Constraint index | `KCL_real[b]` | `KCL_real[b, t]` — bus × time |
| Thermal limit | `line_lim_from[l]` | `line_lim_from[l, t]` |
| Voltage bound | `v_pyo[b]` | `v_constraint[b, t]` |
| Voltage/line params | `Vmax[b]`, `SLmax[l]` (per bus/line) | same — time-independent |
| Bus/load/gen data | fixed scalars | loaded from SimBench profiles per `t` |
| Flexible devices | `add_storage()` on Basemodel | separate `*_multi_period` objects via `get_all(model)` |
| Time set | — | `model.T` = `range(fromT, toT)` |

> **Profile data** — `Basemodel_multi_period` reads SimBench time-series profiles
> and passes them to flexibility objects.  Each device module stamps its own
> time-indexed parameters (e.g. `PD[d, t]`, `P_PV[pv, t]`) onto `model`.

---

## Flexible Device Composition

Flexible devices in multi-period models are **not** class parents; they are
separate objects that attach their own Pyomo components to an existing model
instance via `get_all(model)`.

### Composition diagram

```mermaid
classDiagram
    direction LR

    class ACOPF_multi_period {
        +model : ConcreteModel
        +flexibilities : list
    }

    class Flexibility_multi_period {
        <<base>>
        +net  baseMVA
        +buses_excl_extGrids
        +PD_data  QD_data  GB_data
        +get_sets(model)
        +make_to_dict(...)
    }

    class Battery_multi_period {
        +bat_power  bat_cap : float
        +bat_efficiency  bat_percentage
        +random_indexes : ndarray
        ── Sets ──
        +BAT  BAT_bus
        ── Parameters ──
        +BAT_Pmax[b]  BAT_SOCmax[b]  BAT_SOCmin[b]
        +BAT_Cap[b]   BAT_Eff[b]     BAT_SOC_init[b]
        ── Variables ──
        +BAT_P[b,t]    BAT_SOC[b,t]
        ── Constraints ──
        +bat_power_con[b,t]
        +bat_soc_con[b,t]
        +bat_soc_update_con[b,t]
        +get_all(model)
    }

    class EVs_multi_period {
        ── Sets: EV  EV_bus ──
        ── Vars: EV_P[e,t]  EV_SOC[e,t] ──
        ── Constrs: ev_power_con  ev_soc_update_con ──
        +get_all(model)
    }

    class HeatPump_multi_period {
        ── Sets: HP  HP_bus ──
        ── Vars: HP_P[h,t] ──
        ── Constrs: hp_power_con[h,t] ──
        +get_all(model)
    }

    class PV_multi_period {
        ── Sets: PV  PV_bus ──
        ── Vars: PV_P[pv,t] ──
        ── Constrs: pv_power_con[pv,t] ──
        +get_all(model)
    }

    class Windpower_multi_period {
        ── Sets: WIND  WIND_bus ──
        ── Vars: WIND_P[w,t] ──
        +get_all(model)
    }

    class Demand_multi_period {
        ── Sets: D (re-declared)  ──
        ── Params: PD[d,t]  QD[d,t] ──
        +get_all(model)
    }

    class Sgens_multi_period {
        ── Sets: sG  sGc ──
        ── Params: PsG[sg,t]  QsG[sg,t] ──
        +get_all(model)
    }

    class Generator_multi_period {
        ── Sets: G  eG  gG ──
        ── Params: PG[g,t]  QG[g,t] ──
        ── Vars: pG[g,t]  qG[g,t] ──
        +get_all(model)
    }

    Flexibility_multi_period <|-- Battery_multi_period
    Flexibility_multi_period <|-- EVs_multi_period
    Flexibility_multi_period <|-- HeatPump_multi_period
    Flexibility_multi_period <|-- PV_multi_period
    Flexibility_multi_period <|-- Windpower_multi_period
    Flexibility_multi_period <|-- Demand_multi_period
    Flexibility_multi_period <|-- Sgens_multi_period
    Flexibility_multi_period <|-- Generator_multi_period
    ACOPF_multi_period "1" *-- "0..*" Flexibility_multi_period : get_all(model)
```

### How composition works

Flexible devices are constructed first (they receive the `net` and number of
time steps), then attached to the model before `add_OPF()` is called:

```python
opf = ACOPF_multi_period(net, toT=96)       # builds model.T, model.L, model.delta…

battery = Battery_multi_period(opf.net, T=96, scenario=1)
battery.get_all(opf.model)                  # stamps BAT, BAT_P, bat_*_con onto model

opf.add_OPF()                               # reads BAT sets for KCL battery injection
opf.add_voltage_deviation_objective()
opf.solve(solver="ipopt")
```

`add_OPF()` calls `get_all_opf(model)` on each device in `flexibilities`, so
the KCL constraints automatically include device injections.

### Available device modules

| Module | Class | Adds to model |
|---|---|---|
| `battery.py` | `Battery_multi_period` | `BAT`, `BAT_P[b,t]`, `BAT_SOC[b,t]`, SOC dynamics |
| `EVs.py` | `EVs_multi_period` | `EV`, `EV_P[e,t]`, `EV_SOC[e,t]`, charging constraints |
| `heatpump.py` | `HeatPump_multi_period` | `HP`, `HP_P[h,t]`, thermal power limits |
| `PV.py` | `PV_multi_period` | `PV`, `PV_P[pv,t]`, irradiance-based upper bound |
| `windpower.py` | `Windpower_multi_period` | `WIND`, `WIND_P[w,t]`, wind-speed-based limit |
| `demand.py` | `Demand_multi_period` | `D`, time-varying `PD[d,t]` `QD[d,t]` from profiles |
| `sgens.py` | `Sgens_multi_period` | `sG` `sGc`, time-varying `PsG[sg,t]` from profiles |
| `generator.py` | `Generator_multi_period` | `G` `eG` `gG`, `pG[g,t]` `qG[g,t]` decision variables |

---

## Model lifecycle

Every model (single- and multi-period) follows the same five-step lifecycle:

```
net ──► __init__()       builds model.B, model.L, power-flow eqs
         │
         ▼
        add_OPF()        stamps voltage/thermal/Q limits and generator bounds
         │
         ▼
        add_*_objective() activates exactly one Pyomo Objective
         │
         ▼
        solve(solver)    calls Pyomo SolverFactory; writes self.results
         │
         ▼
        pyo_to_net()     reads Pyomo solution back into net.res_bus,
                         net.res_line, net.res_sgen, …
```

`solve()` calls `pyo_to_net` automatically when `to_net=True` (the default),
so `net.res_bus.vm_pu` is always populated after a successful solve.

### Storage (single-period)

The single-period `Basemodel` supports an optional storage extension:

```python
model = ACOPF(net)
model.add_storage()      # adds STOR, STOR_Pchg[s], STOR_Pdis[s], STOR_SOC[s]
model.add_OPF()
```

`add_storage()` must be called *before* `add_OPF()` so that the storage
injections are included in the KCL constraints.
