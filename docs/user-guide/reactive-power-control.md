# Reactive-Power Control for PV and Wind

Grid codes (VDE-AR-N 4105 for LV, BDEW for MV) require distributed PV and wind
generators to regulate reactive power as a function of active output or local
voltage.  `potpourri` implements these constraints directly in the OPF
formulation so the solver explores the feasible (P, Q) region automatically.

Seven complementary constraint groups are available:

| Constraint | Symbol | Applies to |
|---|---|---|
| Q(P) characteristic | `PV_QP_pos/neg`, `sG_QP_pos/neg` | PV, wind (single- and multi-period) |
| Q(U) droop | `PV_QU_min/max`, `sG_QU_min/max` | PV, wind (single- and multi-period) |
| Inverter S² circle | `sgen_inverter_s2` | any sgen with `net.sgen.sn_mva` set |
| cos(φ) cone | `sgen_cos_phi_upper/lower` | sgens in `sGinv` (**requires the S² circle**) |
| P(U) curtailment | `sgen_pu_curtail` | any controllable PV sgen (AC only) |
| Fixed cos(φ) equality | `sgen_fixed_cos_phi` | any controllable sgen |
| cos(φ)(P) profile | `sgen_cpp` | any controllable sgen (NLP only) |

## Activation reference

Single- and multi-period models activate these constraints **differently**, and
a missing precondition is silent — no warning is raised, the constraint is
simply never added.  Check this table first when a constraint appears to have
no effect:

| Constraint | Single-period (`ACOPF.add_OPF`) | Multi-period (`ACOPF_multi_period.add_OPF`) |
|---|---|---|
| Q(P) / Q(U) | `pv_q_control="qp"` / `"qu"` / `"both"`, `net.sgen.var_q` set, **and** the sgen's `type` in `sgen_types` | automatic from `net.sgen.var_q` (no type filter) |
| Inverter S² circle | `inverter_s2=True` **and** `net.sgen.sn_mva` present | automatic from `net.sgen.sn_mva` |
| cos(φ) cone | `inverter_s2=True` **and** `sn_mva` **and** a `cos_phi_min` value | automatic from `sn_mva` **and** `net.sgen.cos_phi_min` |
| P(U) curtailment | `pu_curtail=True` | automatic from `net.sgen.pu_curtail` |
| Fixed cos(φ) | `fixed_cos_phi=…` or `net.sgen.fixed_cos_phi` | automatic from `net.sgen.fixed_cos_phi` |
| cos(φ)(P) profile | `cos_phi_p_profile=True` **and** `net.sgen.cos_phi_min` | automatic from `net.sgen.cos_phi_p_profile` |

Two consequences worth calling out:

* **`inverter_s2` defaults to `False`** in the single-period model (it is
  opt-in, so that adding `sn_mva` to a network cannot silently change an
  existing model).  The multi-period model has no such flag — it enables the
  circle as soon as `sn_mva` is present.
* **The cos(φ) cone is nested inside the S² circle.**  Setting
  `net.sgen["cos_phi_min"]` alone does nothing; `inverter_s2=True` and a
  usable `sn_mva` are both required in the single-period model.

---

## Which sgens `pv_q_control` reaches

The single-period model additionally filters by the sgen `type` column.
Matching is **exact**, against `sgen_types` (default
`("PV", "PV_MV")` — the constant `DEFAULT_PV_SGEN_TYPES`).  This matters
because SimBench names rooftop PV in LV grids `PV` but medium-voltage PV
`PV_MV`, and labels the aggregated LV renewables in MV grids `lv_RES`:

| SimBench grid | sgen types present | reached by default |
|---|---|---|
| `1-LV-rural1` | `PV` | 4 |
| `1-MV-rural` | `lv_RES`, `Wind_MV`, `Biomass_MV`, `PV_MV`, `Hydro_MV` | 2 |
| `1-MV-semiurb` | `lv_RES`, `Wind_MV`, `PV_MV`, … | 4 |
| `1-MV-urban` | `lv_RES`, `Hydro_MV` | **0** |
| `1-MV-comm` | `lv_RES`, `PV_MV`, `Wind_MV`, … | 5 |

Widen the selection when a study needs the other categories:

```python
opf.add_OPF(
    pv_q_control="both",
    sgen_types=("PV", "PV_MV", "Wind_MV", "lv_RES"),
)
```

On the MV grids above that raises the reach from 2–5 units to 87–133.
`1-MV-urban` has no `PV_MV` at all, so it needs `lv_RES` listed explicitly
before Q-control applies to anything.

!!! note "Two related paths are filtered differently"
    The **wind** Q-control path (`model.WINDc`) is selected separately and
    still matches `type == "Wind"` exactly, so SimBench's `Wind_MV` units are
    not reached through it.  The **multi-period** model applies no type filter
    at all — it keys purely off `var_q`, so it reaches every annotated sgen
    regardless of category.

---

## Selecting a grid code

The technical connection rules (TAR) are represented as selectable
`GridCode` parameter sets in `potpourri.technologies.q_control`, so a study
can target the rule that applies to its voltage level:

| Grid code | Short name | Voltage level | Status |
|---|---|---|---|
| VDE-AR-N 4105 | `"4105"` | low voltage | normative values |
| VDE-AR-N 4110 | `"4110"` | medium voltage | **provisional — placeholder values** |

A `GridCode` carries the Q(U) voltage breakpoints, the Q/Pn capability
table and its variants, the Q(P) breakpoints, and the P(U) and cos(φ)(P)
thresholds.  Select one model-wide:

```python
# Single-period
opf.add_OPF(pv_q_control="both", grid_code="4110")

# Multi-period (applies to every Q-controlled sgen in the model)
mpopf.add_OPF(grid_code="4110")
```

`grid_code` accepts `None` (VDE-AR-N 4105, the default, so existing models
are unaffected), a short name (`"4105"`, `"4110"`, or the full
`"VDE-AR-N 4110"`), or a `GridCode` instance.  An unknown name raises
`ValueError`.

!!! danger "VDE-AR-N 4110 currently holds placeholder values"
    `VDE_AR_N_4110` is wired into the registry but **reuses the VDE-AR-N
    4105 (low-voltage) parameters as a placeholder**.  Its normative
    medium-voltage figures have not been entered yet, so results obtained
    with `grid_code="4110"` are **not** 4110-compliant.

    Selecting it emits a `ProvisionalGridCodeWarning` rather than failing,
    so exploratory runs work, but do not report such results as
    medium-voltage grid-code compliant.  To complete it, replace the
    `vqu_v_points`, `vqu_q_max` and `qp_*` fields of `VDE_AR_N_4110` in
    `potpourri/technologies/q_control.py` and clear its `provisional` flag.

To add a further rule, construct a `GridCode` and register it:

```python
from potpourri.technologies.q_control import GRID_CODES, GridCode

MY_TAR = GridCode(
    name="my-tar",
    title="Operator TAR",
    voltage_level="medium voltage",
    vqu_v_points=...,   # [[V1, V2], [V3, V4]] in p.u.
    vqu_q_max=...,      # shape (2, n_variants), Q/Pn
    qp_p_high=0.1,
    qp_p_low=0.2,
    vpu_v_curtail=1.06,
    vpu_v_max=1.10,
    cpp_p_threshold_pu=0.2,
)
GRID_CODES[MY_TAR.name] = MY_TAR
```

!!! note "One grid code per model, and one path not yet covered"
    The capability curves are computed once into a single table, so the
    grid code applies model-wide; per-sgen grid codes are not supported.

    The hosting-capacity wind path in
    `potpourri/technologies/windpower.py` keeps its own private copy of the
    VDE-AR-N 4105 table and is **not** driven by the registry, so
    `grid_code` does not affect `HC_ACOPF` runs.

---

## Background: VDE-AR-N 4105 Q-control modes

### Q(P) characteristic

Reactive power is bounded as a linear function of the installed active power
capacity *P*_n:

$$Q_{\min}(P) \;\le\; Q \;\le\; Q_{\max}(P)$$

$$Q_{\max}(P) = b^+_{\rm QP} \cdot P_n + m^+_{\rm QP} \cdot P$$

The slope and intercept depend on the selected variant index `var_q` (0–2).
The coefficients are calibrated so that the envelope passes through the
grid-code capability values at the reference point *P* = 0.2 · *P*_n
(`QP_P_LOW` in `potpourri.technologies.q_control`):

| `var_q` | Q_max / P_n at 0.2·P_n (capacitive) | Q_min / P_n at 0.2·P_n (inductive) |
|---|---|---|
| 0 | +0.48 | −0.23 |
| 1 | +0.41 | −0.33 |
| 2 | +0.33 | −0.41 |

These are the two rows of `VQU_Q_MAX`; variant 0 is the widest capacitive
envelope, variant 2 the widest inductive one.  At the lower breakpoint
*P* = 0.1 · *P*_n (`QP_P_HIGH`) the envelope narrows to ±0.1 · *P*_n for
every variant.

!!! note "The Q(P) bound is not clipped above the reference point"
    `compute_q_curves()` returns a single linear segment, so the bound keeps
    widening for *P* > 0.2 · *P*_n rather than holding at the table value —
    for `var_q=0` it reaches ±3.5 · *P*_n at full output.  The Q(P)
    constraint is therefore only binding at low active power; at high output
    the effective reactive limit comes from the **inverter S² circle** and
    the **cos(φ) cone**.  Enable those (`inverter_s2=True` plus `sn_mva`)
    whenever you need a physically meaningful Q limit across the whole
    operating range.

### Q(U) droop

Reactive power is bounded as a linear function of the per-unit bus voltage *v*:

$$Q_{\min}(v) \;\le\; Q \;\le\; Q_{\max}(v)$$

$$Q_{\max}(v) = b^+_{\rm QU} + m^+_{\rm QU} \cdot v$$

The droop coefficients follow from the VDE-AR-N 4105 characteristic table
stored in `potpourri.technologies.q_control`.

---

## Inverter operating region for PV generators

The (P, Q) feasible region for a PV grid-forming inverter is the intersection
of three constraints — a "pizza slice" in the P-Q plane:

1. **P ≥ 0** — solar panels only produce active power.
2. **S² circle** — apparent power is bounded by the inverter rating
   S_inv = s_n × converter_sizing:
   $$P^2 + Q^2 \;\le\; S_{\rm inv}^2$$
3. **cos(φ) cone** — the power factor stays above a minimum value
   cos(φ_min):
   $$|Q| \;\le\; P \cdot \tan\!\bigl(\arccos(\cos\varphi_{\min})\bigr)$$

The crossover between the cone and the circle occurs at
*P*_cross = S_inv · cos(φ_min).  For *P* > *P*_cross the S² circle
is binding; for *P* < *P*_cross the cone is binding.

**Contrast with battery inverters:** battery storage operates in all four
quadrants (P and Q both positive or negative), so only the S² circle applies.
The circle for batteries is already added automatically via
`Basemodel.add_storage()`.

---

## P(U) active-power curtailment (VDE-AR-N 4105 §8.5)

When bus voltage exceeds a threshold, the inverter reduces active output
linearly to zero:

$$P \;\le\; P_n \cdot \frac{V_{\max} - v[b]}{V_{\max} - V_{\rm curtail}}$$

Default thresholds: *V*_curtail = 1.06 p.u., *V*_max = 1.10 p.u.

The constraint is bilinear in P and v[b] — requires an NLP solver (IPOPT).

---

## Fixed cos(φ) mode

An equality constraint that fixes the reactive-to-active ratio at every
operating point:

$$Q[g] = P[g] \cdot \tan(\arccos(\cos\varphi))$$

Unlike the cos(φ) cone (which is a bound), this is an equality: the
inverter tracks the prescribed power factor exactly.

---

## cos(φ)(P) profile (VDE-AR-N 4105)

A piecewise P-Q curve: no reactive power below a threshold *P*_t, then Q
increases with P up to the target power factor at full output:

$$Q \cdot (P_n - P_t) = \tan(\varphi) \cdot P \cdot (P - P_t)$$

Implemented as a quadratic equality (NLP). Q → 0 at P = 0 or P = P_t;
Q_max = P_n · tan(φ) at P = P_n.

Default threshold: *P*_t = 0.2 · P_n.

---

## Single-period usage

### Enabling Q-control modes

Annotate the network sgens with `var_q` (variant index) and `p_inst_mw`
(installed capacity) before building the model:

```python
import simbench as sb
from potpourri.models.ACOPF_base import ACOPF

net = sb.get_simbench_net("1-LV-rural1--0-sw")

# Mark PV sgens with Q-control variant 0 (Qmax = 0.48 Pn)
net.sgen["var_q"] = None       # object-dtype column; None for non-PV rows
mask = net.sgen["type"] == "PV"
net.sgen.loc[mask, "var_q"] = 0
net.sgen["p_inst_mw"] = net.sgen["p_mw"].abs()   # installed capacity

net.sgen["controllable"] = True
net.bus["max_vm_pu"] = 1.05
net.bus["min_vm_pu"] = 0.95
```

Pass `pv_q_control` to `add_OPF()` to select the control mode:

```python
opf = ACOPF(net)

# Mode options:
#   None      — no Q-control (default)
#   "qp"      — Q(P) characteristic only
#   "qu"      — Q(U) droop only
#   "both"    — Q(P) + Q(U) combined  (also accepts True for backward compat.)
opf.add_OPF(pv_q_control="both")
opf.add_voltage_deviation_objective()
opf.solve(solver="ipopt")
```

### Enabling the inverter S² circle

Set `net.sgen.sn_mva` (apparent-power rating) and optionally
`net.sgen.converter_sizing_pu` (default 1.0) before constructing the model,
then pass `inverter_s2=True` explicitly — the flag is opt-in and defaults to
`False`:

```python
net.sgen["sn_mva"] = net.sgen["p_mw"].abs() / 0.9   # 90 % power factor rating
net.sgen["converter_sizing_pu"] = 1.0                # no derating

opf = ACOPF(net)
opf.add_OPF(inverter_s2=True)
```

### Enabling P(U) curtailment

Set `net.sgen["p_inst_mw"]` and optionally per-sgen voltage thresholds,
then pass `pu_curtail=True`:

```python
net.sgen["p_inst_mw"] = net.sgen["p_mw"].abs()
net.sgen["v_curtail_pu"] = 1.06      # optional; defaults to 1.06
net.sgen["v_max_curtail_pu"] = 1.10  # optional; defaults to 1.10

opf = ACOPF(net)
opf.add_OPF(pu_curtail=True)         # applies to PV-type sgens in sGc
```

### Enabling fixed cos(φ) mode

```python
# Scalar: same power factor for all controllable sgens
opf.add_OPF(fixed_cos_phi=0.95)

# Per-sgen: set net.sgen["fixed_cos_phi"] column (takes precedence)
net.sgen["fixed_cos_phi"] = 0.95
opf.add_OPF()
```

### Enabling the cos(φ)(P) profile

```python
net.sgen["cos_phi_min"] = 0.9        # power factor at full output
net.sgen["p_inst_mw"] = net.sgen["p_mw"].abs()
# net.sgen["cpp_p_threshold_pu"] = 0.2  # optional; default 0.2

opf.add_OPF(cos_phi_p_profile=True)
```

### Enabling the cos(φ) cone

The cone is added *inside* the S² circle block, so `inverter_s2=True` and a
usable `net.sgen.sn_mva` are required in **both** forms below.  Setting
`cos_phi_min` without `inverter_s2=True` adds no constraint and raises no
warning.

```python
net.sgen["sn_mva"] = net.sgen["p_mw"].abs() / 0.9   # required

# Per-sgen (different limits per generator) — column takes precedence
net.sgen["cos_phi_min"] = 0.9
opf.add_OPF(inverter_s2=True)

# Or a scalar for every sgen in sGinv, when no column is present
opf.add_OPF(inverter_s2=True, cos_phi_min=0.9)
```

Both forms create `model.sgen_cos_phi_upper/lower` over set `model.sGpf`.  To
confirm the cone was actually built:

```python
assert hasattr(opf.model, "sgen_cos_phi_upper")
```

In the multi-period model no flag is needed — `sn_mva` plus `cos_phi_min` in
`net.sgen` is sufficient.

---

## OLTC tap optimisation (multi-period)

The multi-period model supports continuous OLTC optimisation via
`add_tap_changer_linear()`, which unfixes the time-indexed `Tap[tr, t]`
variable and adds per-step bounds:

```python
mpopf = ACOPF_multi_period(net, toT=96)
mpopf.add_OPF()

# Continuous OLTC — free within [tap_min, tap_max] each step
mpopf.add_tap_changer_linear()

# With optional inter-step rate limit (Δtap ≤ 0.01 per 15 min)
mpopf.add_tap_changer_linear(max_tap_change_per_step=0.01)

mpopf.add_voltage_deviation_objective()
mpopf.solve(solver="ipopt")
```

For discrete tap positions (MIP solver required):

```python
mpopf.add_tap_changer_discrete()
mpopf.solve(solver="mindtpy", mip_solver="gurobi")
```

---

## Multi-period usage

In `ACOPF_multi_period` the Q-control constraints are added automatically from
the network data — no extra keyword argument is needed.

```python
import copy
import simbench as sb
from potpourri.models_multi_period.ACOPF_multi_period import ACOPF_multi_period

net = sb.get_simbench_net("1-LV-rural1--0-sw")

# Annotate sgens as above (var_q, p_inst_mw, sn_mva, cos_phi_min)
net.sgen["var_q"] = None
mask = net.sgen["type"] == "PV"
net.sgen.loc[mask, "var_q"] = 0
net.sgen["p_inst_mw"] = net.sgen["p_mw"].abs()
net.sgen["sn_mva"] = net.sgen["p_mw"].abs() / 0.9
net.sgen["cos_phi_min"] = 0.9

net.sgen["max_p_mw"] = net.sgen["p_mw"]
net.sgen["min_p_mw"] = 0.0
net.bus["max_vm_pu"] = 1.05
net.bus["min_vm_pu"] = 0.95

mpopf = ACOPF_multi_period(net, toT=13848, fromT=13824)
mpopf.add_OPF()
mpopf.add_voltage_deviation_objective()
mpopf.solve(solver="ipopt")
```

`_calc_opf_parameters()` detects `var_q` in `net.sgen` and calls
`static_generation_q_ctrl_data()` and `static_generation_inverter_data()`
automatically.

### Selecting an inverter controller mode

The three inverter controller modes are **alternatives** — pick one per study.
Each is activated by its own `net.sgen` column:

```python
# (a) P(U) curtailment — boolean column
net.sgen.loc[mask, "pu_curtail"] = True

# (b) Fixed cos(φ) — per-sgen float column
net.sgen.loc[mask, "fixed_cos_phi"] = 0.95

# (c) cos(φ)(P) profile — boolean column + cos_phi_min
net.sgen.loc[mask, "cos_phi_p_profile"] = True
net.sgen.loc[mask, "cos_phi_min"] = 0.9
```

!!! warning "Do not stack the equality modes"
    Fixed cos(φ) (b) and the cos(φ)(P) profile (c) are both **equality**
    constraints on the same `qsG[g, t]`.  Setting both columns on the same
    sgens builds `sgen_fixed_cos_phi` *and* `sgen_cpp` simultaneously, which
    over-determines reactive power and typically leaves the model infeasible
    or numerically degenerate.  Combining either equality mode with the Q(P)
    envelope is likewise over-constrained, since Q is then pinned by the
    equality and can no longer move within the grid-code band.

    `scripts/q_control_opf.py` shows the intended pattern: solve each mode
    against its own copy of the network and compare the results.

Adding P(U) curtailment (a) *alongside* Q-control is fine — it constrains
active rather than reactive power.

---

## Pyomo components reference

### Single-period (ACOPF)

| Component | Type | Description |
|---|---|---|
| `model.PVc` | Set | Controllable PV sgen indices |
| `model.PV_var_q[g]` | Param | Q-control variant (0–2) |
| `model.PV_p_inst[g]` | Param | Installed capacity P_n (p.u.) |
| `model.qPV[g]` | Var | PV reactive dispatch (p.u.) |
| `model.PV_QP_pos[g]` | Constraint | Q ≤ Q_max(P) upper bound |
| `model.PV_QP_neg[g]` | Constraint | Q ≥ Q_min(P) lower bound |
| `model.PV_QU_min[g]` | Constraint | Q ≥ Q_min(v) lower bound |
| `model.PV_QU_max[g]` | Constraint | Q ≤ Q_max(v) upper bound |
| `model.sGinv` | Set | sgen indices with inverter ratings |
| `model.S_inv[g]` | Param | Inverter apparent-power limit S_inv (p.u.) |
| `model.sgen_inverter_s2[g]` | Constraint | P² + Q² ≤ S_inv² |
| `model.sGpf` | Set | sgen indices with cos(φ) limits |
| `model.tan_phi[g]` | Param | tan(arccos(cos_phi_min)) |
| `model.sgen_cos_phi_upper[g]` | Constraint | Q ≤ P · tan(φ) |
| `model.sgen_cos_phi_lower[g]` | Constraint | Q ≥ −P · tan(φ) |
| `model.sGpu` | Set | PV sgens with P(U) curtailment |
| `model.P_inst_pu[g]` | Param | Installed capacity P_n (p.u.) |
| `model.V_curtail[g]` | Param | Curtailment threshold voltage (p.u.) |
| `model.V_max_curtail[g]` | Param | Zero-output voltage (p.u.) |
| `model.sgen_pu_curtail[g]` | Constraint | P · ΔV ≤ P_n · (V_max − v[b]) |
| `model.sGfcf` | Set | sgen indices with fixed cos(φ) |
| `model.fixed_tan_phi[g]` | Param | tan(arccos(fixed_cos_phi)) |
| `model.sgen_fixed_cos_phi[g]` | Constraint | Q == tan_phi · P |
| `model.sGcpp` | Set | sgen indices with cos(φ)(P) profile |
| `model.cpp_tan_phi[g]` | Param | tan_phi at full output |
| `model.cpp_Pn[g]` | Param | Installed capacity P_n (p.u.) |
| `model.cpp_P_thresh[g]` | Param | Threshold P_t = P_thresh_pu · P_n |
| `model.sgen_cpp[g]` | Constraint | Q·(Pn−Pt) == tan_phi·P·(P−Pt) |

### Multi-period (ACOPF_multi_period via Sgens_multi_period)

The time-indexed variants use `[g, t]` indices:

| Component | Description |
|---|---|
| `model.sGqc` | Set of Q-controlled sgen indices |
| `model.sG_QP_pos[g, t]` | Q ≤ Q_max(P) upper bound |
| `model.sG_QP_neg[g, t]` | Q ≥ Q_min(P) lower bound |
| `model.sG_QU_min[g, t]` | Q ≥ Q_min(v) lower bound |
| `model.sG_QU_max[g, t]` | Q ≤ Q_max(v) upper bound |
| `model.sgen_inverter_s2[g, t]` | P² + Q² ≤ S_inv² |
| `model.sgen_cos_phi_upper[g, t]` | Q ≤ P · tan(φ) |
| `model.sgen_cos_phi_lower[g, t]` | Q ≥ −P · tan(φ) |
| `model.sgen_pu_curtail[g, t]` | P · ΔV ≤ P_n · (V_max − v[b,t]) |
| `model.sgen_fixed_cos_phi[g, t]` | Q == tan_phi · P |
| `model.sgen_cpp[g, t]` | Q·(Pn−Pt) == tan_phi·P·(P−Pt) |

---

## Assigning a different strategy per sgen

The grid code applies model-wide, but *which controller each sgen follows* is
per-row.  In the multi-period model every strategy is driven by a `net.sgen`
column, so units on the same feeder can follow different rules:

```python
net.sgen["var_q"] = None
net.sgen["fixed_cos_phi"] = float("nan")
net.sgen["cos_phi_p_profile"] = False
net.sgen["cos_phi_min"] = float("nan")
net.sgen["pu_curtail"] = False

net.sgen.at[0, "var_q"] = 0                  # Q(P) + Q(U) envelope
net.sgen.at[1, "fixed_cos_phi"] = 0.95       # fixed power factor
net.sgen.at[2, "cos_phi_p_profile"] = True   # cos(φ)(P) profile
net.sgen.at[2, "cos_phi_min"] = 0.90
net.sgen.at[3, "pu_curtail"] = True          # P(U) curtailment
net.sgen.at[3, "var_q"] = 0                  # …may combine with a Q rule

mpopf = ACOPF_multi_period(net, toT=24)
mpopf.add_OPF(grid_code="4105")
```

Two rules for combining strategies on one sgen:

* **Fixed cos(φ) and the cos(φ)(P) profile must not share an sgen.** Both are
  equality constraints on the same `qsG`, so together they over-determine
  reactive power and the model is typically infeasible or degenerate. Keep
  them on disjoint sets.
* **P(U) curtailment may be combined with a Q rule**, because it constrains
  active rather than reactive power.

!!! note "Per-row assignment is a multi-period feature"
    In the single-period model `pu_curtail` and `cos_phi_p_profile` are
    model-wide switches on `add_OPF()`, so only `fixed_cos_phi` and `var_q`
    can vary per row there. Use the multi-period model when the study needs
    genuinely mixed controllers.

### What a mixed assignment looks like

A fixed power factor couples Q rigidly to P. When the feeder already sits
above 1.0 p.u. and the objective penalises voltage deviation, the only way
for such a unit to shed reactive power is to shed active power — so it
curtails, while its Q(P)/Q(U) neighbours keep producing. That coupling is the
practical cost of a fixed power factor, and it is the usual reason to prefer
the bound-type rules when active yield matters.

---

## Example scripts

`scripts/grid_code_q_strategies.py` covers this page's two selection
mechanisms:

1. **Grid-code selection** — solves one snapshot under every registered grid
   code and reports the capability envelope, objective and voltage band for
   each, surfacing the `ProvisionalGridCodeWarning` for VDE-AR-N 4110 rather
   than silencing it.
2. **One strategy per sgen** — assigns Q(P)/Q(U), fixed cos(φ), cos(φ)(P) and
   P(U) curtailment to different PV units, reports which constraint blocks
   were built, and prints the resulting per-sgen dispatch.

`scripts/q_control_opf.py` demonstrates both use cases:

1. **Single-period comparison** — solves four modes side-by-side
   (uncontrolled, Q(P) only, Q(U) only, Q(P)+Q(U)) and prints an objective /
   voltage-band summary and a per-sgen Q dispatch table with Q(P) bound
   verification.

2. **Multi-period 24-step run** — shows how `ACOPF_multi_period` picks up
   `var_q` automatically and reports a P/Q time series for the first
   Q-controlled sgen.

See [Scripts reference](../scripts/examples.md) for the full list of
runnable examples.
