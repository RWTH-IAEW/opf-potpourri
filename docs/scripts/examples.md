# Example Scripts Overview

All example scripts are in the `scripts/` folder and are fully self-contained.
Run any script from the repository root:

```bash
conda activate potpourri_env
python scripts/<script_name>.py
```

Optional output (plots, tables) is written to `results/`.

---

## Fundamentals

### `minimal_ac_power_flow.py`
End-to-end AC power flow pipeline.  Loads a SimBench LV network, builds a
Pyomo AC model, solves with IPOPT, and compares voltage magnitudes, angles,
and line losses side-by-side with pandapower Newton-Raphson.

### `pandapower_to_pyomo_inspection.py`
Interactive model inspector.  Prints how each pandapower table (`bus`, `line`,
`sgen`, …) maps to Pyomo sets, parameters, variables, and constraints.  Shows
variable counts (free vs. fixed), degrees of freedom, and admittance values.
No solver required.

### `dc_opf.py`
Linearised DC power flow and DC OPF.  Step 1 compares DC bus angles against
pandapower's `rundcpp`.  Step 2 minimises external grid import subject to a
line-loading constraint — solved with GLPK as a linear programme in milliseconds.

---

## Optimisation studies

### `acopf_loadcase_analysis.py`
Two-step AC OPF on a LV rural network:
1. Reactive-power minimisation at a fixed active dispatch.
2. Voltage-deviation minimisation for the SimBench `lW` and `hL` load cases.

### `generator_capability_curve_demo.py`
Shows how reactive capability limits shape AC-OPF dispatch on a 5-bus
self-contained feeder (no SimBench needed):
- PV inverter: rectangular Q bounds from a power-factor rule (cos φ ≥ pf).
- Wind inverter: tighter rectangular bounds.
- Battery: circular P² + Q² ≤ S² constraint via the built-in
  `stor_inverter_cap` Pyomo constraint.

Compares Scenario A (wide limits, pf ≥ 0.70) against Scenario B (grid-code,
PV pf ≥ 0.90 / wind pf ≥ 0.95).  Reports Q dispatch, bus voltages, ext-grid
reactive import, and the battery operating point on its capability circle.

### `objective_tradeoff_demo.py`
Compares four AC-OPF objectives (voltage deviation, reactive generation, active
import, network losses) on the same network and scenario.  Shows how the choice
of objective changes the optimal reactive dispatch.

### `constraint_activation_demo.py`
Systematically activates each constraint group (voltage bounds, line thermal
limits, reactive power limits) and shows how each one restricts the feasible
space and forces additional PV curtailment.

### `time_series_snapshot_opf.py`
Solves eight independent AC-OPFs for representative seasonal/diurnal snapshots
(winter midday, summer noon, etc.).  Illustrates the snapshot-OPF paradigm and
seasonal patterns in voltage and loading.

### `compute_feasible_operation_region.py`
Traces the boundary of the (P, Q) feasible operation region at the grid
connection point using angle-based sampling.  Useful for studying reactive
flexibility and grid hosting capacity.

---

## Multi-period planning

### `multi_period_acopf.py`
24-hour AC OPF (96 × 15-min steps) without storage.  All time steps are
optimised jointly as one large NLP.  Shows hourly ext-grid dispatch and
voltage band over the full day.

### `battery_multi_period_opf.py`
Adds `Battery_multi_period` storage to a 24-hour AC OPF.  Compares voltage
deviation with and without batteries, displays per-battery SOC trajectories,
and demonstrates explicit battery parameter sizing.

### `hosting_capacity_opf.py`
Hosting capacity analysis with binary wind placement.  Sweeps the minimum
turbine size (`SWmin`) and the wind-vs-loss trade-off parameter (`eps`).
Enforces VDE-AR-N 4105 grid-code Q constraints.

---

## Validation

### `validate_ac_model_against_pandapower.py`
Full quantitative validation on six SimBench networks (HV, MV, LV).  Reports
MAE, RMSE, and maximum absolute error for voltage magnitude, angle, active
losses, and reactive losses.
