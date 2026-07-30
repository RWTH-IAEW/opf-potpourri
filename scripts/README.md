# scripts — potpourri examples

Each script is a self-contained, runnable example demonstrating one feature
area of the `potpourri` package.  All scripts write optional output to
`results/` (created automatically).

Run any script from the repository root after activating the environment:

```bash
conda activate potpourri_env
python scripts/minimal_ac_power_flow.py
```

---

## Script overview

### Fundamentals

| Script | What it demonstrates |
|--------|---------------------|
| `minimal_ac_power_flow.py` | End-to-end AC power flow: pandapower network → Pyomo AC model → IPOPT solve → compare against pandapower Newton-Raphson |
| `pandapower_to_pyomo_inspection.py` | How each pandapower table maps to Pyomo sets, parameters, variables, and constraints; prints a structured model summary |
| `dc_opf.py` | Linearised DC power flow and DC OPF (maximise local generation, line-loading constraint); solved with GLPK in milliseconds |
| `ehv_grid_opf.py` | Feasibility test on the large EHV/HV simbench grid (`1-EHVHV-mixed-all-0-no_sw`); runs both DC-OPF (GLPK) and AC-OPF (IPOPT) and reports solve time, voltage range, and ext-grid dispatch |

### Optimisation studies

| Script | What it demonstrates |
|--------|---------------------|
| `acopf_loadcase_analysis.py` | Two-step AC OPF: (1) reactive-power minimisation at a fixed dispatch, (2) voltage-deviation minimisation per SimBench load case |
| `generator_capability_curve_demo.py` | How PV/wind power-factor limits and battery S² inverter circles constrain reactive dispatch; compares voltage profile and ext-grid Q between a wide-limits and a grid-code scenario |
| `custom_objective_weighted_voltage.py` | How to define a **custom objective** outside the ACOPF class: attaches a per-voltage-level weighted voltage-deviation objective (`C_lv * Σ(v_lv−1)² + C_mv * Σ(v_mv−1)²`) directly to `ac.model`; compares equal-weight, LV-priority, and MV-priority scenarios on the 3-level mixed MVLV feeder `1-MVLV-urban-5.303-0-no_sw` (110 kV / 10 kV / 0.4 kV) |
| `objective_tradeoff_demo.py` | How four different objectives (voltage deviation, reactive generation, active import, network losses) produce different dispatch decisions on the same network |
| `constraint_activation_demo.py` | How activating each constraint group (voltage bounds, line loading, Q limits) restricts the feasible space and changes the optimal PV dispatch |
| `time_series_snapshot_opf.py` | Independent AC OPF for 8 representative seasonal/diurnal snapshots; illustrates the snapshot-OPF paradigm |
| `timeseries_acopf.py` | AC OPF at every step of a full-day pandapower `run_timeseries` simulation; passes `run_acopf` as the `run=` kwarg, drives loads and sgens from simbench profiles via `ConstControl`/`DFData`, and logs results with `OutputWriter` |
| `compute_feasible_operation_region.py` | Traces the (P, Q) feasible operation region at the grid connection point using angle-based boundary sampling |

### Multi-period planning

| Script | What it demonstrates |
|--------|---------------------|
| `multi_period_acopf.py` | 24-hour AC OPF (96 × 15 min) without storage; all time steps solved jointly as one NLP |
| `battery_multi_period_opf.py` | Multi-period AC OPF with battery storage; compares voltage deviation with and without batteries; shows SOC trajectories |
| `hosting_capacity_opf.py` | Hosting capacity analysis with binary wind placement, VDE-AR-N 4105 grid-code Q constraints, and eps/SWmin parameter sweeps |
| `q_control_opf.py` | **VDE-AR-N 4105 Q-control**: annotate PV/wind sgens with `var_q` and solve (1) a single-period AC OPF comparing the Q(P)/Q(U) modes, (2) the PV inverter controller modes (P(U) curtailment, fixed cos(φ), cos(φ)(P)), and (3) a 24-step multi-period AC OPF with automatic Q-control detection |
| `grid_code_q_strategies.py` | **Grid-code selection and per-sgen strategies**: solves one snapshot under each registered grid code (`add_OPF(grid_code=…)`), surfacing the provisional-values warning for VDE-AR-N 4110; then assigns Q(P)/Q(U), fixed cos(φ), cos(φ)(P) and P(U) curtailment to different PV units in one multi-period model and reports which constraint blocks were built |

### Validation and benchmarking

| Script | What it demonstrates |
|--------|---------------------|
| `validate_ac_model_against_pandapower.py` | Full AC model validation on 6 SimBench networks; reports MAE, RMSE, and max error for voltage magnitude, angle, and line losses |
| `compare_solvers.py` | Compares NEOS NLP solver backends (IPOPT, KNITRO, BONMIN, …) across 6 SimBench networks; records solve time and accuracy; saves CSV (requires NEOS email) |
| `performance_test_solver.py` | Perfplot-based scaling benchmark: solve time vs. network size across NEOS solvers; saves a PNG chart (requires `pip install potpourri[performance-test]` and NEOS email) |
| `pglib_benchmark.py` | Validates the AC- and DC-OPF formulations against the IEEE PES PGLib-OPF reference values. Loads MATPOWER `.m` cases via `potpourri.benchmarks.load_pglib_case`, applies the PGLib-compatible OPF flags (`thermal_limit='mva'`, `free_slack_vm=True`, `angle_limits=True`), and reports objective gaps. Writes `results/pglib_benchmark.{csv,md}`. |

---

## Solver requirements

| Script | Solver needed |
|--------|--------------|
| `minimal_ac_power_flow.py` | IPOPT |
| `pandapower_to_pyomo_inspection.py` | none (model inspection only) |
| `dc_opf.py` | GLPK |
| `acopf_loadcase_analysis.py` | IPOPT |
| `objective_tradeoff_demo.py` | IPOPT |
| `constraint_activation_demo.py` | IPOPT |
| `time_series_snapshot_opf.py` | IPOPT |
| `timeseries_acopf.py` | IPOPT |
| `compute_feasible_operation_region.py` | IPOPT |
| `multi_period_acopf.py` | IPOPT |
| `battery_multi_period_opf.py` | IPOPT |
| `hosting_capacity_opf.py` | GLPK (MindtPy); Gurobi recommended for larger runs |
| `q_control_opf.py` | IPOPT |
| `grid_code_q_strategies.py` | IPOPT |
| `ehv_grid_opf.py` | GLPK (DC-OPF) + IPOPT (AC-OPF) |
| `custom_objective_weighted_voltage.py` | IPOPT |
| `validate_ac_model_against_pandapower.py` | IPOPT |
| `compare_solvers.py` | NEOS (requires `NEOS_EMAIL`) |
| `performance_test_solver.py` | NEOS + `perfplot` (see `performance-test` optional dep) |
| `pglib_benchmark.py` | IPOPT; needs the `pglib-opf` benchmark cases cloned to `benchmarks/pglib-opf/` and `pip install matpowercaseframes` for parsing `.m` files |

IPOPT and GLPK are installed automatically via `environment.yaml`.
