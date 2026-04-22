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

### Optimisation studies

| Script | What it demonstrates |
|--------|---------------------|
| `acopf_loadcase_analysis.py` | Two-step AC OPF: (1) reactive-power minimisation at a fixed dispatch, (2) voltage-deviation minimisation per SimBench load case |
| `generator_capability_curve_demo.py` | How PV/wind power-factor limits and battery S² inverter circles constrain reactive dispatch; compares voltage profile and ext-grid Q between a wide-limits and a grid-code scenario |
| `objective_tradeoff_demo.py` | How four different objectives (voltage deviation, reactive generation, active import, network losses) produce different dispatch decisions on the same network |
| `constraint_activation_demo.py` | How activating each constraint group (voltage bounds, line loading, Q limits) restricts the feasible space and changes the optimal PV dispatch |
| `time_series_snapshot_opf.py` | Independent AC OPF for 8 representative seasonal/diurnal snapshots; illustrates the snapshot-OPF paradigm |
| `compute_feasible_operation_region.py` | Traces the (P, Q) feasible operation region at the grid connection point using angle-based boundary sampling |

### Multi-period planning

| Script | What it demonstrates |
|--------|---------------------|
| `multi_period_acopf.py` | 24-hour AC OPF (96 × 15 min) without storage; all time steps solved jointly as one NLP |
| `battery_multi_period_opf.py` | Multi-period AC OPF with battery storage; compares voltage deviation with and without batteries; shows SOC trajectories |
| `hosting_capacity_opf.py` | Hosting capacity analysis with binary wind placement, VDE-AR-N 4105 grid-code Q constraints, and eps/SWmin parameter sweeps |

### Validation

| Script | What it demonstrates |
|--------|---------------------|
| `validate_ac_model_against_pandapower.py` | Full AC model validation on 6 SimBench networks; reports MAE, RMSE, and max error for voltage magnitude, angle, and line losses |

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
| `compute_feasible_operation_region.py` | IPOPT |
| `multi_period_acopf.py` | IPOPT |
| `battery_multi_period_opf.py` | IPOPT |
| `hosting_capacity_opf.py` | GLPK (MindtPy); Gurobi recommended for larger runs |
| `validate_ac_model_against_pandapower.py` | IPOPT |

IPOPT and GLPK are installed automatically via `environment.yaml`.
