# PGLib-OPF Benchmark Validation

**Script:** `scripts/pglib_benchmark.py`

Benchmarks `potpourri` DC-OPF and AC-OPF against the published reference
objective values from the
[IEEE PES Power Grid Library (PGLib-OPF)](https://github.com/power-grid-lib/pglib-opf),
solved by [PowerModels.jl](https://github.com/lanl-ansi/PowerModels.jl) with
IPOPT.

---

## Purpose

PGLib-OPF is the standard benchmark suite for comparing OPF solvers and
formulations. Each case ships as a MATPOWER `.m` file with a published
reference objective value (minimum generation cost in $/h) for DC and AC
OPF under three scenario groups:

| Group | Abbreviation | Description |
|---|---|---|
| Typical Operations | TYP | Near-nominal loading, baseline scenario |
| Congested Operations | API | Elevated loading, active constraints |
| Small Angle Difference | SAD | Tight phase-angle bounds |

The benchmark script runs `potpourri` DC and AC OPF on a configurable subset
of TYP cases and reports the objective gap relative to the PowerModels.jl
reference.

---

## Setup

The PGLib-OPF case files are included as a git submodule. Initialise it once:

```bash
git submodule update --init benchmarks/pglib-opf
```

Or clone with submodules from the start:

```bash
git clone --recurse-submodules https://github.com/RWTH-IAEW/opf-potpourri.git
```

---

## Workflow

```
for each selected PGLib case (.m file)
    ┌─ load via load_pglib_case()
    │    (from_mpc → pandapower net, attach ANGMIN/ANGMAX, rebalance dispatch)
    ├─ DC-OPF: DCOPF + add_poly_cost_objective + angle_limits=True
    ├─ AC-OPF: ACOPF + add_poly_cost_objective
    │    + thermal_limit='mva'
    │    + free_slack_vm=True
    │    + fix_hv_buses=False
    │    + angle_limits=True
    └─ report objective vs BASELINE.md reference, write results/
```

The PGLib-compatible flags are passed directly to `ACOPF.add_OPF`; no
monkey-patching of the model is required.

---

## Usage

```bash
conda activate potpourri_env

# Default: small and medium cases up to 300 buses
python scripts/pglib_benchmark.py

# Include larger cases (up to 1500 buses)
python scripts/pglib_benchmark.py --max-buses 1500

# Specific cases only
python scripts/pglib_benchmark.py --cases case5_pjm case14_ieee case30_ieee

# AC only (skip DC-OPF)
python scripts/pglib_benchmark.py --no-dc
```

Results are written to `results/pglib_benchmark.csv` and
`results/pglib_benchmark.md`.

---

## Example output

```
Running 6 case(s):
  - case5_pjm       (  5 buses)
  - case14_ieee     ( 14 buses)
  - case24_ieee_rts ( 24 buses)
  - case30_ieee     ( 30 buses)
  - case57_ieee     ( 57 buses)
  - case118_ieee    (118 buses)

case5_pjm         (    5 buses)  DC:   17551.89 (ref   17551.89, +0.00%, 0.4s) ✓  AC:   17551.89 (ref   17551.89, +0.00%, 1.2s) ✓
case14_ieee       (   14 buses)  DC:    7642.59 (ref    7642.59, +0.00%, 0.5s) ✓  AC:    8081.52 (ref    8081.52, +0.00%, 1.8s) ✓
case24_ieee_rts   (   24 buses)  DC:   76943.24 (ref   76943.24, +0.00%, 0.6s) ✓  AC:   63352.20 (ref   63352.20, +0.00%, 2.1s) ✓
case30_ieee       (   30 buses)  DC:     576.89 (ref     576.89, +0.00%, 0.5s) ✓  AC:     576.89 (ref     576.89, +0.00%, 2.3s) ✓
case57_ieee       (   57 buses)  DC:   41006.74 (ref   41006.74, +0.00%, 0.8s) ✓  AC:   41737.79 (ref   41737.79, +0.00%, 3.4s) ✓
case118_ieee      (  118 buses)  DC:   97212.83 (ref   97212.83, +0.00%, 1.2s) ✓  AC:  129338.15 (ref  129338.15, +0.00%, 6.1s) ✓
```

---

## Key API

### `potpourri.benchmarks.load_pglib_case`

```python
from potpourri.benchmarks import load_pglib_case

net = load_pglib_case("case14_ieee")
```

Converts a PGLib `.m` file into a `pandapowerNet` with:

- All generators flagged `controllable=True`
- Polynomial cost coefficients in `net.poly_cost`
- Phase-angle limits (`angmin_degree` / `angmax_degree`) on `net.line` and
  `net.trafo`, read from the MATPOWER `ANGMIN` / `ANGMAX` fields
- Initial dispatch rebalanced so `pp.runpp()` converges

### `potpourri.models.cost_objective.add_poly_cost_objective`

```python
from potpourri.models.cost_objective import add_poly_cost_objective

acopf = ACOPF(net)
acopf.add_OPF(thermal_limit="mva", free_slack_vm=True, angle_limits=True)
add_poly_cost_objective(acopf)
acopf.solve(solver="ipopt")
```

Wires `net.poly_cost` coefficients (`c2·P² + c1·P + c0`) into a Pyomo
objective over all controllable generators (`ext_grid`, `gen`, `sgen`).

### Reference baselines

```python
from potpourri.benchmarks import PGLIB_BASELINE_TYP, PGLIB_BASELINE_API, PGLIB_BASELINE_SAD

ref = PGLIB_BASELINE_TYP["pglib_opf_case14_ieee"]
print(ref["dc"], ref["ac"])   # 7642.59, 8081.52
```

---

## Known limitations

- **Large transmission cases with many `net.impedance` rows** (e.g.
  `case200_activ`, `case588_sdet`) show objective gaps of ±8–32 % relative to
  the PowerModels.jl reference. Smaller cases (`case118_ieee`, `case197_snem`)
  match to ≤ 0.01 %. The gap on larger cases is most plausibly explained by
  PowerModels.jl's input-preprocessing pipeline (topology simplification,
  low-impedance branch merging) that `potpourri` does not perform.
- Only the **Typical Operations (TYP)** baseline is used by default. API and
  SAD cases can be loaded via `load_pglib_case` but are not included in the
  automated benchmark table.
- DC-OPF in PGLib uses linear costs only; `add_poly_cost_objective` accepts
  quadratic terms but warns when `c2 > 0` entries are present.
