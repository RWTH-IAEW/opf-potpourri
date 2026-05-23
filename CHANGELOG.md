# Changelog

All notable changes to `potpourri` are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

## [0.3.1] — 2026-05-23

### Added

- **PGLib-OPF git submodule** at `benchmarks/pglib-opf/` — the
  [IEEE PES Power Grid Library](https://github.com/power-grid-lib/pglib-opf)
  benchmark dataset is now included as a submodule. Clone with
  `git clone --recurse-submodules` to enable benchmark runs. The submodule
  is not part of the PyPI package; `pip install opf-potpourri` is unaffected.
- Documentation updates: README "Benchmarking against PGLib-OPF" section,
  `docs/scripts/examples.md` entries for `pglib_benchmark.py` and
  `timeseries_acopf.py`, developer clone command updated to
  `--recurse-submodules`.

## [0.3.0] — 2026-05-23

### Added

- `scripts/timeseries_acopf.py` — pandapower `run_timeseries` demo that
  drives loads and sgens from SimBench profiles via
  `ConstControl` / `DFData`, runs the potpourri AC OPF at every step via
  a custom `run=run_acopf` callback, and logs results with
  `OutputWriter`.
- **PGLib-OPF benchmark integration**
  - New `potpourri.benchmarks` package with `load_pglib_case`, which parses
    a PGLib `.m` file (via `pandapower.converter.matpower.from_mpc`) into a
    pandapower network ready for OPF, attaches MATPOWER `ANGMIN` / `ANGMAX`
    onto `net.line` / `net.trafo`, and rebalances the initial dispatch so
    `pp.runpp` converges on PGLib snapshots.
  - `PGLIB_BASELINE_TYP` / `PGLIB_BASELINE_API` / `PGLIB_BASELINE_SAD` dicts
    parsed from upstream `BASELINE.md`, plus `parse_baseline_md` helper.
  - `scripts/pglib_benchmark.py` runs DC and AC OPF on a configurable
    PGLib subset and writes `results/pglib_benchmark.{csv,md}` with a
    BASELINE.md-style markdown table.
- **Polynomial generator cost objective**
  - `potpourri.models.cost_objective.add_poly_cost_objective` wires
    `net.poly_cost` polynomial coefficients (`c2·P² + c1·P + c0`) into a
    Pyomo objective over `ext_grid + gen + sgen`. Raises a clear error if
    `net.pwl_cost` is non-empty (PWL costs not yet supported).
- **Configurable AC-OPF flags** on `ACOPF.add_OPF`:
  - `thermal_limit="current"|"mva"` — choose between the current-based
    constraint `|S|² ≤ SLmax²·v²` (default, physically meaningful for
    distribution conductors) and the MATPOWER / PowerModels-compatible
    constant-MVA limit `|S|² ≤ SLmax²`.
  - `free_slack_vm=True` (new default) — the slack-bus voltage magnitude
    floats within `[Vmin, Vmax]` for AC OPF, while the reference angle
    remains fixed. Set `False` to reproduce the legacy AC-PF behaviour
    where the slack `vm` is pinned to its base-case value.
  - `fix_hv_buses=False` (new default) / `hv_bus_kv=110.0` — replace the
    hardcoded 110 kV bus pinning with an explicit opt-in flag.
  - `angle_limits=False` — optional enforcement of branch
    phase-angle-difference constraints `angmin ≤ δ_from − δ_to ≤ angmax`
    from `net.line.angmin_degree` / `net.line.angmax_degree`. Available
    on both AC and DC OPF.
- **`net.impedance` branch support** in single- and multi-period AC / DC
  OPF. pandapower's third branch table — used for branches whose
  from/to buses have different nominal voltage but no off-nominal tap —
  is now treated as additional lines (per-unit physics is identical) and
  enters `model.L` with synthetic indices `≥ len(net.line)`. Without
  this, several PGLib cases (e.g. `case118_ieee`, `case89_pegase`,
  `case200_activ`) had electrically isolated buses or under-modelled
  parallel paths.
- **Regression test suite for the formulation audit fixes** in
  `tests/unit_tests/test_audit_fixes.py` — covers D1–D14 plus the
  impedance-branch fix and the multi-period DC OPF construction path.
- **CHANGELOG.md** (this file).

### Changed

- **Slack-bus voltage in AC OPF** is now free within `[Vmin, Vmax]` by
  default; the reference angle stays fixed. The previous pinning is
  available via `free_slack_vm=False`.
- **110 kV bus voltage pinning** is no longer applied by default. Re-enable
  with `fix_hv_buses=True` (and optionally `hv_bus_kv=...`).
- **Storage docstring** in `Basemodel.add_storage` now states explicitly
  the round-trip-efficiency convention (`η·Pchg − Pdis/η`), the convex
  relaxation `Pchg + Pdis ≤ Pmax` used in place of complementarity, and
  the single-period nature of the SOC update.
- **AC branch admittance comment** in `AC.py` now documents the explicit
  assumption that branch shunt conductance is zero (the standard MATPOWER
  convention).
- **`Basemodel.preprocess_grid`** now references a comprehensive list of
  element tables (lines, trafos, trafo3w, impedance, dcline, load, sgen,
  gen, ext_grid, shunt, storage, ward, xward, motor, asymmetric loads
  and sgens) for both the bus-bus switch merge and the orphan-bus prune.
- **PGLib benchmark script** no longer monkey-patches the model;
  PGLib-compatible behaviour is now driven entirely by the new
  `ACOPF.add_OPF` flags (`thermal_limit="mva"`, `free_slack_vm=True`,
  `fix_hv_buses=False`, `angle_limits=True`).
- **Generator-cost helper** writes the polynomial cost in MW units
  (per-unit pyomo variables multiplied by `baseMVA`), matching the
  MATPOWER `mpc.gencost` polynomial convention.

### Fixed

- **`OPF.generation_real_power_limits`** and
  **`ACOPF_base.generation_reactive_power_limits`** now filter source
  tables to in-service rows before broadcasting into the
  `net._gen_order[element]` slice. Out-of-service generators no longer
  trigger `ValueError: could not broadcast …` (e.g. PGLib
  `case200_activ` with 11 OOS gens).
- **Degenerate generator P-range** (`Pmin == Pmax`) — synchronous
  condensers in PGLib (`case30_ieee`, `case60_c`, `case240_pserc`, …)
  no longer drive IPOPT into `TOO_FEW_DEGREES_OF_FREEDOM`. The variable
  is pinned via tight ε-padded bounds while the redundant range
  constraint is skipped; the variable remains in the NL file so IPOPT
  doesn't lose a degree of freedom to Pyomo's fixed-variable elimination.
- **`Basemodel.preprocess_grid`** no longer silently drops buses that
  are referenced only by `shunt` or `storage` (the prior
  `referenced_buses` union omitted those tables, which crashed
  `pp.create_continuous_bus_index` on, e.g., PGLib `case200_activ`).
- **Bus-bus switch merging** now updates `shunt`, `storage`, `trafo`,
  `trafo3w`, `ward`, `xward`, and the asymmetric load/sgen tables when
  remapping a removed bus, not only `line` / `load` / `sgen` / `gen` /
  `ext_grid`.
- **`pyo_to_net._generation_results_to_net`** now uses key-aware
  iteration (`model.G` ∩ `_gen_order[et]` slice) rather than positional
  `range(f, t)` indexing, making it robust to non-contiguous generator
  indices.
- **`pyo_to_net._line_results_to_net`** now writes back line flows via
  `net.line.index`-keyed lookup with `dict.get(idx, 0.0)`. Out-of-service
  lines and impedance-only synthetic indices in `model.L` no longer
  produce length mismatches.
- **DC line-susceptance lookup** (`DC.BL_data`, `DC_multi_period.BL_data`)
  is now a `pd.Series` keyed by the actual `model.L` indices, so cases
  with out-of-service lines (e.g. PGLib `case2000_goc` with 6 OOS lines)
  no longer fail with `KeyError: Index 'k' is not valid for indexed
  component 'BL'`.
- **`DC_multi_period.create_model`** has been rewritten:
  - The dead `if self.T is None: ... else: ...` split has been removed
    (`self.T` is always an int set by `Basemodel_multi_period.__init__`).
  - Constraint indexing now uses the Pyomo Set `self.model.T` rather
    than the bare int `self.T` (the cause of
    `TypeError: Cannot create a Set from data that does not support
    __contains__. ... received 'int'`).
  - Single-period parameters (`BL[l]`, `BLT[l]`, `shift[l]`, `GB[s]`)
    are no longer accessed with a spurious time index.
  - Time-indexed variables (`delta[bus, t]`, `deltaL[l, t]`,
    `deltaLT[l, t]`, `pLfrom/pLto/pThv/pTlv[l, t]`) are correctly
    addressed throughout.
- **`DC_multi_period.__init__`** now accepts a `pf=1` keyword so
  subclasses (`DCOPF_multi_period`) can forward it through the MRO.
- **`DCOPF_multi_period`** rewritten to follow the single-period DCOPF
  pattern: `__init__` no longer triggers a duplicate `create_model()`,
  and an explicit `add_OPF` method attaches the line/transformer
  thermal limits using the correct variable names (`pLfrom`, `pThv`
  instead of the previously non-existent `pL`, `pLT`).
- **`Sgens_multi_period.get_opf_parameters`** no longer references
  `self.QsGmax_tuple` / `self.QsGmin_tuple`, which are populated only
  by the AC-only step `static_generation_reactive_power_limits`.
  Reactive-power parameter setup has moved into a new
  `get_acopf_parameters` hook called only via `get_all_acopf`. Without
  this, any multi-period DC OPF call raised
  `AttributeError: 'Sgens_multi_period' object has no attribute
  'QsGmax_tuple'`.

### Notes / unresolved

- For several large PGLib transmission cases with many `net.impedance`
  rows (`case200_activ`, `case588_sdet`, `case793_goc`, `case2312_goc`),
  the AC and DC objectives match the PGLib reference within a consistent
  ±8 % to ±32 % band rather than to machine precision. case197_snem
  (55 impedance rows) and case118_ieee (2 impedance rows) match to
  ≤ 0.01 %, indicating the impedance-branch fix is mathematically
  correct; the remaining gap on the larger cases is the most plausibly
  explained by PowerModels.jl's input-preprocessing pipeline
  (e.g. topology simplification, low-impedance branch merging) that
  potpourri does not perform. Side-by-side investigation against
  PowerModels.jl is the recommended next step.
- The pandapower import warning `some BR_B of transformers in
  ppc['branch'] are positive, but transformers are inductive` was
  investigated and confirmed to be **benign**: it reflects pandapower's
  internal `i0_percent` interpretation. potpourri reads `BR_B` directly
  from `_ppc['branch']` and treats it as the generic line-charging
  susceptance, matching MATPOWER and PowerModels conventions. The
  underlying data is roundtripped to float epsilon.

## [0.2.1]

### Changed
- README: split installation into a user section (`pip install opf-potpourri`)
  and a developer section (conda environment + editable install).
- README: added PyPI badge and direct links to the Read the Docs sections.
- README: moved solver installation guidance into a dedicated Solvers section
  and noted that PyPI users must install solvers separately.

## [0.2.0] — initial public release

### Added
- Single-period AC OPF (`ACOPF`, `HC_ACOPF`) and DC OPF (`DCOPF`) models.
- Multi-period AC OPF (`ACOPF_multi_period`) with a configurable time horizon
  and SimBench profile integration.
- Flexible-resource mix-in modules: `Battery_multi_period`,
  `HeatPump_multi_period`, `PV_multi_period`, `Windpower_multi_period`,
  `Demand_multi_period`, `Sgens_multi_period`, `Generator_multi_period`.
- Hosting-capacity OPF (`HC_ACOPF`, `HC_ACOPF_multi_period`) with binary
  placement variables and VDE-AR-N 4105 reactive-power constraints.
- Post-processing utilities (`pyo_to_net`, `pyo_to_net_multi_period`) that
  write Pyomo solutions back to `net.res_*` DataFrames.
- Warm-start helpers (`init_pyo_from_pp_res`,
  `init_pyo_from_pp_res_multi_period`) based on pandapower power-flow results.
- Runnable example scripts covering single-period OPF, multi-period planning,
  hosting-capacity analysis, feasible-operation-region computation, and
  solver benchmarking.
- MkDocs documentation with API reference, mathematical formulations, and
  class architecture diagrams.
- GitHub Actions CI: lint (Ruff), unit tests (Python 3.10–3.12), and
  installation smoke test.
- MIT licence.
