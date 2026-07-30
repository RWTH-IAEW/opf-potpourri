# Changelog

All notable changes to `potpourri` are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

### Added

- **Selectable grid codes** (`src/potpourri/technologies/q_control.py`): the
  connection rules are now `GridCode` parameter sets in a `GRID_CODES`
  registry instead of module-level constants, selected model-wide via
  `ACOPF.add_OPF(grid_code=…)` or `ACOPF_multi_period.add_OPF(grid_code=…)`.
  Accepts `None` (VDE-AR-N 4105, the default, so existing models are
  unaffected), a short name (`"4105"`, `"4110"`), or a `GridCode`; unknown
  names raise `ValueError`. A `GridCode` carries the Q(U) breakpoints, the
  Q/Pn capability table and its variants, the Q(P) breakpoints, and the P(U)
  and cos(φ)(P) thresholds.
  - **`VDE_AR_N_4110` (medium voltage) is provisional**: it is registered but
    reuses the VDE-AR-N 4105 values as a placeholder, because its normative
    medium-voltage figures have not been entered yet. Results obtained with
    `grid_code="4110"` are **not** 4110-compliant. Selecting it emits a
    `ProvisionalGridCodeWarning` instead of failing.
  - `ACOPF.static_generation_wind_var_q()` previously reimplemented the
    Q-curve maths inline with hard-coded constants, so the single-period path
    would have ignored the selected code; it now routes through
    `compute_q_curves()`. Numerically identical for VDE-AR-N 4105.
  - The module-level constants (`VQU_Q_MAX`, `QP_P_HIGH`, `VPU_V_CURTAIL`, …)
    are retained as aliases of the default grid code, so existing imports
    keep working.
  - Not yet covered: the hosting-capacity wind path in `windpower.py` keeps a
    private copy of the 4105 table and is not registry-driven, so `grid_code`
    does not affect `HC_ACOPF` runs.
- **`ACOPF.add_OPF(sgen_types=…)`** — the sgen `type` values that
  `pv_q_control` treats as PV, defaulting to `DEFAULT_PV_SGEN_TYPES`
  (`("PV", "PV_MV")`). Previously the filter matched `type == "PV"` exactly,
  which reached **zero** sgens on every SimBench medium-voltage grid (they are
  labelled `PV_MV`, `Wind_MV`, `lv_RES`, …), so `pv_q_control` silently built
  no constraints there and raised no warning. The new default reaches 2–5
  units per MV grid; pass a wider list such as
  `("PV", "PV_MV", "Wind_MV", "lv_RES")` to include the aggregated
  LV-renewable and MV wind units, raising the reach to 87–133. LV grids are
  unaffected. Note that the wind Q-control path is selected separately and
  still matches `type == "Wind"` exactly, and that the multi-period model
  applies no type filter at all.
- **`tests/unit_tests/test_q_control.py`** — 71 solver-free unit tests for the
  Q-control code, which previously had no automated coverage at all. Covers
  the grid-code registry (capability curves against the VDE-AR-N 4105 table,
  variant ordering, alias resolution, unknown-code `ValueError`, the
  provisional-value warning contract), single-period activation for every
  `add_OPF` argument, and multi-period column-driven activation. Includes
  regression tests for the two failure modes that are otherwise silent: the
  cos(φ) cone quietly not being added without `inverter_s2=True`, and the
  `var_q` reactive-bound path, which a refactor once broke while the whole
  suite still passed.
- **`scripts/grid_code_q_strategies.py`** — worked example of both selection
  mechanisms: one snapshot solved under each registered grid code (surfacing
  the provisional-values warning rather than silencing it), then Q(P)/Q(U),
  fixed cos(φ), cos(φ)(P) and P(U) curtailment assigned to different PV units
  within a single multi-period model.
- **Documentation** for per-sgen strategy assignment in
  `docs/user-guide/reactive-power-control.md`, including which strategies may
  share an sgen: fixed cos(φ) and the cos(φ)(P) profile must not (both are
  equalities on the same reactive power), while P(U) curtailment may be
  combined with a Q rule since it constrains active power. Per-row assignment
  is a multi-period feature; in the single-period model `pu_curtail` and
  `cos_phi_p_profile` are model-wide switches.
- **VDE-AR-N 4105 / BDEW reactive-power control for PV and wind sgens**
  (`src/potpourri/technologies/q_control.py`,
  `src/potpourri/technologies/pv.py`,
  `src/potpourri/technologies/sgens.py`):
  - `compute_q_curves()` returns the slope/intercept coefficients for all
    three Q-control variants (`var_q` 0–2).
  - `ACOPF.add_OPF(pv_q_control=…)` accepts `None` (off), `"qp"` (Q(P)
    only), `"qu"` (Q(U) only), or `"both"`. Adds `PV_QP_pos/neg` and
    `PV_QU_min/max` on set `PVc`.
  - `ACOPF_multi_period` detects `var_q` in `net.sgen` automatically and
    adds the time-indexed `sG_QP_pos/neg` / `sG_QU_min/max` constraints.
- **Inverter S² apparent-power circle**: `ACOPF.add_OPF(inverter_s2=True)`
  (opt-in; defaults to `False`) reads `net.sgen.sn_mva` and
  `net.sgen.converter_sizing_pu` and adds
  `psG[g]² + qsG[g]² ≤ S_inv[g]²`. The multi-period model enables the
  time-indexed form automatically when `sn_mva` is present.
- **cos(φ) cone** completing the PV operating region (P ≥ 0 ∩ S² circle ∩
  cone): `|qsG[g]| ≤ psG[g] · tan(arccos(cos_phi_min))`, from
  `net.sgen["cos_phi_min"]` or the `cos_phi_min` keyword. In the
  single-period model the cone sits inside the S² block, so it also
  requires `inverter_s2=True` and a usable `sn_mva`.
- **P(U) active-power curtailment** (VDE-AR-N 4105 §8.5): above a voltage
  threshold the allowed active output falls linearly to zero.
  `add_OPF(pu_curtail=True)`, thresholds from `net.sgen.v_curtail_pu` /
  `v_max_curtail_pu` (defaults 1.06 / 1.10 p.u.). Bilinear — needs an NLP
  solver.
- **Fixed cos(φ) mode**: `qsG[g] = psG[g] · tan(arccos(cos_phi))` as an
  equality, via `add_OPF(fixed_cos_phi=…)` or `net.sgen["fixed_cos_phi"]`.
- **cos(φ)(P) profile**: quadratic equality
  `qsG · (Pn − Pt) = tan_phi · psG · (psG − Pt)`, via
  `add_OPF(cos_phi_p_profile=True)`.
- **Multi-period OLTC tap optimisation**: `add_tap_changer_linear()` and
  `add_tap_changer_discrete()` now use the time-indexed `Tap[tr, t]`
  variable, and `add_tap_changer_linear()` accepts an optional
  `max_tap_change_per_step` rate limit.
- **`scripts/q_control_opf.py`** — compares the Q-control modes in the
  single-period case, then the PV inverter controller modes, then runs a
  24-step multi-period AC OPF with automatic Q-control detection.
- **`docs/user-guide/reactive-power-control.md`** — user-guide page for all
  of the above, including an activation-reference table giving the exact
  precondition for each constraint group in both the single- and
  multi-period model, since a missing precondition is silent (the
  constraint is simply not added, with no warning).
- **Docs CI**: `.github/workflows/ci.yml` gained a `docs` job running
  `mkdocs build --strict`, so broken internal links and invalid navigation
  now fail the build.

### Fixed

- `Basemodel` no longer deletes buses aggressively during network
  preparation, which raised `KeyError` on HV/MV grids.
- `mkdocs.yml` now enables the `admonition` markdown extension. It was
  missing, so every `!!! note` block in the documentation rendered as
  literal text.
- `mkdocs.yml` disables `mkdocs-bibtex` inline citations. Inline (bare)
  `@key` parsing runs on the raw markdown before code fences are handled, so
  the Pyomo decorator `@model.Constraint(...)` in the device-development
  example was read as a citation key and failed a strict docs build.
  Bracketed `[@key]` citations are unaffected.

### Changed

- Example scripts in `scripts/` use a uniform style: tunable parameters sit
  in a configuration block of module-level constants instead of
  command-line arguments. Some previously infeasible example setups were
  corrected.

### Known issues

- **Single-period models are built on the wrong bus set when pandapower's ppc
  conversion adds auxiliary buses.** `Basemodel.__init__` truncates the ppc bus
  table to `len(net.bus)` rows and uses those ppc bus numbers as `model.B`.
  Where the ppc is longer — `1-MV-rural--0-sw` has 103 ppc buses for 97
  pandapower buses, from switch handling — the slice keeps auxiliary buses
  while dropping real pandapower buses along with the in-service branches
  attached to them. It also lets `net._pd2ppc_lookups["bus"]` resolve to a bus
  the model does not contain, raising `KeyError` in `pyo_to_net` when writing
  results back (observed on `1-MV-rural--0-sw`; latent on `1-MV-comm--0-sw`).
  Grids whose ppc has no auxiliary buses — all the LV cases — are unaffected,
  which is why this went unnoticed. Fixing it is not a local change:
  `model.B` is in ppc space while `get_v_limits()` returns arrays in
  pandapower positional order, and twelve sites across the AC and LPAC models
  (single- and multi-period) index the latter by the former. A correct fix
  must also decide which voltage limits apply to auxiliary buses, which
  pandapower does not expose directly. **Medium-voltage results should be
  treated as unreliable until this is resolved.**

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
