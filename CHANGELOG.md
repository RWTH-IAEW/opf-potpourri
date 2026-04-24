# Changelog

All notable changes to `potpourri` are documented here.
Format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [Unreleased]

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
