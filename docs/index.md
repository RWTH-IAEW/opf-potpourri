# POTPOURRI

**Multi-Period Optimal Power Flow for Distribution Grids with Storage Application**

POTPOURRI is a Python library for AC and DC Optimal Power Flow (OPF) in distribution grids. It provides a class-based interface over [Pyomo](https://pyomo.readthedocs.io/) for building and solving power system optimisation problems defined on [pandapower](https://pandapower.readthedocs.io/) networks.

## Key features

- **Single-period AC/DC OPF** — full nonlinear AC and linearised DC formulations
- **Multi-period OPF** — time-indexed optimisation over load/generation profiles from [SimBench](https://simbench.de/en/)
- **Modular flexible devices** — batteries, heat pumps, PV, wind
- **Hosting capacity analysis** — binary wind generator placement with grid-code Q-curve constraints
- **Multiple solvers** — IPOPT (NLP), Gurobi/CBC via MindtPy (MINLP), NEOS remote

## Quick example

```python
import simbench as sb
from potpourri.models.ACOPF_base import ACOPF

net = sb.get_simbench_net("1-LV-rural1--0-sw")
opf = ACOPF(net)
opf.add_OPF()
opf.add_voltage_deviation_objective()
opf.solve(solver="ipopt")

# results written to net.res_bus, net.res_line, net.res_sgen, ...
print(net.res_bus[["vm_pu", "va_degree"]])
```

## Navigation

- [**Getting Started**](getting-started.md) — installation and solver setup
- [**Single-Period OPF**](user-guide/single-period.md) — AC/DC OPF tutorial
- [**Multi-Period OPF**](user-guide/multi-period.md) — time-series OPF tutorial
- [**Flexible Devices**](user-guide/devices.md) — batteries, EVs, heat pumps
- [**API Reference**](api/models.md) — full class and method documentation
