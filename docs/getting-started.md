# Getting Started

## Installation

---

### For users — install from PyPI

If you want to use `potpourri` in your own scripts or notebooks without
modifying the source code, install the released package directly from PyPI.
This is the recommended path for most users.

**Requirements:** Python 3.9–3.12.

```bash
pip install opf-potpourri
```

Or with [uv](https://docs.astral.sh/uv/) (faster):

```bash
uv pip install opf-potpourri
```

Verify the installation:

```bash
python -c "import potpourri; print(potpourri.__version__)"
```

!!! note "Solvers not included"
    `potpourri` does not bundle any solvers. You need to install at least one
    separately — see [Solver Setup](#solver-setup) below. For AC OPF, IPOPT is
    required. For DC OPF or hosting-capacity problems, GLPK or CBC suffice.
    Alternatively, use [NEOS](#option-c-neos-no-local-solver-required) for
    free remote solver access.

---

### For contributors — local development install

Clone the repository and install in editable mode if you plan to modify the
source code, run the full test suite, or build the documentation locally.

**Requirements:** Python 3.9–3.12 and Conda (or Mamba).

#### Option A — Conda / Mamba (recommended)

All Python dependencies and the IPOPT / GLPK solvers are pinned in
`environment.yaml`.

**1. Install Miniconda or Mamba**

Download the installer for your platform from
[https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
(Miniconda) or
[https://github.com/conda-forge/miniforge](https://github.com/conda-forge/miniforge)
(Mamba / Miniforge — faster solver, drop-in replacement).

```bash
# Linux / macOS — run the downloaded installer
bash Miniconda3-latest-Linux-x86_64.sh   # or the Mamba equivalent
```

**2. Clone the repository**

```bash
git clone https://github.com/RWTH-IAEW/opf-potpourri.git
cd opf-potpourri
```

**3. Create the environment**

```bash
conda env create -f environment.yaml   # Miniconda
# or
mamba env create -f environment.yaml   # Mamba (faster)
```

This creates the `potpourri_env` environment with all pinned dependencies
including IPOPT 3.14.19 and GLPK 5.0.

**4. Activate the environment**

```bash
conda activate potpourri_env
```

**5. Install the package in editable mode**

```bash
pip install -e ".[dev]"
```

The `-e` flag (editable install) means changes to the source code under
`src/potpourri/` take effect immediately without reinstalling. The `[dev]`
extra installs `ruff`, `pytest`, `pytest-cov`, and `pre-commit`.

**6. Update an existing environment**

When `environment.yaml` changes (e.g. after pulling new commits):

```bash
conda env update -f environment.yaml --prune
```

---

#### Option B — Docker (recommended for Windows)

The Dockerfile provides a fully self-contained Linux environment with IPOPT
3.14.16 compiled from source, CBC, and SHOT. This is the recommended path on
**Windows** because Windows conda channels do not distribute IPOPT.

**Prerequisites**

Install [Docker Desktop](https://www.docker.com/products/docker-desktop/).
On Windows, Docker Desktop uses the WSL 2 backend by default — no additional
configuration is required.

**1. Build the image**

From the repository root (takes 10–20 minutes on first build; IPOPT is
compiled from source):

```bash
docker build -t potpourri:latest .
```

**2. Run an interactive session**

```bash
docker run -it --rm potpourri:latest bash
```

Inside the container the `potpourri_env` conda environment is already active
and `potpourri` is installed.

**3. Mount your working directory (optional)**

To edit files on your host machine and run them inside the container:

```bash
# Windows (PowerShell)
docker run -it --rm -v ${PWD}:/app potpourri:latest bash

# Linux / macOS
docker run -it --rm -v $(pwd):/app potpourri:latest bash
```

**4. Run a script directly**

```bash
docker run --rm -v $(pwd):/app potpourri:latest \
    conda run -n potpourri_env python scripts/loadcases.py
```

**Solvers available inside the container**

| Solver | Type | Source |
|---|---|---|
| IPOPT 3.14.16 | NLP | compiled from source |
| GLPK | LP / MIP | conda-forge (via environment.yaml) |
| CBC | MIP | `coinor-cbc` apt package |
| SHOT | MINLP | compiled from source |

---

#### Option C — NEOS (no local solver required)

[NEOS](https://neos-server.org/) is a free public optimisation server that
accepts Pyomo models over the internet. This option is available to both
regular users and developers when no local solver can be installed.

**1. Register an e-mail address** at [https://neos-server.org](https://neos-server.org).

**2. Set the environment variable** before running any script:

```bash
export NEOS_EMAIL="your@email.address"   # Linux / macOS
set NEOS_EMAIL=your@email.address        # Windows CMD
$env:NEOS_EMAIL = "your@email.address"  # Windows PowerShell
```

Or set it at the top of your Python script:

```python
import os
os.environ["NEOS_EMAIL"] = "your@email.address"
```

**3. Pass `solver='neos'`** to any `solve()` call:

```python
opf.solve(solver='neos', neos_opt='ipopt')
```

!!! note
    NEOS solves run on shared infrastructure. Jobs are queued and results
    can take longer than a local solve. For repeated or large problems a
    local solver installation is strongly preferred.

---

## Verify installation

```python
import pandapower as pp
import simbench as sb
from potpourri.models.ACOPF_base import ACOPF

net = sb.get_simbench_net("1-LV-rural1--0-sw")
opf = ACOPF(net)
opf.add_OPF()
opf.add_voltage_deviation_objective()
opf.solve(solver="ipopt", print_solver_output=False)

print(opf.net.res_bus[["vm_pu", "va_degree"]].head())
```

If the solver returns an optimal solution and `opf.net.res_bus` is populated,
the installation is working correctly.

---

## Solver setup

`potpourri` does not bundle any solvers. Install at least one of the
following before calling `solve()`:

| Solver | Problem type | Install |
|---|---|---|
| **IPOPT** | NLP — AC OPF | `conda install -c conda-forge ipopt` |
| **GLPK** | LP / MIP — DC OPF | `conda install -c conda-forge glpk` |
| **CBC** | LP / MIP | `conda install -c conda-forge coincbc` |
| **Gurobi** | LP / MIP / NLP | `pip install gurobipy` (licence required) |

IPOPT and GLPK are included automatically in the developer Conda environment
(`environment.yaml`). PyPI users must install solvers separately using the
commands above.

---

## Development workflow

Once you have a local development install (Option A or B above):

```bash
ruff check .                        # lint
ruff format .                       # format (auto-fix)
pytest -m "not integration"         # unit tests — no solver required
pytest                              # all tests — integration tests need IPOPT
```

To build and serve the documentation locally:

```bash
pip install -e ".[docs]"
mkdocs serve                        # available at http://127.0.0.1:8000/
```

Analysis and example scripts are in `scripts/`. See the
[Examples](scripts/examples.md) section of this documentation for details.

---

## Acknowledgements

The development of `potpourri` was inspired in part by modelling concepts from
the `oats` package [@bukhsh2020oats]. Early development built on ideas and code
structures from `oats`; over time, the project evolved into a standalone
package focused on providing a pandapower-compatible interface to Pyomo-based
OPF formulations.

The modular and inheritance-based structure of `potpourri` was also influenced
by PowerModels.jl [@coffrin2018powermodels].

The authors thank the contributors to pandapower [@thurner2018pandapower],
Pyomo [@bynum2021pyomo; @hart2011pyomo], and SimBench [@meinecke2020simbench],
on which `potpourri` depends.

\bibliography
