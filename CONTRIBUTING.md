# Contributing to potpourri

Thank you for your interest in contributing to `potpourri`. All contributions — bug reports, feature suggestions, documentation improvements, and code — are welcome.

## How to report issues

Please open a **GitHub Issue** at <https://github.com/RWTH-IAEW/opf-potpourri/issues>.

When reporting a bug, include:
- a minimal reproducible example,
- the Python version and operating system,
- the solver being used (IPOPT, GLPK, CBC, Gurobi, NEOS), and
- the full traceback.

## How to propose changes

1. **Fork** the repository on GitHub.
2. Create a branch for your change:
   ```bash
   git checkout -b fix/my-fix
   ```
3. Make your changes and commit them (see commit style below).
4. Open a **Pull Request** against the `main` branch.

All pull requests are reviewed before merging. Ensure that the CI checks pass before requesting a review.

## Setting up the development environment

```bash
git clone https://github.com/RWTH-IAEW/opf-potpourri.git
cd opf-potpourri

conda env create -f environment.yaml   # creates potpourri_env
conda activate potpourri_env

pip install -e ".[dev]"                # installs ruff, pytest, pytest-cov, pre-commit
pre-commit install                     # optional: install git hooks
```

Solvers (IPOPT, GLPK, CBC) are included in `environment.yaml` via `conda-forge`.

## Running the tests

```bash
pytest                              # all tests (integration tests require IPOPT)
pytest -m "not integration"        # unit tests only — no solver required
pytest tests/installation_with_pip  # verify public API
```

Integration tests are marked with `@pytest.mark.integration` and require a locally installed IPOPT solver. They are skipped automatically in CI.

## Code style

This project uses [Ruff](https://docs.astral.sh/ruff/) for linting and formatting.

```bash
ruff check .       # lint
ruff format .      # format (auto-fix)
```

CI will reject pull requests that fail either check. Run both locally before pushing.

## Commit messages

Follow [Conventional Commits](https://www.conventionalcommits.org/):

```
feat: add battery discharge efficiency parameter
fix: resolve AttributeError on tap_phase_shifter networks
docs: improve quick-start example in README
```

## Questions

For questions about usage or the mathematics behind the models, open a GitHub Discussion or an Issue with the label `question`.
