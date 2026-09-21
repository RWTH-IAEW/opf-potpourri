# Documenting potpourri

This page is the documentation counterpart to
[Licensing](licensing.md): what a docstring in this repository has to
say, where prose belongs, and how to check it before pushing. For how
the modelling code itself is written, see
[Pyomo Conventions](contributing-pyomo.md).

The audience we write for is a power-system researcher who knows OPF but
has never seen this codebase, and a Python developer who knows neither
our conventions nor the grid codes. Both should be able to open a file
and find out what it does without reading the implementation.

---

## Style

**Google-style docstrings, rendered as Markdown.** Ruff enforces the
convention (`[tool.ruff.lint.pydocstyle] convention = "google"` with the
`D` rules selected), and `mkdocstrings` renders docstring bodies as
Markdown, not reStructuredText.

That last point decides the markup. Use:

| for | write | not |
|-----|-------|-----|
| inline maths | `$v_b$` | `` :math:`v_b` `` |
| display maths | `$$ ... $$` | `.. math::` |
| a cross-reference | `` [`ACOPF`][potpourri.models.ACOPF_base.ACOPF] `` | `` :class:`ACOPF` `` |
| code | `` `model.v[b]` `` | ``` ``model.v[b]`` ``` |
| a table | a Markdown pipe table | an RST grid table |

Maths renders through `pymdownx.arithmatex` and MathJax, which the site
already loads. A docstring containing a backslash must be a raw string
(`r"""`), or ruff's `D301` will say so.

Keep lines within 79 characters, the same limit as the code. SPDX header
lines are exempt; ruff knows.

### Sections

Use the Google sections the renderer understands: `Args:`, `Returns:`,
`Yields:`, `Raises:`, `Attributes:`, `Example:`, `Note:`, `See Also:`.
Do not invent section names and do not leave a section empty.

Open with a **one-line summary**, then a blank line, then the detail
(`D205`). If the summary will not fit on one line it is doing too much:
shorten it and move the qualification into the body.

The constructor contract goes in the **class** docstring, not in
`__init__`. `merge_init_into_class` folds them together when rendering,
so documenting both duplicates the text; `D107` and interrogate's
`ignore-init-method` are both configured for that choice.

---

## What a docstring must establish

Beyond restating the name:

- **Units and sign conventions.** p.u. or MW, radians or degrees,
  injection or consumption. Storage is the usual trap: it follows
  pandapower's *load* convention, so positive means charging.
- **Indexing.** Which set an index comes from, and whether it is a
  pandapower label or an internal **ppc** bus number. These differ
  whenever pandapower inserts auxiliary buses for node-node switches.
- **Side effects.** Anything written to `self.model`, `self.net` or the
  result tables. A method that builds a model and returns `None` should
  say where the components can be reached.
- **Prerequisites.** Which call has to happen first (`add_OPF` before
  an objective, `create_model` before inspecting `model`).
- **Repeat calls.** Whether calling twice is safe, replaces components,
  or is unsupported.
- **What it is not.** The single most useful sentence in many of these
  docstrings says which nearby thing does the job instead.

### Model-building methods and Pyomo rules

A constraint rule deserves a docstring whenever its mathematics is not
obvious from three lines of Pyomo -- which, for this codebase, is
nearly always. State the relation it enforces, in maths, and mention:

- **What it returns.** A Pyomo *expression*, not a `bool`. A ranged
  constraint returns the 3-tuple `(lower, expr, upper)`. Say which.
- **`Constraint.Skip`.** When a rule skips, say under what condition.
- **Side effects inside rules.** Several bound rules call `unfix()`
  before returning; that is the important half of what they do.
- **Approximations.** If a constraint is a relaxation, a linearisation
  or an inner approximation, say so and say in which direction it
  errs. Do not call a formulation convex, exact or globally optimal
  without justification.

Keep the derivation itself in
[Mathematical Modelling](mathematical-modelling.md) and link to it;
keep the short statement of the equation next to the code.

Note that mkdocstrings does **not** render a function nested inside a
method, so a constraint rule's docstring reaches readers of the source
and the `show_source` viewer, not the API page. Write it for someone
reading the code -- that is who needs it -- and put anything the API
page must show in the enclosing method's docstring.

---

## Comments

Comments explain **why**, docstrings explain **what**. Write a comment
for a sign convention, an index mapping, a numerical safeguard, an
upstream compatibility workaround, or a construction order that matters.

Do not write `# loop over buses`. Comment density is deliberately *not*
a quality target here, and no check measures it: a high ratio is as
easily achieved by noise as by insight.

---

## Examples

Examples must use the real API and be executed, not eyeballed. Build a
small network in the example itself rather than loading a file or a
benchmark case.

Be explicit about what an example needs. Constructing a model runs a
pandapower power flow but needs no optimization solver; solving needs
IPOPT or equivalent. A doctest belongs in the first category -- put
solver-dependent examples in a test marked `integration` instead of a
docstring, so the doctest suite stays runnable without a solver.

Assert properties with tolerances, never exact floating-point output.

---

## Running the checks

The same four commands run locally and in CI.

```bash
# 1. Coverage: does a docstring exist at all?
interrogate src/potpourri

# 2. Structure: is it well formed, Google convention?
ruff check src/potpourri

# 3. Rendering: do the generated API pages still build?
mkdocs build --strict

# 4. Examples: do the doctests still match the API?
pytest tests/unit_tests/test_docstrings.py
```

`pre-commit run --all-files` runs 1 and 2. The GitLab `documentation`
job runs 1 and 2, and `docs_build` runs 3.

### What the numbers mean

`interrogate` measures **presence, not correctness**. A file can score
100% and still describe the wrong equation. The threshold in
`[tool.interrogate]` is a ratchet against regression, nothing more --
the scientific content is reviewed by people.

Scope is `src/potpourri`. `scripts/`, `tests/` and `tools/` are
measured separately and are not gated, because an example script and a
regression test have different documentation needs from a library API.
Nothing inside `src/potpourri` is excluded: in particular nested
functions stay counted, because in this codebase they are the Pyomo
constraint rules, and they are exactly the code that most needs
explaining.

---

## When you are unsure

If the code and an existing docstring disagree, **investigate before
editing either**. Do not quietly change an equation to match its prose,
and do not restate a comment you cannot verify. If the behaviour looks
wrong, document what it actually does, add a focused regression test,
and raise the discrepancy separately.
