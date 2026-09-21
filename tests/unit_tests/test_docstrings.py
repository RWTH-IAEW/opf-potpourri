# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""The examples in the docstrings must match the real API.

A docstring example that no longer runs is worse than none: it reads as
authoritative. Every module carrying a doctest is executed here, so a
renamed argument or a changed return value fails the suite rather than
quietly misleading the next reader.

None of this needs an optimization solver. Constructing a model does run
a pandapower power flow internally, which is exactly the behaviour the
`AC` example documents, so that cost is intended. Solver-dependent
examples belong in an `integration`-marked test, not in a docstring.
"""

from __future__ import annotations

import doctest
import importlib
import pkgutil

import pytest

import potpourri

# Modules whose docstrings contain `>>>` examples. Discovered rather
# than listed, so a new example is picked up without editing this file.
# research/ is skipped: it is exploratory code with its own write-ups,
# and some of it needs data files that are not in the repository.
SKIP_PREFIXES = ("potpourri.research.",)


def _modules_with_doctests() -> list[str]:
    """Names of importable potpourri modules containing a doctest.

    Returns:
        Sorted module names. Import failures are reported as names too,
        so a module that cannot be imported fails its own test rather
        than silently dropping out of the collection.
    """
    finder = doctest.DocTestFinder()
    found = []
    for info in pkgutil.walk_packages(potpourri.__path__, prefix="potpourri."):
        if info.name.startswith(SKIP_PREFIXES):
            continue
        try:
            module = importlib.import_module(info.name)
        except Exception:  # noqa: BLE001 - not this test's business
            continue
        # DocTestFinder attributes each example to the module that
        # defines it, so a `>>>` in an imported symbol's docstring does
        # not make this module look like it has examples.
        #
        # It walks module attributes to do that, which trips Pyomo's
        # deferred imports (`pint` is optional and usually absent) --
        # an unavailable optional dependency is not a doctest failure.
        try:
            tests = finder.find(module)
        except Exception:  # noqa: BLE001 - optional dependency probing
            continue
        if any(test.examples for test in tests):
            found.append(info.name)
    return sorted(set(found))


DOCTEST_MODULES = _modules_with_doctests()


def test_at_least_one_module_has_a_doctest():
    """Guard against the discovery silently finding nothing.

    If the walk breaks, every parametrised test below would vanish and
    the suite would still look green.
    """
    assert DOCTEST_MODULES, "no doctests discovered - is discovery broken?"


@pytest.mark.parametrize("module_name", DOCTEST_MODULES)
def test_docstring_examples_run(module_name):
    """Run a module's doctests and require every example to pass."""
    module = importlib.import_module(module_name)
    result = doctest.testmod(
        module,
        verbose=False,
        optionflags=doctest.NORMALIZE_WHITESPACE | doctest.ELLIPSIS,
    )
    assert result.failed == 0, (
        f"{result.failed} of {result.attempted} doctest example(s) "
        f"failed in {module_name}"
    )
    assert result.attempted > 0, f"no examples ran in {module_name}"
