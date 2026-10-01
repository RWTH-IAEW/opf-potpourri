# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""The `neos` fixture turns a NEOS service failure into a skip, not a fail.

On 2026-10-01 NEOS's CPLEX queue accepted every job and returned an empty
solution file, which Pyomo surfaces as an `ActionManagerError`; the one test
that used it turned `test_unit` red on unrelated merge requests. The fixture
in `tests/conftest.py` exists so that an outage of the public server reads
as "skipped: NEOS did not return a solution" and the suite stays green. This
module checks that behaviour without needing the outage: a stand-in model
raises the error the kestrel plugin raises, and the fixture must skip.
"""

import pytest
from pyomo.opt.parallel.manager import ActionManagerError


class _FailingSubmission:
    """Stand-in for a model object whose NEOS submission returns nothing."""

    def solve(self, **kwargs):
        raise ActionManagerError(
            "Problem executing an event.  No results are available."
        )


class _Solved:
    """Stand-in for a model object whose NEOS submission succeeds."""

    def solve(self, **kwargs):
        return {"solver": kwargs.get("solver"), "opt": kwargs.get("neos_opt")}


@pytest.mark.integration
def test_neos_fixture_skips_when_the_server_returns_no_solution(neos):
    with pytest.raises(pytest.skip.Exception, match="did not return"):
        neos(_FailingSubmission(), neos_opt="cplex")


@pytest.mark.integration
def test_neos_fixture_passes_results_through(neos):
    results = neos(_Solved(), neos_opt="cbc")
    assert results == {"solver": "neos", "opt": "cbc"}
