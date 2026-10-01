# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Shared pytest fixtures for the potpourri test suite."""

import os
import socket
import xmlrpc.client

import pytest
import pandapower as pp
import simbench as sb
from pyomo.opt.parallel.manager import ActionManagerError

#: Seconds to wait for the NEOS server to answer a ping before the NEOS
#: tests are skipped as unreachable.
NEOS_PING_TIMEOUT_S = 15


def _neos_reachable() -> bool:
    """Whether the public NEOS server answers a ping within the timeout.

    Returns:
        True when `kestrelAMPL().neos.ping()` reports the server alive;
        False on any failure (no network, DNS, timeout, XML-RPC fault).
    """
    previous = socket.getdefaulttimeout()
    socket.setdefaulttimeout(NEOS_PING_TIMEOUT_S)
    try:
        from pyomo.neos.kestrel import kestrelAMPL

        kestrel = kestrelAMPL()
        return kestrel.neos is not None and "alive" in str(kestrel.neos.ping())
    except Exception:  # noqa: BLE001 - any failure means "not reachable"
        return False
    finally:
        socket.setdefaulttimeout(previous)


@pytest.fixture(scope="session")
def neos():
    """Solve a model on NEOS, skipping the test when the service fails us.

    The NEOS tests check the transport to the public NEOS server, not the
    models -- every model they submit is solved locally elsewhere in the
    suite. A server that is unreachable, or that accepts the job and then
    returns no solution (its reply is then ``ERROR: An error occured with
    your submission`` with an empty ``.sol`` file, which Pyomo surfaces as
    an `ActionManagerError`), is a service problem, not a defect in this
    package, so the test is skipped with the reason instead of failing the
    suite. Seen on 2026-10-01, when NEOS's CPLEX queue rejected every
    submission while its IPOPT and CBC queues kept working.

    Returns:
        A callable ``run(model_obj, **solve_kwargs)`` that calls
        ``model_obj.solve(solver="neos", **solve_kwargs)`` and returns the
        results, or skips the test.
    """
    os.environ.setdefault("NEOS_EMAIL", "test@example.com")
    if not _neos_reachable():
        pytest.skip("NEOS server not reachable")

    def run(model_obj, **solve_kwargs):
        """Submit to NEOS; skip the calling test on a service-side failure.

        Args:
            model_obj: A potpourri model object with a `solve` method.
            **solve_kwargs: Forwarded to `solve` (e.g. `neos_opt`).

        Returns:
            The solver results.
        """
        try:
            return model_obj.solve(solver="neos", **solve_kwargs)
        except (ActionManagerError, xmlrpc.client.Error, OSError) as err:
            pytest.skip(f"NEOS did not return a solution: {err}")

    return run


@pytest.fixture(scope="session")
def four_bus_net():
    """Minimal pandapower four-bus network for fast model-construction tests."""
    return pp.networks.simple_four_bus_system()


@pytest.fixture(scope="session")
def lv_rural_net():
    """SimBench 1-LV-rural1--0-sw network (15 buses, 4 PV sgens)."""
    return sb.get_simbench_net("1-LV-rural1--0-sw")
