"""Shared pytest fixtures for the potpourri test suite."""

import pytest
import pandapower as pp
import simbench as sb


@pytest.fixture(scope="session")
def four_bus_net():
    """Minimal pandapower four-bus network for fast model-construction tests."""
    return pp.networks.simple_four_bus_system()


@pytest.fixture(scope="session")
def lv_rural_net():
    """SimBench 1-LV-rural1--0-sw network (15 buses, 4 PV sgens)."""
    return sb.get_simbench_net("1-LV-rural1--0-sw")
