import pandapower as pp
import pyomo.environ as pyo

from potpourri.models.basemodel import Basemodel
from potpourri.models.ACOPF_base import ACOPF
from potpourri.models.DCOPF import DCOPF


def _four_bus_net():
    return pp.networks.simple_four_bus_system()


def test_basemodel_rejects_non_pandapower_input():
    """Basemodel raises ValueError for non-pandapower input."""
    try:
        Basemodel("not a network")
        assert False, "expected ValueError"
    except ValueError:
        pass


def test_basemodel_deep_copies_network():
    """Basemodel must not mutate the original network."""
    net = _four_bus_net()
    original_bus_count = len(net.bus)
    acopf = ACOPF(net)
    acopf.add_OPF()
    assert len(net.bus) == original_bus_count
    assert net is not acopf.net


def test_basemodel_sets_populated():
    """Basemodel.create_model() populates all expected Pyomo sets."""
    net = _four_bus_net()
    acopf = ACOPF(net)

    required_sets = ['B', 'b0', 'D', 'L', 'TRANSF', 'G', 'sG', 'Dbs', 'Gbs']
    for name in required_sets:
        assert hasattr(acopf.model, name), f"missing Pyomo set '{name}'"

    assert len(list(acopf.model.B)) == len(net.bus)
    assert len(list(acopf.model.D)) == len(net.load[net.load.in_service])


def test_basemodel_demand_fixed():
    """Load variables are fixed to profile values before OPF constraints are added."""
    net = _four_bus_net()
    acopf = ACOPF(net)

    for d in acopf.model.D:
        assert acopf.model.pD[d].is_fixed(), f"demand variable pD[{d}] should be fixed"


def test_acopf_voltage_params_within_bounds():
    """Vmax must be strictly greater than Vmin for every bus."""
    net = _four_bus_net()
    acopf = ACOPF(net)
    acopf.add_OPF()

    for b in acopf.model.B:
        assert pyo.value(acopf.model.Vmax[b]) > pyo.value(acopf.model.Vmin[b]), (
            f"Vmax <= Vmin at bus {b}"
        )


def test_dcopf_line_limits_non_negative():
    """Thermal line limits SLmax must be non-negative."""
    net = _four_bus_net()
    dcopf = DCOPF(net)
    dcopf.add_OPF()

    for l in dcopf.model.L:
        assert pyo.value(dcopf.model.SLmax[l]) >= 0, f"negative SLmax for line {l}"


if __name__ == '__main__':
    test_basemodel_rejects_non_pandapower_input()
    test_basemodel_deep_copies_network()
    test_basemodel_sets_populated()
    test_basemodel_demand_fixed()
    test_acopf_voltage_params_within_bounds()
    test_dcopf_line_limits_non_negative()
