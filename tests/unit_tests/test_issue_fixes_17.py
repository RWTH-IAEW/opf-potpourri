# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""Regression tests for GitLab issue #17.

Device placement used the unseeded global NumPy RNG, so two runs of the same
script equipped different buses and gave different results. These tests pin the
three properties that make it reproducible: a fixed default, independence from
the global RNG, and a caller-supplied seed or generator.
"""

from __future__ import annotations

import copy
import warnings

import numpy as np
import pytest
import simbench as sb

from potpourri.models_multi_period.ACOPF_multi_period import (
    ACOPF_multi_period,
)
from potpourri.technologies.battery import Battery_multi_period
from potpourri.technologies.flexibility import (
    DEFAULT_PLACEMENT_SEED,
    Flexibility_multi_period,
)
from potpourri.technologies.heat_pump import Heatpump_multi_period
from potpourri.technologies.pv import PV_multi_period

warnings.filterwarnings("ignore")

# Every device that places itself randomly, with a penetration high enough that
# the draw has real freedom on a 15-bus network.
PLACING_DEVICES = [
    Battery_multi_period,
    PV_multi_period,
    Heatpump_multi_period,
]
PENETRATION = 30.0


@pytest.fixture(scope="module")
def net():
    return sb.get_simbench_net("1-LV-rural1--0-sw")


@pytest.fixture(scope="module")
def model_net(net):
    """A net that has been through a multi-period constructor.

    The devices read `net.profiles` and `net._ppc`, which the model sets up.
    """
    return ACOPF_multi_period(copy.deepcopy(net), toT=3).net


def _buses(device):
    return sorted(int(b) for b in device.random_indexes)


def _make(cls, model_net, **kw):
    return cls(model_net, T=3, penetration=PENETRATION, **kw)


@pytest.mark.parametrize("cls", PLACING_DEVICES, ids=lambda c: c.__name__)
def test_issue17_placement_is_reproducible_by_default(cls, model_net):
    """Two identical constructions must equip the same buses.

    Before, `np.random.choice` on the global generator gave a different draw
    every time, so no study could be re-run.
    """
    first = _buses(_make(cls, model_net))
    second = _buses(_make(cls, model_net))
    third = _buses(_make(cls, model_net))
    assert first == second == third
    assert first, "penetration too low to place anything"


@pytest.mark.parametrize("cls", PLACING_DEVICES, ids=lambda c: c.__name__)
def test_issue17_placement_ignores_the_global_rng(cls, model_net):
    """The draw must not depend on the process-global NumPy state.

    Seeding `np.random` was the only workaround before, which coupled the
    placement to anything else in the process that happened to draw first.
    """
    np.random.seed(1)
    first = _buses(_make(cls, model_net))
    np.random.seed(999_999)
    second = _buses(_make(cls, model_net))
    assert first == second


@pytest.mark.parametrize("cls", PLACING_DEVICES, ids=lambda c: c.__name__)
def test_issue17_seed_selects_the_placement(cls, model_net):
    """A seed has to actually change the draw, or sampling is impossible."""
    a = _buses(_make(cls, model_net, seed=1))
    b = _buses(_make(cls, model_net, seed=2))
    assert a != b
    # And it is repeatable for a given seed.
    assert a == _buses(_make(cls, model_net, seed=1))


@pytest.mark.parametrize("cls", PLACING_DEVICES, ids=lambda c: c.__name__)
def test_issue17_default_seed_matches_the_documented_constant(cls, model_net):
    explicit = _buses(_make(cls, model_net, seed=DEFAULT_PLACEMENT_SEED))
    assert _buses(_make(cls, model_net)) == explicit


def test_issue17_shared_generator_gives_an_independent_replayable_sweep(
    model_net,
):
    """The Monte Carlo case: successive draws differ, the sweep replays.

    A single seed cannot express this — every device would land on the same
    buses — which is why an `rng` argument exists alongside `seed`.
    """

    def sweep():
        rng = np.random.default_rng(7)
        return [
            _buses(_make(Battery_multi_period, model_net, rng=rng))
            for _ in range(3)
        ]

    first = sweep()
    assert len({tuple(s) for s in first}) == len(first), (
        "scenarios in one sweep must differ from each other"
    )
    assert sweep() == first, "the sweep as a whole must replay"


def test_issue17_rng_takes_precedence_over_seed(model_net):
    from_rng = _buses(
        _make(
            Battery_multi_period,
            model_net,
            seed=12345,
            rng=np.random.default_rng(7),
        )
    )
    expected = _buses(
        _make(Battery_multi_period, model_net, rng=np.random.default_rng(7))
    )
    assert from_rng == expected


def test_issue17_no_device_draws_from_the_global_rng():
    """Guard the whole package, not only the three modules known to have done it.

    A bare `np.random.<fn>` call anywhere in the device layer reintroduces the
    bug, so this fails on the call rather than waiting for a study to disagree
    with itself.
    """
    import ast
    from pathlib import Path

    root = Path(__file__).resolve().parents[2] / "src" / "potpourri"
    offenders = []
    for py in root.rglob("*.py"):
        tree = ast.parse(py.read_text(), filename=str(py))
        for node in ast.walk(tree):
            # Matches np.random.<anything>, but not np.random.default_rng(...)
            # which constructs an independent generator.
            if not isinstance(node, ast.Attribute):
                continue
            value = node.value
            if (
                isinstance(value, ast.Attribute)
                and value.attr == "random"
                and isinstance(value.value, ast.Name)
                and value.value.id == "np"
                and node.attr != "default_rng"
            ):
                offenders.append(
                    f"{py.name}:{node.lineno} np.random.{node.attr}"
                )

    assert offenders == [], (
        "device placement must draw from a local Generator, not the global "
        f"NumPy state: {offenders}"
    )


def test_issue17_base_class_exposes_the_draw_helper():
    """Subclasses should not each re-implement the draw."""
    assert hasattr(Flexibility_multi_period, "draw_placement")
    assert isinstance(DEFAULT_PLACEMENT_SEED, int)
