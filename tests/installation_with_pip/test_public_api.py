# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""The installed distribution must report a coherent version.

Checks that `potpourri.__version__` exists, agrees with the
distribution metadata, and is not the source-tree fallback sentinel --
which would mean the check is reading the checkout instead of the
installed package.
"""

from importlib.metadata import version

from potpourri import __version__


def test_version_exists():
    assert isinstance(__version__, str)
    assert len(__version__) > 0


def test_version_matches_distribution_metadata():
    """`__version__` must be read from the installed distribution.

    A hard-coded literal here drifted two minor releases behind
    pyproject.toml before anyone noticed, because the only assertion was
    that the string was non-empty. Comparing against the distribution
    metadata makes a re-introduced literal fail immediately.
    """
    assert __version__ == version("opf-potpourri")


def test_version_is_not_the_source_tree_sentinel():
    """This suite runs against a pip-installed package, never a bare tree.

    `potpourri.__init__` falls back to "0.0.0.dev0" when the distribution
    metadata is missing. Seeing that here means the import resolved to a
    source checkout instead of the installation under test.
    """
    assert __version__ != "0.0.0.dev0"
