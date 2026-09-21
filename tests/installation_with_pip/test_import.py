# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""The installed distribution must import and carry its files.

Run against an installed potpourri from outside the repository,
so a packaging mistake -- a missing subpackage, a bad entry in
`tool.setuptools.packages.find` -- shows up here rather than at
a user's first import.
"""
#
# tests/installation_with_pip/test_installed_package.py

import importlib
import pathlib


def test_package_can_be_imported():
    pkg = importlib.import_module("potpourri")
    assert pkg is not None


def test_package_has_file():
    pkg = importlib.import_module("potpourri")
    assert pathlib.Path(pkg.__file__).exists()
