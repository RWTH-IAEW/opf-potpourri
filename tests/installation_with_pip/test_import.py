# tests/installation_with_pip/test_installed_package.py
import importlib
import pathlib


def test_package_can_be_imported():
    pkg = importlib.import_module("potpourri")
    assert pkg is not None


def test_package_has_file():
    pkg = importlib.import_module("potpourri")
    assert pathlib.Path(pkg.__file__).exists()