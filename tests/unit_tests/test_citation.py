# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""CITATION.cff has to stay in step with the version being released.

The Zenodo archive job builds its metadata from CITATION.cff while the
PyPI job reads pyproject.toml, and only the latter is checked against the
git tag. A stale `version:` here therefore publishes the release under the
wrong number on the concept record, silently — which is exactly what
happened when 0.6.0 was cut and caught by hand.
"""

import os

import pytest

REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)
)
CITATION = os.path.join(REPO_ROOT, "CITATION.cff")
PYPROJECT = os.path.join(REPO_ROOT, "pyproject.toml")


def _citation():
    yaml = pytest.importorskip("yaml", reason="PyYAML parses CITATION.cff")
    with open(CITATION, encoding="utf-8") as handle:
        return yaml.safe_load(handle)


def _pyproject():
    # tomllib is 3.11+; the project still supports 3.10, where this skips
    # rather than failing the whole suite on a missing stdlib module.
    tomllib = pytest.importorskip("tomllib")
    with open(PYPROJECT, "rb") as handle:
        return tomllib.load(handle)


def test_citation_file_exists_and_parses():
    assert os.path.exists(CITATION), "CITATION.cff is missing"
    cff = _citation()
    assert cff["cff-version"] == "1.2.0"
    assert cff["type"] == "software"


def test_citation_version_matches_pyproject():
    """The two version strings a release reads must agree."""
    assert _citation()["version"] == _pyproject()["project"]["version"]


def test_citation_licence_matches_pyproject():
    assert _citation()["license"] == _pyproject()["project"]["license"]


def test_citation_names_every_package_author():
    """An author on the package is an author on the citation.

    The two lists are maintained in different files and drifted apart
    before; pyproject carries given and family names in one string, so
    compare on the family name.
    """
    cff_families = {a["family-names"] for a in _citation()["authors"]}
    for author in _pyproject()["project"]["authors"]:
        family = author["name"].split()[-1]
        assert family in cff_families, (
            f"{author['name']} is not in CITATION.cff"
        )


def test_citation_carries_the_concept_doi():
    """The DOI cited by the README and the paper has to stay resolvable."""
    dois = [
        i["value"]
        for i in _citation().get("identifiers", [])
        if i.get("type") == "doi"
    ]
    assert "10.5281/zenodo.20357011" in dois
