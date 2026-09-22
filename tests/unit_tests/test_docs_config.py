# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

r"""`mkdocs.yml` must not silently lose a plugin to a duplicate key.

Duplicate mapping keys are legal YAML and the last one wins, without a
warning from either the parser or `mkdocs --strict`. A second `plugins:`
block therefore does not extend the list, it replaces it.

That happened: the API-documentation work appended its own `plugins:`
block, which discarded the `bibtex` plugin above it. The strict build
stayed green, the site kept deploying, and the Scientific Work page
rendered a literal `\bibliography` where the reference list should have
been — found by a reader, not by CI.
"""

import os
import re

import pytest

REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)
)
MKDOCS = os.path.join(REPO_ROOT, "mkdocs.yml")


def _text():
    with open(MKDOCS, encoding="utf-8") as handle:
        return handle.read()


def _top_level_keys(text):
    """Every key at column zero, in order, duplicates included."""
    return re.findall(r"^([A-Za-z_][\w-]*):", text, re.M)


def test_no_duplicate_top_level_keys():
    """The bug that lost the bibliography, in general form."""
    keys = _top_level_keys(_text())
    duplicates = {key for key in keys if keys.count(key) > 1}
    assert not duplicates, (
        f"{sorted(duplicates)} appear more than once at the top level of "
        f"mkdocs.yml. YAML keeps only the last, so the earlier block is "
        f"discarded silently."
    )


def test_every_configured_plugin_survives_parsing():
    """What the file says and what mkdocs sees have to agree."""
    yaml = pytest.importorskip("yaml")

    text = _text()
    written = set(re.findall(r"^  - (\w+):?$", text, re.M))

    class Loader(yaml.SafeLoader):
        """Ignores the Python-object tags mkdocs extensions use."""

    Loader.add_multi_constructor(
        "tag:yaml.org,2002:python/name:", lambda loader, suffix, node: suffix
    )
    config = yaml.load(text, Loader=Loader)

    parsed = set()
    for entry in config.get("plugins", []):
        parsed.add(entry if isinstance(entry, str) else next(iter(entry)))

    plugin_names = written & {"search", "bibtex", "mkdocstrings"}
    missing = plugin_names - parsed
    assert not missing, (
        f"{sorted(missing)} is configured in mkdocs.yml but absent after "
        f"parsing, so mkdocs never loads it."
    )


def test_bibliography_directive_has_a_plugin_to_render_it():
    r"""A page using `\bibliography` needs the bibtex plugin loaded."""
    docs = os.path.join(REPO_ROOT, "docs")
    users = []
    for root, _dirs, files in os.walk(docs):
        for name in files:
            if not name.endswith(".md"):
                continue
            path = os.path.join(root, name)
            with open(path, encoding="utf-8") as handle:
                if "\\bibliography" in handle.read():
                    users.append(os.path.relpath(path, REPO_ROOT))
    if not users:
        pytest.skip("no page uses the \\bibliography directive")

    text = _text()
    assert "- bibtex:" in text, (
        f"{users} use \\bibliography but mkdocs.yml configures no bibtex "
        f"plugin, so the directive renders as literal text."
    )
    keys = _top_level_keys(text)
    assert keys.count("plugins") == 1, (
        "bibtex is configured under a duplicated `plugins:` key, so it may "
        "not survive parsing."
    )
