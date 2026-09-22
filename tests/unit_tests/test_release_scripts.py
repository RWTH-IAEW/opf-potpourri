# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT

"""The release guard has to fail on a bad tag, not skip.

`.github/scripts/check_release_version.py` is the only thing standing between
a mistyped or stale version and an irreversible publish, and it runs exactly
once per release — on the tag push, when nobody is watching. These tests run
it against the real tree and against deliberately wrong tags.

`changelog_section.py` supplies the release notes; a silent failure there
would publish an empty release body.
"""

import os
import subprocess
import sys

import pytest

REPO_ROOT = os.path.abspath(
    os.path.join(os.path.dirname(__file__), os.pardir, os.pardir)
)
SCRIPTS = os.path.join(REPO_ROOT, ".github", "scripts")
CHECK = os.path.join(SCRIPTS, "check_release_version.py")
NOTES = os.path.join(SCRIPTS, "changelog_section.py")


def run(script, tag=None, cwd=REPO_ROOT):
    """Run a release script the way the workflow does, from the repo root."""
    env = dict(os.environ)
    env.pop("GITHUB_REF_NAME", None)
    if tag is None:
        env.pop("RELEASE_TAG", None)
    else:
        env["RELEASE_TAG"] = tag
    return subprocess.run(
        [sys.executable, script],
        cwd=cwd,
        env=env,
        capture_output=True,
        text=True,
    )


def current_version():
    tomllib = pytest.importorskip("tomllib")
    with open(os.path.join(REPO_ROOT, "pyproject.toml"), "rb") as handle:
        return tomllib.load(handle)["project"]["version"]


@pytest.fixture(autouse=True)
def _needs_yaml():
    pytest.importorskip("yaml", reason="the release guard parses CITATION.cff")


def test_guard_accepts_the_version_in_the_tree():
    result = run(CHECK, f"v{current_version()}")
    assert result.returncode == 0, result.stderr
    assert "agrees with the tree" in result.stdout


def test_guard_accepts_a_tag_without_the_v_prefix():
    """`RELEASE_TAG=0.6.0` and `v0.6.0` have to mean the same thing."""
    assert run(CHECK, current_version()).returncode == 0


def test_guard_rejects_a_version_that_is_not_in_the_tree():
    """The failure mode this exists for: it must fail, not skip."""
    result = run(CHECK, "v99.99.99")
    assert result.returncode == 1
    assert "Refusing to publish" in result.stderr
    # every source of truth has to be named, not just the first mismatch
    assert "pyproject.toml" in result.stderr
    assert "CITATION.cff" in result.stderr
    assert "CHANGELOG.md" in result.stderr


def test_guard_needs_a_tag():
    assert run(CHECK, None).returncode == 2


def test_guard_catches_a_stale_citation_file(tmp_path):
    """The 0.6.0 regression: pyproject bumped, CITATION.cff left behind.

    This is the shape that got through before — the PyPI job read the
    correct version while the Zenodo job stamped the archive with the
    previous one, and the workflow went green either way.
    """
    import shutil

    version = current_version()
    for name in ("pyproject.toml", "CITATION.cff", "CHANGELOG.md"):
        shutil.copy(os.path.join(REPO_ROOT, name), tmp_path / name)
    scripts = tmp_path / ".github" / "scripts"
    scripts.mkdir(parents=True)
    shutil.copy(CHECK, scripts / os.path.basename(CHECK))

    cff = tmp_path / "CITATION.cff"
    cff.write_text(
        cff.read_text(encoding="utf-8").replace(
            f'version: "{version}"', 'version: "0.0.1"', 1
        ),
        encoding="utf-8",
    )

    result = run(
        str(scripts / os.path.basename(CHECK)), f"v{version}", cwd=tmp_path
    )
    assert result.returncode == 1, result.stdout
    assert "CITATION.cff" in result.stderr
    assert "0.0.1" in result.stderr
    # pyproject agreed, so it must not be reported as a problem
    assert "pyproject.toml says" not in result.stderr


def test_notes_come_from_the_changelog_section():
    version = current_version()
    result = run(NOTES, f"v{version}")
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip(), "release notes must not be empty"
    # the body stops at the next release heading rather than running on
    assert "\n## [" not in result.stdout
    assert f"blob/v{version}/CHANGELOG.md" in result.stdout


def test_notes_fail_loudly_for_an_unreleased_version():
    result = run(NOTES, "v99.99.99")
    assert result.returncode == 1
    assert "no '## [99.99.99]' section" in result.stderr
