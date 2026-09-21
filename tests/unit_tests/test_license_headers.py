# SPDX-FileCopyrightText: 2023-2026 Institute for High Voltage Equipment and Grids, Digitalization and Energy Economics (IAEW), RWTH Aachen University
#
# SPDX-License-Identifier: MIT
#
# REUSE-IgnoreStart -- below this line SPDX tags appear as data
# (policy constants, test fixtures, help text), not as declarations
# for this file. See docs/licensing.md.

"""The licensing-header gate must be precise about *where* a notice is.

A notice only counts when it sits in the physical Python comment header.
The same text in a docstring, a string literal or a documentation
example must not satisfy the check -- that is the whole point of the
tool, so it is tested from several angles.
"""

from __future__ import annotations

import os
import subprocess
import sys

import pytest

REPO_ROOT = os.path.dirname(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
)
sys.path.insert(0, os.path.join(REPO_ROOT, "tools"))

import license_headers as lh  # noqa: E402

CR = lh.FIRST_PARTY_COPYRIGHT
GOOD = f"# SPDX-FileCopyrightText: {CR}\n#\n# SPDX-License-Identifier: MIT\n"


def compliant(body: str = "x = 1\n") -> str:
    return GOOD + "\n" + body


# ---------------------------------------------------------------- check


def test_plain_module_without_a_header_is_flagged():
    problems = lh.check_text("a.py", '"""Doc."""\n\nx = 1\n')
    assert any("SPDX-FileCopyrightText" in p for p in problems)
    assert any("SPDX-License-Identifier" in p for p in problems)


def test_compliant_module_passes():
    assert lh.check_text("a.py", compliant()) == []


def test_empty_init_is_flagged_then_passes_once_fixed():
    assert lh.check_text("pkg/__init__.py", "") != []
    fixed, note = lh.fix_text("pkg/__init__.py", "")
    assert note == "added"
    assert lh.check_text("pkg/__init__.py", fixed) == []


def test_spdx_inside_a_docstring_does_not_count():
    text = (
        '"""Docs.\n\n'
        f"# SPDX-FileCopyrightText: {CR}\n"
        "# SPDX-License-Identifier: MIT\n"
        '"""\n\nx = 1\n'
    )
    problems = lh.check_text("a.py", text)
    assert any("does not count" in p for p in problems)


def test_spdx_inside_a_string_literal_does_not_count():
    text = 'TEMPLATE = "# SPDX-License-Identifier: MIT"\n'
    problems = lh.check_text("a.py", text)
    assert any("SPDX-License-Identifier" in p for p in problems)
    assert lh.parse_header(text).licenses == []


def test_comment_after_code_is_not_a_header():
    text = "x = 1\n# SPDX-License-Identifier: MIT\n"
    assert lh.parse_header(text).licenses == []


def test_duplicate_identical_licence_lines_are_flagged():
    text = (
        f"# SPDX-FileCopyrightText: {CR}\n"
        "# SPDX-License-Identifier: MIT\n"
        "# SPDX-License-Identifier: MIT\n\nx = 1\n"
    )
    assert any("duplicated" in p for p in lh.check_text("a.py", text))


def test_contradictory_licence_lines_are_flagged():
    text = (
        f"# SPDX-FileCopyrightText: {CR}\n"
        "# SPDX-License-Identifier: MIT\n"
        "# SPDX-License-Identifier: GPL-3.0-or-later\n\nx = 1\n"
    )
    assert any("contradictory" in p for p in lh.check_text("a.py", text))


def test_malformed_expression_is_flagged():
    text = (
        f"# SPDX-FileCopyrightText: {CR}\n"
        "# SPDX-License-Identifier: Not-A-Licence\n\nx = 1\n"
    )
    assert any("invalid SPDX" in p for p in lh.check_text("a.py", text))


def test_wrong_but_valid_licence_is_flagged():
    text = (
        f"# SPDX-FileCopyrightText: {CR}\n"
        "# SPDX-License-Identifier: Apache-2.0\n\nx = 1\n"
    )
    assert any(
        "policy for this file" in p for p in lh.check_text("a.py", text)
    )


def test_empty_copyright_line_is_flagged():
    text = "# SPDX-FileCopyrightText:\n# SPDX-License-Identifier: MIT\n\nx=1\n"
    assert any("empty" in p for p in lh.check_text("a.py", text))


def test_first_party_file_missing_the_canonical_holder_is_flagged():
    text = (
        "# SPDX-FileCopyrightText: 2026 Someone Else\n"
        "# SPDX-License-Identifier: MIT\n\nx = 1\n"
    )
    assert any("canonical notice" in p for p in lh.check_text("a.py", text))


def test_multiple_holders_are_accepted():
    text = (
        f"# SPDX-FileCopyrightText: {CR}\n"
        "# SPDX-FileCopyrightText: 2023 Steffen Kortmann\n"
        "#\n# SPDX-License-Identifier: MIT\n\nx = 1\n"
    )
    assert lh.check_text("a.py", text) == []


def test_personal_copyright_claim_in_a_docstring_is_flagged():
    """Copyright is institutional and belongs in the SPDX header.

    A `(c) YEAR, Name` line in prose is a copyright claim, so it is a
    finding even when the header itself is correct.
    """
    text = GOOD + '\n"""Doc.\n\n(c) 2023, Steffen Kortmann\n"""\n'
    problems = lh.check_text("a.py", text)
    assert any("asserts a copyright claim" in p for p in problems)


def test_authorship_line_in_a_docstring_is_fine():
    """Attribution is not a copyright claim, so it must not be flagged."""
    text = GOOD + '\n"""Doc.\n\nAuthor: Steffen Kortmann (2023)\n"""\n'
    assert lh.check_text("a.py", text) == []


def test_third_party_exception_is_honoured(monkeypatch):
    monkeypatch.setitem(
        lh.LICENSE_EXCEPTIONS, "v/x.py", ("BSD-3-Clause", "upstream copy")
    )
    text = (
        "# SPDX-FileCopyrightText: 1996-2015 PSERC\n"
        "# SPDX-License-Identifier: BSD-3-Clause\n\nx = 1\n"
    )
    assert lh.check_text("v/x.py", text) == []
    # ... and the first-party licence is then the wrong answer for it.
    assert lh.check_text("v/x.py", compliant()) != []


# ------------------------------------------------------------------ fix


@pytest.mark.parametrize(
    "body",
    [
        "x = 1\n",
        '"""Doc."""\n\nx = 1\n',
        "from __future__ import annotations\n\nx = 1\n",
        '"""Doc."""\n\nfrom __future__ import annotations\n\nx = 1\n',
        "",
    ],
)
def test_fix_then_check_passes_and_is_idempotent(body):
    once, note = lh.fix_text("a.py", body)
    assert note == "added"
    assert lh.check_text("a.py", once) == []
    twice, note2 = lh.fix_text("a.py", once)
    assert note2 == ""
    assert twice == once


def test_fix_preserves_shebang_and_encoding_and_bom():
    text = (
        "﻿#!/usr/bin/env python3\n"
        "# -*- coding: utf-8 -*-\n"
        '"""Doc."""\n\nx = 1\n'
    )
    fixed, _ = lh.fix_text("a.py", text)
    lines = fixed.split("\n")
    assert lines[0] == "﻿#!/usr/bin/env python3"
    assert lines[1] == "# -*- coding: utf-8 -*-"
    assert lines[2].startswith("# SPDX-FileCopyrightText:")
    assert lh.check_text("a.py", fixed) == []
    assert fixed.count("﻿") == 1


def test_fix_preserves_docstring_and_future_import_legality():
    import ast

    text = '"""Doc."""\n\nfrom __future__ import annotations\n\nx = 1\n'
    fixed, _ = lh.fix_text("a.py", text)
    tree = ast.parse(fixed)  # a misplaced __future__ import is a SyntaxError
    assert ast.get_docstring(tree) == "Doc."


def test_fix_keeps_other_leading_comments():
    text = "# ruff: noqa: E501\n# type: ignore\n\nx = 1\n"
    fixed, _ = lh.fix_text("a.py", text)
    assert "# ruff: noqa: E501" in fixed
    assert "# type: ignore" in fixed
    assert lh.check_text("a.py", fixed) == []


def test_fix_refuses_a_docstring_copyright_claim():
    """Retiring someone's copyright notice is not a mechanical edit."""
    text = '"""Doc.\n\n(c) 2023, Steffen Kortmann\n"""\n\nx = 1\n'
    out, note = lh.fix_text("a.py", text)
    assert note.startswith("refused")
    assert out == text


def test_fix_leaves_an_authorship_line_alone():
    text = '"""Doc.\n\nAuthor: Steffen Kortmann (2023)\n"""\n\nx = 1\n'
    fixed, note = lh.fix_text("a.py", text)
    assert note == "added"
    assert "Author: Steffen Kortmann (2023)" in fixed
    assert f"# SPDX-FileCopyrightText: {CR}" in fixed
    assert lh.check_text("a.py", fixed) == []
    assert lh.fix_text("a.py", fixed)[1] == ""


def test_fix_refuses_a_foreign_licence():
    text = (
        "# SPDX-FileCopyrightText: 1996-2015 PSERC\n"
        "# SPDX-License-Identifier: BSD-3-Clause\n\nx = 1\n"
    )
    out, note = lh.fix_text("a.py", text)
    assert note.startswith("refused")
    assert out == text


def test_fix_refuses_an_unverifiable_holder():
    text = (
        "# SPDX-FileCopyrightText: 2011 Some Upstream Author\n"
        "# SPDX-License-Identifier: MIT\n\nx = 1\n"
    )
    out, note = lh.fix_text("a.py", text)
    assert note.startswith("refused")
    assert out == text


def test_fix_refuses_declared_third_party_files(monkeypatch):
    monkeypatch.setitem(
        lh.LICENSE_EXCEPTIONS, "v/x.py", ("BSD-3-Clause", "upstream copy")
    )
    out, note = lh.fix_text("v/x.py", "x = 1\n")
    assert note.startswith("refused")
    assert out == "x = 1\n"


# ------------------------------------------------------------ discovery


def test_discovery_finds_paths_containing_spaces(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    odd = tmp_path / "a dir with spaces"
    odd.mkdir()
    (odd / "a module.py").write_text("x = 1\n")
    (tmp_path / "plain.py").write_text("x = 1\n")
    found = lh.python_files(str(tmp_path))
    assert "a dir with spaces/a module.py" in found
    assert "plain.py" in found


def test_discovery_covers_pyi_and_pyw(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    for name in ("a.py", "b.pyi", "c.pyw"):
        (tmp_path / name).write_text("x = 1\n")
    assert lh.python_files(str(tmp_path)) == ["a.py", "b.pyi", "c.pyw"]


@pytest.mark.parametrize(
    "layout",
    [
        ("env/lib/python3.12/site-packages/dep", None),
        ("env/lib/dist-packages/dep", None),
        ("env/lib/dep", "env/pyvenv.cfg"),
        ("env/lib/dep", "env/conda-meta/history"),
    ],
)
def test_discovery_skips_untracked_environments_in_the_tree(tmp_path, layout):
    """CI builds its conda env inside the working tree (see .gitignore).

    An un-ignored environment must not drag thousands of installed
    third-party modules into the audit.
    """
    subdir, marker = layout
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    (tmp_path / "mine.py").write_text("x = 1\n")
    pkg = tmp_path / subdir
    pkg.mkdir(parents=True)
    (pkg / "vendored.py").write_text("x = 1\n")
    if marker:
        path = tmp_path / marker
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("")
    assert lh.python_files(str(tmp_path)) == ["mine.py"]


def test_discovery_still_audits_tracked_files_in_such_a_directory(tmp_path):
    """Skipping applies to untracked files only; a committed file counts."""
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    pkg = tmp_path / "vendor" / "site-packages"
    pkg.mkdir(parents=True)
    (pkg / "kept.py").write_text("x = 1\n")
    subprocess.run(
        ["git", "-C", str(tmp_path), "add", "vendor/site-packages/kept.py"],
        check=True,
    )
    assert lh.python_files(str(tmp_path)) == ["vendor/site-packages/kept.py"]


# ------------------------------------------------------------- encoding


def test_pep263_declared_encoding_is_honoured(tmp_path):
    path = tmp_path / "latin.py"
    path.write_bytes(
        "# -*- coding: latin-1 -*-\n# café ± degrees\nx = 1\n".encode(
            "latin-1"
        )
    )
    text, encoding = lh.read_source(str(path))
    assert encoding == "latin-1"
    assert "café ± degrees" in text


def test_fix_writes_back_in_the_files_own_encoding(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    path = tmp_path / "latin.py"
    path.write_bytes(
        "# -*- coding: latin-1 -*-\n# café ±\nx = 1\n".encode("latin-1")
    )
    lh.fix_repository(str(tmp_path), dry_run=False)
    raw = path.read_bytes()
    assert b"\xe9" in raw  # still latin-1 on disk, not re-encoded to UTF-8
    text, encoding = lh.read_source(str(path))
    assert encoding == "latin-1"
    assert lh.check_text("latin.py", text) == []
    lines = text.split("\n")
    assert lines[0] == "# -*- coding: latin-1 -*-"  # declaration stays first


def test_undecodable_file_is_reported_not_crashed(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    (tmp_path / "good.py").write_text(compliant())
    # Declares UTF-8 but is not: the scan must survive and name it.
    (tmp_path / "bad.py").write_bytes(b"# -*- coding: utf-8 -*-\nx = '\xb1'\n")
    failures = lh.check_repository(str(tmp_path))
    assert "bad.py" in failures
    assert "cannot decode" in failures["bad.py"][0]
    assert "good.py" not in failures


def test_discovery_raises_instead_of_reporting_an_empty_pass(tmp_path):
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    with pytest.raises(lh.DiscoveryError):
        lh.python_files(str(tmp_path))


def test_discovery_raises_outside_a_repository(tmp_path):
    with pytest.raises(lh.DiscoveryError):
        lh.python_files(str(tmp_path / "nope"))


def test_checker_reports_discovery_failure_as_error_not_success(tmp_path):
    sys.path.insert(0, os.path.join(REPO_ROOT, "tools"))
    import check_license_headers as gate

    assert gate.main(str(tmp_path / "missing")) == 2


# --------------------------------------------------------------- policy


def test_repository_is_compliant():
    """The gate must be green on this repository itself."""
    failures = lh.check_repository(REPO_ROOT)
    assert failures == {}, failures


def test_license_texts_do_not_drift():
    """LICENSES/MIT.txt is the REUSE copy of the root LICENSE."""
    root = open(os.path.join(REPO_ROOT, "LICENSE"), encoding="utf-8").read()
    reuse_copy = open(
        os.path.join(REPO_ROOT, "LICENSES", "MIT.txt"), encoding="utf-8"
    ).read()
    assert root == reuse_copy


def test_declared_licence_matches_pyproject():
    # tomllib is 3.11+; the project still supports 3.10, where this skips
    # rather than failing the whole suite on a missing stdlib module.
    tomllib = pytest.importorskip("tomllib")

    with open(os.path.join(REPO_ROOT, "pyproject.toml"), "rb") as handle:
        pyproject = tomllib.load(handle)
    assert pyproject["project"]["license"] == lh.FIRST_PARTY_LICENSE


# REUSE-IgnoreEnd
