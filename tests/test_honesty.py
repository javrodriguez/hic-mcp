"""Invariant sweeps over every tracked file: the claims this repo may not make.

These are absolutes, so they are swept across ALL tracked text rather than the one
place a violation is expected. They exist because the cost of a false claim in a
public repo is not a failing test - it is someone believing it.
"""

import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]

# Claims this project is not entitled to make anywhere a reader will see them.
FORBIDDEN = {
    "production-readiness": r"production[\s-]?(ready|grade)",
    "biographical claim": (
        r"\b(PhD|postdoc|years of experience|my research|author of|"
        r"decade of (research|experience|work)|NYU|Langone)\b"
    ),
}
# Priority claims about this project ("the first MCP server to...", "the only Hi-C tool")
# are allowed only inside an explicitly dated form. Ordinary ordinal English ("the first
# two diagonals") is not a claim and must not trip this.
SUPERLATIVE = (
    r"\b(the first|the only|world'?s first)\b[^.\n]{0,40}"
    r"\b(MCP|server|tool|package|library|implementation)\b(?!'s)"  # not "the first tool's ..."
)
DATE_SCOPED = r"as of \d{4}"


def tracked_text_files() -> list[Path]:
    """Tracked files AND new files not yet committed.

    Untracked-but-not-ignored files are included deliberately: a sweep that only saw
    committed files would pass locally on a brand-new file and only fail later in CI,
    which is exactly how these checks were first caught out.
    """
    out = subprocess.run(
        ["git", "-C", str(REPO), "ls-files", "--cached", "--others", "--exclude-standard"],
        capture_output=True,
        text=True,
        check=True,
    ).stdout.split()
    keep = {".md", ".py", ".toml", ".yml", ".yaml", ".json", ".jsonl", ".cff", ".txt"}
    files = [REPO / p for p in out if Path(p).suffix in keep]
    # this file names the forbidden patterns in order to test for them
    return [f for f in files if f.exists() and f.name != "test_honesty.py"]


def test_tracked_files_were_found():
    """Guard the guard: an empty sweep would pass every test below vacuously."""
    files = tracked_text_files()
    assert len(files) > 10
    assert any(f.name == "README.md" for f in files)


@pytest.mark.parametrize("label,pattern", list(FORBIDDEN.items()))
def test_no_forbidden_claims_anywhere(label, pattern):
    hits = []
    for f in tracked_text_files():
        for i, line in enumerate(f.read_text(errors="ignore").splitlines(), 1):
            if re.search(pattern, line, re.IGNORECASE):
                hits.append(f"{f.relative_to(REPO)}:{i}: {line.strip()}")
    assert not hits, f"{label} found:\n" + "\n".join(hits)


def test_superlative_claims_are_date_scoped():
    """A 'first/only' claim without a date is a claim that quietly expires."""
    hits = []
    for f in tracked_text_files():
        for i, line in enumerate(f.read_text(errors="ignore").splitlines(), 1):
            if re.search(SUPERLATIVE, line, re.IGNORECASE) and not re.search(
                DATE_SCOPED, line, re.IGNORECASE
            ):
                hits.append(f"{f.relative_to(REPO)}:{i}: {line.strip()}")
    assert not hits, "undated superlative claim:\n" + "\n".join(hits)


def test_no_local_machine_paths_leak():
    """Covers the captured session log too - the file most likely to carry one."""
    hits = []
    for f in tracked_text_files():
        for i, line in enumerate(f.read_text(errors="ignore").splitlines(), 1):
            # a real home path, not a pattern that merely mentions one
            if re.search(r"/Users/[a-z]|/home/[a-z]", line, re.IGNORECASE):
                hits.append(f"{f.relative_to(REPO)}:{i}: {line.strip()[:160]}")
    assert not hits, "machine-local path in a tracked file:\n" + "\n".join(hits)


def test_captured_session_is_a_real_tool_driven_session():
    """The demo log must show this server actually being called, not a model reciting."""
    import json

    log = REPO / "demo" / "raw-session.jsonl"
    records = [json.loads(line) for line in log.read_text().splitlines() if line.strip()]
    calls = [
        c["name"]
        for r in records
        if r.get("type") == "assistant"
        for c in r["message"]["content"]
        if c.get("type") == "tool_use" and str(c.get("name", "")).startswith("mcp__hic-mcp__")
    ]
    assert len(calls) >= 2, f"session made too few hic-mcp calls: {calls}"
    assert any(r.get("mcp_servers") for r in records), "no server connection recorded"


def test_readme_structure_and_config_parity():
    """The client config is early, and the README and example file cannot drift apart."""
    result = subprocess.run(
        [sys.executable, str(REPO / "scripts" / "check_readme.py")],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr


def test_license_scopes_itself_to_code():
    """The data is not ours to license; the LICENSE must say so."""
    text = " ".join((REPO / "LICENSE").read_text().split())  # normalize line wrapping
    assert "source code in this repository only" in text
    assert "4DN Data Release and Use Policy" in text


def test_every_tool_is_documented_in_the_readme_table():
    from hic_mcp import analysis

    readme = (REPO / "README.md").read_text()
    for tool in (
        "matrix_summary",
        "contacts_at_locus",
        "insulation_tads",
        "compartments",
        "virtual_4c",
        "expected_observed",
    ):
        assert hasattr(analysis, tool), f"{tool} is documented but does not exist"
        assert f"`{tool}`" in readme, f"{tool} exists but is not in the README table"
