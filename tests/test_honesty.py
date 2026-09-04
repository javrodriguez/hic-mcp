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


def test_readme_sample_output_is_what_the_code_actually_returns():
    """A block labelled "Real output" must stay real as the code changes.

    It drifted once already: a rounding fix changed every value by a digit and the
    README kept the old numbers, which is exactly the small dishonesty this repo
    claims to design against.
    """
    import json
    import re

    from hic_mcp.analysis import insulation_tads

    readme = (REPO / "README.md").read_text()
    block = re.search(r"```json\n(\{\n  \"resolution_used\".*?)```", readme, re.DOTALL)
    assert block, "README no longer contains the sample-output block"
    claimed = json.loads(block.group(1))
    live = insulation_tads(region="chr17:65,000,000-67,000,000", top_n=2)
    assert claimed["resolution_used"] == live["resolution_used"]
    assert claimed["windows_bp"] == live["windows_bp"]
    assert claimed["ranked_by"] == live["ranked_by"]
    assert claimed.get("view") == live["view"]
    for shown, actual in zip(claimed["top_boundaries"], live["top_boundaries"], strict=True):
        # EVERY field, not a chosen subset: log2_insulation once drifted while
        # strength stayed put, because the guard checked strength and not log2
        assert set(shown) == set(actual), "README block and live output have different fields"
        for key, shown_value in shown.items():
            actual_value = actual[key]
            if isinstance(shown_value, float):
                assert shown_value == pytest.approx(actual_value, rel=1e-3), key
            else:
                assert shown_value == actual_value, key


def test_documented_landmarks_are_re_derived_from_live_code():
    """Every measured figure in PROVENANCE must still be what the code produces.

    Three "measured" artifacts once went stale in a single commit because a fix
    changed the computation and nobody re-derived them. This test re-derives them,
    so that class of drift fails the build instead of shipping.
    """
    import re

    import numpy as np
    import pandas as pd
    from cooltools import eigs_cis

    from hic_mcp.analysis import compartments, insulation_tads
    from hic_mcp.data import load_arms_view, load_gc_track, open_matrix, resolve_input_path

    doc = (REPO / "data" / "PROVENANCE.md").read_text()

    ins = insulation_tads(region="chr17:65,000,000-67,000,000", top_n=2)
    top = ins["top_boundaries"][0]
    line = re.search(r"- Strongest insulation boundary: (.+)", doc).group(1)
    assert top["locus"] in line
    assert f"strength {top['strength']}" in line
    assert f"log2 insulation {top['log2_insulation']}" in line

    clr = open_matrix(resolve_input_path(None), 100_000, default=100_000)
    gc = load_gc_track()
    _, vecs = eigs_cis(clr, phasing_track=gc, view_df=load_arms_view(), n_eigs=1)
    merged = pd.merge(vecs, gc, on=["chrom", "start", "end"]).dropna(subset=["E1", "GC"])
    r_live = float(np.corrcoef(merged["E1"], merged["GC"])[0, 1])
    r_doc = float(re.search(r"Pearson r = ([-0-9.]+) at 100 kb", doc).group(1))
    assert r_live == pytest.approx(r_doc, abs=5e-4)

    a = compartments(region="chr17:50,100,000-51,100,000")["region_mean_E1"]
    b = compartments(region="chr17:51,400,000-52,400,000")["region_mean_E1"]
    flip = re.search(r"- Adjacent compartment flip: (.+)", doc).group(1)
    assert f"mean E1 {a}" in flip and f"mean E1 {b}" in flip


def test_transcript_quotes_are_verbatim_from_the_raw_log():
    """Every quoted line in TRANSCRIPT.md must appear in the log word for word.

    A paraphrase presented as a quote is the exact failure this document exists to
    rule out, and one slipped through a hand check ("Now check the" for "Now checking
    its"). Prose lines are compared; table rows and the prompt are not, since the
    prompt is the user's and tables are explicitly marked as abridged.
    """
    import json

    log = REPO / "demo" / "raw-session.jsonl"
    said = []
    for line in log.read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        if rec.get("type") == "assistant":
            said += [
                c["text"] for c in rec["message"]["content"] if c.get("type") == "text"
            ]
    blob = " ".join(said).replace("**", "").replace("\u2013", "-").replace("\u2014", "-")

    prompt_marker = "Using the hic-mcp tools"
    quoted = []
    for raw in (REPO / "demo" / "TRANSCRIPT.md").read_text().splitlines():
        if not raw.startswith("> "):
            continue
        text = raw[2:].strip()
        if not text or text.startswith("|") or text.startswith("#") or text.startswith("-"):
            continue
        if prompt_marker in text:  # the prompt is the user's line, not the agent's
            continue
        quoted.append(text.replace("**", "").replace("\u2013", "-").replace("\u2014", "-"))

    assert quoted, "no prose quotes found - the check would pass vacuously"
    missing = [q for q in quoted if q not in blob]
    assert not missing, "quoted but not in the raw log:\n" + "\n".join(missing[:5])
