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
    """The data is not ours to license, and the repository must say so.

    The scope note lives in README.md and data/PROVENANCE.md rather than inside LICENSE,
    so that automated licence detection reads the source licence as plain MIT. The
    guarantee is unchanged and this test is stricter than the LICENSE-only version it
    replaces: the statement must be present, AND LICENSE must stay clean so the badge
    a reader sees matches the claim made about it.
    """
    readme = " ".join((REPO / "README.md").read_text().split())  # normalize line wrapping
    assert "source code only" in readme
    assert "4DN Data Release and Use Policy" in readme
    assert "not covered by the MIT licence" in readme

    provenance = " ".join((REPO / "data" / "PROVENANCE.md").read_text().split())
    assert "4DN" in provenance

    # LICENSE stays unmodified MIT: anything appended here breaks GitHub's detection and
    # makes the repo page contradict the README. This is the half that regressed once.
    license_text = " ".join((REPO / "LICENSE").read_text().split())
    assert license_text.startswith("MIT License")
    assert "Scope note" not in license_text
    assert license_text.rstrip().endswith("DEALINGS IN THE SOFTWARE.")


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
    block = re.search(r"```json\n(\{\n  \"region\".*?)```", readme, re.DOTALL)
    assert block, "README no longer contains the sample-output block"
    claimed = json.loads(block.group(1))
    live = insulation_tads(region="chr17:65,000,000-67,000,000", top_n=2)
    # the whole response, top level included: the block once dropped `region` and
    # `boundary_counts_per_window` while every per-boundary field matched, because
    # this check only ever looked inside the boundaries
    # the client sees the response MODEL, nulls included - comparing against the raw
    # dict with nulls stripped is what let the block omit scale_note while claiming
    # to be "the entire response, nothing removed"
    from hic_mcp import models

    live_shown = json.loads(models.InsulationTads(**live).model_dump_json())
    assert set(claimed) == set(live_shown), (
        f"README block and live response differ at the top level: "
        f"missing {sorted(set(live_shown) - set(claimed))}, "
        f"extra {sorted(set(claimed) - set(live_shown))}"
    )
    for key in set(claimed) - {"top_boundaries"}:
        assert claimed[key] == live_shown[key], key
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


def test_transcript_quotes_the_conclusion_completely_and_unaltered():
    """The conclusion block must equal the log's final answer exactly - no cuts at all.

    An earlier version declared its cuts in prose and asked a test to police the claim,
    but a test can verify that quoted lines EXIST in the log while being structurally
    blind to lines silently left out. Quoting the answer whole removes the class: the
    check below is equality, not containment.
    """
    import json

    said = []
    for line in (REPO / "demo" / "raw-session.jsonl").read_text().splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        if rec.get("type") == "assistant":
            said += [c["text"] for c in rec["message"]["content"] if c.get("type") == "text"]
    final = said[-1].strip()

    transcript = (REPO / "demo" / "TRANSCRIPT.md").read_text()
    marker = "**The agent's conclusion, in full"
    assert marker in transcript, "the conclusion block is no longer quoted in full"
    after = transcript[transcript.index(marker) :]
    body = after[after.index("\n\n") + 2 :]
    # bound on a real section separator - a line that is exactly '---'. Slicing on the
    # first '---' once stopped inside a markdown table rule ('---|---|---|'), leaving
    # eight duplicated lines below it unchecked while this test passed.
    cut = None
    offset = 0
    for raw_line in body.splitlines(keepends=True):
        if raw_line.strip() == "---":
            cut = offset
            break
        offset += len(raw_line)
    assert cut is not None, "no section separator after the conclusion block"
    block = body[:cut]
    stray = [ln for ln in block.splitlines() if ln.strip() and not ln.startswith(">")]
    assert not stray, f"unquoted lines inside the conclusion block: {stray[:3]}"

    unquoted = "\n".join(
        ln[2:] if ln.startswith("> ") else ("" if ln.strip() == ">" else ln)
        for ln in block.splitlines()
    ).strip()
    assert unquoted == final, (
        "the quoted conclusion differs from the raw log - re-run scripts/capture_demo.py "
        "and re-quote it whole rather than editing either one"
    )


def _assert_same(recorded, live, path: str) -> None:
    """Deep equality with a float tolerance, so no value can drift unchecked."""
    assert type(recorded) is type(live) or (
        isinstance(recorded, (int, float)) and isinstance(live, (int, float))
    ), f"{path}: type changed ({type(recorded).__name__} -> {type(live).__name__})"
    if isinstance(recorded, dict):
        assert set(recorded) == set(live), (
            f"{path}: keys differ (log-only {sorted(set(recorded) - set(live))}, "
            f"live-only {sorted(set(live) - set(recorded))})"
        )
        for k in recorded:
            _assert_same(recorded[k], live[k], f"{path}.{k}")
    elif isinstance(recorded, list):
        assert len(recorded) == len(live), f"{path}: length {len(recorded)} -> {len(live)}"
        for i, (r, v) in enumerate(zip(recorded, live, strict=True)):
            _assert_same(r, v, f"{path}[{i}]")
    elif isinstance(recorded, bool) or recorded is None:
        assert recorded == live, f"{path}: {recorded!r} -> {live!r}"
    elif isinstance(recorded, (int, float)):
        assert recorded == pytest.approx(live, rel=1e-6), f"{path}: {recorded} -> {live}"
    else:
        assert recorded == live, f"{path}: {recorded!r} -> {live!r}"


def test_committed_demo_results_replay_against_live_code():
    """Replay every recorded tool result: the log must not drift from the server.

    The log went stale twice - once carrying pre-fix numbers, once carrying field
    names the code no longer emits - while TRANSCRIPT.md claimed a test would catch
    exactly that. No test did: the landmark guard re-derives PROVENANCE, not the log.
    This one calls the tools with the recorded arguments and compares the responses.
    """
    import json

    from hic_mcp import analysis, models

    records = [
        json.loads(line)
        for line in (REPO / "demo" / "raw-session.jsonl").read_text().splitlines()
        if line.strip()
    ]
    calls = {}
    for rec in records:
        if rec.get("type") == "assistant":
            for c in rec["message"]["content"]:
                if c.get("type") == "tool_use" and str(c.get("name", "")).startswith(
                    "mcp__hic-mcp__"
                ):
                    calls[c["id"]] = (c["name"].rsplit("__", 1)[-1], c["input"])
    results = {}
    for rec in records:
        msg = rec.get("message", {})
        if rec.get("type") == "user" and isinstance(msg.get("content"), list):
            for c in msg["content"]:
                if c.get("type") == "tool_result" and c.get("tool_use_id") in calls:
                    body = c.get("content")
                    text = (
                        " ".join(x.get("text", "") for x in body if isinstance(x, dict))
                        if isinstance(body, list)
                        else str(body)
                    )
                    results[c["tool_use_id"]] = text

    assert results, "no recorded tool results found in the demo log"
    for call_id, recorded_text in results.items():
        name, kwargs = calls[call_id]
        try:
            recorded = json.loads(recorded_text)
        except json.JSONDecodeError:
            continue  # not a JSON payload; nothing to compare
        # wrap through the same response model the server returns, so optional
        # fields serialise identically - otherwise this compares a raw dict against
        # an MCP payload and reports serialisation as drift
        model = getattr(models, "".join(w.title() for w in name.split("_")), None)
        live_raw = getattr(analysis, name)(**kwargs)
        live = json.loads(model(**live_raw).model_dump_json()) if model else live_raw
        assert set(recorded) == set(live), (
            f"{name}: the committed log and the live server return different fields "
            f"(log-only {sorted(set(recorded) - set(live))}, "
            f"live-only {sorted(set(live) - set(recorded))}). Re-run "
            "scripts/capture_demo.py rather than editing the log."
        )
        # EVERY value, to any depth. This once compared only str/int/bool/None, which
        # skipped region_mean_E1, top_boundaries, boundary_counts_per_window, eigenvalues
        # and E1_track - i.e. every scientific number in the log. An evaluator mutated
        # those and the whole honesty module still passed.
        _assert_same(recorded, live, name)
