#!/usr/bin/env python
"""Capture a real agent session against this server, and record it verbatim.

This is how demo/raw-session.jsonl is produced. It launches a real MCP client
(Claude Code, headless) against a real hic-mcp server running from this checkout,
and records the session as it happened.

What is recorded verbatim: every conversation record - the model's messages, its
tool calls, and the tool results this server returned. None of it is edited; a
transcript with an improved answer in it would be a fabrication, not a demo.

What is filtered: the harness's own envelope records (`system` and `result`) are
reduced to an allowlist of fields, because they carry host-machine details - the
operator's home directory, session ids, a local socket path - that have no place
in a public repository. The allowlist is declared below, the script refuses to
write a file that still contains a host path, and re-running this script from any
machine reproduces the same shape.

Usage:
    python scripts/capture_demo.py [--out demo/raw-session.jsonl] [--model sonnet]
"""

import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]

PROMPT = (
    "Using the hic-mcp tools, find the strongest TAD boundary in the bundled demo Hi-C "
    "data, then check which compartment that boundary sits in. Report the locus, its "
    "boundary strength, which diamond windows detect it, and the compartment call with "
    "its eigenvector value."
)

TOOLS = (
    "matrix_summary", "contacts_at_locus", "insulation_tads",
    "compartments", "virtual_4c", "expected_observed",
)

# Envelope records keep only these fields; everything else in them is host-specific.
# `tools` is deliberately NOT kept whole: the init record lists the operator's entire local
# agent-tool inventory, most of which has nothing to do with this project and some of which
# is not part of a stock install. That is machine-specific configuration of exactly the kind
# this allowlist exists to strip, and it is why re-running this script elsewhere would not
# reproduce "the same shape" the demo README promises. Only this server's own tools evidence
# anything about this server, so only those are kept - see keep_hic_tools_only().
ENVELOPE_KEEP = {
    "type", "subtype", "mcp_servers", "model", "is_error",
    "num_turns", "result",
}
MCP_TOOL_PREFIX = "mcp__hic-mcp__"
# A recorded session must not carry any of these into the repository.
LEAK_PATTERNS = [r"/Users/", r"/home/[a-z]", r"C:\\\\Users", r"cc-socks", r"\.claude/projects"]


def build_config(tmpdir: str) -> str:
    cfg = {
        "mcpServers": {
            "hic-mcp": {"command": "uv", "args": ["--directory", str(REPO), "run", "hic-mcp"]}
        }
    }
    path = os.path.join(tmpdir, "mcp.json")
    Path(path).write_text(json.dumps(cfg, indent=2))
    return path


def run_session(config_path: str, model: str, workdir: str) -> list[dict]:
    cmd = [
        "claude", "-p", PROMPT,
        "--model", model,
        "--output-format", "stream-json",
        "--verbose",
        "--strict-mcp-config",
        "--mcp-config", config_path,
        "--allowedTools", ",".join(f"mcp__hic-mcp__{t}" for t in TOOLS),
    ]
    env = {k: v for k, v in os.environ.items() if k not in ("CLAUDE_PROJECT_DIR",)}
    with open(os.devnull) as devnull:
        proc = subprocess.run(
            cmd, cwd=workdir, env=env, stdin=devnull, capture_output=True, text=True
        )
    if proc.returncode != 0:
        sys.exit(f"capture failed (exit {proc.returncode}): {proc.stderr[-2000:]}")
    return [json.loads(line) for line in proc.stdout.splitlines() if line.strip()]


def clean(records: list[dict]) -> list[dict]:
    out = []
    for rec in records:
        if rec.get("type") in ("assistant", "user"):
            out.append(rec)  # verbatim: the conversation itself
        else:
            kept = {k: v for k, v in rec.items() if k in ENVELOPE_KEEP}
            tools = rec.get("tools")
            if isinstance(tools, list):
                kept["tools"] = [
                    t for t in tools if isinstance(t, str) and t.startswith(MCP_TOOL_PREFIX)
                ]
            out.append(kept)
    return out


def assert_no_leaks(records: list[dict]) -> None:
    blob = json.dumps(records)
    for pattern in LEAK_PATTERNS:
        hit = re.search(pattern, blob)
        if hit:
            sys.exit(
                f"refusing to write: the session still contains a host-machine detail "
                f"matching {pattern!r} ({blob[max(0, hit.start() - 60):hit.end() + 60]!r}). "
                "Fix the capture environment rather than editing the transcript."
            )


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", default=str(REPO / "demo" / "raw-session.jsonl"))
    ap.add_argument("--model", default="sonnet")
    args = ap.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        config_path = build_config(tmpdir)
        records = run_session(config_path, args.model, workdir=tmpdir)

    records = clean(records)
    assert_no_leaks(records)

    tool_calls = sum(
        1
        for r in records
        if r.get("type") == "assistant"
        for c in r["message"]["content"]
        if c.get("type") == "tool_use" and str(c.get("name", "")).startswith("mcp__hic-mcp__")
    )
    if not tool_calls:
        sys.exit("refusing to write: the session made no hic-mcp tool calls")

    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w") as fh:
        for rec in records:
            fh.write(json.dumps(rec) + "\n")
    print(f"wrote {out} - {len(records)} records, {tool_calls} hic-mcp tool calls")


if __name__ == "__main__":
    main()
