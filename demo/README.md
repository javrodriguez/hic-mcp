# The demo session

[`TRANSCRIPT.md`](TRANSCRIPT.md) — the readable version.
[`raw-session.jsonl`](raw-session.jsonl) — the machine-captured session it quotes from.

## Reproduce it

```bash
python scripts/capture_demo.py
```

This launches a real MCP client (Claude Code, headless) against a real `hic-mcp` server running from your checkout, and rewrites `raw-session.jsonl` from that session.
It needs `claude` on your PATH and a logged-in Claude Code; nothing else.
A different model or a different day will phrase the answer differently — the tool calls and the numbers they return are what should reproduce, because those come from the data.

## What the raw log contains, exactly

- **Conversation records, verbatim** (`assistant` and `user` types): the model's messages, its tool calls, and the tool results this server returned. These are recorded exactly as they happened. They are never edited — a transcript improved after the fact would be a fabrication, and the whole point of shipping the raw log is that you can check the pretty one against it.
- **Envelope records, filtered** (`system` and `result` types): reduced to an allowlist of fields (declared as `ENVELOPE_KEEP` in `scripts/capture_demo.py`). The harness puts host-machine details in these records — the operator's home directory, session identifiers, a local socket path — and those have no place in a public repository. The capture script refuses to write a file that still matches any host-path pattern, so this is enforced rather than promised.

The script also refuses to write a session that made no `hic-mcp` tool calls, so this file cannot quietly become a transcript of a model answering from memory.
