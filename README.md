# HiC-MCP

Hi-C / 3D-chromatin analysis for AI agents: an [MCP](https://modelcontextprotocol.io) server exposing the [open2c](https://open2c.github.io/) stack ([cooler](https://github.com/open2c/cooler), [cooltools](https://github.com/open2c/cooltools)) as tools over local `.mcool`/`.cool` contact-matrix files.
Every tool runs the real computation on real data — nothing is mocked.

> **Status: under construction.** The analysis tools land in the next commits; this section is removed when they do.

## Quickstart

```bash
git clone https://github.com/javrodriguez/hic-mcp
cd hic-mcp
uv sync
```

Paste into `claude_desktop_config.json` (Claude Desktop):

```json
{
  "mcpServers": {
    "hic-mcp": {
      "command": "uv",
      "args": ["--directory", "/ABS/PATH/TO/hic-mcp", "run", "hic-mcp"]
    }
  }
}
```

Or with Claude Code:

```bash
claude mcp add hic-mcp -- uv --directory /ABS/PATH/TO/hic-mcp run hic-mcp
```

Claude Desktop launches without a login-shell `PATH` on some systems — if the server does not appear, use the absolute path to `uv` (`which uv`) as `"command"`.

## What this is

A demonstration system for agent-driven 3D-genome analysis, built on the official MCP Python SDK.
It is not a substitute for a full interactive analysis environment.

## License

MIT (source code only — bundled data carries its own terms; see `LICENSE`).
