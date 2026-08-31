"""MCP integration: every tool round-trips through a real client.

The in-memory tests drive the actual server object over MCP; the stdio test
launches the packaged console script as a real subprocess, exactly as an MCP
client (Claude Desktop, Claude Code) would.
"""

import sys
from pathlib import Path

from mcp import StdioServerParameters
from mcp.client.client import Client

from hic_mcp.server import server

REPO_ROOT = Path(__file__).resolve().parents[1]
TOOLS = {
    "matrix_summary",
    "contacts_at_locus",
    "insulation_tads",
    "compartments",
    "virtual_4c",
    "expected_observed",
}


async def test_list_tools_names_and_schemas():
    async with Client(server, raise_exceptions=True) as client:
        tools = (await client.list_tools()).tools
        assert {t.name for t in tools} == TOOLS
        for t in tools:
            assert t.description, f"{t.name} has no description"
            assert t.output_schema is not None, f"{t.name} has no output schema"


async def test_matrix_summary_roundtrip():
    async with Client(server, raise_exceptions=True) as client:
        r = await client.call_tool("matrix_summary", {})
        assert not r.is_error
        sc = r.structured_content
        assert sc["is_bundled_demo"] is True
        assert sc["chromosomes"] == {"chr17": 83_257_441}
        assert {x["resolution_bp"] for x in sc["resolutions"]} == {10_000, 100_000, 1_000_000}


async def test_each_analysis_tool_roundtrip():
    calls = {
        "contacts_at_locus": {"region": "chr17:50,000,000-50,500,000", "resolution": 100_000},
        "insulation_tads": {"region": "chr17:65,000,000-67,000,000"},
        "compartments": {"region": "chr17:50,100,000-51,100,000"},
        "virtual_4c": {"viewpoint": "chr17:63,000,000-63,100,000", "resolution": 100_000},
        "expected_observed": {"region": "chr17:50,000,000-52,500,000"},
    }
    async with Client(server, raise_exceptions=True) as client:
        for name, args in calls.items():
            r = await client.call_tool(name, args)
            assert not r.is_error, f"{name} errored: {r.content}"
            sc = r.structured_content
            assert sc["method"], name
            assert sc["resolution_used"] > 0, name
    # spot ground truths through the MCP layer
    async with Client(server, raise_exceptions=True) as client:
        comp = (await client.call_tool("compartments", calls["compartments"])).structured_content
        assert comp["region_call"] == "A"
        ins_r = await client.call_tool("insulation_tads", calls["insulation_tads"])
        boundaries = ins_r.structured_content["top_boundaries"]
        assert any(b["locus"] == "chr17:66,180,000-66,190,000" for b in boundaries)


async def test_anticipated_errors_reach_the_agent():
    async with Client(server) as client:
        r = await client.call_tool(
            "virtual_4c", {"viewpoint": "chr17:45,500,000-45,510,000", "resolution": 10_000}
        )
        assert r.is_error
        assert "ICE-filtered" in r.content[0].text
        r2 = await client.call_tool("contacts_at_locus", {"region": "chr9:1-1000"})
        assert r2.is_error
        assert "Chromosomes here" in r2.content[0].text


async def test_no_bare_error_ever_reaches_the_agent():
    """Every reachable failure names a reason the model can act on."""
    calls = [
        ("expected_observed", {"region": "chr17:20,000,000-25,000,000"}),  # spans the arms
        ("contacts_at_locus", {"region": "chr17:50,000,000-50,000,000"}),  # zero width
        ("virtual_4c", {"viewpoint": "chr17:24,000,000-24,010,000"}),  # filtered
        ("insulation_tads", {"windows_bp": [1000]}),  # window under 3 bins
    ]
    async with Client(server) as client:
        for name, args in calls:
            r = await client.call_tool(name, args)
            assert r.is_error, f"{name} should have failed"
            text = r.content[0].text
            assert text.strip() != f"Error executing tool {name}", f"{name} gave a bare error"
            assert len(text) > 60, f"{name} error too terse to act on: {text}"


async def test_unexpected_failures_are_wrapped_not_bare(monkeypatch):
    """The backstop turns a genuine bug into something the agent can still read."""
    from hic_mcp import analysis

    def boom(**kwargs):
        raise RuntimeError("simulated internal defect")

    monkeypatch.setattr(analysis, "matrix_summary", boom)
    async with Client(server) as client:
        r = await client.call_tool("matrix_summary", {})
        assert r.is_error
        assert "simulated internal defect" in r.content[0].text
        assert "bug in hic-mcp" in r.content[0].text


async def test_real_stdio_subprocess_roundtrip():
    """The packaged entry point answers over stdio, as a real MCP client launches it."""
    params = StdioServerParameters(
        command=sys.executable,
        args=["-m", "hic_mcp"],
        cwd=str(REPO_ROOT),
    )
    async with Client(params, raise_exceptions=True) as client:
        tools = (await client.list_tools()).tools
        assert {t.name for t in tools} == TOOLS
        r = await client.call_tool("matrix_summary", {})
        assert r.structured_content["is_bundled_demo"] is True
