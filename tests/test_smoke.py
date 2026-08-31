"""Smoke: the server constructs and a real in-memory client handshake round-trips."""

from mcp.client.client import Client

from hic_mcp.server import server


async def test_server_handshake_and_list_tools() -> None:
    async with Client(server, raise_exceptions=True) as client:
        result = await client.list_tools()
        assert isinstance(result.tools, list)
