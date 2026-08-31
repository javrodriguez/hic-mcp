"""The MCP layer: server construction, tool registration, stdio entrypoint."""

from mcp.server.mcpserver import MCPServer

server = MCPServer(
    name="hic-mcp",
    version="0.1.0",
    instructions=(
        "Hi-C / 3D-chromatin analysis over local .mcool/.cool contact-matrix files. "
        "Every tool runs a real open2c (cooler/cooltools) computation; nothing is mocked. "
        "A small real demo dataset (human chr17, HFF Micro-C, hg38) is bundled - tools "
        "default to it when no file path is given."
    ),
)


def main() -> None:
    """Console-script entry point (stdio transport)."""
    server.run()
