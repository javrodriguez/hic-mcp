"""The MCP layer: server construction, the six tools, stdio entrypoint.

Thin by design - each tool validates through the shared helpers, delegates to
the pure functions in analysis.py, and wraps the result in a typed model.
Anticipated failures surface as ToolError so the calling agent can read the
reason and correct course; nothing here prints to stdout (stdio transport).
"""

from typing import Annotated

from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError
from pydantic import Field

from hic_mcp import analysis, models
from hic_mcp.analysis import AnalysisError
from hic_mcp.data import DataError

server = MCPServer(
    name="hic-mcp",
    version="0.1.0",
    instructions=(
        "Hi-C / 3D-chromatin analysis over local .mcool/.cool contact-matrix files. "
        "Every tool runs a real open2c (cooler/cooltools) computation; nothing is mocked. "
        "A small real demo dataset (human chr17, HFF Micro-C, hg38) is bundled - tools "
        "default to it when no file path is given. Start with matrix_summary to see what "
        "a file contains. Regions are UCSC-style, e.g. 'chr17:50,000,000-52,000,000'."
    ),
)

_FILE = Annotated[
    str | None,
    Field(description="Path to a local .cool/.mcool file; omit to use the bundled demo data"),
]
_RESOLUTION = Annotated[
    int | None,
    Field(description="Bin size in bp; omit to auto-select from the file's resolutions"),
]


def _run(fn, /, **kwargs):
    try:
        return fn(**kwargs)
    except (DataError, AnalysisError) as e:
        raise ToolError(str(e)) from e
    except Exception as e:  # backstop: an agent must never receive a reason-less failure
        args = ", ".join(f"{k}={v!r}" for k, v in kwargs.items() if v is not None)
        raise ToolError(
            f"{fn.__name__} failed unexpectedly on ({args}): {type(e).__name__}: {e}. "
            "This is a bug in hic-mcp, not in your request - try a different region or "
            "resolution, and please report it."
        ) from e


@server.tool()
def matrix_summary(file: _FILE = None) -> models.MatrixSummary:
    """Summarize a Hi-C file: chromosomes, resolutions, contact totals, balancing status."""
    return models.MatrixSummary(**_run(analysis.matrix_summary, file=file))


@server.tool()
def contacts_at_locus(
    region: Annotated[
        str, Field(description="UCSC-style region, e.g. 'chr17:50,000,000-52,000,000'")
    ],
    region2: Annotated[
        str | None, Field(description="Optional second region for an off-diagonal block")
    ] = None,
    file: _FILE = None,
    resolution: _RESOLUTION = None,
    balanced: Annotated[
        bool, Field(description="Report ICE-balanced values when available")
    ] = True,
) -> models.ContactsAtLocus:
    """Contact-matrix statistics at a locus (plus the matrix itself for small windows)."""
    return models.ContactsAtLocus(
        **_run(
            analysis.contacts_at_locus,
            file=file,
            region=region,
            region2=region2,
            resolution=resolution,
            balanced=balanced,
        )
    )


@server.tool()
def insulation_tads(
    region: Annotated[
        str | None, Field(description="Restrict reported boundaries to this region; omit for all")
    ] = None,
    file: _FILE = None,
    resolution: _RESOLUTION = None,
    windows_bp: Annotated[
        list[int] | None,
        Field(description="Diamond window sizes in bp; omit to scale with the bin size"),
    ] = None,
    top_n: Annotated[int, Field(description="How many strongest boundaries to report")] = 10,
) -> models.InsulationTads:
    """Insulation score and TAD-boundary calls (diamond insulation, Crane et al. 2015)."""
    return models.InsulationTads(
        **_run(
            analysis.insulation_tads,
            file=file,
            region=region,
            resolution=resolution,
            windows_bp=windows_bp,
            top_n=top_n,
        )
    )


@server.tool()
def compartments(
    region: Annotated[
        str | None,
        Field(description="Report the A/B call and E1 for this region; omit for a genome summary"),
    ] = None,
    file: _FILE = None,
    resolution: _RESOLUTION = None,
) -> models.Compartments:
    """A/B compartment eigenvector from the cis observed/expected matrix (eigs_cis)."""
    return models.Compartments(
        **_run(analysis.compartments, file=file, region=region, resolution=resolution)
    )


@server.tool()
def virtual_4c(
    viewpoint: Annotated[
        str, Field(description="Viewpoint region, e.g. 'chr17:63,000,000-63,100,000'")
    ],
    file: _FILE = None,
    resolution: _RESOLUTION = None,
) -> models.Virtual4C:
    """Virtual-4C profile: balanced contact frequency of one viewpoint with its chromosome."""
    return models.Virtual4C(
        **_run(analysis.virtual_4c, file=file, viewpoint=viewpoint, resolution=resolution)
    )


@server.tool()
def expected_observed(
    region: Annotated[str, Field(description="UCSC-style region for the O/E matrix")],
    file: _FILE = None,
    resolution: _RESOLUTION = None,
) -> models.ExpectedObserved:
    """Distance-expected contact curve, P(s) slope, and the observed/expected matrix."""
    return models.ExpectedObserved(
        **_run(analysis.expected_observed, file=file, region=region, resolution=resolution)
    )


def main() -> None:
    """Console-script entry point (stdio transport)."""
    server.run()
