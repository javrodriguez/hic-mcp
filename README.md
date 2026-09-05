# HiC-MCP

[![CI](https://github.com/javrodriguez/hic-mcp/actions/workflows/ci.yml/badge.svg)](https://github.com/javrodriguez/hic-mcp/actions/workflows/ci.yml)

**Hi-C / 3D-chromatin analysis for AI agents.** An [MCP](https://modelcontextprotocol.io) server that exposes the [open2c](https://open2c.github.io/) stack ([cooler](https://github.com/open2c/cooler), [cooltools](https://github.com/open2c/cooltools)) as tools over local `.mcool`/`.cool` contact matrices — TAD boundaries, A/B compartments, virtual 4C, observed/expected, and more.

**A real 6.8 MB Hi-C dataset ships with it**, so the quickstart works offline with no API keys, no accounts, and no data of your own. Every tool runs the real computation; nothing is mocked.

## Quickstart

```bash
git clone https://github.com/javrodriguez/hic-mcp
cd hic-mcp
uv sync
```

Then add it to your MCP client. **Claude Desktop** — paste into `claude_desktop_config.json`:

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

**Claude Code** — one line:

```bash
claude mcp add hic-mcp -- uv --directory /ABS/PATH/TO/hic-mcp run hic-mcp
```

Now ask your agent: *"Find the strongest TAD boundary in the demo Hi-C data, then tell me which compartment it sits in."*

## What that actually returns

Real output from the bundled data — the entire response to `insulation_tads(region="chr17:65,000,000-67,000,000", top_n=2)`, nothing removed:

```json
{
  "region": "chr17:65,000,000-67,000,000",
  "resolution_used": 10000,
  "view": "chr17 p/q arms (bundled)",
  "windows_bp": [
    100000,
    250000,
    500000
  ],
  "ranked_by": "boundary_strength at the 250000 bp window - top_boundaries lists ONLY boundaries called at that window (capped at top_n), so it is a subset of boundary_counts_per_window, whose other entries count different populations",
  "boundary_counts_per_window": {
    "100000": 5,
    "250000": 4,
    "500000": 3
  },
  "top_boundaries": [
    {
      "locus": "chr17:66,180,000-66,190,000",
      "strength": 2.76292,
      "log2_insulation": -1.63746,
      "windows_detected": [
        100000,
        250000,
        500000
      ]
    },
    {
      "locus": "chr17:66,680,000-66,690,000",
      "strength": 1.72627,
      "log2_insulation": -1.1018,
      "windows_detected": [
        100000,
        250000
      ]
    }
  ],
  "balanced": true,
  "scale_note": null,
  "method": "cooltools.insulation (diamond insulation score; Li threshold boundary calls). Returns the called boundaries and their scores, not the full per-bin insulation track."
}
```

Every response names the method it used, the resolution it ran at, and whether values are ICE-balanced — so an agent can quote the result without overstating it.

## The tools

| Tool | What it computes | Method |
|---|---|---|
| `matrix_summary` | Chromosomes, resolutions, contact totals, balancing status, provenance | `cooler.Cooler.info` across every resolution |
| `contacts_at_locus` | Contact statistics at a locus or between two loci; the balanced matrix for small windows | `cooler.Cooler.matrix().fetch` |
| `insulation_tads` | Insulation score and TAD-boundary calls, across several diamond windows | `cooltools.insulation` (Crane et al. 2015; Li-threshold calls) |
| `compartments` | A/B compartment eigenvector, GC-phased so the sign is meaningful | `cooltools.eigs_cis` on the cis observed/expected map |
| `virtual_4c` | Contact profile of one viewpoint against its chromosome | `cooltools.virtual4c` |
| `expected_observed` | Distance-expected contact curve, P(s) slope, and the O/E matrix | `cooltools.expected_cis` |

Point any tool at your own file with `file="/path/to/yours.mcool"`; omit it to use the bundled demo.

**Bringing your own file.** Two optional arguments make the analyses as meaningful on your data as they are on the demo:

- `view` — a tab-separated BED-like file (`chrom`/`start`/`end`[/`name`]) partitioning the genome into regions analysed independently. Hi-C statistics should not be normalised across a centromere, so pass chromosome arms here; without it the tools fall back to whole chromosomes and say so in the `view` field they return.
- `phasing_track` — a `chrom`/`start`/`end`/`value` file (GC fraction, gene density) that orients the compartment eigenvector. Without it the sign of E1 is arbitrary, so `compartments` reports `unphased` and refuses to label A or B rather than guessing.

The bundled demo ships both, which is why it gets phased A/B calls and arm-partitioned scores out of the box.

**Two things worth knowing before you point it at a large file.** The `region` argument filters what is *reported*, not what is *computed* — insulation, compartments and expected curves are calculated across the whole chromosome (or arm) either way, so a genome-wide 10 kb `.mcool` will take minutes per call, not seconds. And the bundled demo dataset lives in this repository, not in the built package: install from a clone to get it, or pass your own file.

## Example prompts

- *"What's in the demo Hi-C file — which chromosomes and resolutions?"*
- *"Call TAD boundaries on chr17 between 65 and 67 Mb and rank them by strength."*
- *"Is chr17:50.1-51.1 Mb in the A or B compartment, and how confident is that call?"*
- *"Plot me the contact-decay curve — what's the P(s) slope?"*
- *"Run a virtual 4C from chr17:63 Mb and describe how contacts fall off with distance."*

More in [`examples/PROMPTS.md`](examples/PROMPTS.md). The same client config is in [`examples/claude_desktop_config.json`](examples/claude_desktop_config.json) (Claude Desktop) and [`examples/mcp.json`](examples/mcp.json) (project-scoped clients).

## A real session

[`demo/TRANSCRIPT.md`](demo/TRANSCRIPT.md) walks through an actual headless agent session answering exactly that question — chaining `insulation_tads` into `compartments` to find the strongest boundary and place it in a compartment. The machine-captured log it quotes from is committed beside it ([`demo/raw-session.jsonl`](demo/raw-session.jsonl)), and `python scripts/capture_demo.py` reproduces it against your own checkout.

## What this is, and what it isn't

This is a **demonstration system**: a small, honest, end-to-end example of exposing a real scientific analysis stack to agents over MCP. It runs genuine open2c computations on real published data, and its outputs are the library's own — but it is not a replacement for an interactive analysis environment, and it deliberately ships one chromosome rather than a genome.

Where a result can't be trusted, the tools say so rather than returning a number: a viewpoint in an ICE-filtered region raises a clear error instead of a null profile, the first two diagonals of an O/E matrix come back null because `cooltools` does not measure expected there, an unphased eigenvector reports a sign-neutral fraction instead of claiming an "A compartment" share, and a bin size too coarse for TADs says so in the response.

## Data

Derived subset (chr17, hg38) of HFFc6 Micro-C generated by the Dekker and Rando labs (UMass Chan) for the NIH 4D Nucleome Network (1U54DK107980-01), accession [`4DNESWST3UBH`](https://data.4dnucleome.org/experiment-set-replicates/4DNESWST3UBH/), obtained via the Open2C [cooltools](https://github.com/open2c/cooltools) test-data registry ([osf.io/3h9js](https://osf.io/3h9js/)).
Cite Krietenstein et al., *Mol Cell* 78:554-565 (2020), doi:[10.1016/j.molcel.2020.03.003](https://doi.org/10.1016/j.molcel.2020.03.003); 4DN White Paper doi:[10.1038/nature23884](https://doi.org/10.1038/nature23884); 4DN Portal doi:[10.1038/s41467-022-29697-4](https://doi.org/10.1038/s41467-022-29697-4).
Redistributed under the [4DN Data Release and Use Policy](https://github.com/4dn-dcic/4dn-policies/blob/master/4dn-data-release-and-use-policy.md).
Full provenance, the rebuild script, and the measured landmarks the test suite asserts against: [`data/PROVENANCE.md`](data/PROVENANCE.md).

## Related work

Surveyed 31 August 2026. Hi-C analysis over MCP is thinly covered, and neighbours differ in scope rather than quality:

- [`zhouhufeng/IGVFagent`](https://github.com/zhouhufeng/IGVFagent) — a 161-tool IGVF/ENCODE agent whose MCP surface includes three Hi-C tools (contact heatmap, Crane insulation, loop QC) built on cooler and hic-straw; it does not expose the cooltools analysis surface (compartments, expected curves, pileups).
- [`aidenlab/juicebox-mcp`](https://github.com/aidenlab/juicebox-mcp) — Hi-C contact-matrix visualization over MCP.
- [`ammawla/encode-toolkit`](https://github.com/ammawla/encode-toolkit) — discovery and QC of ENCODE Hi-C experiments; no contact-matrix operations.
- [`GPTomics/bioSkills`](https://github.com/GPTomics/bioSkills) — a broad Hi-C skill set (including cooltools-based compartment and TAD analysis) as Claude Agent Skills rather than an MCP server.

## Development

```bash
uv run pytest -q      # includes ground-truth assertions against the bundled data
uv run ruff check .
```

The test suite asserts known measured results — a specific boundary locus, compartment signs on named blocks, the P(s) slope, exact contact totals — so a regression that changes the science fails the build rather than passing quietly. Tests run with sockets blocked, proving the demo needs no network.

## License

MIT — **source code only**. The bundled dataset carries its own terms; see [`LICENSE`](LICENSE) and [`data/PROVENANCE.md`](data/PROVENANCE.md).
