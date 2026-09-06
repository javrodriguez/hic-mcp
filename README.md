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

Both roads assume [`uv`](https://docs.astral.sh/uv/) is installed and on your PATH. GUI clients
launch without a login shell, so a `uv` in `~/.local/bin` is often invisible to them — if your
client reports it cannot spawn `uv`, replace `"command": "uv"` with the absolute path that
`which uv` prints.

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

The bundled demo ships both. The arm view applies at any resolution; the GC track is binned at 100 kb, so out-of-the-box phased A/B calls come at `resolution=100000` (which is `compartments`' default). At another resolution the response says so and names the fix, rather than telling you to go and find a track you already have.

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

Surveyed 5 September 2026. **The cooltools surface is already reachable over MCP** — this server
is not the first thing to do that, and the honest distinction is *how*, not *whether*:

- [`coala-info/coala`](https://github.com/coala-info/coala) — a general CWL→MCP adapter. Its
  [Hi-C example](https://coala.info/use-cases/Hi-C.html) (Feb 2026) stands up an MCP server over
  `cooler_dump`, `cooltools_insulation`, `cooltools_eigs_cis`, `cooltools_expected_cis` and
  `cooltools_saddle`, and its companion
  [`coala-repo`](https://github.com/coala-info/coala-repo) carries CWL definitions for the whole
  cooltools CLI, `virtual4c` included. **Anything hic-mcp computes, coala can already run.** The
  difference is that coala wraps the *command line* generically, so a tool's output is whatever
  the CLI emits; hic-mcp wraps the *Python API* in six typed responses that each name their
  method, their scope, and the cases where the number should not be trusted.
- [`zhouhufeng/IGVFagent`](https://github.com/zhouhufeng/IGVFagent) — a large IGVF/ENCODE agent
  whose MCP surface now includes ~12 Hi-C tools (contact matrix, Crane insulation, loops, and a
  `spatial_hic_*` set added Sep 2026 covering A/B compartments). Built on cooler and hic-straw
  with its own reimplementations rather than cooltools; no virtual 4C.
- [`aidenlab/juicebox-mcp`](https://github.com/aidenlab/juicebox-mcp) — Hi-C contact-matrix
  visualization over MCP (Dec 2025), and the earliest Hi-C MCP server I found.
- [`BIsnake2001/ChromSkills`](https://github.com/BIsnake2001/ChromSkills) — the broadest Hi-C
  analysis surface of any of these (compartment shifts, nested and differential TADs, loop
  annotation), as Claude Skills in Docker rather than MCP.
- [`ammawla/encode-toolkit`](https://github.com/ammawla/encode-toolkit) — discovery and QC of
  ENCODE Hi-C experiments; no contact-matrix operations over MCP.
- [`GPTomics/bioSkills`](https://github.com/GPTomics/bioSkills) — a broad Hi-C skill set as
  Claude Agent Skills; archived and read-only since Aug 2026.

**So this repository claims no priority of any kind.** As of 5 September 2026, "the first Hi-C
MCP server" would be false — juicebox-mcp precedes it by nine months — and "the only one exposing
cooltools" would be false too, because coala already does. What it offers instead is a purpose-built, typed surface
over the cooltools Python API, with a real dataset bundled so the whole thing runs offline, and
a response contract that says what each number means and when it should not be believed.

## Development

```bash
uv run pytest -q      # 1-2 min: every assertion runs the real computation on the bundled data
uv run ruff check .
```

`-q` prints nothing until it finishes, so the wait is expected rather than a hang — and if you
only want a first look, `uv run pytest -q tests/test_data.py tests/test_server.py` covers the
plumbing and the MCP round-trips in well under a minute. There is no faster honest version of the
full run: the ground-truth assertions decompose the real matrix, which is the point of them.

The test suite asserts known measured results — a specific boundary locus, compartment signs on named blocks, the P(s) slope, exact contact totals — so a regression that changes the science fails the build rather than passing quietly. The in-process tests run with sockets blocked, proving the demo needs no network. (The one test that launches the packaged server as a child process is outside that in-process patch — it is the round-trip a real client makes, and it reaches nothing but your own file.)

## How this was built, and what the commit messages mean

Many commit bodies in this repository refer to "rounds", "evaluators" and "findings". They are not
noise, and this is what they mean.

Every change here was graded by **independent fresh-context evaluators**, each reading a clean
clone blind, each given a byte-identical prompt, each writing its report before anything was
fixed. Thirteen rounds ran: eight against the first design with **three evaluators per round**,
then five more against a frozen surface with **two per round** — dropped to two on the evidence
that every material finding in the last three rounds of the first design had been raised by at
least two of the three. Across them, **120 unique findings** were drawn
from 159 raw evaluator reports — the gap is the same defect found independently by two or three
evaluators. Every one was reproduced and fixed at root. No finding was overridden or downgraded.

**The stopping rule was never met, and that is worth saying plainly.** The first design required
two consecutive rounds where every evaluator returned clean; the second relaxed the *shape* but
not the bar — two consecutive rounds with no MATERIAL finding on a frozen surface, minors still
fixed. Neither was ever reached. The loop was stopped by a judgement call, not by convergence, so
this repository is *not* certified by its own standard.

That does not mean defects are outstanding — there is no known unfixed defect. It means each round
kept finding a new instance of one class: **a road the test corpus never walks**. Region order,
then region overlap, then bin-grid alignment. What decided the stop was round 5's largest finding
being a defect inside round 4's own fix.

What the loop actually caught, since that is the fairer measure:

- an observed/expected matrix reading **936×** where it should read ~1
- a P(s) slope computed file-wide but labelled with whatever region you asked for
- a confidence figure computed over a single bin, which made the demo read as more certain than the
  data supported

Where it stands now: **191 tests**, ruff and mypy clean, CI green on Python 3.12 and 3.13, and a
clean-clone walk that installs and passes in one to two minutes.

## License

MIT — **source code only**.

The bundled demo dataset under `data/` is not covered by the MIT licence. It is a derived subset of
published, publicly released 4D Nucleome data, redistributed under the
[4DN Data Release and Use Policy](https://github.com/4dn-dcic/4dn-policies/blob/master/4dn-data-release-and-use-policy.md).
Full provenance, citations and terms are in [`data/PROVENANCE.md`](data/PROVENANCE.md).

_(This scope note lives here rather than inside `LICENSE`, so that automated licence detection reads the
source licence as plain MIT. The scope itself is unchanged.)_
