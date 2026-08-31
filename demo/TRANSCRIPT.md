# A real agent session

This is a curated reading of [`raw-session.jsonl`](raw-session.jsonl) — an actual headless Claude Code session driving this server over MCP.
Every quotation below is copied from that file; the raw log is the evidence, this page is only the readable version of it.

**Reproduce it:** `python scripts/capture_demo.py` (see [`README.md`](README.md) in this folder for what that does).

---

**Prompt given to the agent:**

> Using the hic-mcp tools, find the strongest TAD boundary in the bundled demo Hi-C data, then check which compartment that boundary sits in. Report the locus, its boundary strength, which diamond windows detect it, and the compartment call with its eigenvector value.

**The server connects:**

```json
{"type": "system", "subtype": "init", "mcp_servers": [{"name": "hic-mcp", "status": "connected"}]}
```

**Turn 1 — the agent calls `insulation_tads`:**

```json
{"name": "mcp__hic-mcp__insulation_tads", "input": {"top_n": 5}}
```

and this server answers with a real cooltools computation over the bundled data:

```json
{"resolution_used": 10000,
 "windows_bp": [100000, 250000, 500000],
 "ranked_by": "boundary_strength at the 250000 bp window",
 "boundary_counts_per_window": {"100000": 270, "250000": 174, "500000": 109},
 "top_boundaries": [{"locus": "chr17:66,180,000-66,190,000", "strength": 2.762923,
                     "log2_insulation": -1.658919,
                     "windows_detected": [100000, 250000, 500000]}, ...]}
```

> The strongest boundary is at chr17:66,180,000-66,190,000. Now checking its compartment call.

**Turn 2 — the agent follows up with `compartments` at that exact locus:**

```json
{"name": "mcp__hic-mcp__compartments", "input": {"region": "chr17:66,180,000-66,190,000"}}
```

```json
{"resolution_used": 100000,
 "view": "chr17 p/q arms (bundled)",
 "sign_convention": "oriented by the bundled GC track (positive E1 = A = gene-dense/GC-rich)",
 "region_mean_E1": -0.243223, "region_call": "B", "region_sign_consistency": 1.0}
```

**The agent's conclusion, verbatim:**

> **Locus:** chr17:66,180,000–66,190,000
> **Boundary strength:** 2.763 (log2 insulation minimum: −1.659)
> **Detected by diamond windows:** 100 kb, 250 kb, and 500 kb (consistent across all three scales — a robust, not scale-specific, boundary)
> **Ranking basis:** strength at the 250 kb window (out of 174 boundaries called at that scale)
>
> **Compartment call:** B (heterochromatic/gene-poor)
> **Eigenvector value (E1):** −0.243 (sign consistency 1.0 — unambiguous)

---

## Reading this honestly

**What it demonstrates.** The agent was given a two-part scientific question and answered it by chaining two real computations, using the second tool's region argument to interrogate the first tool's result. The numbers in its answer are the ones this server computed from the bundled contact matrix — the boundary locus and strength are the same values the test suite asserts against (see [`data/PROVENANCE.md`](../data/PROVENANCE.md)), and the sign convention it quotes is the server telling it how the eigenvector was oriented.

**One overstatement, left uncorrected.** The model's closing line calls this "the strongest TAD boundary genome-wide". The bundled demo is chr17 only, so the correct scope is "the strongest boundary in this dataset". The transcript is not edited to fix it — a demo transcript that has been improved after the fact is not evidence of anything. It is noted here instead, and it is a fair illustration of why the tools state their own scope in every response: `matrix_summary` reports exactly one chromosome, which is what a careful agent would have checked first.
