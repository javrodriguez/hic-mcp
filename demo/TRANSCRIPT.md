# A real agent session

This is a curated reading of [`raw-session.jsonl`](raw-session.jsonl) — an actual headless Claude Code session driving this server over MCP.
Every quotation below is copied from that file word for word, and the agent's conclusion is quoted in full rather than abridged — so there are no cuts to declare. The raw log is the evidence; this page is only the readable version of it.

**Reproduce it:** `python scripts/capture_demo.py` (see [`README.md`](README.md) in this folder for what that does).

---

**Prompt given to the agent:**

> Using the hic-mcp tools, find the strongest TAD boundary in the bundled demo Hi-C data, then check which compartment that boundary sits in. Report the locus, its boundary strength, which diamond windows detect it, and the compartment call with its eigenvector value.

**Turn 1 — the agent orients itself with `matrix_summary`, then calls `insulation_tads`:**

```json
{"name": "mcp__hic-mcp__insulation_tads", "input": {"top_n": 5}}
```

**Turn 2 — it follows up with `compartments` around the strongest boundary:**

```json
{"name": "mcp__hic-mcp__compartments", "input": {"region": "chr17:66,150,000-66,220,000"}}
```

**The agent's conclusion, in full — every line, unaltered.** Nothing is cut, so nothing has to be declared; you can diff this block against [`raw-session.jsonl`](raw-session.jsonl) and a test does exactly that on every run:

> ## Strongest TAD boundary — bundled Hi-C demo (chr17, HFF Micro-C, hg38)
>
> **Locus:** chr17:66,180,000–66,190,000
> **Boundary strength:** 2.763 (log2 insulation minimum: −1.637)
> **Detected by:** all three diamond windows tested — 100 kb, 250 kb, and 500 kb
>
> This was the top-ranked boundary out of 174 called at the 250 kb window (270 at 100 kb, 109 at 500 kb), and it's one of only a few boundaries robust across all window sizes simultaneously, which is why it ranks #1 overall.
>
> ## Compartment context
>
> This boundary sits almost exactly on an **A/B compartment transition**, not comfortably inside one compartment:
>
> | Bin start | E1 |
> |---|---|
> | 66,000,000 | −0.313 |
> | 66,100,000 | **−0.243** ← boundary's bin (B) |
> | 66,200,000 | **−0.002** ← next bin, essentially zero |
> | 66,300,000 | +0.502 (A) |
> | 66,400,000 | +0.540 (A) |
>
> - The 100 kb bin containing the boundary (66,100,000–66,200,000) has E1 = −0.243, calling it weakly **B**.
> - The immediately adjacent bin (66,200,000–66,300,000) has E1 ≈ 0, and the region flips sharply to strongly **A** (+0.50 to +0.54) just downstream.
> - The tool's own transition flag confirmed this: querying the region returned `region_call: "B"` (mean E1 = −0.123) but explicitly warned *"E1 changes sign within the flanking window: this locus is at or near a compartment transition, so a single-region A/B label understates it."*
>
> **Bottom line:** the strongest TAD boundary in the demo data isn't just a local insulation minimum — it marks (or lies right at the edge of) a B→A compartment switch on chr17q, consistent with the general pattern that strong, multi-scale TAD boundaries often coincide with compartment domain edges.

---

## Why this log has been re-captured four times

**Once because the tools were wrong.** An early capture of this same question concluded the opposite — that the boundary sat "inside a B-compartment domain … rather than a compartment-transition artifact". The fault was this server's: `compartments` reported `region_sign_consistency: 1.0` for a query covering a **single** 100 kb bin. A single bin is trivially consistent with itself, so the figure carried no information, but it reads as confidence and the agent reasonably treated it as such. It never saw that E1 flips sign one bin away. The tools now return `bins_used`, withhold a consistency figure below three bins, and attach a `transition_note` with the flanking track — which is why the session above reads the flip correctly.

**Once because the tools changed.** A later fix gave the analyses a chromosome-arm view, which shifts `log2_insulation` by a per-arm constant. The committed log still held pre-fix numbers and field names the code no longer emits. Rather than edit it, the session was run again.

**Once more, for the same reason — and that is the point.** Another round of fixes changed a response field again, and the committed log went stale again, while this page claimed a test would catch exactly that. No test did: the guard re-derived the documented landmarks, not the log. One was added that replays every recorded tool call against the live server.

**And once because that replay guard was itself hollow.** An external evaluator mutated this log — boundary strength to 42.0, an eigenvector value to −99.9, the boundary counts to 1/1/1 — and the entire honesty suite still passed. The guard compared only strings, integers and booleans, so it skipped `region_mean_E1`, `top_boundaries`, `boundary_counts_per_window`, `eigenvalues` and `E1_track`: every scientific number in the file. It now compares every value to any depth, floats included. The same mutation now fails the build, which is the only reason this page's claim is worth reading.

Every time, the rule was the same: **the transcript is never edited to match the code — the session is re-run.** That is the whole point of committing the raw log beside this page.

External reviewers found all three, on this repository, by reading the demo it was meant to showcase — twice catching a claim this page made about its own integrity that no test was actually enforcing.
