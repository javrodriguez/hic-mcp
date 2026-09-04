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

> The strongest boundary is at **chr17:66,180,000-66,190,000** (strength 2.76292). Now checking its compartment context.

**Turn 2 — it follows up with `compartments` around that locus:**

```json
{"name": "mcp__hic-mcp__compartments", "input": {"region": "chr17:66,150,000-66,220,000"}}
```

**The agent's conclusion, in full — every line, unaltered.** Nothing is cut, so nothing has to be declared; you can diff this block against [`raw-session.jsonl`](raw-session.jsonl) and a test does exactly that on every run:

> ## Strongest TAD boundary — bundled demo Hi-C (HFF Micro-C, chr17, hg38)
>
> **Locus:** chr17:66,180,000–66,190,000
> **Boundary strength:** 2.763 (log2 insulation = −1.637)
> **Detected by:** all three diamond windows — 100 kb, 250 kb, and 500 kb
>
> This is the single strongest of 174 boundaries called at the 250 kb window (270 at 100 kb, 109 at 500 kb), and the fact that it's picked up consistently across all three window sizes indicates a robust, multi-scale boundary rather than an artifact of one window setting.
>
> ### Compartment context
>
> Querying `compartments` for the flanking region (chr17:66,150,000–66,220,000, 100 kb resolution):
>
> - **Region call:** B (region mean E1 = −0.1228)
> - **But flagged as a compartment transition** — the E1 track crosses zero right at this locus:
>
> | bin start | E1 | call |
> |---|---|---|
> | 65,800,000 | −0.422 | B |
> | 65,900,000 | −0.409 | B |
> | 66,000,000 | −0.313 | B |
> | 66,100,000 | −0.243 | B |
> | **66,200,000** | **−0.002** | **~boundary** |
> | 66,300,000 | +0.502 | A |
> | 66,400,000 | +0.540 | A |
> | 66,500,000 | +0.362 | A |
>
> **Interpretation:** the strongest TAD boundary in this dataset sits almost exactly at a B→A compartment switch — its immediate 100 kb bin (E1 ≈ −0.12, nominally "B") is really the tail end of the B domain right before it flips to A one bin downstream. So rather than a boundary sitting cleanly inside one compartment, this is a case where insulation and compartmentalization line up: the strongest structural (TAD) boundary coincides with the strongest compartment transition nearby.

---

## Why this log has been re-captured twice

**Once because the tools were wrong.** An early capture of this same question concluded the opposite — that the boundary sat "inside a B-compartment domain … rather than a compartment-transition artifact". The fault was this server's: `compartments` reported `region_sign_consistency: 1.0` for a query covering a **single** 100 kb bin. A single bin is trivially consistent with itself, so the figure carried no information, but it reads as confidence and the agent reasonably treated it as such. It never saw that E1 flips sign one bin away. The tools now return `bins_used`, withhold a consistency figure below three bins, and attach a `transition_note` with the flanking track — which is why the session above reads the flip correctly.

**Once because the tools changed.** A later fix gave the analyses a chromosome-arm view, which shifts `log2_insulation` by a per-arm constant. The committed log still held pre-fix numbers and field names the code no longer emits. Rather than edit it, the session was run again.

Both times the rule was the same: **the transcript is never edited to match the code — the session is re-run.** That is the whole point of committing the raw log beside this page. A test now re-derives every documented figure from the live code, so this particular drift fails the build instead of shipping.

An external reviewer found both problems, on this repository, by reading the demo it was meant to showcase.
