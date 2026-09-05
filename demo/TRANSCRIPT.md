# A real agent session

This is a curated reading of [`raw-session.jsonl`](raw-session.jsonl) — an actual headless Claude Code session driving this server over MCP.
Every quotation below is copied from that file word for word, and the agent's conclusion is quoted in full rather than abridged — so there are no cuts to declare. The raw log is the evidence; this page is only the readable version of it.

**Reproduce it:** `python scripts/capture_demo.py` (see [`README.md`](README.md) in this folder for what that does).

---

**Prompt given to the agent:**

> Using the hic-mcp tools, find the strongest TAD boundary in the bundled demo Hi-C data, then check which compartment that boundary sits in. Report the locus, its boundary strength, which diamond windows detect it, and the compartment call with its eigenvector value.

**Turn 1 — the agent orients itself with `matrix_summary`, then ranks boundaries with `insulation_tads`:**

```json
{"name": "mcp__hic-mcp__matrix_summary", "input": {}}
{"name": "mcp__hic-mcp__insulation_tads", "input": {"top_n": 5}}
```

**Turn 2 — it calls `compartments` twice, first on the boundary itself and then on the 100 kb bin containing it:**

```json
{"name": "mcp__hic-mcp__compartments", "input": {"region": "chr17:66,150,000-66,220,000"}}
{"name": "mcp__hic-mcp__compartments", "input": {"region": "chr17:66,100,000-66,200,000"}}
```

**The agent's conclusion, in full — every line, unaltered.** Nothing is cut, so nothing has to be declared; you can diff this block against [`raw-session.jsonl`](raw-session.jsonl) and a test does exactly that on every run:

> ## Strongest TAD boundary
>
> **Locus:** chr17:66,180,000–66,190,000 (10 kb bin, hg38, chr17 q-arm)
>
> **Boundary strength:** 2.763 (log2 insulation = −1.637) — the single strongest call genome-wide (chr17) in this dataset, ranked by the 250 kb diamond window.
>
> **Diamond windows detecting it:** all three tested window sizes — 100 kb, 250 kb, and 500 kb — flagged this bin as a boundary, indicating it's a robust, multi-scale insulator rather than an artifact of one window choice.
>
> ## Compartment context
>
> - **Containing 100 kb bin (66,100,000–66,200,000): E1 = −0.243 → compartment B**
> - However, this boundary sits almost exactly on an **A/B transition**: E1 climbs steadily through the neighborhood (−0.42 at 65.8 Mb → −0.31 at 66.0 Mb → −0.24 at the boundary bin → −0.002 at 66.2 Mb → +0.50 at 66.3 Mb, flipping firmly into A just one bin downstream).
>
> So the strongest insulation boundary in the demo data calls as **B compartment** (E1 ≈ −0.24), but it's positioned right at the shoulder of a B→A compartment switch — consistent with the common pattern where strong, robust TAD boundaries coincide with compartment domain edges.

---

## Why this log keeps being re-captured

**It is re-captured, never edited, every time the tools change what they return.** That rule is the reason the page is worth reading: if a number here disagreed with the server, the honest fix is a new session, not a new sentence. Six re-captures so far, each for a reason worth stating:

1. **The tools were wrong.** An early capture concluded the opposite of the truth — that the strongest boundary sat safely inside a B compartment. `compartments` had reported `region_sign_consistency: 1.0` for a query covering a **single** 100 kb bin, which is trivially true of one bin and reads as certainty. The agent believed it and never saw that E1 flips sign one bin away. The tools now report `bins_used`, withhold a consistency figure below three bins, and attach a `transition_note` with the flanking track.
2. **A fix changed the numbers.** Giving the analyses a chromosome-arm view shifts `log2_insulation` by a per-arm constant, and the committed log still held pre-fix values.
3. **The same again**, after another round of fixes — while this page claimed a test would catch exactly that. None did: the guard re-derived the documented landmarks, not the log. One was added that replays every recorded call against the live server.
4. **That replay guard turned out to be hollow.** An external evaluator mutated this file — boundary strength to 42.0, an eigenvector value to −99.9, the counts to 1/1/1 — and the whole honesty suite still passed, because the comparison accepted only strings, integers and booleans and skipped every scientific number in the file. It now compares every value to any depth, and the attack was reproduced to watch it fail.
5. **Wording changed in two response fields** (`ranked_by` and `method`), which the now-strict guard caught immediately.
6. **Three separate defects changed what the tools return**, all found in one round: `virtual_4c` raised a bare library error for any anchor not sitting exactly on the bin grid, `contacts_at_locus` reported a bin count that contradicted the matrix beside it, and `compartments` produced a confident A/B call from a phasing track that correlated with nothing. The responses now carry a Pearson r for the orientation, disclose a snapped anchor, and count bins the way the file does.

External reviewers found the first four and the sixth, on this repository, by reading the demo it was meant to showcase — twice catching a claim this page made about its own integrity that no test was actually enforcing.
