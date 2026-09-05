# A real agent session

This is a curated reading of [`raw-session.jsonl`](raw-session.jsonl) — an actual headless Claude Code session driving this server over MCP.
Every quotation below is copied from that file word for word, and the agent's conclusion is quoted in full rather than abridged — so there are no cuts to declare. The raw log is the evidence; this page is only the readable version of it.

**Reproduce it:** `python scripts/capture_demo.py` (see [`README.md`](README.md) in this folder for what that does).

---

**Prompt given to the agent:**

> Using the hic-mcp tools, find the strongest TAD boundary in the bundled demo Hi-C data, then check which compartment that boundary sits in. Report the locus, its boundary strength, which diamond windows detect it, and the compartment call with its eigenvector value.

**Turn 1 — the agent orients itself with `matrix_summary`, then calls `insulation_tads`:**

```json
{"name": "mcp__hic-mcp__insulation_tads", "input": {"top_n": 1}}
```

**Turn 2 — it follows up with `compartments` at the boundary it found:**

```json
{"name": "mcp__hic-mcp__compartments", "input": {"region": "chr17:66,180,000-66,190,000"}}
```

**The agent's conclusion, in full — every line, unaltered.** Nothing is cut, so nothing has to be declared; you can diff this block against [`raw-session.jsonl`](raw-session.jsonl) and a test does exactly that on every run:

> ## Strongest TAD boundary — chr17:66,180,000–66,190,000 (hg38, HFF Micro-C demo)
>
> **TAD boundary**
> - Boundary strength: **2.76292** (log2 insulation dip: −1.63746)
> - Detected consistently across **all three diamond windows**: 100 kb, 250 kb, and 500 kb — the most robust boundary in the dataset, and the only one in the top hit list flagged at every window.
>
> **Compartment context**
> - E1 at this bin (100 kb res): **−0.243223** → nominal call **B** (oriented so positive = A, from the bundled GC-content phasing track)
> - But this boundary sits almost exactly on a **compartment transition**: the flanking E1 track flips sign right after it —
>
> | bin start | E1 |
> |---|---|
> | 65,800,000 | −0.422 |
> | 65,900,000 | −0.409 |
> | 66,000,000 | −0.313 |
> | **66,100,000** | **−0.243** (boundary bin) |
> | 66,200,000 | −0.002 |
> | 66,300,000 | +0.502 |
> | 66,400,000 | +0.540 |
>
> So the locus is technically called **B**, but it's the B-side bin immediately adjacent to a sharp B→A switch about 100–200 kb downstream — a classic case of the strongest insulation boundary in the region coinciding with an A/B compartment transition, rather than sitting deep inside one compartment.

---

## Why this log keeps being re-captured

**It is re-captured, never edited, every time the tools change what they return.** That rule is the reason the page is worth reading: if a number here disagreed with the server, the honest fix is a new session, not a new sentence. Five re-captures so far, each for a reason worth stating:

1. **The tools were wrong.** An early capture concluded the opposite of the truth — that the strongest boundary sat safely inside a B compartment. `compartments` had reported `region_sign_consistency: 1.0` for a query covering a **single** 100 kb bin, which is trivially true of one bin and reads as certainty. The agent believed it and never saw that E1 flips sign one bin away. The tools now report `bins_used`, withhold a consistency figure below three bins, and attach a `transition_note` with the flanking track.
2. **A fix changed the numbers.** Giving the analyses a chromosome-arm view shifts `log2_insulation` by a per-arm constant, and the committed log still held pre-fix values.
3. **The same again**, after another round of fixes — while this page claimed a test would catch exactly that. None did: the guard re-derived the documented landmarks, not the log. One was added that replays every recorded call against the live server.
4. **That replay guard turned out to be hollow.** An external evaluator mutated this file — boundary strength to 42.0, an eigenvector value to −99.9, the counts to 1/1/1 — and the whole honesty suite still passed, because the comparison accepted only strings, integers and booleans and skipped every scientific number in the file. It now compares every value to any depth, and the attack was reproduced to watch it fail.
5. **Wording changed in two response fields** (`ranked_by` and `method`), which the now-strict guard caught immediately.

External reviewers found the first four, on this repository, by reading the demo it was meant to showcase — twice catching a claim this page made about its own integrity that no test was actually enforcing.
