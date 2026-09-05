# A real agent session

This is a curated reading of [`raw-session.jsonl`](raw-session.jsonl) — an actual headless Claude Code session driving this server over MCP.
Every quotation below is copied from that file word for word, and the agent's conclusion is quoted in full rather than abridged — so there are no cuts to declare. The raw log is the evidence; this page is only the readable version of it.

**Reproduce it:** `python scripts/capture_demo.py` (see [`README.md`](README.md) in this folder for what that does).

---

**Prompt given to the agent:**

> Using the hic-mcp tools, find the strongest TAD boundary in the bundled demo Hi-C data, then check which compartment that boundary sits in. Report the locus, its boundary strength, which diamond windows detect it, and the compartment call with its eigenvector value.

**Turn 1 — the agent orients itself, then ranks boundaries:**

```json
{"name": "mcp__hic-mcp__matrix_summary", "input": {}}
{"name": "mcp__hic-mcp__insulation_tads", "input": {"top_n": 5}}
```

**Turn 2 — it follows up on the boundary it found:**

```json
{"name": "mcp__hic-mcp__compartments", "input": {"region": "chr17:66,100,000-66,270,000"}}
{"name": "mcp__hic-mcp__compartments", "input": {"region": "chr17:66,100,000-66,200,000"}}
```

**The agent's conclusion, in full — every line, unaltered.** Nothing is cut, so nothing has to be declared; you can diff this block against [`raw-session.jsonl`](raw-session.jsonl) and a test does exactly that on every run:

> ## Strongest TAD boundary
>
> **Locus:** chr17:66,180,000–66,190,000 (hg38, chr17q)
> **Boundary strength:** 2.763 (log2 insulation dip: −1.637) — the highest of 109–270 boundaries called across the chromosome depending on window size
> **Diamond windows detecting it:** all three tested — 100 kb, 250 kb, and 500 kb — meaning the insulation dip is consistent across scales, a hallmark of a robust structural boundary (not a small-window artifact)
>
> ## Compartment context
>
> This boundary sits almost exactly on an **A/B compartment transition**:
>
> | 100 kb bin | E1 |
> |---|---|
> | 66,000,000–66,100,000 | −0.313 |
> | **66,100,000–66,200,000** (boundary's bin) | **−0.243 → B** |
> | 66,200,000–66,300,000 | −0.002 (essentially 0) |
> | 66,300,000–66,400,000 | +0.502 → A |
>
> The bin containing the boundary (66,100,000–66,200,000) is formally called **compartment B** with E1 = −0.243, but the tool explicitly flags this as a transition zone — E1 climbs steadily from −0.42 five bins upstream through near-zero just downstream, crossing into strongly positive (A, E1 ≈ +0.50) one bin later at 66,300,000–66,400,000.
>
> **Interpretation:** the strongest insulation boundary in the demo chr17 dataset coincides with a genuine A→B compartment switch, which is biologically expected — TAD boundaries often demarcate compartment domains, and the strongest ones frequently sit right at these A/B transitions rather than inside a homogeneous compartment.

---

## Why this log keeps being re-captured

**It is re-captured, never edited, every time the tools change what they return.** That rule is the reason the page is worth reading: if a number here disagreed with the server, the honest fix is a new session, not a new sentence. Seven re-captures so far, each for a reason worth stating:

1. **The tools were wrong.** An early capture concluded the opposite of the truth — that the strongest boundary sat safely inside a B compartment. `compartments` had reported `region_sign_consistency: 1.0` for a query covering a **single** 100 kb bin, which is trivially true of one bin and reads as certainty. The agent believed it and never saw that E1 flips sign one bin away. The tools now report `bins_used`, withhold a consistency figure below three bins, and attach a `transition_note` with the flanking track.
2. **A fix changed the numbers.** Giving the analyses a chromosome-arm view shifts `log2_insulation` by a per-arm constant, and the committed log still held pre-fix values.
3. **The same again**, after another round of fixes — while this page claimed a test would catch exactly that. None did: the guard re-derived the documented landmarks, not the log. One was added that replays every recorded call against the live server.
4. **That replay guard turned out to be hollow.** An external evaluator mutated this file — boundary strength to 42.0, an eigenvector value to −99.9, the counts to 1/1/1 — and the whole honesty suite still passed, because the comparison accepted only strings, integers and booleans and skipped every scientific number in the file. It now compares every value to any depth, and the attack was reproduced to watch it fail.
5. **Wording changed in two response fields** (`ranked_by` and `method`), which the now-strict guard caught immediately.
6. **Three separate defects changed what the tools return**, all found in one round: `virtual_4c` raised a bare library error for any anchor not sitting exactly on the bin grid, `contacts_at_locus` reported a bin count that contradicted the matrix beside it, and `compartments` produced a confident A/B call from a phasing track that correlated with nothing. The responses now carry a Pearson r for the orientation, disclose a snapped anchor, and count bins the way the file does.
7. **The round-4 fixes were themselves defective, and the next round found it.** The orientation check measured one correlation pooled over the whole file, while the eigenvector's sign is flipped independently inside each view region — so a track that oriented one chromosome arm and not the other produced confident A/B calls on the arm it never oriented, flipping with the noise seed. The correlation is now measured, reported and gated per region. Two more response fields were saying things that were not so: a `coverage_note` claiming 10 Mb over a profile carrying 20 Mb, and a `profile_note` reporting a bin count far larger than the list beside it.

External reviewers found the first four and the sixth, on this repository, by reading the demo it was meant to showcase — twice catching a claim this page made about its own integrity that no test was actually enforcing.
