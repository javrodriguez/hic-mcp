# A real agent session

This is a curated reading of [`raw-session.jsonl`](raw-session.jsonl) — an actual headless Claude Code session driving this server over MCP.
Every quotation below is copied from that file word for word, and the agent's conclusion is quoted in full rather than abridged — so there are no cuts to declare. The raw log is the evidence; this page is only the readable version of it.

**Reproduce it:** `python scripts/capture_demo.py` (see [`README.md`](README.md) in this folder for what that does).

---

**Prompt given to the agent:**

> Using the hic-mcp tools, find the strongest TAD boundary in the bundled demo Hi-C data, then check which compartment that boundary sits in. Report the locus, its boundary strength, which diamond windows detect it, and the compartment call with its eigenvector value.

**Turn 1 — the agent calls `insulation_tads`:**

```json
{"name": "mcp__hic-mcp__insulation_tads", "input": {"top_n": 5}}
```

> The strongest boundary is at chr17:66,180,000-66,190,000. Now checking its compartment call.

**Turn 2 — it follows up with `compartments` at that exact locus:**

```json
{"name": "mcp__hic-mcp__compartments", "input": {"region": "chr17:66,180,000-66,190,000"}}
```

**The agent's conclusion, in full — every line, unaltered.** Nothing is cut, so nothing has to be declared; you can diff this block against [`raw-session.jsonl`](raw-session.jsonl) and a test does exactly that on every run:

> ## Strongest TAD boundary — chr17:66,180,000-66,190,000 (chr17q, hg38)
>
> - **Boundary strength:** 2.763 (top-ranked of 174 boundaries at the 250 kb window; log2 insulation = -1.637)
> - **Detected by all three diamond windows:** 100 kb, 250 kb, and 500 kb — a robust, scale-consistent call
> - **Compartment call:** **B**, with E1 = **-0.243** (100 kb resolution, oriented by the bundled GC-content track)
>
> One caveat worth flagging: this boundary sits right at a compartment transition. The E1 track around it goes from clearly negative (-0.42 at 65.8 Mb) through the boundary bin (-0.24 at 66.1 Mb) and flips positive just downstream (+0.50 at 66.3 Mb). So the boundary is functionally at (or just inside) a B→A switch — a biologically sensible place to find the strongest insulation, since TAD boundaries often coincide with compartment transitions.

---

## Why this log has been re-captured three times

**Once because the tools were wrong.** An early capture of this same question concluded the opposite — that the boundary sat "inside a B-compartment domain … rather than a compartment-transition artifact". The fault was this server's: `compartments` reported `region_sign_consistency: 1.0` for a query covering a **single** 100 kb bin. A single bin is trivially consistent with itself, so the figure carried no information, but it reads as confidence and the agent reasonably treated it as such. It never saw that E1 flips sign one bin away. The tools now return `bins_used`, withhold a consistency figure below three bins, and attach a `transition_note` with the flanking track — which is why the session above reads the flip correctly.

**Once because the tools changed.** A later fix gave the analyses a chromosome-arm view, which shifts `log2_insulation` by a per-arm constant. The committed log still held pre-fix numbers and field names the code no longer emits. Rather than edit it, the session was run again.

**Once more, for the same reason — and that is the point.** Another round of fixes changed a response field again, and the committed log went stale again, while this page claimed a test would catch exactly that. No test did: the guard re-derived the documented landmarks, not the log. There is now one that replays every recorded tool call against the live server and compares the responses, so this page's claim is enforced rather than asserted.

Every time, the rule was the same: **the transcript is never edited to match the code — the session is re-run.** That is the whole point of committing the raw log beside this page.

External reviewers found all three, on this repository, by reading the demo it was meant to showcase — twice catching a claim this page made about its own integrity that no test was actually enforcing.
