# A real agent session

This is a curated reading of [`raw-session.jsonl`](raw-session.jsonl) — an actual headless Claude Code session driving this server over MCP.
Every quotation below is copied from that file; the raw log is the evidence, this page is only the readable version of it.

**Reproduce it:** `python scripts/capture_demo.py` (see [`README.md`](README.md) in this folder for what that does).

---

**Prompt given to the agent:**

> Using the hic-mcp tools, find the strongest TAD boundary in the bundled demo Hi-C data, then check which compartment that boundary sits in. Report the locus, its boundary strength, which diamond windows detect it, and the compartment call with its eigenvector value.

**Turn 1 — the agent orients itself with `matrix_summary`, then calls `insulation_tads`:**

```json
{"name": "mcp__hic-mcp__insulation_tads", "input": {"top_n": 5}}
```

> The strongest boundary is at **chr17:66,180,000-66,190,000** (strength 2.76292, detected across all three diamond windows). Now let's check the compartment call for that locus.

**Turn 2 — it follows up with `compartments` around that locus:**

```json
{"name": "mcp__hic-mcp__compartments", "input": {"region": "chr17:66,150,000-66,220,000"}}
```

**The agent's conclusion, verbatim:**

> **Locus:** chr17:66,180,000–66,190,000 (10 kb bin)
> **Boundary strength:** 2.763 (log2 insulation minimum: −1.659)
> **Detected by diamond windows:** all three — 100 kb, 250 kb, and 500 kb
>
> | bin start | E1 |
> |---|---|
> | 66,100,000 | −0.243 |
> | 66,200,000 | **−0.002** |
> | 66,300,000 | +0.502 |
>
> The boundary bin … falls in the 66,100,000–66,200,000 bin, where **E1 = −0.243** — a weak **B compartment** call, but E1 is climbing toward zero and flips positive (A) just one bin downstream. The tool explicitly flagged this as a compartment transition rather than clean B territory.
>
> **Bottom line:** the strongest TAD boundary in the demo dataset sits almost exactly at an A/B compartment transition — weakly B on the boundary bin itself, flipping to strongly A immediately downstream. This is a biologically sensible result: TAD boundaries frequently coincide with compartment switches.

---

## Why this transcript was recaptured

An earlier capture of this same question reached the **opposite** conclusion — that the boundary sat "inside a B-compartment domain … a robust structural boundary rather than a compartment-transition artifact." That was wrong, and the fault was this server's, not the model's.

The old `compartments` response reported `region_sign_consistency: 1.0` for a query covering a **single** 100 kb bin. A single bin is trivially "100% consistent with itself", so the number carried no information — but it reads as confidence, and the agent reasonably treated it as such. It never saw that E1 flips sign one bin away.

Three things changed in response, and the difference above is the result:

- `region_sign_consistency` is now `null` when fewer than three bins support it, with a note saying why.
- `bins_used` is always reported, so "1 bin" is visible rather than implied.
- Every region query returns the **flanking** E1 track, and a `transition_note` when the sign changes nearby.

The transcript was not edited to fix the conclusion — the tools were fixed and the session was run again. That distinction is the whole point of committing the raw log: you can check this page against it, and you can re-run the capture yourself.

An external reviewer found this defect, on this repository, by reading the demo it was meant to showcase.
