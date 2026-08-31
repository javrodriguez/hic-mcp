# Example prompts

Everything here runs against the bundled demo dataset (human chr17, HFF Micro-C, hg38) with no setup beyond the quickstart.
Add `using my file /path/to/data.mcool` to any of them to work on your own data.

## Getting oriented

- "What's in the demo Hi-C file — which chromosomes, what resolutions, how many contacts?"
- "Is this data ICE-balanced, and where did it come from?"

## TADs and insulation

- "Call TAD boundaries on chr17 between 65 and 67 Mb and rank them by strength."
- "Find the strongest TAD boundary in the demo data, then tell me which compartment it sits in."
- "Which boundaries are called at all three diamond windows? Those are the robust ones."
- "Run insulation at 100 kb instead of 10 kb — does the strongest boundary survive the coarser view?"

## Compartments

- "Is chr17:50.1-51.1 Mb in the A or B compartment? How consistent is the eigenvector across that block?"
- "Compare chr17:50.1-51.1 Mb with chr17:51.4-52.4 Mb — is there a compartment flip between them?"
- "What fraction of chr17 is in the A compartment, and how is the sign oriented?"

## Contact structure

- "Show me the observed/expected matrix for chr17:50-52.5 Mb. Is there enrichment off the diagonal?"
- "What's the P(s) slope for this dataset, and over what range was it fitted?"
- "Run a virtual 4C from chr17:63 Mb and describe how contacts decay with distance."
- "Compare contacts within chr17:50-50.5 Mb against contacts between that region and chr17:60-60.5 Mb."

## Probing the honest edges

These are worth trying because the answer is a clear explanation rather than a number:

- "Run a virtual 4C from the middle of the chr17 centromere." → refuses, naming the ICE-filtered region
- "Give me the O/E matrix for chr17:20-25 Mb." → refuses, because that span crosses the centromeric gap between two arms
- "Call TADs at 1 Mb resolution." → runs, and says those are megabase-scale features rather than TADs
