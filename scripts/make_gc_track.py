#!/usr/bin/env python
"""Rebuild (and verify) the bundled GC phasing track from public reference sequence.

The track orients the compartment eigenvector: without it the sign of E1 is
arbitrary, so this file decides which side is called "A". That makes it as
load-bearing as the contact matrix beside it, and it gets the same treatment -
a re-runnable derivation and a verification mode.

Source: the UCSC hg38 chr17 FASTA (public, no account), binned to match the
demo cooler and reduced to GC fraction per bin.

Usage:
    uv run python scripts/make_gc_track.py --out /tmp/gc.tsv \
        [--verify-against data/gc_100kb.tsv]
"""

import argparse
import gzip
import sys
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[1]
FASTA_URL = "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr17.fa.gz"
BINSIZE = 100_000
CHROM = "chr17"


def fetch_sequence(cache_dir: Path) -> str:
    cache_dir.mkdir(parents=True, exist_ok=True)
    local = cache_dir / "chr17.fa.gz"
    if not local.exists():
        print(f"downloading {FASTA_URL} ...")
        urllib.request.urlretrieve(FASTA_URL, local)  # noqa: S310 - fixed public URL
    with gzip.open(local, "rt") as fh:
        lines = [ln.strip() for ln in fh if not ln.startswith(">")]
    return "".join(lines).upper()


def gc_per_bin(seq: str, binsize: int) -> pd.DataFrame:
    rows = []
    for start in range(0, len(seq), binsize):
        chunk = seq[start : start + binsize]
        acgt = sum(chunk.count(b) for b in "ACGT")
        gc = sum(chunk.count(b) for b in "GC")
        rows.append(
            {
                "chrom": CHROM,
                "start": start,
                "end": min(start + binsize, len(seq)),
                # N-masked bins carry no usable base composition
                "GC": (gc / acgt) if acgt else np.nan,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True)
    ap.add_argument("--cache-dir", default="/tmp/hic-mcp-reference")
    ap.add_argument("--verify-against")
    args = ap.parse_args()

    seq = fetch_sequence(Path(args.cache_dir))
    track = gc_per_bin(seq, BINSIZE)
    track.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {args.out} - {len(track)} bins of {BINSIZE:,} bp")

    if args.verify_against:
        shipped = pd.read_csv(args.verify_against, sep="\t")
        if len(shipped) != len(track):
            sys.exit(f"MISMATCH: shipped has {len(shipped)} rows, rebuild has {len(track)}")
        if not np.array_equal(shipped["start"].to_numpy(), track["start"].to_numpy()):
            sys.exit("MISMATCH: bin boundaries differ")
        a, b = shipped["GC"].to_numpy(), track["GC"].to_numpy()
        if not np.array_equal(np.isnan(a), np.isnan(b)):
            sys.exit("MISMATCH: different bins are N-masked")
        if not np.allclose(a, b, rtol=1e-6, equal_nan=True):
            worst = np.nanmax(np.abs(a - b))
            sys.exit(f"MISMATCH: GC values differ (max {worst:.2e})")
        print("VERIFIED: rebuild matches the shipped GC track")


if __name__ == "__main__":
    main()
