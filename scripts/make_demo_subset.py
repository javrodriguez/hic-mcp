#!/usr/bin/env python
"""Rebuild (and verify) the bundled demo subset from the public source file.

Source: HFF Micro-C, hg38 (Krietenstein et al. 2020, Mol Cell 78:554-565),
4DN accession 4DNESWST3UBH, obtained via the Open2C cooltools test-data
registry (https://osf.io/3h9js/download; anonymous, no account needed).

Road: download and md5-verify the source .mcool; extract the chr17 cis pixel
block at the 10 kb base resolution with the cooler Python API; write a fresh
single-resolution .cool with its own rebased bin table; `cooler zoomify` to
100 kb and 1 Mb with ICE balancing recomputed at every level; stamp assembly
and provenance metadata on every resolution.

Balancing note: the source weights are a genome-wide ICE solution over
{chr2, chr17}. Once chr2 is dropped that solution no longer satisfies the
marginals of the remaining matrix, so the weight column is deliberately NOT
carried over - it is recomputed from scratch by `cooler zoomify --balance`.

The rebuilt file is NOT byte-identical to the shipped one (cooler writes a
creation date), so verification compares the pixel tables and bin tables of
all three resolutions instead: `--verify-against data/hff_microc_chr17_hg38.mcool`.

Usage:
    python scripts/make_demo_subset.py --out /tmp/rebuild.mcool \
        [--source /path/to/hff_microc_hg38.mcool] [--cache-dir /tmp/hic-mcp-src] \
        [--verify-against data/hff_microc_chr17_hg38.mcool]
"""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import urllib.request

import cooler
import h5py
import numpy as np

SOURCE_URL = "https://osf.io/3h9js/download"
SOURCE_MD5 = "e4a0fc25c8dc3d38e9065fd74c565dd1"
CHROM = "chr17"
BASE_RES = 10000
RESOLUTIONS = [10000, 100000, 1000000]


def md5_of(path: str) -> str:
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fetch_source(cache_dir: str) -> str:
    """Download the public source .mcool into cache_dir, verifying its md5."""
    os.makedirs(cache_dir, exist_ok=True)
    dest = os.path.join(cache_dir, "hff_microc_hg38.mcool")
    if os.path.exists(dest) and md5_of(dest) == SOURCE_MD5:
        print(f"source already cached: {dest}")
        return dest
    print(f"downloading {SOURCE_URL} (~152 MB) ...")
    tmp = dest + ".part"
    urllib.request.urlretrieve(SOURCE_URL, tmp)  # noqa: S310 - fixed public URL
    got = md5_of(tmp)
    if got != SOURCE_MD5:
        os.remove(tmp)
        raise SystemExit(f"source md5 mismatch: got {got}, expected {SOURCE_MD5}")
    os.replace(tmp, dest)
    print(f"downloaded and md5-verified: {dest}")
    return dest


def build(source: str, out: str) -> None:
    src = cooler.Cooler(f"{source}::/resolutions/{BASE_RES}")

    # bins: chr17 only, reindexed from 0, weight column dropped
    bins = src.bins().fetch(CHROM)
    offset = int(bins.index[0])
    bins = bins[["chrom", "start", "end"]].reset_index(drop=True)
    bins["chrom"] = bins["chrom"].astype(str)

    # pixels: chr17 x chr17 cis block, bin ids rebased
    pix = src.matrix(as_pixels=True, balance=False).fetch(CHROM)
    pix = pix[["bin1_id", "bin2_id", "count"]].copy()
    pix["bin1_id"] -= offset
    pix["bin2_id"] -= offset

    assert pix["bin1_id"].min() >= 0
    assert pix["bin2_id"].max() < len(bins)
    assert (pix["bin2_id"] >= pix["bin1_id"]).all(), "expected symmetric-upper"
    print(f"bins={len(bins):,}  pixels={len(pix):,}  counts={int(pix['count'].sum()):,}")

    with tempfile.TemporaryDirectory() as tmpdir:
        base = os.path.join(tmpdir, f"subset_base_{BASE_RES}.cool")
        cooler.create_cooler(
            base,
            bins,
            pix,
            dtypes={"count": np.int32},
            assembly="hg38",  # verified: chromsizes match UCSC hg38 exactly
            ordered=True,
            symmetric_upper=True,
            metadata={
                "source": "HFF Micro-C (Krietenstein et al. 2020, Mol Cell 78:554-565)",
                "source_file": SOURCE_URL,
                "source_md5": SOURCE_MD5,
                "subset": f"{CHROM} cis only, base {BASE_RES} bp",
                "balancing": "recomputed with cooler balance after subsetting",
            },
        )
        if os.path.exists(out):
            os.remove(out)
        cmd = [
            sys.executable, "-m", "cooler", "zoomify",
            "--resolutions", ",".join(str(r) for r in RESOLUTIONS),
            "--balance",
            "--balance-args", "--cis-only",
            "--out", out,
            base,
        ]
        print("+", " ".join(cmd))
        subprocess.run(cmd, check=True)

    # zoomify does not propagate assembly/metadata to the coarsened levels;
    # stamp every resolution so the shipped file is self-describing.
    meta = {
        "source": "HFF Micro-C (Krietenstein et al. 2020, Mol Cell 78:554-565)",
        "source_4dn_accession": "4DNESWST3UBH",
        "source_file": SOURCE_URL,
        "source_md5": SOURCE_MD5,
        "subset": f"{CHROM} cis only; base {BASE_RES} bp; coarsened to 100000, 1000000",
        "balancing": "ICE recomputed after subsetting (cooler balance --cis-only)",
    }
    with h5py.File(out, "r+") as f:
        for r in f["resolutions"]:
            f["resolutions"][r].attrs["genome-assembly"] = "hg38"
            f["resolutions"][r].attrs["metadata"] = json.dumps(meta)
    print(f"wrote {out} ({os.path.getsize(out) / 1e6:.2f} MB)")


def verify_against(rebuilt: str, shipped: str) -> None:
    """Compare bin and pixel tables of every resolution; exit non-zero on any drift."""
    for res in RESOLUTIONS:
        a = cooler.Cooler(f"{rebuilt}::/resolutions/{res}")
        b = cooler.Cooler(f"{shipped}::/resolutions/{res}")
        cols = ["chrom", "start", "end"]
        for name, col_a, col_b in (
            ("bins", a.bins()[:][cols], b.bins()[:][cols]),
            ("pixels", a.pixels()[:], b.pixels()[:]),
        ):
            cols_a = col_a.select_dtypes(include=[np.number]).to_numpy()
            cols_b = col_b.select_dtypes(include=[np.number]).to_numpy()
            if not np.array_equal(cols_a, cols_b):
                raise SystemExit(f"MISMATCH at {res} bp: {name} tables differ")
        print(f"{res} bp: bins + pixels identical")
    print("VERIFIED: rebuild matches the shipped subset at every resolution")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", required=True, help="output .mcool path")
    ap.add_argument("--source", help="existing source .mcool (skips download; md5 still checked)")
    ap.add_argument("--cache-dir", default=os.path.join(tempfile.gettempdir(), "hic-mcp-source"),
                    help="where the 152 MB source download is cached")
    ap.add_argument("--verify-against", help="shipped .mcool to compare the rebuild against")
    args = ap.parse_args()

    if args.source:
        got = md5_of(args.source)
        if got != SOURCE_MD5:
            raise SystemExit(f"--source md5 mismatch: got {got}, expected {SOURCE_MD5}")
        source = args.source
    else:
        source = fetch_source(args.cache_dir)

    build(source, args.out)
    if args.verify_against:
        verify_against(args.out, args.verify_against)


if __name__ == "__main__":
    main()
