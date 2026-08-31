"""File location, region parsing, and resolution selection for Hi-C inputs.

Pure helpers - no MCP imports. Anticipated failures raise DataError with a
message written for the calling agent to read and act on.
"""

import os
from pathlib import Path

import cooler
import pandas as pd

DEMO_FILENAME = "hff_microc_chr17_hg38.mcool"
GC_TRACK_FILENAME = "gc_100kb.tsv"
ARMS_FILENAME = "chr17_arms.bed"


class DataError(ValueError):
    """A problem with the input file, region, or resolution - message is agent-facing."""


def data_dir() -> Path:
    """The bundled data directory (repo `data/`; override with HIC_MCP_DATA_DIR)."""
    env = os.environ.get("HIC_MCP_DATA_DIR")
    if env:
        return Path(env)
    return Path(__file__).resolve().parents[2] / "data"


def demo_path() -> Path:
    return data_dir() / DEMO_FILENAME


def is_demo(path: Path) -> bool:
    return path.name == DEMO_FILENAME


def resolve_input_path(file: str | None) -> Path:
    """Resolve the user-supplied path, or fall back to the bundled demo file."""
    if file is None:
        p = demo_path()
        if not p.exists():
            raise DataError(
                f"The bundled demo dataset is missing (expected at {p}). "
                "Run from a full checkout of the repository, or pass an explicit "
                "path to a .mcool/.cool file."
            )
        return p
    p = Path(file).expanduser()
    if not p.exists():
        raise DataError(f"No such file: {p}. Pass a path to a local .cool or .mcool file.")
    if p.suffix not in (".cool", ".mcool"):
        raise DataError(f"{p.name} is not a .cool or .mcool file.")
    return p


def list_resolutions(path: Path) -> list[int]:
    """Resolutions available in the file (a plain .cool has exactly one)."""
    uris = cooler.fileops.list_coolers(str(path))
    if uris == ["/"]:
        return [int(cooler.Cooler(str(path)).binsize)]
    res = []
    for uri in uris:
        res.append(int(cooler.Cooler(f"{path}::{uri}").binsize))
    return sorted(res)


def open_matrix(path: Path, resolution: int | None, default: int) -> cooler.Cooler:
    """Open the file at the requested resolution (or the nearest sensible default)."""
    available = list_resolutions(path)
    if resolution is None:
        resolution = default if default in available else available[0]
    if resolution not in available:
        raise DataError(
            f"Resolution {resolution} bp is not in {path.name}. "
            f"Available: {', '.join(str(r) for r in available)} bp."
        )
    uris = cooler.fileops.list_coolers(str(path))
    if uris == ["/"]:
        return cooler.Cooler(str(path))
    return cooler.Cooler(f"{path}::/resolutions/{resolution}")


def parse_region_checked(clr: cooler.Cooler, region: str) -> tuple[str, int, int]:
    """Parse a UCSC-style region string against the file's own chromosomes."""
    try:
        chrom, start, end = cooler.util.parse_region(region, clr.chromsizes)
    except ValueError as e:
        chroms = ", ".join(clr.chromnames[:25])
        raise DataError(
            f"Could not interpret region {region!r} against this file ({e}). "
            f"Use e.g. 'chr17' or 'chr17:1,000,000-2,000,000'. Chromosomes here: {chroms}."
        ) from e
    return str(chrom), int(start), int(end)


def load_gc_track() -> pd.DataFrame:
    """The bundled GC phasing track for the demo file (chrom/start/end/GC at 100 kb)."""
    p = data_dir() / GC_TRACK_FILENAME
    if not p.exists():
        raise DataError(f"Bundled GC track missing (expected at {p}).")
    return pd.read_csv(p, sep="\t")


def load_arms_view() -> pd.DataFrame:
    """The bundled chr17 p/q arm view for the demo file."""
    p = data_dir() / ARMS_FILENAME
    if not p.exists():
        raise DataError(f"Bundled arm view missing (expected at {p}).")
    return pd.read_csv(p, sep="\t", names=["chrom", "start", "end", "name"])
