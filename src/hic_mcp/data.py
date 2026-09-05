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


def _unreadable(path: Path, exc: Exception) -> DataError:
    """A file that will not open is a fixable input, not an internal defect."""
    head = b""
    try:
        head = path.open("rb").read(64)
    except OSError:
        pass
    if head.startswith(b"version https://git-lfs"):
        why = (
            "it is a Git-LFS pointer, not the data itself - run `git lfs pull` in the "
            "repository that provided it"
        )
    elif not head.startswith(b"\x89HDF"):
        why = "it is not an HDF5 file at all (a .cool/.mcool is HDF5)"
    else:
        why = (
            "it is HDF5 but not a readable cooler - most often a truncated or "
            "partially-downloaded file"
        )
    return DataError(f"Could not open {path.name}: {why}. ({type(exc).__name__}: {exc})")


def _resolution_uris(path: Path) -> dict[int, str]:
    """Every cooler inside the file, keyed by bin size -> its actual URI.

    Multi-resolution files are not always laid out as /resolutions/<binsize>; opening
    the URI that was really found is what keeps a custom layout from raising a KeyError.
    """
    try:
        uris = cooler.fileops.list_coolers(str(path))
        if uris == ["/"]:
            return {int(cooler.Cooler(str(path)).binsize): "/"}
        if not uris:
            raise IndexError("no coolers in file")
        return {int(cooler.Cooler(f"{path}::{u}").binsize): u for u in uris}
    except DataError:
        raise
    except Exception as e:  # noqa: BLE001 - every failure here is about the file
        raise _unreadable(path, e) from e


def list_resolutions(path: Path) -> list[int]:
    """Resolutions available in the file (a plain .cool has exactly one)."""
    return sorted(_resolution_uris(path))


def open_matrix(path: Path, resolution: int | None, default: int) -> cooler.Cooler:
    """Open the file at the requested resolution (or the nearest sensible default)."""
    by_res = _resolution_uris(path)
    available = sorted(by_res)
    if resolution is None:
        if default in available:
            resolution = default
        else:
            # nearest at-or-coarser-than the default, else the coarsest there is: an
            # analysis that asks for 100 kb wants 100 kb-ish, not the finest level
            at_or_above = [r for r in available if r >= default]
            resolution = at_or_above[0] if at_or_above else available[-1]
    if resolution not in available:
        raise DataError(
            f"Resolution {resolution} bp is not in {path.name}. "
            f"Available: {', '.join(str(r) for r in available)} bp."
        )
    uri = by_res[resolution]
    return cooler.Cooler(str(path) if uri == "/" else f"{path}::{uri}")


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
    if int(end) <= int(start):
        raise DataError(
            f"Region {region!r} has zero width (start {int(start):,} is not before end "
            f"{int(end):,}). Ask for a span of at least one bin."
        )
    return str(chrom), int(start), int(end)


def load_gc_track() -> pd.DataFrame:
    """The bundled GC phasing track for the demo file (chrom/start/end/GC at 100 kb)."""
    p = data_dir() / GC_TRACK_FILENAME
    if not p.exists():
        raise DataError(f"Bundled GC track missing (expected at {p}).")
    return pd.read_csv(p, sep="\t")


def _header_kwargs(p: Path) -> dict:
    """Decide how to parse a BED-like file: header row or not, comments or not.

    Detected by SHAPE, not by a literal prefix - a header is any first row whose
    start/end fields are not numeric. A column called "Chromosome" is still a header.
    A '#chrom' line is BED convention, but treating '#' as a comment marker at the same
    time blanks it and silently eats the first real data row, so the two never combine.
    """
    lines = [ln for ln in p.read_text(errors="ignore").splitlines() if ln.strip()]
    # a '#' line with fewer than three tab fields is prose, not a header row
    skip = 0
    for ln in lines:
        if ln.lstrip().startswith("#") and len(ln.split("\t")) < 3:
            skip += 1
        else:
            break
    if skip >= len(lines):
        return {"sep": "\t", "header": None, "comment": "#", "skiprows": skip}
    fields = lines[skip].split("\t")
    looks_like_header = True
    if len(fields) >= 3:
        try:
            int(float(fields[1]))
            int(float(fields[2]))
            looks_like_header = False
        except ValueError:
            looks_like_header = True
    if looks_like_header:
        # the '#' is part of the header row, so it must not also be a comment marker
        return {"sep": "\t", "header": 0, "skiprows": skip}
    return {"sep": "\t", "header": None, "comment": "#", "skiprows": skip}


def load_view_file(path: str) -> pd.DataFrame:
    """A user-supplied region view: BED-like, chrom/start/end[/name], tab-separated."""
    p = Path(path).expanduser()
    if not p.exists():
        raise DataError(f"No such view file: {p}. Expected a tab-separated BED-like file.")
    if not p.is_file():
        raise DataError(f"{p} is a directory, not a view file. Pass the BED-like file itself.")
    kwargs = _header_kwargs(p)
    try:
        df = pd.read_csv(p, **kwargs)
    except Exception as e:  # noqa: BLE001 - the message is what the agent acts on
        raise DataError(f"Could not read view file {p.name}: {e}") from e
    if df.empty:
        raise DataError(
            f"View file {p.name} has a header but no regions. Add at least one "
            "tab-separated chrom/start/end row."
        )
    if df.shape[1] < 3:
        raise DataError(
            f"View file {p.name} has {df.shape[1]} column(s); expected at least "
            "chrom, start, end (tab-separated)."
        )
    df = df.iloc[:, :4] if df.shape[1] >= 4 else df.iloc[:, :3]
    df.columns = ["chrom", "start", "end", "name"][: df.shape[1]]
    if "name" not in df.columns:
        df["name"] = df["chrom"].astype(str) + ":" + df["start"].astype(str)
    df["chrom"] = df["chrom"].astype(str)
    for col in ("start", "end"):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if df[["start", "end"]].isna().any().any():
        raise DataError(
            f"View file {p.name} has non-numeric start/end values. Expected tab-separated "
            "chrom, start, end (an optional name column is allowed, and a header row is fine)."
        )
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def load_track_file(path: str) -> pd.DataFrame:
    """A user-supplied phasing track: chrom/start/end/value, tab-separated.

    The fourth column orients the compartment eigenvector (higher = A by convention),
    e.g. GC fraction or gene density.
    """
    p = Path(path).expanduser()
    if not p.exists():
        raise DataError(f"No such phasing track: {p}. Expected a tab-separated bedGraph.")
    if not p.is_file():
        raise DataError(f"{p} is a directory, not a phasing track. Pass the file itself.")
    kwargs = _header_kwargs(p)
    try:
        df = pd.read_csv(p, **kwargs)
    except Exception as e:  # noqa: BLE001
        raise DataError(f"Could not read phasing track {p.name}: {e}") from e
    if df.empty:
        raise DataError(f"Phasing track {p.name} has a header but no rows.")
    if df.shape[1] < 4:
        raise DataError(
            f"Phasing track {p.name} has {df.shape[1]} column(s); expected chrom, start, "
            "end and a value column (e.g. GC fraction)."
        )
    df = df.iloc[:, :4]
    value_name = str(df.columns[3]) if kwargs.get("header") == 0 else "value"
    df.columns = ["chrom", "start", "end", value_name]
    df["chrom"] = df["chrom"].astype(str)
    for col in ("start", "end", value_name):
        df[col] = pd.to_numeric(df[col], errors="coerce")
    if df[["start", "end"]].isna().any().any():
        raise DataError(
            f"Phasing track {p.name} has non-numeric start/end values. Expected "
            "tab-separated chrom, start, end and a numeric value column."
        )
    if df[value_name].isna().all():
        raise DataError(
            f"Phasing track {p.name}: its fourth column ({value_name!r}) holds no numbers. "
            "A phasing track needs a numeric value per bin (GC fraction, gene density...) - "
            "a BED file whose fourth column is a name or '.' is a region list, not a track."
        )
    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    return df


def load_arms_view() -> pd.DataFrame:
    """The bundled chr17 p/q arm view for the demo file."""
    p = data_dir() / ARMS_FILENAME
    if not p.exists():
        raise DataError(f"Bundled arm view missing (expected at {p}).")
    return pd.read_csv(p, sep="\t", names=["chrom", "start", "end", "name"])
