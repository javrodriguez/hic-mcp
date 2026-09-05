"""Input FILES a human would never construct, driven through the tools that read them.

`test_awkward_coordinates.py` generalised one class - a road the test corpus never walks -
along the axis it had just been bitten on: coordinates. The next round found the same class
one axis over. A phasing track shifted 50 kb off the bin grid, or mixing two interval widths,
reached cooltools untouched and came back to the agent as:

    "This is a bug in hic-mcp, not in your request - try a different region or resolution"

which was wrong three times over: not a bug, it IS the request, and the remedy named cannot
fix it. Neither track is exotic - a bedGraph from another pipeline routinely starts at a
different offset, and a track lifted from a 50 kb annotation is the ordinary case.

So the class is not "coordinates". It is **any input the author would not naturally build**,
and this file covers the other axis: the files the tools accept. The property is the same one
the coordinate sweep asserts - a result, or a refusal that names what is wrong with the input
and what to do about it, but never a bare library error dressed as an internal defect.
"""

import numpy as np
import pandas as pd
import pytest

from hic_mcp.analysis import AnalysisError, compartments
from hic_mcp.data import DataError, data_dir


def _bundled() -> pd.DataFrame:
    return pd.read_csv(data_dir() / "gc_100kb.tsv", sep="\t")


def _write(df: pd.DataFrame, path) -> str:
    df.to_csv(path, sep="\t", index=False)
    return str(path)


def _mangle(kind: str, df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    col = df.columns[3]
    if kind == "off_grid":                       # right width, wrong offset
        df["start"] += 50_000
        df["end"] = df["start"] + 100_000
        return df[(df["end"] - df["start"]) == 100_000]
    if kind == "mixed_widths":                   # two binnings in one file
        df.loc[df.index[5:], "end"] = df.loc[df.index[5:], "start"] + 50_000
        return df
    if kind == "wrong_binning":                  # uniform, but not this matrix's
        df["end"] = df["start"] + 250_000
        return df
    if kind == "all_nan":                        # structurally fine, no values
        df[col] = np.nan
        return df
    if kind == "constant":                       # no variance to correlate against
        df[col] = 0.42
        return df
    if kind == "one_row":
        return df.iloc[:1]
    if kind == "unknown_chrom":
        df["chrom"] = "chrZZ"
        return df
    raise AssertionError(kind)


AWKWARD_TRACKS = [
    "off_grid", "mixed_widths", "wrong_binning", "all_nan",
    "constant", "one_row", "unknown_chrom",
]


@pytest.mark.parametrize("kind", AWKWARD_TRACKS)
def test_a_malformed_track_is_the_callers_problem_stated_plainly(kind, tmp_path):
    """Never a bare library error, and never blamed on this server."""
    path = _write(_mangle(kind, _bundled()), tmp_path / f"{kind}.tsv")
    try:
        out = compartments(region="chr17:50,100,000-51,100,000", resolution=100_000,
                           phasing_track=path)
    except (AnalysisError, DataError) as e:
        msg = str(e)
        assert msg.strip(), kind
        assert "bug in hic-mcp" not in msg, kind
        # a refusal has to leave the caller somewhere to go
        assert any(w in msg.lower() for w in
                   ("pass", "supply", "re-bin", "extend", "narrow", "check", "expected",
                    "needs", "omit")), \
            f"{kind}: refusal names no way forward -> {msg}"
        return
    except Exception as e:  # noqa: BLE001
        pytest.fail(f"{kind} leaked {type(e).__name__}: {e}")
    # if it did answer, it must not have silently claimed an orientation it cannot justify
    assert "sign_convention" in out, kind


def test_the_off_grid_track_is_diagnosed_as_off_grid(tmp_path):
    """The specific diagnosis, not merely a readable one.

    This track keeps 100 kb widths and full chromosome coverage, so every guard that existed
    before passed it: only the offset is wrong, and only the offset explains the failure.
    """
    path = _write(_mangle("off_grid", _bundled()), tmp_path / "off.tsv")
    with pytest.raises(AnalysisError, match="do not start on this matrix's bin grid"):
        compartments(region="chr17:50,100,000-51,100,000", resolution=100_000,
                     phasing_track=path)


def test_a_final_partial_interval_is_normal_and_not_refused():
    """The bundled track's last chr17 bin is 57,441 bp; a width check must allow that.

    Written because the first version of the grid check refused the repo's own demo data -
    a guard strict enough to reject the file it ships with is a new defect, not a fix.
    """
    widths = (_bundled()["end"] - _bundled()["start"]).unique()
    assert len(widths) > 1, "the bundled track no longer has a partial final bin"
    out = compartments(region="chr17:50,100,000-51,100,000", resolution=100_000)
    assert out["region_call"] in {"A", "B"}


def test_citation_version_matches_the_package():
    """CITATION.cff asks to be cited; a citation naming the wrong version is worse than none."""
    import tomllib
    from pathlib import Path

    repo = Path(__file__).resolve().parents[1]
    pyproject = tomllib.loads((repo / "pyproject.toml").read_text())
    cff = (repo / "CITATION.cff").read_text()
    version = pyproject["project"]["version"]
    assert f'version: "{version}"' in cff, f"CITATION.cff does not name version {version}"
    assert "date-released:" in cff
