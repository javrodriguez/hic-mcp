"""Coordinates a human would never type, driven through every tool.

Three consecutive evaluation rounds found the same SHAPE of defect: a road the test corpus
never walks. Round 2 found region ORDER (descending region2 returned a false zero), round 3
found region OVERLAP (16% undercount), round 4 found bin-grid ALIGNMENT (shape_bins
contradicting the matrix beside it, and virtual_4c raising a bare ValueError for any anchor
not on a bin edge). Each was fixed as an instance; each time the next round found the next
instance, because every coordinate in this repo's tests, README, examples and demo is a round
multiple of the bin size, ascending and disjoint.

So this file does not add a fourth special case. It generates the awkward cases - unaligned,
sub-bin, single-bp, descending, overlapping, chromosome-edge - and drives every tool with
them, asserting the properties that must hold whatever the caller typed:

  1. no bare library exception ever escapes (the agent must never be told its own request is
     a bug in this server);
  2. a reported bin count describes the object actually returned;
  3. a fraction stays a fraction;
  4. what the response SAYS was analysed is what was analysed.

A real genomic anchor - a TSS, an enhancer, a coordinate pasted from a genome browser - is a
multiple of 10,000 essentially never. These are the ordinary inputs, not the exotic ones.
"""

import numpy as np
import pytest

from hic_mcp.analysis import (
    AnalysisError,
    compartments,
    contacts_at_locus,
    expected_observed,
    insulation_tads,
    virtual_4c,
)
from hic_mcp.data import DataError, demo_path, open_matrix

# unaligned starts, unaligned ends, sub-bin spans, 1 bp, and the chromosome's last bin
AWKWARD = [
    "chr17:50,005,000-50,205,000",   # unaligned start and end, 10 kb grid
    "chr17:50,000,001-50,200,000",   # 1 bp off the grid
    "chr17:50,005,000-50,006,000",   # sub-bin, wholly inside one bin
    "chr17:50,000,000-50,000,001",   # 1 bp
    "chr17:50,005,000-50,505,000",   # unaligned, wide enough to cross the sparse-road cap
    "chr17:83,250,001-83,257,441",   # the chromosome's final, partial bin
]


def _extent(region: str) -> int:
    clr = open_matrix(demo_path(), 10_000, default=10_000)
    lo, hi = clr.extent(region.replace(",", ""))
    return int(hi - lo)


@pytest.mark.parametrize("region", AWKWARD)
def test_contacts_reports_the_object_it_returns(region):
    """shape_bins, the matrix, and cooler's own extent must describe one thing."""
    out = contacts_at_locus(region=region, resolution=10_000)
    assert out["shape_bins"][0] == _extent(region), region
    matrix = out.get("balanced_matrix")
    if matrix is not None:
        assert len(matrix) == out["shape_bins"][0], region
        assert len(matrix[0]) == out["shape_bins"][1], region
    assert 0.0 <= out["nonzero_fraction"] <= 1.0, f"{region} -> {out['nonzero_fraction']}"


@pytest.mark.parametrize("region", AWKWARD)
def test_the_sparse_road_agrees_on_awkward_coordinates_too(region):
    """The parity that round 3 established, on coordinates round 3 never tried."""
    import hic_mcp.analysis as an

    dense = contacts_at_locus(region=region, resolution=10_000)
    saved = an.DENSE_FETCH_CAP_GB
    try:
        an.DENSE_FETCH_CAP_GB = 0.0
        sparse = contacts_at_locus(region=region, resolution=10_000)
    finally:
        an.DENSE_FETCH_CAP_GB = saved
    for key in ("raw_contacts_sum", "raw_contacts_max", "nonzero_fraction",
                "balanced_mean", "balanced_max"):
        d, sp = dense.get(key), sparse.get(key)
        if isinstance(d, float) and isinstance(sp, float):
            assert sp == pytest.approx(d, rel=1e-6), f"{region} {key}"
        else:
            assert sp == d, f"{region} {key}"
    assert 0.0 <= sparse["nonzero_fraction"] <= 1.0, region


@pytest.mark.parametrize("region2", ["chr17:50,003,000-50,203,000", "chr17:49,995,000-50,105,000"])
def test_two_awkward_regions_together(region2):
    """Unaligned AND overlapping AND, in the second case, starting before region 1."""
    r1 = "chr17:50,005,000-50,205,000"
    out = contacts_at_locus(region=r1, region2=region2, resolution=10_000)
    assert out["shape_bins"] == [_extent(r1), _extent(region2)]
    assert 0.0 <= out["nonzero_fraction"] <= 1.0


# a virtual-4C anchor is a small locus by contract, so the awkward set for it is the
# narrow members plus off-grid anchors of its own - the wide ones belong to the refusal
# test below, where an over-wide anchor is correct product behaviour, not a defect
AWKWARD_ANCHORS = [
    "chr17:50,005,000-50,006,000",
    "chr17:50,000,000-50,000,001",
    "chr17:63,000,001-63,000,002",
    "chr17:63,005,000-63,025,000",
    "chr17:83,250,001-83,257,441",
]


@pytest.mark.parametrize("viewpoint", AWKWARD_ANCHORS)
def test_virtual_4c_takes_any_anchor_and_says_what_it_used(viewpoint):
    """An anchor off the bin grid is ordinary input, not a server bug."""
    try:
        out = virtual_4c(viewpoint=viewpoint)
    except AnalysisError as e:
        # a refusal is a legitimate answer (chr17's final bin is ICE-filtered); a CRASH is
        # not, and the message must be about the data, never about this server
        assert "bug in hic-mcp" not in str(e), viewpoint
        assert str(e).strip(), viewpoint
        return
    assert out["profile_points"], viewpoint
    vp_start = int(out["viewpoint"].split(":")[1].split("-")[0].replace(",", ""))
    assert vp_start % 10_000 == 0, f"{viewpoint} -> {out['viewpoint']} is not on the grid"
    raw_start = int(viewpoint.split(":")[1].split("-")[0].replace(",", ""))
    if raw_start % 10_000:
        # the caller must be told the anchor moved, not silently handed a different locus
        assert "does not sit on" in out["profile_note"], viewpoint


@pytest.mark.parametrize("region", AWKWARD)
def test_expected_observed_caps_on_the_matrix_it_returns(region):
    """The bin cap must count the same bins the response carries."""
    out = expected_observed(region=region, resolution=100_000)
    oe = out.get("oe_matrix")
    if oe is not None:
        assert len(oe) <= 50, f"{region} returned {len(oe)} rows past the cap"


@pytest.mark.parametrize("region", AWKWARD)
@pytest.mark.parametrize("tool", [insulation_tads, compartments])
def test_no_awkward_coordinate_is_reported_as_a_server_bug(tool, region):
    """Whatever a caller types, the answer is a result or an agent-readable refusal."""
    try:
        tool(region=region, resolution=100_000 if tool is compartments else 10_000)
    except (AnalysisError, DataError) as e:
        assert str(e).strip(), region
        assert "bug in hic-mcp" not in str(e), region
    except Exception as e:  # noqa: BLE001
        pytest.fail(f"{tool.__name__}({region!r}) leaked {type(e).__name__}: {e}")


def test_the_awkward_set_is_actually_awkward():
    """A guard that quietly stopped being adversarial would pass forever.

    Every region above must be off the 10 kb grid at one end or the other; if someone
    'tidies' these coordinates into round numbers this file stops testing anything, exactly
    the way the aligned corpus it exists to compensate for did.
    """
    off_grid = [r for r in AWKWARD
                if int(r.split(":")[1].split("-")[0].replace(",", "")) % 10_000
                or int(r.split("-")[1].replace(",", "")) % 10_000]
    assert len(off_grid) == len(AWKWARD), f"these are on the grid: {set(AWKWARD) - set(off_grid)}"
    assert np.all([r.startswith("chr17:") for r in AWKWARD])


@pytest.mark.parametrize(
    "viewpoint", ["chr17:50,005,000-50,205,000", "chr17:50,000,001-50,200,000"]
)
def test_an_over_wide_anchor_is_refused_readably_not_crashed(viewpoint):
    """The width cap is the contract; an unaligned request must still meet it by refusal."""
    with pytest.raises(AnalysisError, match="viewpoint is a small anchor"):
        virtual_4c(viewpoint=viewpoint)
