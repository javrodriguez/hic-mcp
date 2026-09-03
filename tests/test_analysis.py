"""Ground-truth tests: every tool asserts a known, measured result on the demo data.

Loci and values come from data/PROVENANCE.md (measured on the shipped file).
Numeric assertions carry tolerances because CI runs a different BLAS/platform
than the machine the values were measured on; integer raw-count sums are exact.
"""

import numpy as np
import pandas as pd
import pytest

from hic_mcp.analysis import (
    AnalysisError,
    compartments,
    contacts_at_locus,
    expected_observed,
    insulation_tads,
    matrix_summary,
    virtual_4c,
)

BOUNDARY_LOCUS = "chr17:66,180,000-66,190,000"
A_BLOCK = "chr17:50,100,000-51,100,000"
B_BLOCK = "chr17:51,400,000-52,400,000"


# --- matrix_summary --------------------------------------------------------


def test_matrix_summary_ground_truth():
    out = matrix_summary()
    assert out["is_bundled_demo"] is True
    assert out["assembly"] == "hg38"
    assert list(out["chromosomes"]) == ["chr17"]
    assert out["chromosomes"]["chr17"] == 83_257_441
    res = {r["resolution_bp"]: r for r in out["resolutions"]}
    assert set(res) == {10_000, 100_000, 1_000_000}
    for r in res.values():
        assert r["total_contacts"] == 129_434_772  # exact: integer sum
        assert r["balanced"] is True
    assert out["provenance"]["source_4dn_accession"] == "4DNESWST3UBH"


# --- contacts_at_locus -----------------------------------------------------


def test_contacts_ground_truth_exact_raw_sum():
    out = contacts_at_locus(region="chr17:45,000,000-46,000,000", resolution=10_000)
    assert out["raw_contacts_sum"] == 2_323_926  # exact: measured integer sum
    assert out["raw_contacts_max"] == 28_614
    assert out["resolution_used"] == 10_000
    assert out["balanced"] is True
    assert out["balanced_matrix"] is None  # 100 bins > cap -> stats only
    assert "cap" in out["note"]


def test_contacts_near_diagonal_enrichment():
    near = contacts_at_locus(
        region="chr17:50,000,000-50,500,000", resolution=10_000
    )["raw_contacts_sum"]
    far = contacts_at_locus(
        region="chr17:50,000,000-50,500,000",
        region2="chr17:60,000,000-60,500,000",
        resolution=10_000,
    )["raw_contacts_sum"]
    assert near > 5 * far


def test_contacts_small_window_returns_matrix_and_autopicks_resolution():
    out = contacts_at_locus(region="chr17:50,000,000-52,000,000")  # 2 Mb, no resolution
    assert out["resolution_used"] == 100_000  # finest fitting the 50-bin cap
    m = out["balanced_matrix"]
    assert m is not None and len(m) == 20
    flat = [v for row in m for v in row if v is not None]
    assert all(v >= 0 for v in flat) and len(flat) > 0


def test_contacts_centromeric_region_is_nearly_empty():
    """ICE-filtered bins keep a handful of raw reads; balanced values are all NaN."""
    out = contacts_at_locus(region="chr17:23,000,000-23,500,000", resolution=10_000)
    assert out["raw_contacts_sum"] < 100  # measured: 17 stray reads in 0.5 Mb
    assert out["balanced_mean"] is None  # every bin here is ICE-filtered


# --- insulation_tads -------------------------------------------------------


def test_insulation_strongest_boundary_ground_truth():
    out = insulation_tads()
    assert out["resolution_used"] == 10_000
    assert out["ranked_by"] == "boundary_strength at the 250000 bp window"
    top = out["top_boundaries"][0]
    assert top["locus"] == BOUNDARY_LOCUS  # measured: strongest at the TAD-scale window
    assert top["strength"] == pytest.approx(2.7629, rel=0.05)
    assert top["windows_detected"] == [100_000, 250_000, 500_000]
    assert out["boundary_counts_per_window"]["100000"] > 100


def test_insulation_region_filter_and_window_guard():
    out = insulation_tads(region="chr17:65,000,000-67,000,000", top_n=3)
    assert any(b["locus"] == BOUNDARY_LOCUS for b in out["top_boundaries"])
    with pytest.raises(AnalysisError, match="at least 3 bins"):
        insulation_tads(windows_bp=[20_000])


def test_insulation_runs_at_every_advertised_resolution():
    """Every resolution is usable, and TAD-scale windows are kept wherever bins allow."""
    expected_windows = {
        10_000: [100_000, 250_000, 500_000],  # classic TAD scale
        100_000: [300_000, 500_000, 1_000_000],  # still TAD scale: 3 bins is 300 kb
        1_000_000: [3_000_000, 5_000_000, 10_000_000],  # too coarse - must warn
    }
    for res, windows in expected_windows.items():
        out = insulation_tads(resolution=res, top_n=1)
        assert out["windows_bp"] == windows, res
        assert out["top_boundaries"], res
        assert (out.get("scale_note") is not None) == (min(windows) > 500_000), res


def test_insulation_at_100kb_stays_at_tad_scale():
    """At 100 kb bins the windows must stay TAD-scale, not drift to compartment scale.

    Ranking is deliberately not asserted here: chr17:51.2 Mb (the documented A/B flip)
    is a genuinely strong insulation feature and competes for the top slot within 0.3%,
    so pinning an order would be a brittle test of real biology rather than of the code.
    """
    out = insulation_tads(resolution=100_000, top_n=3)
    assert max(out["windows_bp"]) <= 1_000_000
    assert out.get("scale_note") is None  # TAD scale is achievable here, so no warning
    coarse = insulation_tads(resolution=1_000_000, top_n=1)
    assert coarse["scale_note"] is not None  # but not here


def test_insulation_needs_weights(tiny_unbalanced_cool):
    with pytest.raises(AnalysisError, match="no ICE weights"):
        insulation_tads(file=tiny_unbalanced_cool)


# --- compartments ----------------------------------------------------------


def test_compartment_signs_ground_truth():
    a = compartments(region=A_BLOCK)
    b = compartments(region=B_BLOCK)
    assert a["region_call"] == "A" and a["region_mean_E1"] > 0
    assert b["region_call"] == "B" and b["region_mean_E1"] < 0
    assert a["region_mean_E1"] == pytest.approx(1.18, rel=0.15)
    assert b["region_mean_E1"] == pytest.approx(-1.17, rel=0.15)
    assert a["region_sign_consistency"] == 1.0
    assert b["region_sign_consistency"] == 1.0
    assert "GC" in a["sign_convention"]


def test_compartment_eigenvector_correlates_with_gc():
    """Data-contract guard: the bundled phasing track really orients this file's E1."""
    from cooltools import eigs_cis

    from hic_mcp.data import load_arms_view, load_gc_track, open_matrix, resolve_input_path

    clr = open_matrix(resolve_input_path(None), 100_000, default=100_000)
    gc = load_gc_track()
    _, vecs = eigs_cis(clr, phasing_track=gc, view_df=load_arms_view(), n_eigs=1)
    merged = pd.merge(vecs, gc, on=["chrom", "start", "end"]).dropna(subset=["E1", "GC"])
    r = np.corrcoef(merged["E1"], merged["GC"])[0, 1]
    assert r > 0.3  # measured ~ +0.48


def test_compartments_centromere_region_is_a_clear_error():
    with pytest.raises(AnalysisError, match="ICE-filtered"):
        compartments(region="chr17:23,000,000-26,000,000")


def test_compartments_genome_summary():
    out = compartments()
    assert 0.2 < out["genome_A_fraction"] < 0.8
    assert out["bins_used"] > 500


def test_compartments_never_claims_an_A_fraction_when_unphased():
    """An arbitrary sign cannot support an 'A fraction' - the field must stay null."""
    out = compartments(resolution=1_000_000)  # phasing track exists only at 100 kb
    assert "UNPHASED" in out["sign_convention"]
    assert out["genome_A_fraction"] is None
    assert 0.0 < out["positive_E1_fraction"] < 1.0


# --- virtual_4c ------------------------------------------------------------


def test_virtual4c_ground_truth_shape():
    out = virtual_4c(viewpoint="chr17:63,000,000-63,100,000")  # default 10 kb
    bands = out["distance_band_means"]
    assert set(bands) == {"0kb-100kb", "100kb-1Mb", "1Mb-5Mb", "5Mb-10Mb"}
    vals = list(bands.values())
    assert all(v is not None for v in vals)
    assert vals == sorted(vals, reverse=True)  # decaying distance-band means
    assert all(p["balanced"] is not None for p in out["profile_points"])
    assert "genuine zeros" in out["profile_note"]
    # every band must average over real bins, and mappable zero-contact bins count
    assert out["distance_band_bins"]["5Mb-10Mb"] == 1000
    assert out["distance_band_means"]["5Mb-10Mb"] == pytest.approx(0.000146, rel=0.05)


def test_virtual4c_counts_genuine_zeros_not_just_contacted_bins():
    """Dropping mappable zero-contact bins would flatten the decay this tool describes."""
    out = virtual_4c(viewpoint="chr17:63,000,000-63,100,000")
    far = out["distance_band_means"]["5Mb-10Mb"]
    near = out["distance_band_means"]["1Mb-5Mb"]
    # with zeros dropped these collapse together; with zeros counted they stay apart
    assert near > 1.5 * far


def test_virtual4c_coarse_resolution_drops_empty_band():
    out = virtual_4c(viewpoint="chr17:63,000,000-63,100,000", resolution=100_000)
    # only the viewpoint's own (deliberately masked) bin falls in that band here
    assert "0kb-100kb" not in out["distance_band_means"]
    vals = list(out["distance_band_means"].values())
    assert vals == sorted(vals, reverse=True)


def test_virtual4c_viewpoint_at_the_chromosome_end():
    """Padding a sub-bin viewpoint must not push the request past the chromosome.

    chr17's final bins are ICE-filtered (telomeric), so the honest answer there is the
    clear filtered-bin message - never an out-of-bounds crash from the padding itself.
    """
    with pytest.raises(AnalysisError, match="ICE-filtered"):
        virtual_4c(viewpoint="chr17:83,250,000-83,257,441", resolution=10_000)
    # the last usable bin, one bin short of the end, still answers
    out = virtual_4c(viewpoint="chr17:83,200,000-83,210,000", resolution=10_000)
    assert out["profile_points"]


def test_virtual4c_filtered_viewpoint_is_a_clear_error():
    with pytest.raises(AnalysisError, match="ICE-filtered"):
        virtual_4c(viewpoint="chr17:45,500,000-45,510,000", resolution=10_000)
    with pytest.raises(AnalysisError, match="ICE-filtered"):
        virtual_4c(viewpoint="chr17:24,000,000-24,010,000", resolution=10_000)


# --- expected_observed -----------------------------------------------------


def test_expected_ps_slope_ground_truth():
    """The slope is fitted only over diagonals cooltools actually measured."""
    out = expected_observed(region=A_BLOCK)  # default 100 kb
    assert out["ignored_diagonals"] == 2
    assert out["ps_fit_range_bp"] == [200_000, 10_000_000]  # starts past the masked head
    assert out["ps_slope"] == pytest.approx(-1.2254, rel=0.05)
    assert len(out["expected_curve_points"]) > 20
    # the masked head must not appear in the reported curve
    assert min(p["dist_bp"] for p in out["expected_curve_points"]) >= 200_000


def test_expected_ps_slope_is_negative_at_every_resolution():
    """Contact frequency falls with distance; a positive slope means a contaminated fit."""
    for res, expected in ((10_000, -1.2416), (100_000, -1.2254), (1_000_000, -1.0547)):
        out = expected_observed(region=A_BLOCK, resolution=res)
        assert out["ps_slope"] == pytest.approx(expected, rel=0.08), res
        assert out["ps_fit_range_bp"][0] >= 2 * res


def test_expected_oe_matrix_is_centred_on_one_per_diagonal():
    """O/E is ~1 by construction: assert it diagonal by diagonal, not on a pooled median."""
    out = expected_observed(region="chr17:50,000,000-52,500,000")
    m = out["oe_matrix"]
    assert m is not None and len(m) == 25
    arr = np.array([[np.nan if v is None else v for v in row] for row in m])
    # the unmeasured head is null, not an invented ratio
    for d in range(out["ignored_diagonals"]):
        assert np.isnan(np.diagonal(arr, offset=d)).all(), f"diagonal {d} should be null"
    assert "null by construction" in out["note"]
    for d in range(out["ignored_diagonals"], 6):
        vals = np.diagonal(arr, offset=d)
        vals = vals[np.isfinite(vals)]
        assert vals.size > 5
        assert 0.3 < float(np.median(vals)) < 3.0, f"diagonal {d} median off 1"
    assert np.nanmax(arr) < 10.0  # no order-of-magnitude artifact anywhere


def test_expected_large_region_stats_only():
    out = expected_observed(region="chr17:27,100,000-83,257,441", resolution=100_000)
    assert out["oe_matrix"] is None and "cap" in out["note"]


def test_expected_curve_is_scoped_to_the_region_not_the_whole_file():
    """A slope labelled with the caller's region must describe that region's arm."""
    p_arm = expected_observed(region="chr17:5,000,000-7,500,000")
    q_arm = expected_observed(region="chr17:50,000,000-52,500,000")
    assert p_arm["view"].startswith("chr17p")
    assert q_arm["view"].startswith("chr17q")
    assert p_arm["ps_slope"] != q_arm["ps_slope"]
    assert p_arm["ps_slope"] == pytest.approx(-1.1625, rel=0.05)
    assert q_arm["ps_slope"] == pytest.approx(-1.2407, rel=0.05)
    assert "not only the requested region" in p_arm["curve_scope"]


def test_expected_curve_tail_is_not_flattened_by_rounding():
    """Six-decimal rounding collapsed a decade of decay into one repeated constant."""
    pts = expected_observed(region=A_BLOCK, resolution=10_000)["expected_curve_points"]
    tail = [p["expected"] for p in pts[-10:] if p["expected"] is not None]
    assert len(set(tail)) == len(tail), f"curve tail is flat: {tail}"
    assert tail == sorted(tail, reverse=True) or tail[0] > tail[-1]


def test_expected_region_straddling_the_centromere_is_a_clear_error():
    with pytest.raises(AnalysisError, match="not contained in a single chromosome arm"):
        expected_observed(region="chr17:20,000,000-25,000,000")


# --- every tool names its method (honest-output invariant) ------------------


def test_every_tool_names_method_resolution_and_balancing():
    outs = [
        matrix_summary(),
        contacts_at_locus(region="chr17:50,000,000-50,500,000", resolution=100_000),
        insulation_tads(region="chr17:65,000,000-67,000,000"),
        compartments(region=A_BLOCK),
        virtual_4c(viewpoint="chr17:63,000,000-63,100,000", resolution=100_000),
        expected_observed(region=A_BLOCK),
    ]
    for out in outs:
        assert "method" in out and ("cooler" in out["method"] or "cooltools" in out["method"])
        assert "balanced" in out
    assert all("resolution_used" in o for o in outs[1:])


# --- the bring-your-own-file path (two chromosomes, not the demo's one) -------


def test_every_tool_works_on_a_multi_chromosome_file(two_chromosome_cool):
    """A single-chromosome demo cannot catch single-chromosome assumptions."""
    f = two_chromosome_cool
    summary = matrix_summary(file=f)
    assert set(summary["chromosomes"]) == {"cA", "cB"}
    assert contacts_at_locus(file=f, region="cA:0-200,000")["raw_contacts_sum"] > 0
    assert insulation_tads(file=f, region="cA:0-500,000", top_n=1)["resolution_used"] == 10_000
    v = virtual_4c(file=f, viewpoint="cA:200,000-210,000")
    assert v["distance_band_means"], "virtual_4c returned no bands"
    assert all(x is not None for x in v["distance_band_means"].values())


def test_virtual4c_profile_stays_on_the_viewpoint_chromosome(two_chromosome_cool):
    """The profile, its weights and its masks must all describe the same rows."""
    v = virtual_4c(file=two_chromosome_cool, viewpoint="cA:200,000-210,000")
    starts = [p["start"] for p in v["profile_points"]]
    assert max(starts) < 60 * 10_000  # cA is 60 bins; cB values would exceed this


def test_insulation_refuses_an_unmeasurable_region_rather_than_reporting_zero():
    """"0 boundaries" is an answer; a fully ICE-filtered region has no measurement.

    The sibling tools already refuse on this exact region, so answering here would be
    the one inconsistent - and quietly misleading - surface.
    """
    with pytest.raises(AnalysisError, match="no measurement"):
        insulation_tads(region="chr17:23,000,000-25,000,000")


def test_compartment_view_does_not_change_with_resolution():
    """The arm view is resolution-independent; only the GC phasing track is 100 kb-bound."""
    for res in (10_000, 100_000, 1_000_000):
        out = compartments(region="chr17:50,100,000-51,100,000", resolution=res)
        assert "arms" in out["view"], res
    # phasing, and therefore the A/B call, is only claimed where the track applies
    assert compartments(region=A_BLOCK, resolution=100_000)["region_call"] == "A"
    assert compartments(region=A_BLOCK, resolution=10_000)["region_call"] == "unphased"


# --- round-4 regressions: fix classes, not instances ---------------------------


def test_insulation_uses_the_same_arm_view_as_its_siblings():
    """Without the view cooltools normalises across the centromeric gap the other tools
    refuse to compute across, shifting every score by a per-arm offset."""
    out = insulation_tads(resolution=100_000, top_n=1)
    assert "arms" in out["view"]
    # and the 10 kb ground truth survives the view (measured: it does)
    top = insulation_tads()["top_boundaries"][0]
    assert top["locus"] == BOUNDARY_LOCUS
    assert top["strength"] == pytest.approx(2.7629, rel=0.05)


def test_virtual4c_refuses_a_region_sized_viewpoint():
    """A viewpoint is an anchor; a whole chromosome must be refused, never profiled as empty."""
    with pytest.raises(AnalysisError, match="at most"):
        virtual_4c(viewpoint="chr17")
    with pytest.raises(AnalysisError, match="at most"):
        virtual_4c(viewpoint="chr17:63,000,000-83,000,000")


def test_virtual4c_never_returns_an_empty_profile_silently():
    """Every successful call carries at least one measured point."""
    out = virtual_4c(viewpoint="chr17:63,000,000-63,100,000")
    assert out["profile_points"]
    assert "every 1th" not in out["profile_note"]
