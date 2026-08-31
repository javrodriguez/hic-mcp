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
    with pytest.raises(AnalysisError, match="too small"):
        insulation_tads(windows_bp=[20_000])


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


# --- virtual_4c ------------------------------------------------------------


def test_virtual4c_ground_truth_shape():
    out = virtual_4c(viewpoint="chr17:63,000,000-63,100,000")  # default 10 kb
    bands = out["distance_band_means"]
    assert set(bands) == {"0-100kb", "100kb-1Mb", "1000kb-5Mb", "5000kb-10Mb"}
    vals = list(bands.values())
    assert all(v is not None for v in vals)
    assert vals == sorted(vals, reverse=True)  # decaying distance-band means
    assert all(p["balanced"] is not None for p in out["profile_points"])
    assert "NaN by construction" in out["profile_note"]


def test_virtual4c_coarse_resolution_drops_empty_band():
    out = virtual_4c(viewpoint="chr17:63,000,000-63,100,000", resolution=100_000)
    assert "0-100kb" not in out["distance_band_means"]  # narrower than one bin
    vals = list(out["distance_band_means"].values())
    assert vals == sorted(vals, reverse=True)


def test_virtual4c_filtered_viewpoint_is_a_clear_error():
    with pytest.raises(AnalysisError, match="ICE-filtered"):
        virtual_4c(viewpoint="chr17:45,500,000-45,510,000", resolution=10_000)
    with pytest.raises(AnalysisError, match="ICE-filtered"):
        virtual_4c(viewpoint="chr17:24,000,000-24,010,000", resolution=10_000)


# --- expected_observed -----------------------------------------------------


def test_expected_ps_slope_ground_truth():
    out = expected_observed(region=A_BLOCK)
    assert out["ps_slope_100kb_10Mb"] is not None
    assert -1.6 <= out["ps_slope_100kb_10Mb"] <= -0.9  # measured ~ -1.245
    assert len(out["expected_curve_points"]) > 20


def test_expected_oe_matrix_for_small_region():
    out = expected_observed(region="chr17:50,000,000-52,500,000")
    m = out["oe_matrix"]
    assert m is not None and len(m) == 25
    flat = [v for row in m for v in row if v is not None]
    assert len(flat) > 100 and all(v >= 0 for v in flat)
    med = float(np.median(flat))
    assert 0.2 < med < 5.0  # O/E is centred near 1 by construction


def test_expected_large_region_stats_only():
    out = expected_observed(region="chr17", resolution=100_000)
    assert out["oe_matrix"] is None and "cap" in out["note"]


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
