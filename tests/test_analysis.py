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
# P(s) slopes at 100 kb, per arm (measured on the shipped file; one constant per quantity)
P_ARM_SLOPE_100KB = -1.16246
Q_ARM_SLOPE_100KB = -1.24074
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
    # exact integer sum, each contact once (upper triangle incl. diagonal); the whole
    # symmetric square summed instead would read 2,323,926
    assert out["raw_contacts_sum"] == 1_890_939
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
    # the ranking window is the ground truth here; the rest of the sentence discloses
    # that top_boundaries is a subset, and is asserted by its own test
    assert out["ranked_by"].startswith("boundary_strength at the 250000 bp window")
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
    """The centromere is BOTH outside every arm and unmappable - say both."""
    with pytest.raises(AnalysisError) as e:
        compartments(region="chr17:23,000,000-26,000,000")
    assert "not contained in any single region of the view" in str(e.value)
    assert "ICE-filtered" in str(e.value)


def test_compartments_genome_summary():
    out = compartments()
    assert 0.2 < out["A_fraction_of_file"] < 0.8
    assert "chr17" in out["fraction_scope"]
    assert out["bins_used"] > 500


def test_compartments_never_claims_an_A_fraction_when_unphased():
    """An arbitrary sign cannot support an 'A fraction' - the field must stay null."""
    out = compartments(resolution=1_000_000)  # phasing track exists only at 100 kb
    assert "UNPHASED" in out["sign_convention"]
    assert out["A_fraction_of_file"] is None
    assert 0.0 < out["positive_E1_fraction_of_file"] < 1.0


# --- virtual_4c ------------------------------------------------------------


def test_virtual4c_ground_truth_shape():
    out = virtual_4c(viewpoint="chr17:63,000,000-63,100,000")  # default 10 kb
    bands = out["distance_band_means"]
    assert set(bands) == {"0kb-100kb", "100kb-1Mb", "1Mb-5Mb", "5Mb-10Mb"}
    vals = list(bands.values())
    assert all(v is not None for v in vals)
    assert vals == sorted(vals, reverse=True)  # decaying distance-band means
    assert all(p["balanced"] is not None for p in out["profile_points"])
    assert "genuine zero" in out["profile_note"]  # one convention, stated once
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
    # chr17q at 100 kb, per-arm expected: the ONE value for this quantity, pinned tightly
    # (a wider band once hid a 0.015 drift between two tests of the same number)
    assert out["ps_slope"] == pytest.approx(Q_ARM_SLOPE_100KB, rel=0.01)
    assert len(out["expected_curve_points"]) > 20
    # the masked head must not appear in the reported curve
    assert min(p["dist_bp"] for p in out["expected_curve_points"]) >= 200_000


def test_expected_ps_slope_is_negative_at_every_resolution():
    """Contact frequency falls with distance; a positive slope means a contaminated fit."""
    for res, expected in ((10_000, -1.25937), (100_000, Q_ARM_SLOPE_100KB), (1_000_000, -1.07962)):
        out = expected_observed(region=A_BLOCK, resolution=res)
        assert out["ps_slope"] == pytest.approx(expected, rel=0.02), res
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
    assert p_arm["ps_slope"] == pytest.approx(P_ARM_SLOPE_100KB, rel=0.01)
    assert q_arm["ps_slope"] == pytest.approx(Q_ARM_SLOPE_100KB, rel=0.01)
    assert "not only the requested region" in p_arm["curve_scope"]


def test_expected_curve_tail_is_not_flattened_by_rounding():
    """Six-decimal rounding collapsed a decade of decay into one repeated constant."""
    pts = expected_observed(region=A_BLOCK, resolution=10_000)["expected_curve_points"]
    tail = [p["expected"] for p in pts[-10:] if p["expected"] is not None]
    assert len(set(tail)) == len(tail), f"curve tail is flat: {tail}"
    assert tail == sorted(tail, reverse=True) or tail[0] > tail[-1]


def test_expected_region_straddling_the_centromere_is_a_clear_error():
    """Names the arms it is not inside, and offers them as the road forward."""
    with pytest.raises(AnalysisError) as e:
        expected_observed(region="chr17:20,000,000-25,000,000")
    assert "not contained in any single region of the view" in str(e.value)
    assert "chr17p" in str(e.value) and "chr17q" in str(e.value)


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
    c = compartments(file=f, region="cA:0-300,000")
    assert c["bins_used"] > 0 and "UNPHASED" in c["sign_convention"]
    e = expected_observed(file=f, region="cA:0-300,000")
    assert e["view"].startswith("cA") and e["expected_curve_points"]


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
    with pytest.raises(AnalysisError) as e:
        insulation_tads(region="chr17:23,000,000-25,000,000")
    # the centromere lies outside both arms of the bundled view, and that is the more
    # precise diagnosis than "unmappable" - it names the view and the road forward
    assert "lies outside every region of the view" in str(e.value)
    assert "chr17p" in str(e.value) and "chr17q" in str(e.value)


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


def test_contact_counts_agree_with_the_file_total():
    """One counting convention: a whole-chromosome region must sum to matrix_summary's
    total, which counts each contact once. A symmetric square summed whole would not."""
    total = matrix_summary()["resolutions"][0]["total_contacts"]
    whole = contacts_at_locus(region="chr17", resolution=1_000_000)
    assert whole["raw_contacts_sum"] == total
    assert "counted once" in whole["counting"]
    # two overlapping regions say so, rather than double-counting in silence
    ov = contacts_at_locus(
        region="chr17:50,000,000-50,500,000",
        region2="chr17:50,200,000-50,700,000",
        resolution=10_000,
    )
    assert "overlap" in ov["counting"]


# --- round-5: the advertised road for a user's own file must actually exist -----


def test_user_can_supply_a_view_and_a_phasing_track(tmp_path, two_chromosome_cool):
    """The unphased message names phasing_track; that parameter has to work."""
    import numpy as np
    import pandas as pd

    f = two_chromosome_cool
    view = tmp_path / "view.bed"
    view.write_text("cA\t0\t600000\tcA_all\ncB\t0\t400000\tcB_all\n")
    out = insulation_tads(file=f, view=str(view), region="cA:0-500,000", top_n=1)
    assert "supplied view" in out["view"]

    # This track used to be uniform random noise, and the tool returned a confident A/B
    # call from it - so the test asserted the very defect a later round found: a track that
    # correlates with nothing still oriented the eigenvector. A phasing track earns its
    # orientation by tracking compartment identity, so the one here does: it is built from
    # the file's own unphased E1, which is what a real GC or gene-density track approximates.
    unphased = compartments(file=f, view=str(view))
    assert "UNPHASED" in unphased["sign_convention"]
    e1 = compartments(file=f, view=str(view), region="cA:0-600,000")["E1_track"]
    by_start = {int(pt["start"]): pt["E1"] for pt in e1 or []}

    rng = np.random.default_rng(0)
    frames = []
    for chrom, n in (("cA", 60), ("cB", 40)):  # the track must cover every view region
        starts = np.arange(n) * 10_000
        if chrom == "cA":
            base = np.array([by_start.get(int(x)) or 0.0 for x in starts], dtype=float)
        else:
            base = np.zeros(n)
        value = 0.45 + 0.1 * base + rng.normal(0, 0.005, n)
        frames.append(
            pd.DataFrame(
                {"chrom": [chrom] * n, "start": starts,
                 "end": starts + 10_000, "GC": value}
            )
        )
    track = tmp_path / "gc.tsv"
    pd.concat(frames).to_csv(track, sep="\t", index=False)
    phased = compartments(file=f, view=str(view), phasing_track=str(track), region="cA:0-300,000")
    assert phased["region_call"] in {"A", "B"}, phased["sign_convention"]
    assert "UNPHASED" not in phased["sign_convention"]
    assert "Pearson r" in phased["sign_convention"]

    # and the counterpart the fix exists for: noise orients nothing, so nothing is called
    noise = tmp_path / "noise.tsv"
    pd.concat(
        [
            pd.DataFrame(
                {"chrom": [c] * n, "start": np.arange(n) * 10_000,
                 "end": (np.arange(n) + 1) * 10_000,
                 "GC": np.random.default_rng(1).uniform(0.3, 0.6, n)}
            )
            for c, n in (("cA", 60), ("cB", 40))
        ]
    ).to_csv(noise, sep="\t", index=False)
    guessed = compartments(file=f, view=str(view), phasing_track=str(noise), region="cA:0-300,000")
    assert guessed["region_call"] == "unphased", guessed["sign_convention"]


def test_phasing_track_not_covering_the_view_is_an_agent_readable_error(
    tmp_path, two_chromosome_cool
):
    """cooltools raises a raw error here; the caller needs to know it is fixable input."""
    import numpy as np
    import pandas as pd

    view = tmp_path / "view.bed"
    view.write_text("cA\t0\t600000\tcA_all\ncB\t0\t400000\tcB_all\n")
    track = tmp_path / "partial.tsv"
    pd.DataFrame(
        {
            "chrom": ["cA"] * 60,
            "start": np.arange(60) * 10_000,
            "end": (np.arange(60) + 1) * 10_000,
            "GC": np.random.default_rng(1).uniform(0.3, 0.6, 60),
        }
    ).to_csv(track, sep="\t", index=False)
    # the coverage check now names the uncovered region specifically
    with pytest.raises(AnalysisError, match="covers less than half of cB_all"):
        compartments(
            file=two_chromosome_cool,
            view=str(view),
            phasing_track=str(track),
            region="cA:0-300,000",
        )


def test_bad_view_and_track_files_are_agent_readable(tmp_path, two_chromosome_cool):
    bad = tmp_path / "bad.bed"
    bad.write_text("cA\t0\n")
    with pytest.raises(Exception, match="at least"):
        insulation_tads(file=two_chromosome_cool, view=str(bad), top_n=1)
    with pytest.raises(Exception, match="No such view file"):
        insulation_tads(file=two_chromosome_cool, view=str(tmp_path / "nope.bed"), top_n=1)


def test_compartments_refuses_to_average_two_eigendecompositions():
    """A region spanning both arms gets one A/B call from two independent decompositions."""
    with pytest.raises(AnalysisError) as e:
        compartments(region="chr17:20,000,000-30,000,000")
    # caught by the view-containment check, which names both arms
    assert "chr17p" in str(e.value) and "chr17q" in str(e.value)


def test_assembly_never_leaks_a_pandas_series_name(tiny_unbalanced_cool):
    """clr.chromsizes.name is 'length'; it must never be reported as an assembly."""
    out = matrix_summary(file=tiny_unbalanced_cool)
    assert out["assembly"] == "unknown"


def test_virtual4c_window_narrows_the_profile():
    wide = virtual_4c(viewpoint="chr17:63,000,000-63,100,000")
    near = virtual_4c(viewpoint="chr17:63,000,000-63,100,000", window_bp=1_000_000)
    assert near["window_bp"] == 1_000_000
    assert len(near["profile_points"]) > 0
    assert max(p["start"] for p in near["profile_points"]) < max(
        p["start"] for p in wide["profile_points"]
    )


def test_expected_curve_discloses_its_downsampling():
    out = expected_observed(region=A_BLOCK, resolution=10_000)
    assert "every" in out["curve_note"] or "all" in out["curve_note"]


# --- round-6: the bring-your-own-file road, where every finding landed ----------


def test_view_file_with_the_documented_header_is_accepted(tmp_path, two_chromosome_cool):
    """The README names chrom/start/end/name; a file in that exact shape must load."""
    view = tmp_path / "view.bed"
    view.write_text("chrom\tstart\tend\tname\ncA\t0\t600000\tcA_all\ncB\t0\t400000\tcB_all\n")
    out = insulation_tads(file=two_chromosome_cool, view=str(view), region="cA:0-500,000", top_n=1)
    assert "supplied view" in out["view"]


def test_view_naming_chromosomes_the_file_lacks_is_named_as_such(tmp_path, two_chromosome_cool):
    """'17' and 'chr17' are different names, and the caller can fix that."""
    view = tmp_path / "bad.bed"
    view.write_text("17\t0\t600000\tx\n")
    with pytest.raises(AnalysisError, match="does not"):
        insulation_tads(file=two_chromosome_cool, view=str(view), top_n=1)


def test_region_outside_the_view_is_not_blamed_on_the_data(tmp_path, two_chromosome_cool):
    """Blaming ICE filtering for mappable bins sends the caller hunting a phantom."""
    view = tmp_path / "part.bed"
    view.write_text("cA\t0\t300000\tcA_part\n")
    outside = "cA:400,000-500,000"
    f = two_chromosome_cool
    v = str(view)
    # the aggregating tools need one curve / one decomposition, so they say the region
    # is not contained in a single view region
    for tool in (compartments, expected_observed):
        with pytest.raises(AnalysisError, match="not contained in any single region of the view"):
            tool(file=f, view=v, region=outside)
    # insulation is per-bin, so it is not gated on containment - but it must still name
    # the view rather than blaming ICE filtering for bins it simply never computed
    with pytest.raises(AnalysisError, match="lies outside every region of the view"):
        insulation_tads(file=f, view=v, region=outside)


def test_phasing_coverage_is_checked_without_a_view_too(tmp_path, two_chromosome_cool):
    """The no-view case IS the default for a user's own file; it must be guarded."""
    import numpy as np
    import pandas as pd

    track = tmp_path / "partial.tsv"
    pd.DataFrame(
        {
            "chrom": ["cA"] * 60,
            "start": np.arange(60) * 10_000,
            "end": (np.arange(60) + 1) * 10_000,
            "GC": np.random.default_rng(2).uniform(0.3, 0.6, 60),
        }
    ).to_csv(track, sep="\t", index=False)
    with pytest.raises(AnalysisError, match="covers less than half of cB"):
        compartments(file=two_chromosome_cool, phasing_track=str(track), region="cA:0-300,000")


def test_virtual4c_window_excludes_separations_from_the_bands_too():
    """Band means must not summarise separations the caller removed."""
    out = virtual_4c(viewpoint="chr17:63,000,000-63,100,000", window_bp=1_000_000)
    assert set(out["distance_band_means"]) == {"0kb-100kb", "100kb-1Mb"}
    assert "1Mb-5Mb" not in out["distance_band_bins"]  # absent, not reported as zero
    assert "limited to 1,000,000 bp" in out["coverage_note"]
    assert "whole chromosome" not in out["coverage_note"]


def test_all_null_balanced_matrix_explains_itself():
    """A wall of nulls beside balanced:true reads like a result; it must not."""
    out = contacts_at_locus(region="chr17:23,000,000-23,300,000", resolution=10_000)
    assert out["balanced_matrix"] is not None
    assert all(v is None for row in out["balanced_matrix"] for v in row)
    assert "ICE-filtered" in out["note"]


# --- round-7: malformed inputs are the caller's to fix, and credit goes where due ---


@pytest.mark.parametrize(
    "label,body,expected",
    [
        (
            "duplicate names",
            "cA\t0\t300000\tdup\ncA\t300000\t600000\tdup\n",
            "reuses the region name",
        ),
        ("overlapping", "cA\t0\t400000\ta\ncA\t300000\t600000\tb\n", "overlap"),
        ("past the end", "cA\t0\t9000000\tlong\n", "runs past the end"),
        ("zero width", "cA\t300000\t300000\tz\n", "zero or negative width"),
    ],
)
def test_malformed_views_are_named_not_reported_as_internal_bugs(
    tmp_path, two_chromosome_cool, label, body, expected
):
    """cooltools rejects these with one opaque sentence; each is a fixable input."""
    view = tmp_path / "v.bed"
    view.write_text(body)
    with pytest.raises(AnalysisError, match=expected):
        insulation_tads(file=two_chromosome_cool, view=str(view), top_n=1)


def test_an_unsorted_view_is_normalised_not_refused(tmp_path, two_chromosome_cool):
    """Row order carries no meaning, so sorting is a normalisation, not a refusal."""
    view = tmp_path / "unsorted.bed"
    view.write_text("cA\t300000\t600000\tb\ncA\t0\t300000\ta\n")
    out = insulation_tads(file=two_chromosome_cool, view=str(view), top_n=1)
    assert "supplied view (2 regions)" in out["view"]


def test_phasing_track_binsize_mismatch_is_named():
    """The repo's own 100 kb track at 10 kb once raised a raw library error."""
    from hic_mcp.data import data_dir

    with pytest.raises(AnalysisError, match="binned at 100,000 bp"):
        compartments(
            resolution=10_000,
            phasing_track=str(data_dir() / "gc_100kb.tsv"),
            region=A_BLOCK,
        )


def test_sign_convention_credits_the_track_actually_used(tmp_path):
    """A user's track must not be described as the bundled one - it changes the answer."""
    import pandas as pd

    from hic_mcp.data import data_dir

    gc = pd.read_csv(data_dir() / "gc_100kb.tsv", sep="\t")
    gc["GC"] = -gc["GC"]  # a legitimate, sign-inverted orientation track
    inverted = tmp_path / "inverted.tsv"
    gc.to_csv(inverted, sep="\t", index=False)

    bundled = compartments(region=A_BLOCK)
    supplied = compartments(region=A_BLOCK, phasing_track=str(inverted))
    assert "bundled GC track" in bundled["sign_convention"]
    assert "track you supplied" in supplied["sign_convention"]
    # and the supplied track really is in force: the call flips with the sign
    assert bundled["region_call"] == "A"
    assert supplied["region_call"] == "B"


def test_curve_note_distinguishes_measured_from_fitted():
    """Saying N 'were fitted' when the slope used far fewer overstates the fit."""
    out = expected_observed(region=A_BLOCK, resolution=10_000)
    note = out["curve_note"]
    assert "measured separations" in note and "fitted range" in note


def test_a_narrow_region_is_not_called_unmappable():
    """A 1-bin region is all masked diagonals, not bad data - and the note must say so."""
    out = expected_observed(region="chr17:50,100,000-50,200,000", resolution=100_000)
    assert out["oe_matrix"] is not None
    assert all(v is None for row in out["oe_matrix"] for v in row)
    assert "The data here is fine" in out["note"]
    assert "ICE-filtered" not in out["note"]


def test_band_coverage_reflects_the_window():
    out = virtual_4c(viewpoint="chr17:63,000,000-63,100,000", window_bp=1_000_000)
    assert out["distance_bands_cover_bp"] == [0, 1_000_000]


# --- round-8 ------------------------------------------------------------------


@pytest.mark.parametrize(
    "label,body,expected_regions",
    [
        (
            "#chrom header",
            "#chrom\tstart\tend\tname\ncA\t0\t200000\tr1\ncA\t200000\t400000\tr2\n",
            2,
        ),
        (
            "plain header",
            "chrom\tstart\tend\tname\ncA\t0\t200000\tr1\ncA\t200000\t400000\tr2\n",
            2,
        ),
        ("comment line", "# a note\ncA\t0\t200000\tr1\ncA\t200000\t400000\tr2\n", 2),
        ("bare rows", "cA\t0\t200000\tr1\ncA\t200000\t400000\tr2\n", 2),
    ],
)
def test_every_view_header_shape_keeps_every_region(tmp_path, label, body, expected_regions):
    """A '#chrom' header once ate the first region silently - wrong science, not an error."""
    from hic_mcp.data import load_view_file

    f = tmp_path / "v.bed"
    f.write_text(body)
    assert len(load_view_file(str(f))) == expected_regions, label


def test_trans_block_is_not_called_balanced_under_a_cis_only_solution(two_chromosome_cool):
    """Weights fitted per chromosome were never fitted for trans."""
    out = contacts_at_locus(
        file=two_chromosome_cool,
        region="cA:0-200,000",
        region2="cB:0-200,000",
        resolution=10_000,
    )
    assert out["balanced"] is False
    assert out["balanced_mean"] is None
    assert "cis-only" in out["note"]
    # the cis path on the same file is still balanced
    cis = contacts_at_locus(file=two_chromosome_cool, region="cA:0-200,000", resolution=10_000)
    assert cis["balanced"] is True


def test_compartments_refuses_before_an_out_of_memory_kill(monkeypatch):
    """An OOM kill hands the agent a dead transport; a refusal hands it a reason."""
    import hic_mcp.analysis as an

    monkeypatch.setattr(an, "COMPARTMENT_MEMORY_CAP_GB", 0.001)
    with pytest.raises(AnalysisError, match="Refusing to run"):
        an.compartments(region=A_BLOCK)


def test_the_real_memory_cap_lets_the_bundled_demo_through():
    """The guard must protect against a genome, not obstruct the demo."""
    assert compartments(region=A_BLOCK)["region_call"] == "A"


def test_unreadable_files_are_diagnosed_not_blamed_on_the_server(tmp_path):
    """A Git-LFS pointer or truncated download is the file's problem, and fixable."""
    from hic_mcp.data import DataError

    lfs = tmp_path / "pointer.cool"
    lfs.write_text("version https://git-lfs.github.com/spec/v1\noid sha256:deadbeef\n")
    with pytest.raises(DataError, match="Git-LFS pointer"):
        matrix_summary(file=str(lfs))

    truncated = tmp_path / "truncated.mcool"
    truncated.write_bytes(b"\x89HDF\r\n\x1a\n" + b"\x00" * 64)
    with pytest.raises(DataError, match="not a readable cooler"):
        matrix_summary(file=str(truncated))


def test_a_thin_phasing_track_is_refused_with_a_reason(tmp_path):
    """Naming the chromosome is not covering it; cooltools needs nearly every bin."""
    track = tmp_path / "thin.bed"
    track.write_text("chr17\t0\t100000\t0.45\nchr17\t100000\t200000\t0.5\n")
    with pytest.raises(AnalysisError, match="covers less than half"):
        compartments(phasing_track=str(track), region=A_BLOCK)


def test_a_header_is_detected_by_shape_not_by_its_first_word(tmp_path):
    """'Chromosome' is a header too; only the numeric shape of row one can tell."""
    from hic_mcp.data import load_track_file, load_view_file

    odd = tmp_path / "odd.tsv"
    odd.write_text("Chromosome\tstart\tend\tGC\nchr17\t0\t100000\t0.45\n")
    assert len(load_track_file(str(odd))) == 1  # the header is not read as data

    bare = tmp_path / "bare.bed"
    bare.write_text("cA\t0\t200000\tr1\n")
    assert len(load_view_file(str(bare))) == 1


def test_insulation_accepts_a_whole_chromosome_region():
    """`region` filters which bins are reported; it is not a containment claim.

    Requiring containment refused region='chr17' while the same call with no region
    returned exactly those boundaries - a refusal the tool contradicted itself on.
    """
    whole = insulation_tads(region="chr17", top_n=1)
    unfiltered = insulation_tads(top_n=1)
    assert whole["top_boundaries"][0]["locus"] == unfiltered["top_boundaries"][0]["locus"]


# --- loop v2, round 1: the resolutions the demo never exercised ---------------


def _balanced_cool(tmp_path, binsize, n=80):
    import cooler

    p = str(tmp_path / f"f_{binsize}.cool")
    bins = pd.DataFrame(
        {"chrom": "cA", "start": np.arange(n) * binsize, "end": (np.arange(n) + 1) * binsize}
    )
    i, j = np.triu_indices(n)
    c = np.random.default_rng(1).poisson(150.0 / (1 + (j - i)) ** 0.8).astype(int)
    k = c > 0
    cooler.create_cooler(p, bins, pd.DataFrame({"bin1_id": i[k], "bin2_id": j[k], "count": c[k]}))
    cooler.balance_cooler(cooler.Cooler(p), store=True, cis_only=True)
    return p


@pytest.mark.parametrize("binsize", [4_000, 8_000, 15_000, 20_000, 30_000])
def test_auto_windows_work_at_every_ordinary_resolution(tmp_path, binsize):
    """Zero user arguments must never crash: windows are snapped to whole bins."""
    out = insulation_tads(file=_balanced_cool(tmp_path, binsize), top_n=1)
    assert all(w % binsize == 0 for w in out["windows_bp"])
    assert min(out["windows_bp"]) >= 3 * binsize


def test_user_windows_off_the_bin_grid_are_named(tmp_path):
    with pytest.raises(AnalysisError, match="not multiples of the 8,000 bp bin size"):
        insulation_tads(file=_balanced_cool(tmp_path, 8_000), windows_bp=[100_000], top_n=1)


def test_phasing_track_with_a_name_column_is_refused_as_a_track(tmp_path):
    """A BED4 whose fourth column is a name is a region list, not a phasing track."""
    from hic_mcp.data import DataError

    p = _balanced_cool(tmp_path, 10_000)
    track = tmp_path / "names.bed"
    track.write_text("".join(f"cA\t{i * 10000}\t{(i + 1) * 10000}\tregion{i}\n" for i in range(80)))
    with pytest.raises(DataError, match="holds no numbers"):
        compartments(file=p, phasing_track=str(track), region="cA:0-300,000")


def test_top_n_below_one_is_refused():
    with pytest.raises(AnalysisError, match="top_n must be at least 1"):
        insulation_tads(top_n=0)


def test_contacts_answers_a_whole_chromosome_without_densifying():
    """A whole chromosome at 10 kb once densified 1.5 GB to return a few scalars."""
    out = contacts_at_locus(region="chr17", resolution=10_000)
    assert out["shape_bins"] == [8326, 8326]
    assert out["raw_contacts_sum"] == 129_434_772  # same answer, sparse road
    assert out["balanced_matrix"] is None
    assert "sparse pixel table" in out["note"]


def test_contacts_sparse_and_dense_roads_agree():
    """The sparse road must give the dense road's numbers, not merely finish."""
    dense = contacts_at_locus(region="chr17:50,000,000-50,500,000", resolution=10_000)
    assert "sparse" not in (dense.get("note") or "")
    import hic_mcp.analysis as an

    saved = an.DENSE_FETCH_CAP_GB
    try:
        an.DENSE_FETCH_CAP_GB = 0.0  # force the sparse road for the same region
        sparse = contacts_at_locus(region="chr17:50,000,000-50,500,000", resolution=10_000)
    finally:
        an.DENSE_FETCH_CAP_GB = saved
    # EVERY statistic, not a chosen few: this guard once compared sum/max/balanced_max -
    # the three that happened to agree - while nonzero_fraction was out by 2x (upper
    # triangle counted against the full square) and balanced_mean used another denominator
    for key in ("raw_contacts_sum", "raw_contacts_max", "nonzero_fraction",
                "balanced_mean", "balanced_max", "counting", "shape_bins"):
        d, sp = dense.get(key), sparse.get(key)
        if isinstance(d, float) and isinstance(sp, float):
            assert sp == pytest.approx(d, rel=1e-6), key
        else:
            assert sp == d, key


@pytest.mark.parametrize(
    "region,region2",
    [
        # Every shape a caller can ask for. Three of these were added after a fix that
        # was checked across all six shapes but on one field only, while this sweep
        # checked all five fields but on three shapes - so each guard missed the other's
        # blind spot and a broken refactor went green in both.
        ("chr17:50,000,000-50,500,000", None),          # symmetric square
        ("chr17:50,000,000-50,300,000", "chr17:60,000,000-60,300,000"),  # disjoint, ascending
        ("chr17:60,000,000-60,300,000", "chr17:50,000,000-50,300,000"),  # disjoint, DESCENDING
        ("chr17:50,000,000-50,400,000", "chr17:50,200,000-50,600,000"),  # OVERLAPPING
        ("chr17:50,000,000-50,300,000", "chr17:50,000,000-50,300,000"),  # IDENTICAL
        ("chr17:50,000,000-50,300,000", "chr17:50,300,000-50,600,000"),  # adjacent
        ("chr17:23,000,000-23,300,000", None),          # fully ICE-filtered
    ],
)
def test_sparse_and_dense_agree_on_every_region_shape(region, region2):
    """cooler stores one triangle: a square's cells and sums must be mirrored, a block's not."""
    import hic_mcp.analysis as an

    kwargs = {"region": region, "resolution": 10_000}
    if region2:
        kwargs["region2"] = region2
    dense = contacts_at_locus(**kwargs)
    saved = an.DENSE_FETCH_CAP_GB
    try:
        an.DENSE_FETCH_CAP_GB = 0.0
        sparse = contacts_at_locus(**kwargs)
    finally:
        an.DENSE_FETCH_CAP_GB = saved
    for key in ("raw_contacts_sum", "raw_contacts_max", "nonzero_fraction",
                "balanced_mean", "balanced_max"):
        d, sp = dense.get(key), sparse.get(key)
        if isinstance(d, float) and isinstance(sp, float):
            assert sp == pytest.approx(d, rel=1e-6), f"{region} {key}"
        else:
            assert sp == d, f"{region} {key}"


def test_virtual4c_points_and_bands_share_one_zero_convention():
    """A mappable zero-contact bin must be a zero in the points, not silently dropped."""
    out = virtual_4c(viewpoint="chr17:63,000,000-63,100,000")
    vals = [p["balanced"] for p in out["profile_points"]]
    assert any(v == 0.0 for v in vals), "no genuine zeros reported in the profile points"
    assert "One convention throughout" in out["profile_note"]


def test_a_directory_is_not_a_view_or_a_track(tmp_path, two_chromosome_cool):
    from hic_mcp.data import DataError

    with pytest.raises(DataError, match="is a directory"):
        insulation_tads(file=two_chromosome_cool, view=str(tmp_path), top_n=1)
    with pytest.raises(DataError, match="is a directory"):
        compartments(file=two_chromosome_cool, phasing_track=str(tmp_path), region="cA:0-300,000")


# --- v2 round 2: the axis the parity tests never varied -----------------------


@pytest.mark.parametrize("resolution", [10_000, 100_000])
def test_region_order_does_not_change_the_answer(resolution):
    """Same two loci, reversed. cooler stores one triangle, so a block entirely below
    the diagonal fetched as pixels returns NOTHING - which once read as "zero contacts,
    balanced: true". The parity tests varied region shape but never region order."""
    upstream = "chr17:50,000,000-53,000,000"
    downstream = "chr17:60,000,000-63,000,000"
    asc = contacts_at_locus(region=upstream, region2=downstream, resolution=resolution)
    desc = contacts_at_locus(region=downstream, region2=upstream, resolution=resolution)
    assert asc["raw_contacts_sum"] > 0, "fixture should have contacts between these loci"
    for key in ("raw_contacts_sum", "raw_contacts_max", "nonzero_fraction",
                "balanced_mean", "balanced_max"):
        d, u = asc.get(key), desc.get(key)
        if isinstance(d, float) and isinstance(u, float):
            assert u == pytest.approx(d, rel=1e-6), key
        else:
            assert u == d, key


def test_sparse_note_states_the_reason_that_actually_fired():
    """It once said "too large to hold as a dense matrix (0.0 GB)" - visibly false."""
    out = contacts_at_locus(region="chr17:50,000,000-53,000,000", resolution=10_000)
    assert "bins on a side" in out["note"]
    assert "0.0 GB" not in out["note"]


def test_ranked_by_discloses_that_it_filters():
    """boundary_counts_per_window and top_boundaries describe different populations."""
    out = insulation_tads(region="chr17:65,000,000-67,000,000", top_n=10)
    assert "ONLY" in out["ranked_by"] and "subset" in out["ranked_by"]
    listed = len(out["top_boundaries"])
    assert listed <= out["boundary_counts_per_window"][str(out["windows_bp"][1])]


# --- v2 round 3: the shape axis the parity tests never covered ----------------


@pytest.mark.parametrize(
    "label,region,region2",
    [
        ("single region", "chr17:50,000,000-53,000,000", None),
        ("overlapping", "chr17:50,000,000-53,000,000", "chr17:52,000,000-55,000,000"),
        ("identical pair", "chr17:50,000,000-53,000,000", "chr17:50,000,000-53,000,000"),
        ("disjoint ascending", "chr17:50,000,000-53,000,000", "chr17:60,000,000-63,000,000"),
        ("disjoint descending", "chr17:60,000,000-63,000,000", "chr17:50,000,000-53,000,000"),
        ("adjacent", "chr17:50,000,000-53,000,000", "chr17:53,000,000-56,000,000"),
    ],
)
def test_sparse_road_equals_cooler_dense_truth(label, region, region2):
    """Checked against cooler itself, not against our own dense branch.

    cooler stores ONE triangle and a stored pixel fills two dense cells. Fetching one
    orientation undercounted overlapping regions by 16% and identical ones by 20%,
    while the response's own `counting` field promised the opposite convention.
    """
    import numpy as np

    from hic_mcp.data import open_matrix, resolve_input_path

    clr = open_matrix(resolve_input_path(None), 10_000, default=10_000)
    dense = np.nan_to_num(clr.matrix(balance=False).fetch(region, region2 or region))
    truth = int(np.triu(dense).sum()) if region2 is None else int(dense.sum())

    kwargs = {"region": region, "resolution": 10_000}
    if region2:
        kwargs["region2"] = region2
    assert contacts_at_locus(**kwargs)["raw_contacts_sum"] == truth, label


def test_downsampling_ordinals_read_as_english():
    """Agent-facing text once said "every 2th of the curve"."""
    note = expected_observed(region="chr17:5,000,000-7,000,000", resolution=100_000)["curve_note"]
    assert "2th" not in note and "3th" not in note and "1th" not in note


def test_coverage_note_never_claims_more_than_the_bands_cover():
    """It once reported a 20 Mb window while the bands stop at 10 Mb."""
    out = virtual_4c(viewpoint="chr17:63,000,000-63,100,000", window_bp=20_000_000)
    assert "10,000,000 bp" in out["coverage_note"]
    assert out["distance_bands_cover_bp"] == [0, 10_000_000]


def test_raw_only_request_says_why_there_is_no_matrix():
    """The tool promises a matrix for small windows; silence reads as a missing result."""
    out = contacts_at_locus(
        region="chr17:50,000,000-50,200,000", resolution=10_000, balanced=False
    )
    assert out["balanced_matrix"] is None
    assert "balanced=false" in out["note"].lower()
    # and a normal balanced call still returns the matrix with no spurious note
    ok = contacts_at_locus(region="chr17:50,000,000-50,200,000", resolution=10_000)
    assert ok["balanced_matrix"] is not None and ok.get("note") is None


def test_a_track_that_orients_nothing_does_not_produce_an_ab_call(tmp_path):
    """cooltools flips E1 to a positive correlation with the track WHATEVER that is.

    Binning, coverage and density guards all pass on a track whose values have been shuffled
    - same bins, same numbers, relationship destroyed - and the tool used to return the
    repo's own documented A block as B, with sign consistency 1.0 and no caveat anywhere.
    A tool that refuses to guess without a track must not guess with a useless one.
    """
    import numpy as np
    import pandas as pd

    from hic_mcp.data import data_dir

    src = pd.read_csv(data_dir() / "gc_100kb.tsv", sep="\t")
    shuffled = src.copy()
    col = shuffled.columns[3]
    shuffled[col] = np.random.default_rng(0).permutation(shuffled[col].to_numpy())
    track = tmp_path / "shuffled_gc.tsv"
    shuffled.to_csv(track, sep="\t", index=False)

    region = "chr17:50,100,000-51,100,000"  # the documented A block
    good = compartments(region=region, resolution=100_000)
    assert good["region_call"] == "A"
    assert "Pearson r" in good["sign_convention"]

    bad = compartments(region=region, resolution=100_000, phasing_track=str(track))
    assert bad["region_call"] == "unphased", bad["sign_convention"]
    assert "UNPHASED" in bad["sign_convention"]
    assert "Pearson r" in bad["sign_convention"]  # the number that justified the refusal

    # and the whole-file survey: an A fraction is a claim the sign has to earn
    whole_bad = compartments(resolution=100_000, phasing_track=str(track))
    assert whole_bad["A_fraction_of_file"] is None
    assert whole_bad["positive_E1_fraction_of_file"] is not None
    whole_good = compartments(resolution=100_000)
    assert whole_good["A_fraction_of_file"] is not None


@pytest.mark.parametrize("balanced", [True, False])
def test_the_sparse_road_names_itself_whichever_flag_the_caller_passed(balanced):
    """This server discloses which road produced a number - it used to go silent on one.

    With balanced=False the response attributed the missing matrix entirely to the caller's
    flag, when the bin cap was the operative reason and the statistics had come from the
    sparse pixel table rather than a dense fetch.
    """
    out = contacts_at_locus(region="chr17:50,000,000-53,000,000", resolution=10_000,
                            balanced=balanced)
    assert out["balanced_matrix"] is None
    assert "sparse pixel table" in out["note"], out["note"]
    if not balanced:
        assert "you passed balanced=false" in out["note"]
