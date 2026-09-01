"""The six analyses, as pure functions over local .cool/.mcool files.

Every function runs the real open2c computation and returns a plain dict that
names its method, the resolution used, and whether values are ICE-balanced.
No MCP imports here; anticipated failures raise AnalysisError (or DataError
from the data helpers) with agent-facing messages.
"""

import json

import numpy as np
from cooler import Cooler
from cooltools import eigs_cis, expected_cis, insulation, virtual4c

from hic_mcp.data import (
    is_demo,
    list_resolutions,
    load_arms_view,
    load_gc_track,
    open_matrix,
    parse_region_checked,
    resolve_input_path,
)


class AnalysisError(ValueError):
    """A problem with the requested computation - message is agent-facing."""


MATRIX_BIN_CAP = 50  # a full matrix is returned only at or under this many bins per side
PROFILE_POINT_CAP = 400
IGNORE_DIAGS = 2  # cooltools' expected default: the first 2 diagonals carry no measurement
CLASSIC_TAD_WINDOWS = (100_000, 250_000, 500_000)  # mammalian TAD scale
COARSE_WINDOW_BINS = (3, 5, 10)  # fallback multipliers when bins are too big for the above
MIN_BINS_FOR_CONSISTENCY = 3  # below this a sign-consistency figure is vacuous


def _finite(x: float) -> float | None:
    """Round to 6 SIGNIFICANT figures, not 6 decimals.

    Contact frequencies span decades; fixed-decimal rounding would flatten the tail of
    a P(s) curve into a run of identical constants.
    """
    return float(f"{float(x):.6g}") if np.isfinite(x) else None


def _weights_present(clr: Cooler) -> bool:
    return "weight" in clr.bins().columns


def _check_viewpoint_not_filtered(clr: Cooler, chrom: str, start: int, end: int) -> None:
    """A viewpoint inside ICE-filtered (NaN-weight) bins has no balanced signal."""
    w = clr.bins().fetch(f"{chrom}:{start}-{end}")["weight"]
    if w.isna().all():
        raise AnalysisError(
            f"The viewpoint {chrom}:{start:,}-{end:,} falls entirely in ICE-filtered bins "
            "(low-mappability or centromeric - no balanced signal exists there). "
            "Choose a viewpoint outside filtered regions; matrix_summary and "
            "insulation_tads can help locate usable regions."
        )


def matrix_summary(file: str | None = None) -> dict:
    """What is in this Hi-C file: chromosomes, resolutions, contacts, balancing."""
    path = resolve_input_path(file)
    resolutions = list_resolutions(path)
    per_res = []
    meta: dict = {}
    for res in resolutions:
        clr = open_matrix(path, res, res)
        info = clr.info
        per_res.append(
            {
                "resolution_bp": res,
                "bins": int(info["nbins"]),
                "nonzero_pixels": int(info["nnz"]),
                "total_contacts": int(info["sum"]),
                "balanced": _weights_present(clr),
            }
        )
        if not meta:
            raw = info.get("metadata")
            if isinstance(raw, str):
                try:
                    meta = json.loads(raw)
                except json.JSONDecodeError:
                    meta = {"note": raw}
            elif isinstance(raw, dict):
                meta = raw
    clr = open_matrix(path, resolutions[0], resolutions[0])
    return {
        "file": path.name,
        "is_bundled_demo": is_demo(path),
        "assembly": clr.info.get("genome-assembly") or clr.chromsizes.name or "unknown",
        "chromosomes": {str(c): int(s) for c, s in clr.chromsizes.items()},
        "resolutions": per_res,
        "balanced": all(r["balanced"] for r in per_res),
        "provenance": meta,
        "method": "cooler.Cooler.info over every resolution in the file",
    }


def contacts_at_locus(
    file: str | None = None,
    region: str = "chr17:50,000,000-52,500,000",
    region2: str | None = None,
    resolution: int | None = None,
    balanced: bool = True,
) -> dict:
    """Contact-matrix statistics (and, for small windows, the matrix itself) at a locus."""
    path = resolve_input_path(file)
    available = list_resolutions(path)
    clr0 = open_matrix(path, available[0], available[0])
    chrom, start, end = parse_region_checked(clr0, region)
    if region2 is not None:
        c2, s2, e2 = parse_region_checked(clr0, region2)
    else:
        c2, s2, e2 = chrom, start, end
    span = max(end - start, e2 - s2)
    if resolution is None:
        fitting = [r for r in available if span / r <= MATRIX_BIN_CAP]
        resolution = fitting[0] if fitting else available[-1]
    clr = open_matrix(path, resolution, resolution)
    use_weights = balanced and _weights_present(clr)
    r1 = f"{chrom}:{start}-{end}"
    r2 = f"{c2}:{s2}-{e2}"
    raw = clr.matrix(balance=False).fetch(r1, r2)
    out: dict = {
        "region": r1,
        "region2": r2 if region2 is not None else None,
        "resolution_used": resolution,
        "shape_bins": list(raw.shape),
        "raw_contacts_sum": int(np.nansum(raw)),
        "raw_contacts_max": int(np.nanmax(raw)) if raw.size else 0,
        "nonzero_fraction": round(float((raw > 0).mean()), 4) if raw.size else 0.0,
        "balanced": use_weights,
        "method": "cooler.Cooler.matrix().fetch (raw counts; ICE-balanced values when available)",
    }
    if use_weights and raw.size:
        bal = clr.matrix(balance=True).fetch(r1, r2)
        finite_bal = bal[np.isfinite(bal)]
        out["balanced_mean"] = _finite(finite_bal.mean()) if finite_bal.size else None
        out["balanced_max"] = _finite(finite_bal.max()) if finite_bal.size else None
        if max(bal.shape) <= MATRIX_BIN_CAP:
            out["balanced_matrix"] = [[_finite(v) for v in row] for row in bal]
        else:
            out["balanced_matrix"] = None
            out["note"] = (
                f"Matrix is {bal.shape[0]}x{bal.shape[1]} bins (cap {MATRIX_BIN_CAP} per side); "
                "statistics only. Narrow the region or pass a coarser resolution "
                f"(available: {', '.join(str(r) for r in available)})."
            )
    return out


def insulation_tads(
    file: str | None = None,
    region: str | None = None,
    resolution: int | None = None,
    windows_bp: list[int] | None = None,
    top_n: int = 10,
) -> dict:
    """Diamond insulation score and TAD-boundary calls (Crane et al. 2015 method)."""
    path = resolve_input_path(file)
    clr = open_matrix(path, resolution, default=10000)
    if not _weights_present(clr):
        raise AnalysisError(
            "This file has no ICE weights ('weight' column); insulation needs a balanced "
            "matrix. Balance it first with `cooler balance`."
        )
    binsize = int(clr.binsize)
    # Prefer the classic TAD-scale windows whenever the bin size can carry them (a
    # diamond needs >= 3 bins); only fall back to bin-scaled windows on coarse data,
    # so a 100 kb file still gets TAD-scale calls rather than compartment-scale ones.
    if windows_bp:
        windows = windows_bp
    elif 3 * binsize <= min(CLASSIC_TAD_WINDOWS):
        windows = list(CLASSIC_TAD_WINDOWS)
    else:
        windows = [n * binsize for n in COARSE_WINDOW_BINS]
    for w in windows:
        if w < 3 * binsize:
            raise AnalysisError(
                f"Window {w} bp is too small for {binsize} bp bins (a diamond needs at "
                f"least 3 bins, i.e. {3 * binsize:,} bp here). Pass larger windows_bp, or "
                "omit resolution to use this file's finest level."
            )
    ins = insulation(clr, windows, verbose=False)
    if region is not None:
        chrom, start, end = parse_region_checked(clr, region)
        ins = ins[(ins["chrom"] == chrom) & (ins["end"] > start) & (ins["start"] < end)]
    # rank at the middle window - the classic TAD scale for mammalian 10 kb data;
    # a boundary that is also called at the flanking windows is the robust kind
    rank_w = sorted(windows)[len(windows) // 2]
    strength_col = f"boundary_strength_{rank_w}"
    bound_cols = {w: f"is_boundary_{w}" for w in windows}
    counts = {str(w): int(ins[c].sum()) for w, c in bound_cols.items()}
    called = ins[ins[bound_cols[rank_w]]].copy()
    called = called.sort_values(strength_col, ascending=False).head(top_n)
    boundaries = [
        {
            "locus": f"{r.chrom}:{int(r.start):,}-{int(r.end):,}",
            "strength": _finite(getattr(r, strength_col)),
            "log2_insulation": _finite(getattr(r, f"log2_insulation_score_{rank_w}")),
            "windows_detected": [w for w in windows if bool(getattr(r, bound_cols[w]))],
        }
        for r in called.itertuples()
    ]
    out = {
        "region": region,
        "resolution_used": binsize,
        "windows_bp": windows,
        "ranked_by": f"boundary_strength at the {rank_w} bp window",
        "boundary_counts_per_window": counts,
        "top_boundaries": boundaries,
        "balanced": True,
        "method": "cooltools.insulation (diamond insulation score; Li threshold boundary calls)",
    }
    if min(windows) > max(CLASSIC_TAD_WINDOWS):
        out["scale_note"] = (
            f"The smallest window here is {min(windows):,} bp, above the ~100-500 kb scale "
            "of mammalian TADs: these are large-scale insulation features (often "
            "compartment boundaries), not TAD boundaries. Use a finer resolution, or pass "
            "windows_bp, for TAD calls."
        )
    return out


def compartments(
    file: str | None = None,
    region: str | None = None,
    resolution: int | None = None,
) -> dict:
    """A/B compartment eigenvector (cis eigendecomposition of the observed/expected map)."""
    path = resolve_input_path(file)
    clr = open_matrix(path, resolution, default=100_000)
    if not _weights_present(clr):
        raise AnalysisError(
            "This file has no ICE weights; compartments need a balanced matrix. "
            "Balance it first with `cooler balance`."
        )
    demo = is_demo(path)
    phasing = view = None
    if demo and int(clr.binsize) == 100_000:
        phasing = load_gc_track()
        view = load_arms_view()
    eigvals, eigvecs = eigs_cis(clr, phasing_track=phasing, view_df=view, n_eigs=3)
    sign_convention = (
        "oriented by the bundled GC track (positive E1 = A = gene-dense/GC-rich)"
        if phasing is not None
        else "UNPHASED - the sign of E1 is mathematically arbitrary; supply your own "
        "orientation (e.g. GC or gene density) before calling A vs B"
    )
    out: dict = {
        "resolution_used": int(clr.binsize),
        "view": "chr17 p/q arms (bundled)" if view is not None else "whole chromosomes",
        "sign_convention": sign_convention,
        "eigenvalues": [
            {
                "region": str(r["name"] if "name" in r else r["chrom"]),
                "eigval1": _finite(r["eigval1"]),
                # raw eigenvalues scale with region size, so a bigger arm always looks
                # "stronger"; the share is the comparable number
                "variance_share": _finite(
                    abs(r["eigval1"])
                    / sum(abs(r[f"eigval{i}"]) for i in (1, 2, 3))
                ),
            }
            for _, r in eigvals.iterrows()
        ],
        "eigenvalue_note": (
            "eigval1 is unnormalised and grows with region size, so a longer arm always "
            "looks 'stronger'. variance_share is eigval1 as a fraction of the top three "
            "eigenvalues' absolute weight - use that to compare regions."
        ),
        "balanced": True,
        "method": "cooltools.eigs_cis (eigendecomposition of the cis observed/expected matrix)",
    }
    vec = eigvecs.dropna(subset=["E1"])
    if region is not None:
        chrom, start, end = parse_region_checked(clr, region)
        sub = vec[(vec["chrom"] == chrom) & (vec["end"] > start) & (vec["start"] < end)]
        if sub.empty:
            raise AnalysisError(
                f"No usable E1 bins in {region} (the region may be entirely ICE-filtered, "
                "e.g. centromeric)."
            )
        mean_e1 = float(sub["E1"].mean())
        out["region"] = region
        out["region_mean_E1"] = _finite(mean_e1)
        out["region_call"] = ("A" if mean_e1 > 0 else "B") if phasing is not None else "unphased"
        out["bins_used"] = int(len(sub))
        # a single bin is trivially "100% consistent with itself"; reporting that as
        # confidence invites an agent to call a knife-edge locus unambiguous
        if len(sub) >= MIN_BINS_FOR_CONSISTENCY:
            consistency = float((np.sign(sub["E1"]) == np.sign(mean_e1)).mean())
            out["region_sign_consistency"] = round(consistency, 3)
        else:
            out["region_sign_consistency"] = None
            out["confidence_note"] = (
                f"This region covers {len(sub)} bin(s) at {int(clr.binsize):,} bp, too few "
                "for a sign-consistency figure. The flanking track below shows whether the "
                "call is stable or sits on a compartment transition - read it before "
                "describing this locus as clearly A or B."
            )
        # always give the neighbourhood: a compartment call means little without the
        # bins either side of it
        flank = max(3, MIN_BINS_FOR_CONSISTENCY)
        lo = start - flank * int(clr.binsize)
        hi = end + flank * int(clr.binsize)
        ctx = vec[(vec["chrom"] == chrom) & (vec["end"] > lo) & (vec["start"] < hi)]
        if len(ctx) <= 200:
            out["E1_track"] = [
                {"start": int(r.start), "E1": _finite(r.E1)} for r in ctx.itertuples()
            ]
            signs = np.sign(ctx["E1"].to_numpy())
            if len(set(signs[signs != 0])) > 1:
                out["transition_note"] = (
                    "E1 changes sign within the flanking window: this locus is at or near "
                    "a compartment transition, so a single-region A/B label understates it."
                )
    else:
        positive_fraction = round(float((vec["E1"] > 0).mean()), 3)
        if phasing is not None:
            out["genome_A_fraction"] = positive_fraction
        else:
            # unphased: the sign is arbitrary, so "A fraction" would be an unfounded claim
            out["genome_A_fraction"] = None
            out["positive_E1_fraction"] = positive_fraction
        out["bins_used"] = int(len(vec))
    return out


def virtual_4c(
    file: str | None = None,
    viewpoint: str = "chr17:63,000,000-63,100,000",
    resolution: int | None = None,
) -> dict:
    """A virtual-4C profile: balanced contact frequency of one viewpoint with everything else."""
    path = resolve_input_path(file)
    clr = open_matrix(path, resolution, default=10_000)
    if not _weights_present(clr):
        raise AnalysisError(
            "This file has no ICE weights; virtual 4C reports balanced contact frequencies. "
            "Balance it first with `cooler balance`."
        )
    chrom, start, end = parse_region_checked(clr, viewpoint)
    if end - start < int(clr.binsize):
        # pad to a whole bin, but never past the chromosome's own end
        chrom_end = int(clr.chromsizes[chrom])
        end = min(start + int(clr.binsize), chrom_end)
        start = max(0, min(start, end - 1))
    _check_viewpoint_not_filtered(clr, chrom, start, end)
    prof = virtual4c(clr, f"{chrom}:{start}-{end}")
    prof = prof[prof["chrom"] == chrom].reset_index(drop=True)
    vals = prof["balanced"].to_numpy()
    pos = prof["start"].to_numpy()
    center = (start + end) // 2
    dist = np.abs(pos + int(clr.binsize) // 2 - center)
    bands = [(0, 100_000), (100_000, 1_000_000), (1_000_000, 5_000_000), (5_000_000, 10_000_000)]

    def _bound(v: int) -> str:
        return f"{v // 1_000_000}Mb" if v >= 1_000_000 else f"{v // 1000}kb"

    def _band_label(lo: int, hi: int) -> str:
        return f"{_bound(lo)}-{_bound(hi)}"

    # A NaN here means one of two different things, and averaging them together is
    # wrong: an ICE-filtered bin carries no measurement (exclude it), while a mappable
    # bin with no contacts is a genuine zero (include it, or the decay looks flatter
    # than it is).
    all_bins = clr.bins()[:]
    weights = all_bins["weight"].to_numpy()
    mappable = ~np.isnan(weights)
    # the viewpoint's own bins are masked by cooltools deliberately (self-contact), so
    # they are neither a measurement nor a genuine zero - exclude them entirely
    own = (
        (all_bins["chrom"].astype(str).to_numpy() == chrom)
        & (all_bins["end"].to_numpy() > start)
        & (all_bins["start"].to_numpy() < end)
    )
    usable = np.where(np.isnan(vals) & mappable, 0.0, vals)

    band_means = {}
    band_bins = {}
    for lo, hi in bands:
        in_band = (dist >= lo) & (dist < hi) & mappable & ~own
        sel = usable[in_band]
        sel = sel[np.isfinite(sel)]
        if sel.size:
            label = _band_label(lo, hi)
            band_means[label] = _finite(sel.mean())
            band_bins[label] = int(sel.size)
    finite = np.isfinite(vals)
    n = int(finite.sum())
    stride = max(1, n // PROFILE_POINT_CAP)
    idx = np.where(finite)[0][::stride]
    return {
        "viewpoint": f"{chrom}:{start:,}-{end:,}",
        "resolution_used": int(clr.binsize),
        "profile_points": [
            {"start": int(pos[i]), "balanced": _finite(vals[i])} for i in idx
        ],
        "profile_note": (
            f"{n} bins carry a measured contact value; downsampled by taking every "
            f"{stride}th point for transport. cooltools masks the viewpoint's own bin, so "
            "it reads null. Bins that are mappable but share no contacts with the "
            "viewpoint are genuine zeros and are counted as zero in the band means; "
            "ICE-filtered bins carry no measurement and are excluded from them."
        ),
        "distance_band_means": band_means,
        "distance_band_bins": band_bins,
        "distance_bands_cover_bp": [0, bands[-1][1]],
        "coverage_note": (
            f"Band means summarise separations up to {bands[-1][1] // 1_000_000} Mb; the "
            f"profile itself spans the whole chromosome "
            f"({int(clr.chromsizes[chrom]) // 1_000_000} Mb here)."
        ),
        "balanced": True,
        "method": "cooltools.virtual4c (balanced row extraction at the viewpoint)",
    }


def expected_observed(
    file: str | None = None,
    region: str = "chr17:50,000,000-52,500,000",
    resolution: int | None = None,
) -> dict:
    """Distance-expected contact curve and the observed/expected matrix for a region."""
    path = resolve_input_path(file)
    clr = open_matrix(path, resolution, default=100_000)
    if not _weights_present(clr):
        raise AnalysisError(
            "This file has no ICE weights; expected/observed uses balanced values. "
            "Balance it first with `cooler balance`."
        )
    chrom, start, end = parse_region_checked(clr, region)
    binsize = int(clr.binsize)
    # the arm view is resolution-independent, so it applies at every resolution - the
    # region model must not change silently with bin size
    view = load_arms_view() if is_demo(path) else None

    scope_name = chrom
    if view is not None:
        row = view[(view["chrom"] == chrom) & (view["start"] <= start) & (view["end"] >= end)]
        if not len(row):
            arms = "; ".join(
                f"{r['name']} {int(r['start']):,}-{int(r['end']):,}"
                for _, r in view[view["chrom"] == chrom].iterrows()
            )
            raise AnalysisError(
                f"{chrom}:{start:,}-{end:,} is not contained in a single chromosome arm, "
                "so no single expected curve applies (the arms are separated by the "
                f"ICE-filtered centromeric gap). Arms here: {arms}. Request a region "
                "inside one arm."
            )
        scope_name = str(row.iloc[0]["name"])

    exp = expected_cis(clr, view_df=view, ignore_diags=IGNORE_DIAGS, nproc=1)
    exp = exp[exp["region1"] == exp["region2"]]
    # cooltools masks the first `ignore_diags` diagonals: `balanced.avg` is NaN there, and
    # the SMOOTHED column carries a smoother extrapolation rather than a measurement.
    measured = exp[exp["balanced.avg"].notna()]
    # the curve belongs to THIS region's scope, never the whole file - the per-region
    # column, not the cross-region `.agg` aggregate
    scoped = measured[measured["region1"] == scope_name]
    curve = (
        scoped[["dist_bp", "balanced.avg.smoothed"]]
        .rename(columns={"balanced.avg.smoothed": "expected"})
        .dropna()
        .sort_values("dist_bp")
    )

    fit_lo = max(100_000, IGNORE_DIAGS * binsize)
    sl = curve[(curve["dist_bp"] >= fit_lo) & (curve["dist_bp"] <= 10_000_000)]
    slope: float | None = None
    fit_range: list[int] | None = None
    if len(sl) > 3:
        slope = _finite(np.polyfit(np.log10(sl["dist_bp"]), np.log10(sl["expected"]), 1)[0])
        fit_range = [int(sl["dist_bp"].min()), int(sl["dist_bp"].max())]

    n_bins = (end - start) // binsize
    out: dict = {
        "region": f"{chrom}:{start:,}-{end:,}",
        "resolution_used": binsize,
        "view": (
            f"{scope_name} (bundled arm view)"
            if view is not None
            else f"{scope_name} (whole chromosome)"
        ),
        "curve_scope": (
            f"the expected curve and slope describe all of {scope_name}, not only the "
            "requested region"
        ),
        "ps_slope": slope,
        "ps_fit_range_bp": fit_range,
        "ignored_diagonals": IGNORE_DIAGS,
        "expected_curve_points": [
            {"dist_bp": int(r.dist_bp), "expected": _finite(r.expected)}
            for r in curve.iloc[:: max(1, len(curve) // 100)].itertuples()
        ],
        "balanced": True,
        "method": (
            "cooltools.expected_cis (smoothed distance-expected, first "
            f"{IGNORE_DIAGS} diagonals unmeasured); O/E = balanced observed over expected "
            "at each separation"
        ),
    }
    if n_bins <= MATRIX_BIN_CAP:
        obs = clr.matrix(balance=True).fetch(f"{chrom}:{start}-{end}")
        exp_by_diag = scoped.groupby("dist")["balanced.avg.smoothed"].mean()
        oe = np.full_like(obs, np.nan)
        for d in range(obs.shape[0]):
            e = exp_by_diag.get(d, np.nan)
            if np.isfinite(e) and e > 0:
                diag = np.diagonal(obs, offset=d) / e
                for k, v in enumerate(diag):
                    oe[k, k + d] = v
                    oe[k + d, k] = v
        out["oe_matrix"] = [[_finite(v) for v in rowv] for rowv in oe]
        if not np.isfinite(oe).any():
            out["note"] = (
                f"Every cell is null: {chrom}:{start:,}-{end:,} has no balanced signal at "
                f"{binsize:,} bp (its bins are ICE-filtered, e.g. centromeric or "
                "low-mappability). Ask for a region with mappable bins."
            )
        else:
            out["note"] = (
                f"The first {IGNORE_DIAGS} diagonals read null by construction: cooltools "
                "does not measure expected there, so no honest observed/expected ratio "
                f"exists for separations under {IGNORE_DIAGS * binsize:,} bp."
            )
    else:
        out["oe_matrix"] = None
        out["note"] = (
            f"Region spans {n_bins} bins (cap {MATRIX_BIN_CAP} per side for the O/E matrix); "
            "curve and slope only. Narrow the region or use a coarser resolution."
        )
    return out
