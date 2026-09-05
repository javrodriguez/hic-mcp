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
    load_track_file,
    load_view_file,
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
VIEWPOINT_MAX_BINS = 10  # a virtual-4C viewpoint is an anchor, not a region
COMPARTMENT_MEMORY_CAP_GB = 4.0  # eigs_cis densifies per region; refuse before the OOM
DENSE_FETCH_CAP_GB = 0.5  # above this, contacts_at_locus answers from sparse pixels instead


def _ordinal(n: int) -> str:
    """2 -> '2nd'. The naive f"{n}th" produced agent-facing text reading 'every 2th'."""
    if 10 <= n % 100 <= 20:
        return f"{n}th"
    return f"{n}{ {1: 'st', 2: 'nd', 3: 'rd'}.get(n % 10, 'th') }".replace(" ", "")


def _finite(x: float) -> float | None:
    """Round to 6 SIGNIFICANT figures, not 6 decimals.

    Contact frequencies span decades; fixed-decimal rounding would flatten the tail of
    a P(s) curve into a run of identical constants.
    """
    return float(f"{float(x):.6g}") if np.isfinite(x) else None


def _dense_cells_from_pixels(clr: Cooler, balance: bool, r1: str, r2: str, symmetric: bool):
    """The pixels that fill the dense r1 x r2 rectangle, expanded from the stored triangle.

    cooler stores one triangle, and a stored pixel (i,j) fills two dense cells: (i,j) and
    (j,i). Fetching one orientation therefore under-reports any rectangle not wholly above
    the diagonal - it returned nothing for a descending order and undercounted overlapping
    regions by 16%. For a SYMMETRIC request the caller wants the upper triangle once for
    the contact total, but the full square for cell counts and means, so both are returned.
    """
    import pandas as pd

    forward = clr.matrix(balance=balance, as_pixels=True).fetch(r1, r2)
    if symmetric:
        upper = forward[forward["bin1_id"] <= forward["bin2_id"]]
        return upper, upper[upper["bin1_id"] == upper["bin2_id"]], None
    reverse = clr.matrix(balance=balance, as_pixels=True).fetch(r2, r1)
    shared = set(
        zip(
            reverse.loc[reverse["bin1_id"] == reverse["bin2_id"], "bin1_id"],
            reverse.loc[reverse["bin1_id"] == reverse["bin2_id"], "bin2_id"],
            strict=True,
        )
    )
    if shared and len(forward):
        dup = np.array(
            [
                (i, j) in shared
                for i, j in zip(forward["bin1_id"], forward["bin2_id"], strict=True)
            ],
            dtype=bool,
        )
        forward = forward[~dup]
    return pd.concat([forward, reverse], ignore_index=True), None, True


def _resolve_view(path, view: str | None, clr: Cooler | None = None):
    """A user-supplied view wins; the bundled arm view applies only to the demo file."""
    if view is not None:
        view_df = load_view_file(view)
        if clr is not None:
            _validate_view(view_df, clr)
        # cooltools needs a sorted viewframe; region order carries no meaning here
        # (each region is analysed independently), so sorting is a normalisation, not
        # a change of answer - and refusing over row order would be pedantry.
        return view_df.sort_values(["chrom", "start"]).reset_index(drop=True)
    return load_arms_view() if is_demo(path) else None


def _validate_view(view_df, clr: Cooler) -> None:
    """Reject a malformed view with a message the caller can act on.

    cooltools requires a proper viewframe - unique names, no overlaps, in bounds - and
    rejects anything else with one opaque sentence. Every failure here is a fixable
    input, so each is named specifically rather than surfacing as an internal bug.
    """
    sizes = {str(c): int(s) for c, s in clr.chromsizes.items()}
    unknown = sorted(set(view_df["chrom"]) - set(sizes))
    if unknown:
        raise AnalysisError(
            f"The view names {', '.join(unknown)}, which this file does not contain "
            f"(it has: {', '.join(list(sizes)[:10])}). Chromosome names must match the "
            "file exactly - '17' and 'chr17' are different names."
        )
    dupes = sorted(view_df["name"][view_df["name"].duplicated()].unique())
    if dupes:
        raise AnalysisError(
            f"The view reuses the region name(s) {', '.join(map(str, dupes))}. Every "
            "region needs its own name, since results are reported per region."
        )
    for _, r in view_df.iterrows():
        if int(r["start"]) < 0 or int(r["end"]) > sizes[str(r["chrom"])]:
            raise AnalysisError(
                f"View region {r['name']} ({r['chrom']}:{int(r['start']):,}-"
                f"{int(r['end']):,}) runs past the end of {r['chrom']}, which is "
                f"{sizes[str(r['chrom'])]:,} bp in this file."
            )
        if int(r["end"]) <= int(r["start"]):
            raise AnalysisError(
                f"View region {r['name']} has zero or negative width "
                f"({int(r['start']):,}-{int(r['end']):,})."
            )
    ordered = view_df.sort_values(["chrom", "start"])
    for chrom, grp in ordered.groupby("chrom"):
        ends = grp["end"].to_numpy()
        starts = grp["start"].to_numpy()
        if len(grp) > 1 and (starts[1:] < ends[:-1]).any():
            raise AnalysisError(
                f"View regions on {chrom} overlap. Each position must belong to at most "
                "one region, or a result would be attributed to two of them."
            )


def _check_region_in_view(view_df, chrom: str, start: int, end: int, clr=None) -> None:
    """A region outside the supplied view is a fixable request, not bad data.

    Without this the tools blame ICE filtering for bins that carry perfectly good
    weights, which sends the caller looking for a data problem that does not exist.
    """
    if view_df is None:
        return
    covering = view_df[
        (view_df["chrom"] == chrom) & (view_df["start"] <= start) & (view_df["end"] >= end)
    ]
    if len(covering):
        return
    same_chrom = view_df[view_df["chrom"] == chrom]
    if len(same_chrom):
        spans = "; ".join(
            f"{r['name']} {int(r['start']):,}-{int(r['end']):,}" for _, r in same_chrom.iterrows()
        )
        where = f"Regions covering {chrom}: {spans}."
    else:
        where = f"The view covers no part of {chrom} at all."
    # if the gap it falls in is also unmappable, say that too - for an arm view the
    # space between regions IS the centromere, and the caller deserves both facts
    why = ""
    if clr is not None and _weights_present(clr):
        try:
            w = clr.bins().fetch(f"{chrom}:{start}-{end}")["weight"]
            if w.isna().all():
                why = (
                    " Its bins are also entirely ICE-filtered (no balanced signal exists "
                    "there) - for an arm view, the space between regions is the "
                    "centromeric gap."
                )
        except ValueError:
            pass
    raise AnalysisError(
        f"{chrom}:{start:,}-{end:,} is not contained in any single region of the view "
        f"being used, so no result applies to it.{why} {where} Ask for a region inside "
        "one of them, or supply a view that covers this one."
    )


def _view_label(view_df, view_arg: str | None, path) -> str:
    if view_df is None:
        return "whole chromosomes (no view supplied; pass view= to partition, e.g. by arm)"
    if view_arg is not None:
        return f"supplied view ({len(view_df)} regions)"
    return "chr17 p/q arms (bundled)"


def _balanced_cis_only(clr: Cooler) -> bool:
    """True when the ICE solution was solved per-chromosome (cooler's cis_only flag).

    Trans values under a cis-only solution are normalised by weights that were never
    fitted for them, so reporting them as "balanced" would overstate what was computed.
    """
    try:
        import h5py

        uri = str(clr.filename)
        group = clr.root.rstrip("/") if getattr(clr, "root", None) else ""
        with h5py.File(uri, "r") as f:
            node = f[f"{group}/bins/weight"] if group else f["bins/weight"]
            return bool(node.attrs.get("cis_only", False))
    except Exception:  # noqa: BLE001 - absence of the flag is simply "not known"
        return False


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
        # NOT clr.chromsizes.name - that is the pandas Series name ("length"), never an assembly
        "assembly": clr.info.get("genome-assembly") or "unknown",
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
    if use_weights and c2 != chrom and _balanced_cis_only(clr):
        # weights fitted per chromosome were never fitted for this trans block
        use_weights = False
        out_trans_note = (
            "This file was ICE-balanced cis-only (per chromosome), so no balanced values "
            f"exist for the {chrom} x {c2} block - its weights were never fitted for trans "
            "contacts. Raw counts are reported instead."
        )
    else:
        out_trans_note = None
    # Project the dense cost BEFORE fetching anything: every sibling that densifies guards
    # first, and this one fetched a whole chromosome at 10 kb (1.5 GB) to return a few
    # scalars. Above the cap the same statistics come from the sparse pixel table instead,
    # so the tool answers rather than refuses.
    n1 = -(-(end - start) // resolution)
    n2 = -(-(e2 - s2) // resolution)
    dense_gb = 8 * n1 * n2 / 1e9
    sparse_mode = dense_gb > DENSE_FETCH_CAP_GB or max(n1, n2) > MATRIX_BIN_CAP * 4
    if region2 is None:
        counting = "each contact counted once (upper triangle incl. diagonal)"
    else:
        counting = "each pixel of the region x region2 block counted once"
        if c2 == chrom and s2 < end and start < e2:
            counting += " - the regions overlap, so contacts inside the overlap appear twice"
    if sparse_mode:
        # Same STATISTICS as the dense road, computed from stored pixels. cooler stores
        # only the upper triangle, so a symmetric square's cell counts and sums must be
        # mirrored - counting stored pixels against the full square would halve
        # nonzero_fraction, and averaging over stored pixels alone would inflate the mean.
        # cooler stores ONE triangle: a stored pixel (i,j) fills two dense cells, (i,j)
        # and (j,i). A single as_pixels fetch therefore under-reports any rectangle that
        # is not wholly above the diagonal - it returned nothing at all for a descending
        # order, and undercounted overlapping regions by 16%. Both orientations are
        # fetched and combined, minus the diagonal they share, which is exact for every
        # shape: disjoint, overlapping, identical, and either order.
        symmetric = region2 is None
        pix, diag, _ = _dense_cells_from_pixels(clr, False, r1, r2, symmetric)
        raw_max = int(pix["count"].max()) if len(pix) else 0
        if symmetric:
            # contact total: the upper triangle once - the file's own convention.
            # cell count: the full square, so mirror the off-diagonal pixels.
            raw_sum = int(pix["count"].sum()) if len(pix) else 0
            nonzero_cells = 2 * len(pix) - len(diag)
        else:
            raw_sum = int(pix["count"].sum()) if len(pix) else 0
            nonzero_cells = len(pix)
        nonzero_fraction = round(nonzero_cells / max(1, n1 * n2), 4)
        raw = None
    else:
        raw = clr.matrix(balance=False).fetch(r1, r2)
        raw_once = np.triu(np.nan_to_num(raw)) if region2 is None else np.nan_to_num(raw)
        raw_sum = int(raw_once.sum())
        raw_max = int(np.nanmax(raw)) if raw.size else 0
        nonzero_fraction = round(float((raw > 0).mean()), 4) if raw.size else 0.0
    out: dict = {
        "region": r1,
        "region2": r2 if region2 is not None else None,
        "resolution_used": resolution,
        "shape_bins": [n1, n2],
        "raw_contacts_sum": raw_sum,
        "counting": counting,
        "raw_contacts_max": raw_max,
        "nonzero_fraction": nonzero_fraction,
        "balanced": use_weights,
        "method": "cooler.Cooler.matrix().fetch (raw counts; ICE-balanced values when available)",
    }
    if out_trans_note:
        out["note"] = out_trans_note
    if not use_weights and out_trans_note is None:
        # the tool description promises a matrix for small windows; when the caller asked
        # for raw counts (or the file has no weights) there is simply nothing balanced to
        # return, and silence reads as a missing result rather than a deliberate one
        out["note"] = (
            "No balanced values were requested or available, so balanced_matrix is null "
            "and only raw counts are reported"
            + ("" if balanced else " (you passed balanced=false)")
            + "."
            if not _weights_present(clr) or not balanced
            else out.get("note")
        )
    # always present, null when not computed - an absent key reads as "not applicable"
    # in one client and "missing" in another
    out.setdefault("balanced_mean", None)
    out.setdefault("balanced_max", None)
    out.setdefault("balanced_matrix", None)
    if use_weights and sparse_mode:
        bpix, bdiag, _ = _dense_cells_from_pixels(clr, True, r1, r2, symmetric)
        bvals = bpix["balanced"].to_numpy()
        finite = np.isfinite(bvals)
        if symmetric:
            dvals = bdiag["balanced"].to_numpy()
            total = 2 * float(bvals[finite].sum()) - float(dvals[np.isfinite(dvals)].sum())
        else:
            total = float(bvals[finite].sum())
        # the dense road averages over every cell that HAS a balanced value, i.e. cells
        # whose two bins are both mappable - not merely over the stored pixels
        mappable1 = int(clr.bins().fetch(r1)["weight"].notna().sum())
        mappable2 = int(clr.bins().fetch(r2)["weight"].notna().sum())
        cells = mappable1 * mappable2
        out["balanced_mean"] = _finite(total / cells) if cells else None
        out["balanced_max"] = _finite(bvals[finite].max()) if finite.any() else None
        out["balanced_matrix"] = None
        reason = (
            f"would need {dense_gb:.1f} GB as a dense matrix"
            if dense_gb > DENSE_FETCH_CAP_GB
            else f"exceeds {MATRIX_BIN_CAP * 4} bins on a side"
        )
        out["note"] = (
            f"{n1}x{n2} bins {reason}, so the statistics were computed from the sparse "
            "pixel table and the matrix itself is not returned. Narrow the region or use "
            "a coarser resolution."
        )
    elif use_weights and raw is not None and raw.size:
        bal = clr.matrix(balance=True).fetch(r1, r2)
        finite_bal = bal[np.isfinite(bal)]
        out["balanced_mean"] = _finite(finite_bal.mean()) if finite_bal.size else None
        out["balanced_max"] = _finite(finite_bal.max()) if finite_bal.size else None
        if max(bal.shape) <= MATRIX_BIN_CAP:
            out["balanced_matrix"] = [[_finite(v) for v in row] for row in bal]
            if not np.isfinite(bal).any():
                # a wall of nulls with balanced:true reads like a result; it is not
                out["note"] = (
                    f"Every balanced value is null: the bins in {r1} are ICE-filtered "
                    "(centromeric or low-mappability), so no balanced contact frequency "
                    "exists there. The raw counts above are unnormalised and not "
                    "comparable with balanced values elsewhere."
                )
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
    view: str | None = None,
    top_n: int = 10,
) -> dict:
    """Diamond insulation score and TAD-boundary calls (Crane et al. 2015 method)."""
    if top_n < 1:
        raise AnalysisError(f"top_n must be at least 1 (got {top_n}).")
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
        off = [w for w in windows if w % binsize]
        if off:
            raise AnalysisError(
                f"Window(s) {', '.join(f'{w:,}' for w in off)} bp are not multiples of the "
                f"{binsize:,} bp bin size; cooltools needs whole bins. Nearest multiples: "
                f"{', '.join(f'{max(3, round(w / binsize)) * binsize:,}' for w in off)} bp."
            )
    elif 3 * binsize <= min(CLASSIC_TAD_WINDOWS):
        # the classic windows, snapped to whole bins: 100/250/500 kb are exact at 10 kb
        # but not at 4, 8, 15, 20 or 30 kb - all ordinary Hi-C resolutions - and cooltools
        # rejects a window that is not a multiple of the bin size
        windows = sorted({max(3, round(w / binsize)) * binsize for w in CLASSIC_TAD_WINDOWS})
    else:
        windows = [n * binsize for n in COARSE_WINDOW_BINS]
    for w in windows:
        if w < 3 * binsize:
            raise AnalysisError(
                f"Window {w} bp is too small for {binsize} bp bins (a diamond needs at "
                f"least 3 bins, i.e. {3 * binsize:,} bp here). Pass larger windows_bp, or "
                "omit resolution to use this file's finest level."
            )
    # the same arm view the sibling tools use: without it cooltools median-normalises
    # the diamond score and measures boundary prominence ACROSS the centromeric gap,
    # which shifts every score by a per-arm offset and can reorder the top boundaries
    view_df = _resolve_view(path, view, clr)
    ins = insulation(clr, windows, view_df=view_df, verbose=False)
    if region is not None:
        chrom, start, end = parse_region_checked(clr, region)
        # NO single-region containment check here: insulation is per bin, so `region`
        # only filters which rows are reported. Requiring containment refused
        # region="chr17" while the same call with no region happily returned those very
        # boundaries - the gate belongs on the tools that need one curve or one
        # eigendecomposition, not on this one.
        ins = ins[(ins["chrom"] == chrom) & (ins["end"] > start) & (ins["start"] < end)]
        # "0 boundaries" is an ANSWER; a region with no insulation score at all is a
        # refusal, and the sibling tools already refuse on exactly this input
        score_col = f"log2_insulation_score_{sorted(windows)[len(windows) // 2]}"
        if (ins.empty or ins[score_col].isna().all()) and view_df is not None:
            overlapping = view_df[
                (view_df["chrom"] == chrom)
                & (view_df["end"] > start)
                & (view_df["start"] < end)
            ]
            if not len(overlapping):
                spans = "; ".join(
                    f"{r['name']} {int(r['start']):,}-{int(r['end']):,}"
                    for _, r in view_df[view_df["chrom"] == chrom].iterrows()
                ) or f"none on {chrom}"
                raise AnalysisError(
                    f"{chrom}:{start:,}-{end:,} lies outside every region of the view "
                    f"in use, so nothing was computed there. Regions on {chrom}: "
                    f"{spans}. Ask inside one of them, or drop the view."
                )
        if ins.empty or ins[score_col].isna().all():
            raise AnalysisError(
                f"No insulation score exists anywhere in {chrom}:{start:,}-{end:,} at "
                f"{binsize:,} bp - its bins are ICE-filtered (centromeric or "
                "low-mappability), so this is not a region with zero boundaries, it is a "
                "region with no measurement. Choose a region with mappable bins."
            )
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
        "view": _view_label(view_df, view, path),
        "windows_bp": windows,
        "ranked_by": (
            f"boundary_strength at the {rank_w} bp window - top_boundaries lists ONLY "
            f"boundaries called at that window (capped at top_n), so it is a subset of "
            f"boundary_counts_per_window, whose other entries count different populations"
        ),
        "boundary_counts_per_window": counts,
        "top_boundaries": boundaries,
        "balanced": True,
        "method": (
            "cooltools.insulation (diamond insulation score; Li threshold boundary calls). "
            "Returns the called boundaries and their scores, not the full per-bin "
            "insulation track."
        ),
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
    view: str | None = None,
    phasing_track: str | None = None,
) -> dict:
    """A/B compartment eigenvector (cis eigendecomposition of the observed/expected map)."""
    path = resolve_input_path(file)
    clr = open_matrix(path, resolution, default=100_000)
    if not _weights_present(clr):
        raise AnalysisError(
            "This file has no ICE weights; compartments need a balanced matrix. "
            "Balance it first with `cooler balance`."
        )
    view_arg_used = view
    view_df = _resolve_view(path, view_arg_used, clr)
    # eigs_cis densifies each view region: an n-bin region costs ~8n^2 bytes, and an
    # OOM kill gives the calling agent a dead transport rather than a readable refusal
    biggest = 0
    if view_df is not None:
        for _, r in view_df.iterrows():
            biggest = max(biggest, (int(r["end"]) - int(r["start"])) // int(clr.binsize) + 1)
    else:
        biggest = max(
            (int(size) // int(clr.binsize) + 1 for size in clr.chromsizes), default=0
        )
    projected_gb = 8 * biggest * biggest / 1e9
    if projected_gb > COMPARTMENT_MEMORY_CAP_GB:
        raise AnalysisError(
            f"Refusing to run: the largest region to decompose is {biggest:,} bins at "
            f"{int(clr.binsize):,} bp, which needs roughly {projected_gb:.1f} GB of dense "
            f"matrix (cap {COMPARTMENT_MEMORY_CAP_GB} GB). Use a coarser resolution, or "
            "pass a view of smaller regions - compartments are usually called at 100 kb "
            "or coarser anyway."
        )
    if phasing_track is not None:
        phasing = load_track_file(phasing_track)
    elif is_demo(path) and int(clr.binsize) == 100_000:
        # the bundled GC track is tied to its own 100 kb binning
        phasing = load_gc_track()
    else:
        phasing = None
    coverage_warning = None
    if phasing is not None:
        spans = (phasing["end"] - phasing["start"]).unique()
        if len(spans) and int(spans[0]) != int(clr.binsize):
            raise AnalysisError(
                f"The phasing track is binned at {int(spans[0]):,} bp but this analysis "
                f"runs at {int(clr.binsize):,} bp; cooltools needs them to match. Pass "
                f"resolution={int(spans[0])} to use the track's own binning, or supply a "
                f"track binned at {int(clr.binsize):,} bp."
            )
        # cooltools requires the track to cover EVERY region it decomposes. With no view
        # that means every chromosome in the FILE - which is the default for a user's own
        # file, so gating this check on a view left the common case unguarded.
        covered = set(phasing["chrom"].astype(str))
        if view_df is not None:
            needed = {str(c) for c in view_df["chrom"]}
            scope = "the view includes those regions"
        else:
            needed = {str(c) for c in clr.chromnames}
            scope = "with no view supplied, every chromosome in the file is decomposed"
        # a track that covers the regions but only sparsely still yields a weak
        # orientation; say so rather than letting a thin track pass silently
        # the bins actually decomposed - the view's regions when one is in force, not
        # the whole file, or a track covering its scope completely is called sparse
        if view_df is not None:
            bins_needed = 0
            for _, r in view_df.iterrows():
                bins_needed += max(
                    1, (int(r["end"]) - int(r["start"]) + int(clr.binsize) - 1) // int(clr.binsize)
                )
        else:
            bins_needed = int(clr.info["nbins"])
        if bins_needed and len(phasing) / bins_needed < 0.5:
            coverage_warning = (
                f"The phasing track has {len(phasing):,} rows against {bins_needed:,} bins "
                "in this file (under half). The orientation it gives is correspondingly "
                "weak - treat A/B calls as provisional."
            )
        else:
            coverage_warning = None
        # the track must supply values for the BINS decomposed, not merely mention the
        # chromosome: a two-row BED4 names chr17 and covers almost none of it
        span_by_chrom = phasing.groupby("chrom").apply(
            lambda g: int((g["end"] - g["start"]).sum()), include_groups=False
        )
        thin = []
        if view_df is not None:
            for _, r in view_df.iterrows():
                width = int(r["end"]) - int(r["start"])
                if span_by_chrom.get(str(r["chrom"]), 0) < 0.5 * width:
                    thin.append(str(r["name"]))
        else:
            for chrom_name, size in clr.chromsizes.items():
                if span_by_chrom.get(str(chrom_name), 0) < 0.5 * int(size):
                    thin.append(str(chrom_name))
        if thin:
            raise AnalysisError(
                f"The phasing track covers less than half of {', '.join(thin[:5])}. "
                "cooltools needs a value for essentially every bin it decomposes - "
                "supply a track binned across the whole region, or narrow the view to "
                "what the track actually covers."
            )
        missing = sorted(needed - covered)
        if missing:
            raise AnalysisError(
                f"The phasing track has no values for {', '.join(missing)}, and "
                f"{scope} - every one must be covered. Extend the track, or pass a view "
                "limited to the regions it covers."
            )
    eigvals, eigvecs = eigs_cis(clr, phasing_track=phasing, view_df=view_df, n_eigs=3)
    phasing_source = (
        "the track you supplied"
        if phasing_track is not None
        else "the bundled GC track"
    )
    sign_convention = (
        f"oriented by {phasing_source} (positive E1 = A, i.e. higher track value = A)"
        if phasing is not None
        else "UNPHASED - the sign of E1 is mathematically arbitrary. Pass phasing_track "
        "(a tab-separated chrom/start/end/value file, e.g. GC fraction) to orient it "
        "before calling A vs B"
    )
    out: dict = {
        "resolution_used": int(clr.binsize),
        "view": _view_label(view_df, view_arg_used, path),
        "sign_convention": sign_convention,
        "phasing_coverage_note": coverage_warning if phasing is not None else None,
        "eigenvalues": [
            {
                "region": str(r["name"] if "name" in r else r["chrom"]),
                "eigval1": _finite(r["eigval1"]),
                # raw eigenvalues scale with region size, so a bigger arm always looks
                # "stronger"; the share is the comparable number
                "eigval1_share_of_top3": _finite(
                    abs(r["eigval1"])
                    / sum(abs(r[f"eigval{i}"]) for i in (1, 2, 3))
                ),
            }
            for _, r in eigvals.iterrows()
        ],
        "eigenvalue_note": (
            "eigval1 is unnormalised and grows with region size, so a longer arm always "
            "looks 'stronger'. eigval1_share_of_top3 is eigval1 as a fraction of the three "
            "computed eigenvalues' absolute weight - a scale-free way to compare regions, "
            "not a share of total variance."
        ),
        "balanced": True,
        "method": "cooltools.eigs_cis (eigendecomposition of the cis observed/expected matrix)",
    }
    vec = eigvecs.dropna(subset=["E1"])
    if region is not None:
        chrom, start, end = parse_region_checked(clr, region)
        _check_region_in_view(view_df, chrom, start, end, clr)
        sub = vec[(vec["chrom"] == chrom) & (vec["end"] > start) & (vec["start"] < end)]
        if sub.empty:
            raise AnalysisError(
                f"No usable E1 bins in {region} (the region may be entirely ICE-filtered, "
                "e.g. centromeric)."
            )
        if view_df is not None:
            spanned = view_df[
                (view_df["chrom"] == chrom) & (view_df["end"] > start) & (view_df["start"] < end)
            ]
            if len(spanned) > 1:
                names = ", ".join(str(n) for n in spanned["name"])
                raise AnalysisError(
                    f"{chrom}:{start:,}-{end:,} spans {len(spanned)} view regions ({names}). "
                    "Each is eigendecomposed independently, so averaging their E1 into one "
                    "A/B call would be meaningless. Ask for a region inside one of them."
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
        # the fraction covers whatever the FILE holds - one chromosome for the demo - so
        # its name says "file", never "genome"
        if phasing is not None:
            out["A_fraction_of_file"] = positive_fraction
        else:
            # unphased: the sign is arbitrary, so an "A fraction" would be an unfounded claim
            out["A_fraction_of_file"] = None
            out["positive_E1_fraction_of_file"] = positive_fraction
        out["fraction_scope"] = (
            f"all usable bins in the file ({', '.join(str(c) for c in clr.chromnames)})"
        )
        out["bins_used"] = int(len(vec))
    return out


def virtual_4c(
    file: str | None = None,
    viewpoint: str = "chr17:63,000,000-63,100,000",
    resolution: int | None = None,
    window_bp: int | None = None,
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
    max_width = max(VIEWPOINT_MAX_BINS * int(clr.binsize), 100_000)
    if end - start > max_width:
        raise AnalysisError(
            f"Viewpoint {chrom}:{start:,}-{end:,} spans {end - start:,} bp; a virtual-4C "
            f"viewpoint is a small anchor, at most {max_width:,} bp at this resolution "
            f"({VIEWPOINT_MAX_BINS} bins). Name a locus, e.g. '{chrom}:{start:,}-"
            f"{start + int(clr.binsize):,}', or use contacts_at_locus for a wide region."
        )
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
    # the profile was filtered to the viewpoint's chromosome, so every mask below must
    # be built from the SAME rows - indexing against the whole-genome bin table breaks
    # on any multi-chromosome file
    if window_bp is not None:
        # applied BEFORE the band means, or they would summarise separations the caller
        # explicitly excluded and that appear nowhere in the returned profile
        if window_bp < int(clr.binsize):
            raise AnalysisError(
                f"window_bp {window_bp:,} is smaller than one {int(clr.binsize):,} bp bin."
            )
        in_window = dist <= window_bp
        vals = np.where(in_window, vals, np.nan)
    else:
        in_window = np.ones(vals.shape, dtype=bool)
    all_bins = clr.bins()[:]
    all_bins = all_bins[all_bins["chrom"].astype(str) == chrom].reset_index(drop=True)
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
        # a band the window excluded is ABSENT, not zero: reporting 0.0 for separations
        # the caller removed would be inventing a measurement
        in_band = (dist >= lo) & (dist < hi) & mappable & ~own & in_window
        sel = usable[in_band]
        sel = sel[np.isfinite(sel)]
        if sel.size:
            label = _band_label(lo, hi)
            band_means[label] = _finite(sel.mean())
            band_bins[label] = int(sel.size)
    # the profile points and the band means must share one convention: a mappable bin
    # with no contacts is a zero in both, never counted in one and dropped from the other
    # ...and the window applies here too: zero-filling happens before this point, so
    # without `in_window` a windowed-out mappable bin would come back as a 0.0 point
    finite = np.isfinite(usable) & mappable & ~own & in_window
    n = int(finite.sum())
    if n == 0:
        raise AnalysisError(
            f"No measured contacts anywhere on {chrom} from viewpoint {chrom}:{start:,}-"
            f"{end:,} - nothing to profile. Check the viewpoint sits in mappable bins."
        )
    stride = max(1, n // PROFILE_POINT_CAP)
    idx = np.where(finite)[0][::stride]
    return {
        "viewpoint": f"{chrom}:{start:,}-{end:,}",
        "resolution_used": int(clr.binsize),
        "window_bp": window_bp,
        "profile_points": [
            {"start": int(pos[i]), "balanced": _finite(usable[i])} for i in idx
        ],
        "profile_note": (
            f"{n} mappable bins are reported"
            + (
                f"; downsampled to every {_ordinal(stride)} point for transport"
                if stride > 1
                else ""
            )
            + ". One convention throughout: a mappable bin with no contacts is a genuine "
            "zero, in the points and in the band means alike; ICE-filtered bins and the "
            "viewpoint's own bin (masked by cooltools) carry no measurement and appear in "
            "neither."
        ),
        "distance_band_means": band_means,
        "distance_band_bins": band_bins,
        "distance_bands_cover_bp": (
            [0, min(window_bp, bands[-1][1])] if window_bp is not None else [0, bands[-1][1]]
        ),
        "coverage_note": (
            (
                f"Both the profile and the band means are limited to "
                f"{min(window_bp, bands[-1][1]):,} bp around the viewpoint"
                + (
                    f" (you asked for {window_bp:,} bp; the bands themselves stop at "
                    f"{bands[-1][1]:,} bp)"
                    if window_bp > bands[-1][1]
                    else ", as requested"
                )
                + ". A band straddling the limit keeps its full label but covers only "
                "the part inside it."
            )
            if window_bp is not None
            else (
                f"Band means summarise separations up to {bands[-1][1] // 1_000_000} Mb; "
                "the profile itself spans the whole chromosome "
                f"({int(clr.chromsizes[chrom]) // 1_000_000} Mb here)."
            )
        ),
        "balanced": True,
        "method": "cooltools.virtual4c (balanced row extraction at the viewpoint)",
    }


def expected_observed(
    file: str | None = None,
    region: str = "chr17:50,000,000-52,500,000",
    resolution: int | None = None,
    view: str | None = None,
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
    view_arg = view
    view_df = _resolve_view(path, view_arg, clr)

    scope_name = chrom
    if view_df is not None:
        # raises with the view's own regions named - never invents a centromere
        _check_region_in_view(view_df, chrom, start, end, clr)
        row = view_df[
            (view_df["chrom"] == chrom) & (view_df["start"] <= start) & (view_df["end"] >= end)
        ]
        scope_name = str(row.iloc[0]["name"])

    exp = expected_cis(clr, view_df=view_df, ignore_diags=IGNORE_DIAGS, nproc=1)
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
            f"{scope_name} ({'supplied view' if view_arg else 'bundled arm view'})"
            if view_df is not None
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
        "curve_note": (
            f"The curve covers {len(curve)} measured separations, of which {len(sl)} fall "
            f"in the fitted range; the points listed above are every "
            f"{_ordinal(max(1, len(curve) // 100))} of the curve, for transport."
            if len(curve) > 100
            else f"All {len(curve)} measured separations are listed; {len(sl)} were fitted."
        ),
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
            if obs.shape[0] <= IGNORE_DIAGS:
                # not a data problem: every separation inside the region is one of the
                # diagonals cooltools does not measure, so no ratio can exist here
                out["note"] = (
                    f"Every cell is null because {chrom}:{start:,}-{end:,} is only "
                    f"{obs.shape[0]} bin(s) wide at {binsize:,} bp, and the first "
                    f"{IGNORE_DIAGS} diagonals carry no expected measurement. The data "
                    f"here is fine - ask for a region of at least {IGNORE_DIAGS + 2} bins "
                    f"({(IGNORE_DIAGS + 2) * binsize:,} bp) to see an observed/expected "
                    "ratio, or use a finer resolution."
                )
            else:
                out["note"] = (
                    f"Every cell is null: {chrom}:{start:,}-{end:,} has no balanced signal "
                    f"at {binsize:,} bp (its bins are ICE-filtered, e.g. centromeric or "
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
