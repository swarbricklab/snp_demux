# ===================
# present_absent_screen.py
# -------------------
# Per-capture, clustering-INDEPENDENT "is this donor present or absent?" screen.
# Productionises a validated experiment (private-site present/absent test). The
# reference implementation correctly calls a genuine-but-missed donor PRESENT and
# clean dropouts ABSENT; this script reuses its EXACT private-site definition and
# verdict thresholds so results match.
#
# Method (per target donor d in the capture pool):
#   1. PRIVATE sites of d = biallelic panel sites where d is an alt-carrier
#      (het 0/1 or hom-alt 1/1) AND EVERY OTHER pool donor is a confident 0/0
#      (hom-ref). At such a site an ALT read can only come from d, or from
#      ambient/error. Missing (./.) in any other donor -> NOT confident hom-ref
#      -> the site is excluded (conservative). Multi-allelic sites are excluded.
#   2. Intersect d's private sites with the cellSNP base sites -> matrix rows.
#   3. Per cell: n_hit = # covered private sites carrying an alt read (AD>=1).
#   4. Ambient null p = sum(AD)/sum(DP) over d's private rows, restricted to cells
#      the demux did NOT assign to d. Pooling over ALL cells (including d's own)
#      lets a DOMINANT donor's reads inflate its own background so its null tail
#      matches its observed tail and it self-masks -- the bug this fixes. Using
#      only not-d cells gives the true ambient/error rate. Falls back to the
#      pooled estimate when the assignment table is missing/empty.
#   5. Null model: in a cell that does NOT contain d, each covered private site
#      carries an alt with prob p, so n_hit ~ Binom(n_covered_in_cell, p).
#      null_exp_ge3 = sum over cells of P(Binom >= 3). Compare the observed tail
#      (cells with n_hit >= 3 / >= 5) against this null.
#
# Verdict (thresholds are module constants below, overridable via params):
#   * present (assigned)    : the demux confidently assigns >= present_min_assigned
#                             cells to d -> present outright, regardless of the tail
#                             (dominant donors need no private-site significance).
#   Otherwise the corrected-ambient private-site test decides:
#   * untestable            : too few of d's private sites covered to decide.
#   * absent                : no tail -- zero cells reach >=5 private-alt sites and
#                             the >=3 tail is at/near the ambient null.
#   * present               : a clear tail well above null -- cells_ge5 above a
#                             small floor and the >=3 tail exceeds the null.
#   * present_underdetermined : testable, but neither cleanly present nor absent.
#
# NOTE: snakemake-driven script (reads the `snakemake` object), mirroring
# report_distinctiveness.py / hierarchical_consensus.py. Deterministic. Robust to
# missing/empty inputs -> emits an empty-but-valid TSV rather than crashing.
#
# Env: envs/embedding.yaml (numpy + pandas + scipy). scipy IS used, for both the
# MatrixMarket read (scipy.io.mmread) and the exact binomial upper tail
# (scipy.stats.binom) -- matching the validated reference implementation.
#
# cellSNP matrix orientation: sites x cells (rows = SNPs from base.vcf, cols =
# barcodes) -- same convention as report_distinctiveness.py / the reference.
# ===================

import os
import gzip
import logging as log

import numpy as np
import pandas as pd
from scipy.io import mmread
from scipy.stats import binom

log.basicConfig(
    filename=snakemake.log[0],
    level=log.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s",
)

# ---- thresholds (module constants; documented; overridable via params) -------
# Reproduce the validated reference calls:
#   DONOR_D -> present (cells_ge5=26, cells_ge3=69 vs null 29.25)
#   DONOR_J -> absent  (cells_ge5=0,  cells_ge3=1  vs null 0.07)
#   DONOR_P -> present (cells_ge5=5,  cells_ge3=32 vs null 0.21)
MIN_COVERED_SITES = 200    # < this many covered private sites -> untestable
PRESENT_GE5_FLOOR = 5      # cells with >=5 private-alt sites must reach this to call present
PRESENT_TAIL_RATIO = 2.0   # AND cells_ge3 must exceed this * null_exp_ge3
ABSENT_TAIL_RATIO = 3.0    # absent needs cells_ge3 <= max(ABSENT_GE3_ABS, this*null)
ABSENT_GE3_ABS = 3         # absolute floor on the >=3 tail tolerated for "absent"
PRESENT_MIN_ASSIGNED = 15  # demux cells assigned to d at/above which d is present outright

# Demux statuses that count as "this cell IS donor d" (vs ambient/other).
ASSIGNED_STATUSES = {"singlet", "recovered"}

OUT_COLS = [
    "capture_id", "donor", "verdict", "n_assigned",
    "n_private", "n_covered", "ambient_p", "ambient_p_pooled",
    "cells_ge1", "cells_ge3", "cells_ge5", "null_exp_ge3",
]

# Non-`present` verdicts sort FIRST so problems surface at the top of the table.
_VERDICT_RANK = {
    "absent": 0,
    "present_underdetermined": 1,
    "untestable": 2,
    "present": 3,
}


def _p(key, default):
    """Read a params.present_absent override, falling back to the module const."""
    try:
        pa = snakemake.params["present_absent"]
    except (KeyError, TypeError):
        return default
    if isinstance(pa, dict) and key in pa and pa[key] is not None:
        return pa[key]
    return default


MIN_COVERED_SITES = int(_p("min_covered_sites", MIN_COVERED_SITES))
PRESENT_GE5_FLOOR = int(_p("present_ge5_floor", PRESENT_GE5_FLOOR))
PRESENT_TAIL_RATIO = float(_p("present_tail_ratio", PRESENT_TAIL_RATIO))
ABSENT_TAIL_RATIO = float(_p("absent_tail_ratio", ABSENT_TAIL_RATIO))
ABSENT_GE3_ABS = int(_p("absent_ge3_abs", ABSENT_GE3_ABS))
PRESENT_MIN_ASSIGNED = int(_p("present_min_assigned", PRESENT_MIN_ASSIGNED))


def _in(key):
    """Optional named input -> path or None (robust when a key is absent)."""
    try:
        return snakemake.input[key]
    except (KeyError, TypeError):
        return None


f_vcf = snakemake.input["known_genotypes_vcf"]
f_base = snakemake.input["base_vcf"]
f_ad = snakemake.input["ad_mtx"]
f_dp = snakemake.input["dp_mtx"]
f_samples = _in("samples")
f_hier = _in("hier")
f_tsv = snakemake.output["tsv"]
capture = snakemake.wildcards["capture_id"]


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _exists(path):
    return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0


def _norm(chrom, pos):
    return (str(chrom).lower().replace("chr", ""), str(pos))


def _write(df):
    df.to_csv(f_tsv, sep="\t", index=False)
    log.info("wrote %s (%d donor rows)", f_tsv, len(df))


def _empty(note):
    """Emit an empty-but-valid TSV rather than crashing."""
    log.warning("present_absent: %s", note)
    _write(pd.DataFrame(columns=OUT_COLS))


# ---- 1. parse the per-donor genotypes into private-site sets ----------------
def classify(gt):
    """'homref' (all 0), 'alt' (any 1), or 'missing' (any '.')."""
    g = gt.split(":", 1)[0].replace("|", "/")
    a = g.split("/")
    if any(x == "." for x in a) or len(a) < 2:
        return "missing"
    if all(x == "0" for x in a):
        return "homref"
    if any(x == "1" for x in a):
        return "alt"
    return "missing"  # multi-allelic leftover -> treat cautiously


def parse_private(path):
    """Return (samples, {donor: set((chrom,pos))}) of biallelic private sites."""
    samples = []
    private = {}
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                private = {s: set() for s in samples}
                continue
            if not samples:
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 + len(samples):
                continue
            ref, alt = f[3], f[4]
            # biallelic only: single REF base, single ALT allele
            if "," in alt or len(ref) != 1 or len(alt) != 1:
                continue
            cls = [classify(g) for g in f[9:9 + len(samples)]]
            for i, s in enumerate(samples):
                if cls[i] != "alt":
                    continue
                if all(cls[j] == "homref" for j in range(len(samples)) if j != i):
                    private[s].add(_norm(f[0], f[1]))
    return samples, private


# ---- 2. cellSNP base coords -> row index ------------------------------------
def base_index(path):
    idx = {}
    with _open(path) as fh:
        i = 0
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            idx[_norm(f[0], f[1])] = i
            i += 1
    return idx


# ---- 3. demux per-cell assignment -> per-donor assigned barcodes ------------
def read_barcodes(path):
    """cellSNP.samples.tsv: one barcode per line, in AD/DP matrix COLUMN order."""
    bcs = []
    with _open(path) as fh:
        for line in fh:
            b = line.strip()
            if b:
                bcs.append(b)
    return bcs


def read_assigned(path):
    """{donor: set(barcodes)} the demux confidently assigned (status in
    ASSIGNED_STATUSES). Empty dict if the table is missing/empty/malformed."""
    if not _exists(path):
        return {}
    try:
        h = pd.read_csv(path, sep="\t", dtype=str)
    except Exception as e:  # noqa: BLE001 - stay robust, never crash the screen
        log.warning("could not read hierarchical consensus %s: %s", path, e)
        return {}
    if not {"barcode", "donor", "status"}.issubset(h.columns):
        log.warning("hier table missing barcode/donor/status columns: %s", list(h.columns))
        return {}
    keep = h[h["status"].str.lower().isin(ASSIGNED_STATUSES)]
    assigned = {}
    for donor, grp in keep.groupby("donor"):
        d = str(donor)
        if d and d.lower() not in ("nan", "doublet", "unassigned"):
            assigned[d] = set(grp["barcode"].astype(str))
    return assigned


# ---- 4. verdict --------------------------------------------------------------
def decide(n_assigned, n_covered, cells_ge3, cells_ge5, null_exp_ge3):
    # A donor the demux confidently assigns is present, tail or no tail.
    if n_assigned >= PRESENT_MIN_ASSIGNED:
        return "present"
    if n_covered < MIN_COVERED_SITES:
        return "untestable"
    if cells_ge5 >= PRESENT_GE5_FLOOR and cells_ge3 > PRESENT_TAIL_RATIO * null_exp_ge3:
        return "present"
    if cells_ge5 == 0 and cells_ge3 <= max(ABSENT_GE3_ABS, ABSENT_TAIL_RATIO * null_exp_ge3):
        return "absent"
    return "present_underdetermined"


def main():
    if not _exists(f_vcf):
        _empty("known_genotypes.vcf missing/empty; no donors screened")
        return
    samples, private = parse_private(f_vcf)
    if not samples:
        _empty("no donors in known_genotypes.vcf")
        return

    have_mtx = _exists(f_ad) and _exists(f_dp) and _exists(f_base)
    if have_mtx:
        idx = base_index(f_base)
        AD = mmread(f_ad).tocsr()
        DP = mmread(f_dp).tocsr()
        ncells = AD.shape[1]
        log.info("pileup: %d base sites, AD shape %s", len(idx), AD.shape)
    else:
        idx, AD, DP, ncells = {}, None, None, 0
        log.warning("cellSNP AD/DP/base missing; every donor -> untestable")

    # Align the demux per-cell assignment onto the AD/DP matrix COLUMNS via the
    # cellSNP barcode order. When barcodes/hier are missing, col_barcodes is []
    # and every donor's "not-assigned" mask is all-True -> corrected == pooled.
    col_barcodes = read_barcodes(f_samples) if _exists(f_samples) else []
    if col_barcodes and have_mtx and len(col_barcodes) != ncells:
        log.warning("cellSNP.samples.tsv (%d) != AD columns (%d); ambient not corrected",
                    len(col_barcodes), ncells)
        col_barcodes = []
    assigned = read_assigned(f_hier)
    col_bc_arr = np.array(col_barcodes, dtype=object) if col_barcodes else None
    log.info("demux-assigned donors: %s",
             {d: len(v) for d, v in assigned.items()} or "none")

    rows_out = []
    for donor in samples:
        psites = private.get(donor, set())
        n_private = len(psites)
        rows = [idx[k] for k in psites if k in idx] if idx else []
        n_covered = len(rows)
        assigned_set = assigned.get(donor, set())
        n_assigned = len(assigned_set)

        ambient_p = ambient_p_pooled = 0.0
        cells_ge1 = cells_ge3 = cells_ge5 = 0
        null_exp_ge3 = 0.0
        if have_mtx and rows:
            ad = AD[rows, :]
            dp = DP[rows, :]
            n_hit = np.asarray((ad >= 1).sum(axis=0)).ravel()      # per cell: #private-alt sites
            n_cov_cell = np.asarray((dp >= 1).sum(axis=0)).ravel()  # per cell: #private sites covered
            cells_ge1 = int((n_hit >= 1).sum())
            cells_ge3 = int((n_hit >= 3).sum())
            cells_ge5 = int((n_hit >= 5).sum())

            # Pooled ambient (all cells) -- reported for transparency.
            tot_ad = float(ad.sum())
            tot_dp = float(dp.sum())
            ambient_p_pooled = (tot_ad / tot_dp) if tot_dp > 0 else 0.0

            # Corrected ambient: restrict to cells the demux did NOT assign to d,
            # so d's own reads cannot inflate its own background. Fall back to the
            # pooled estimate if we cannot resolve not-d cells (no barcodes/hier,
            # or no coverage among not-d cells).
            ambient_p = ambient_p_pooled
            if col_bc_arr is not None and assigned_set:
                not_d_cols = np.nonzero(~np.isin(col_bc_arr, list(assigned_set)))[0]
                if not_d_cols.size:
                    dp_nd = dp[:, not_d_cols]
                    tot_dp_nd = float(dp_nd.sum())
                    if tot_dp_nd > 0:
                        ambient_p = float(ad[:, not_d_cols].sum()) / tot_dp_nd

            if ambient_p > 0:
                n = n_cov_cell.astype(int)
                pv = 1.0 - binom.cdf(2, n, ambient_p)   # P(Binom(n, p) >= 3)
                pv = np.where(n >= 3, pv, 0.0)
                null_exp_ge3 = float(np.nansum(pv))

        verdict = decide(n_assigned, n_covered, cells_ge3, cells_ge5, null_exp_ge3)
        rows_out.append({
            "capture_id": capture,
            "donor": donor,
            "verdict": verdict,
            "n_assigned": n_assigned,
            "n_private": n_private,
            "n_covered": n_covered,
            "ambient_p": round(ambient_p, 5),
            "ambient_p_pooled": round(ambient_p_pooled, 5),
            "cells_ge1": cells_ge1,
            "cells_ge3": cells_ge3,
            "cells_ge5": cells_ge5,
            "null_exp_ge3": round(null_exp_ge3, 3),
        })
        log.info(
            "donor=%s verdict=%s n_assigned=%d n_private=%d n_covered=%d "
            "ambient_p=%.5f (pooled=%.5f) cells_ge1=%d cells_ge3=%d cells_ge5=%d null_exp_ge3=%.3f",
            donor, verdict, n_assigned, n_private, n_covered, ambient_p,
            ambient_p_pooled, cells_ge1, cells_ge3, cells_ge5, null_exp_ge3,
        )

    out = pd.DataFrame(rows_out, columns=OUT_COLS)
    # Non-`present` verdicts first (problems surface), then by donor for stability.
    out["_rank"] = out["verdict"].map(_VERDICT_RANK).fillna(9).astype(int)
    out = out.sort_values(["_rank", "donor"], kind="mergesort").drop(columns="_rank")
    _write(out)


if __name__ == "__main__":
    main()
