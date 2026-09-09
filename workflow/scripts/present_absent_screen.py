# ===================
# present_absent_screen.py
# -------------------
# Per-capture, clustering-INDEPENDENT PRE-FLIGHT "is this donor present or absent?"
# coarse QC gate. It reads ONLY the pool's known genotypes and the sc_pileup
# cellSNP pileup -- no demux dependency -- so it can run early/in parallel and give
# an at-a-glance "who is here / who dropped out" verdict the report surfaces at the
# top of every page. It is a FLAG only: nothing downstream prunes the panel or sets
# souporcell k from it (both tested and rejected).
#
# Method (per target donor d in the capture pool):
#   1. PRIVATE sites of d = biallelic panel sites where d is an alt-carrier
#      (het 0/1 or hom-alt 1/1) AND EVERY OTHER pool donor is a confident 0/0
#      (hom-ref). At such a site an ALT read can only come from d, or from
#      ambient/error. Missing (./.) in any other donor -> NOT confident hom-ref
#      -> the site is excluded (conservative). Multi-allelic sites are excluded.
#   2. Intersect d's private sites with the cellSNP base sites -> matrix rows.
#   3. Per cell: n_hit = # covered private sites carrying an alt read (AD>=1).
#   4. n_evidence = # cells with n_hit >= K_SITES (a cell showing several of d's
#      private alleles at once is unlikely from ambient alone -> evidence for d).
#
# Verdict (coarse; thresholds are module constants, overridable via params):
#   * untestable : too few of d's private sites covered to decide (n_covered < MIN).
#   * present    : n_evidence >= MIN_CELLS.
#   * absent     : otherwise (testable but no multi-site evidence).
#
# ambient_p (= pooled sum(AD)/sum(DP) over d's private rows across all cells) is
# reported for context only -- it is NOT used in the verdict.
#
# ROSTER-AWARENESS (params.donor_ids = the full declared pool roster):
#   * A donor omitted from the panel (no genotype in any source) is UNGENOTYPED.
#     By-elimination can recover at most ONE per pool, so >1 ungenotyped -> the
#     screen HARD-FAILS the capture at pre-flight (a clear early error; the deep
#     annotate.py raise stays as a backstop). The single ungenotyped donor gets a
#     `verdict=ungenotyped` row so the report names it (it has no private sites and
#     cannot be positively screened -- it is validated post-hoc instead).
#
# SCORE MATRIX (output score_matrix): per-barcode x genotyped-donor private-alt
#   scores -- columns `<donor>__hit` (# private sites with an ALT read) and
#   `<donor>__cov` (# with any read). The post-hoc negative-test diagnostic reads
#   this to check whether a by-elimination donor's assigned barcodes look like a
#   genuinely distinct donor (baseline hits on every known donor) or like
#   mislabeled known-donor / doublet cells (own-donor-level hits on one/multiple).
#
# NOTE: snakemake-driven script (reads the `snakemake` object), mirroring
# report_distinctiveness.py / hierarchical_consensus.py. Deterministic. Robust to
# missing/empty inputs -> emits an empty-but-valid TSV rather than crashing.
#
# Env: envs/embedding.yaml (numpy + pandas + scipy). scipy is used only for the
# MatrixMarket read (scipy.io.mmread).
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

log.basicConfig(
    filename=snakemake.log[0],
    level=log.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s",
)

# ---- thresholds (module constants; documented; overridable via params) -------
MIN_COVERED_SITES = 200  # < this many covered private sites -> "untestable"
K_SITES = 3              # a cell needs >= this many private-alt sites to be "evidence"
MIN_CELLS = 3            # >= this many evidence cells -> "present"

OUT_COLS = [
    "capture_id", "donor", "verdict", "n_evidence",
    "n_private", "n_covered", "ambient_p",
    "cells_ge1", "cells_ge2", "cells_ge3", "cells_ge5", "cells_ge10", "max_hits",
]

# Non-`present` verdicts sort FIRST so problems surface at the top of the table.
# `ungenotyped` (the single recover-by-elimination donor) sorts at the very top.
_VERDICT_RANK = {"ungenotyped": -1, "absent": 0, "untestable": 1, "present": 2}


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
K_SITES = int(_p("k_sites", K_SITES))
MIN_CELLS = int(_p("min_cells", MIN_CELLS))

f_vcf = snakemake.input["known_genotypes_vcf"]
f_base = snakemake.input["base_vcf"]
f_ad = snakemake.input["ad_mtx"]
f_dp = snakemake.input["dp_mtx"]
f_tsv = snakemake.output["tsv"]
capture = snakemake.wildcards["capture_id"]


def _opt_input(key):
    try:
        return snakemake.input[key]
    except (KeyError, TypeError):
        return None


def _opt_output(key):
    try:
        return snakemake.output[key]
    except (KeyError, TypeError):
        return None


f_barcodes = _opt_input("barcodes")     # cellSNP.samples.tsv (one barcode/line)
f_matrix = _opt_output("score_matrix")  # per-barcode x donor private-alt score

# Full declared roster (comma/space separated); genotyped donors are read from the
# panel below. Ungenotyped = roster - genotyped (donors omitted from the panel by
# build_harmonized_panel.sh because they have no genotype in any source).
try:
    _roster_raw = snakemake.params["donor_ids"]
except (KeyError, TypeError):
    _roster_raw = ""
roster = [d for d in str(_roster_raw).replace(",", " ").split() if d]


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _exists(path):
    return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0


def _norm(chrom, pos):
    return (str(chrom).lower().replace("chr", ""), str(pos))


def _write(df):
    df.to_csv(f_tsv, sep="\t", index=False)
    log.info("wrote %s (%d donor rows)", f_tsv, len(df))


def _write_matrix(df):
    """Persist the per-barcode x donor score matrix (or an empty stub). Declared as
    a rule output, so it must exist on every non-failing exit."""
    if not f_matrix:
        return
    if df is None:
        df = pd.DataFrame(columns=["barcode"])
    df.to_csv(f_matrix, sep="\t", index=False, compression="gzip")
    log.info("wrote %s (%d barcode rows, %d cols)", f_matrix, len(df), df.shape[1])


def _empty(note):
    """Emit an empty-but-valid TSV (+ matrix stub) rather than crashing."""
    log.warning("present_absent: %s", note)
    _write(pd.DataFrame(columns=OUT_COLS))
    _write_matrix(None)


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


# ---- 3. verdict --------------------------------------------------------------
def decide(n_covered, n_evidence):
    if n_covered < MIN_COVERED_SITES:
        return "untestable"
    if n_evidence >= MIN_CELLS:
        return "present"
    return "absent"


# ---- 3b. per-barcode x donor score matrix -----------------------------------
def _read_barcodes(n_expected):
    """Barcodes from cellSNP.samples.tsv (col 1), aligned to the AD/DP columns.
    Returns None if unavailable or the count disagrees with the matrix width."""
    if not (f_barcodes and _exists(f_barcodes)):
        return None
    with _open(f_barcodes) as fh:
        bc = [ln.strip().split("\t")[0] for ln in fh if ln.strip()]
    if n_expected and len(bc) != n_expected:
        log.warning("barcodes count %d != matrix cols %d; using positional ids",
                    len(bc), n_expected)
        return None
    return bc


def _emit_matrix(donors, hit_by_donor, cov_by_donor, n_cells):
    """Write barcode x genotyped-donor private-alt scores: for each donor d two
    columns `<d>__hit` (# of d's private sites with an ALT read in the barcode) and
    `<d>__cov` (# with any read). The post-hoc negative test compares a
    by-elimination donor's barcodes against these known-donor columns."""
    if not f_matrix:
        return
    if not n_cells:
        _write_matrix(None)
        return
    bc = _read_barcodes(n_cells)
    if bc is None:
        bc = ["cell_%d" % i for i in range(n_cells)]
    data = {"barcode": bc}
    for d in donors:
        data["%s__hit" % d] = hit_by_donor[d]
        data["%s__cov" % d] = cov_by_donor[d]
    _write_matrix(pd.DataFrame(data))


def main():
    if not _exists(f_vcf):
        _empty("known_genotypes.vcf missing/empty; no donors screened")
        return
    samples, private = parse_private(f_vcf)
    if not samples:
        _empty("no donors in known_genotypes.vcf")
        return

    # ---- roster-awareness + pre-flight ungenotyped-donor gate ----------------
    # `samples` = genotyped donors (panel columns). Ungenotyped = declared roster
    # donors omitted from the panel (no genotype in any source). By-elimination can
    # recover AT MOST ONE; >1 is unrecoverable -> hard-fail the capture up front
    # (clearer than the deep annotate.py raise, which remains a backstop).
    genotyped_set = set(samples)
    ungenotyped = [d for d in roster if d not in genotyped_set]
    if len(ungenotyped) > 1:
        msg = (
            "PRE-FLIGHT FAIL [%s]: %d pool donors have NO known genotype (%s). "
            "Recovery-by-elimination resolves at most ONE ungenotyped donor per "
            "pool; provide genotypes for all but one (or drop the extras from the "
            "roster) before demuxing this capture."
            % (capture, len(ungenotyped), ", ".join(ungenotyped))
        )
        log.error(msg)
        raise RuntimeError(msg)
    if ungenotyped:
        log.info(
            "donor=%s is UNGENOTYPED (recover-by-elimination candidate; no private "
            "sites, so it cannot be positively screened -- validate post-hoc)",
            ungenotyped[0],
        )

    have_mtx = _exists(f_ad) and _exists(f_dp) and _exists(f_base)
    if have_mtx:
        idx = base_index(f_base)
        AD = mmread(f_ad).tocsr()
        DP = mmread(f_dp).tocsr()
        log.info("pileup: %d base sites, AD shape %s", len(idx), AD.shape)
    else:
        idx, AD, DP = {}, None, None
        log.warning("cellSNP AD/DP/base missing; every donor -> untestable")

    n_cells = AD.shape[1] if have_mtx else 0
    hit_by_donor = {}  # donor -> per-cell #private-alt sites hit (for the matrix)
    cov_by_donor = {}  # donor -> per-cell #private sites covered (any read)

    rows_out = []
    for donor in samples:
        psites = private.get(donor, set())
        n_private = len(psites)
        rows = [idx[k] for k in psites if k in idx] if idx else []
        n_covered = len(rows)

        ambient_p = 0.0
        cells_ge1 = cells_ge2 = cells_ge3 = cells_ge5 = cells_ge10 = 0
        max_hits = 0
        n_evidence = 0
        n_hit = np.zeros(n_cells, dtype=int)
        n_cov_cell = np.zeros(n_cells, dtype=int)
        if have_mtx and rows:
            ad = AD[rows, :]
            dp = DP[rows, :]
            n_hit = np.asarray((ad >= 1).sum(axis=0)).ravel()  # per cell: #private-alt sites
            n_cov_cell = np.asarray((dp >= 1).sum(axis=0)).ravel()  # per cell: #covered
            cells_ge1 = int((n_hit >= 1).sum())
            cells_ge2 = int((n_hit >= 2).sum())
            cells_ge3 = int((n_hit >= 3).sum())
            cells_ge5 = int((n_hit >= 5).sum())
            cells_ge10 = int((n_hit >= 10).sum())
            max_hits = int(n_hit.max()) if n_hit.size else 0
            n_evidence = int((n_hit >= K_SITES).sum())
            # Descriptive pooled ambient alt-fraction over d's private rows (report-only).
            tot_ad = float(ad.sum())
            tot_dp = float(dp.sum())
            ambient_p = (tot_ad / tot_dp) if tot_dp > 0 else 0.0
        hit_by_donor[donor] = n_hit
        cov_by_donor[donor] = n_cov_cell

        verdict = decide(n_covered, n_evidence)
        rows_out.append({
            "capture_id": capture,
            "donor": donor,
            "verdict": verdict,
            "n_evidence": n_evidence,
            "n_private": n_private,
            "n_covered": n_covered,
            "ambient_p": round(ambient_p, 5),
            "cells_ge1": cells_ge1,
            "cells_ge2": cells_ge2,
            "cells_ge3": cells_ge3,
            "cells_ge5": cells_ge5,
            "cells_ge10": cells_ge10,
            "max_hits": max_hits,
        })
        log.info(
            "donor=%s verdict=%s n_evidence=%d n_private=%d n_covered=%d ambient_p=%.5f "
            "cells_ge1=%d cells_ge2=%d cells_ge3=%d cells_ge5=%d cells_ge10=%d max_hits=%d",
            donor, verdict, n_evidence, n_private, n_covered, ambient_p,
            cells_ge1, cells_ge2, cells_ge3, cells_ge5, cells_ge10, max_hits,
        )

    # The single ungenotyped donor: a placeholder row so the report names it and
    # flags it explicitly (no private sites -> nothing to positively screen).
    if ungenotyped:
        rows_out.append({
            "capture_id": capture, "donor": ungenotyped[0], "verdict": "ungenotyped",
            "n_evidence": 0, "n_private": 0, "n_covered": 0, "ambient_p": 0.0,
            "cells_ge1": 0, "cells_ge2": 0, "cells_ge3": 0, "cells_ge5": 0,
            "cells_ge10": 0, "max_hits": 0,
        })

    _emit_matrix(samples, hit_by_donor, cov_by_donor, n_cells)

    out = pd.DataFrame(rows_out, columns=OUT_COLS)
    # Non-`present` verdicts first (problems surface), then by donor for stability.
    out["_rank"] = out["verdict"].map(_VERDICT_RANK).fillna(9).astype(int)
    out = out.sort_values(["_rank", "donor"], kind="mergesort").drop(columns="_rank")
    _write(out)


if __name__ == "__main__":
    main()
