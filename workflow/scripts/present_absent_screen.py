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
_VERDICT_RANK = {"absent": 0, "untestable": 1, "present": 2}


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


# ---- 3. verdict --------------------------------------------------------------
def decide(n_covered, n_evidence):
    if n_covered < MIN_COVERED_SITES:
        return "untestable"
    if n_evidence >= MIN_CELLS:
        return "present"
    return "absent"


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
        log.info("pileup: %d base sites, AD shape %s", len(idx), AD.shape)
    else:
        idx, AD, DP = {}, None, None
        log.warning("cellSNP AD/DP/base missing; every donor -> untestable")

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
        if have_mtx and rows:
            ad = AD[rows, :]
            dp = DP[rows, :]
            n_hit = np.asarray((ad >= 1).sum(axis=0)).ravel()  # per cell: #private-alt sites
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

    out = pd.DataFrame(rows_out, columns=OUT_COLS)
    # Non-`present` verdicts first (problems surface), then by donor for stability.
    out["_rank"] = out["verdict"].map(_VERDICT_RANK).fillna(9).astype(int)
    out = out.sort_values(["_rank", "donor"], kind="mergesort").drop(columns="_rank")
    _write(out)


if __name__ == "__main__":
    main()
