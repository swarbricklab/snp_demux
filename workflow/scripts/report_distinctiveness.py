# ===================
# report_distinctiveness.py
# -------------------
# Per-capture donor DISTINCTIVENESS diagnostic (design doc §D). Two parts:
#
#  1. From the real per-donor genotypes (known_genotypes.vcf, bare donor ids,
#     the pool-subset informative panel) count each donor's UNIQUE sites: sites
#     where that donor is the SOLE alt-carrier (het 0/1 or hom-alt 1/1) among the
#     pool's donors. Missing (./.) is treated as a NON-carrier (and counted/noted
#     separately) so a donor is "sole carrier" iff exactly one donor carries alt.
#
#  2. From the cellSNP pileup (AD/DP matrices + base.vcf coords) compute the
#     pooled allele fraction per site (sum AD / sum DP across cells) and, for each
#     donor's unique sites, the fraction SURVIVING pooled-AF thresholds
#     [0, 0.05, 0.1, 0.2] -- the MAF-survival curve. This quantifies the design's
#     failure mode 1 (minMAF culling of minority-donor private sites).
#
# Outputs a TSV (donor, n_unique, keep@0.05, keep@0.1, keep@0.2) + a barplot PNG.
#
# NOTE: snakemake-driven script (reads the `snakemake` object), mirroring
# hierarchical_consensus.py / annotate.py. Env: envs/vireosnp.yaml (pandas +
# numpy + matplotlib, NO scipy -> MatrixMarket is parsed by hand). Deterministic.
#
# cellSNP matrix orientation: sites x cells (rows = SNPs from base.vcf, cols =
# barcodes from cellSNP.samples.tsv) -- verified against a real pileup.
# ===================

import os
import gzip
import logging as log

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

log.basicConfig(
    filename=snakemake.log[0],
    level=log.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s",
)

f_vcf = snakemake.input["known_genotypes_vcf"]
f_ad = snakemake.input["ad_mtx"]
f_dp = snakemake.input["dp_mtx"]
f_base = snakemake.input["base_vcf"]
f_tsv = snakemake.output["tsv"]
f_png = snakemake.output["png"]

THRESHOLDS = [0.0, 0.05, 0.1, 0.2]
OUT_COLS = ["donor", "n_unique", "keep@0.05", "keep@0.1", "keep@0.2"]


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def _exists(path):
    return path and os.path.exists(path) and os.path.getsize(path) > 0


def _norm(chrom, pos):
    return (str(chrom).lower().replace("chr", ""), str(pos))


def _empty(note):
    """Emit an empty-but-valid TSV and a placeholder PNG rather than crashing."""
    log.warning("distinctiveness: %s", note)
    pd.DataFrame(columns=OUT_COLS).to_csv(f_tsv, sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.text(0.5, 0.5, note, ha="center", va="center", wrap=True, fontsize=9)
    ax.axis("off")
    fig.savefig(f_png, dpi=120, bbox_inches="tight")
    plt.close(fig)


# ---- 1. per-donor unique sites from the real genotypes ---------------------
def parse_gt_vcf(path):
    samples = []
    sites = []  # (chrom, pos, [carriers], n_missing)
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                continue
            if not samples:
                continue
            f = line.rstrip("\n").split("\t")
            gts = f[9:]
            carriers = []
            n_missing = 0
            for s, cell in zip(samples, gts):
                gt = cell.split(":", 1)[0].replace("|", "/")
                al = gt.split("/")
                if len(al) != 2 or "." in al:
                    n_missing += 1
                    continue  # ./. -> non-carrier (noted)
                if al[0] != "0" or al[1] != "0":  # het or hom-alt
                    carriers.append(s)
            sites.append((f[0], f[1], carriers, n_missing))
    return samples, sites


if not _exists(f_vcf):
    _empty("known_genotypes.vcf missing/empty; no distinctiveness computed")
    raise SystemExit(0)

samples, sites = parse_gt_vcf(f_vcf)
if not samples or not sites:
    _empty("no donors/sites in known_genotypes.vcf")
    raise SystemExit(0)

n_unique = {s: 0 for s in samples}
unique_sites = {s: [] for s in samples}
n_missing_total = 0
for chrom, pos, carriers, nm in sites:
    n_missing_total += nm
    if len(carriers) == 1:
        d = carriers[0]
        n_unique[d] += 1
        unique_sites[d].append(_norm(chrom, pos))
log.info("donors=%d sites=%d missing-GT-calls=%d unique-per-donor=%s",
         len(samples), len(sites), n_missing_total, n_unique)


# ---- 2. pooled AF per cellSNP site -----------------------------------------
def read_mtx_rowsums(path):
    """Sum each row of a MatrixMarket coordinate file. Returns (nrows, rowsum)."""
    rowsum = None
    with _open(path) as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            parts = line.split()
            if rowsum is None:
                rowsum = np.zeros(int(parts[0]), dtype=np.float64)
                continue
            rowsum[int(parts[0]) - 1] += float(parts[2])
    return (0 if rowsum is None else len(rowsum)), (rowsum if rowsum is not None else np.zeros(0))


def read_base_coords(path):
    coords = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split("\t")
            coords.append(_norm(f[0], f[1]))
    return coords


site_af = {}
if _exists(f_ad) and _exists(f_dp) and _exists(f_base):
    _, ad_rs = read_mtx_rowsums(f_ad)
    _, dp_rs = read_mtx_rowsums(f_dp)
    coords = read_base_coords(f_base)
    for i, key in enumerate(coords):
        dp = dp_rs[i] if i < len(dp_rs) else 0.0
        ad = ad_rs[i] if i < len(ad_rs) else 0.0
        site_af[key] = (ad / dp) if dp > 0 else 0.0
    log.info("pileup sites with AF: %d", len(site_af))
else:
    log.warning("cellSNP AD/DP/base missing; AF treated as 0 (no site survives)")


# ---- 3. MAF-survival per donor + output ------------------------------------
rows = []
for d in samples:
    us = unique_sites[d]
    n = len(us)
    rec = {"donor": d, "n_unique": n}
    for t in THRESHOLDS:
        if t == 0.0:
            continue  # keep@0 is trivially 1.0; not reported
        if n == 0:
            rec["keep@%g" % t] = float("nan")
        else:
            surv = sum(1 for key in us if site_af.get(key, 0.0) >= t)
            rec["keep@%g" % t] = surv / n
    rows.append(rec)

out = pd.DataFrame(rows, columns=OUT_COLS)
out.to_csv(f_tsv, sep="\t", index=False)
log.info("wrote %s", f_tsv)

# barplot: grouped keep@thresholds per donor (annotated with n_unique)
fig, ax = plt.subplots(figsize=(max(6, 1.1 * len(samples)), 4.2))
donors = out["donor"].tolist()
x = np.arange(len(donors))
keep_cols = ["keep@0.05", "keep@0.1", "keep@0.2"]
w = 0.8 / len(keep_cols)
colors = ["#2c6fa8", "#5aa3d6", "#a9cbe8"]
for j, col in enumerate(keep_cols):
    vals = out[col].fillna(0.0).to_numpy()
    ax.bar(x + (j - 1) * w, vals, width=w, label=col, color=colors[j])
for i, nu in enumerate(out["n_unique"].tolist()):
    ax.text(x[i], 1.02, "n=%d" % nu, ha="center", va="bottom", fontsize=7, rotation=90)
ax.set_xticks(x)
ax.set_xticklabels([str(d) for d in donors], rotation=45, ha="right")
ax.set_ylim(0, 1.25)
ax.set_ylabel("fraction of unique sites surviving")
ax.set_title("Donor distinctiveness: MAF survival of unique sites")
ax.legend(fontsize=8, loc="upper right")
fig.tight_layout()
fig.savefig(f_png, dpi=120)
plt.close(fig)
log.info("wrote %s", f_png)
