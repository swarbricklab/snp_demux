# ===================
# report_variant_embedding.py
# -------------------
# Per-capture variant-space EMBEDDING diagnostic (design doc §D). Mirrors the
# investigation's UMAP analysis: build a cells x sites alt-fraction matrix from
# the cellSNP pileup, reduce with TruncatedSVD then UMAP, colour points by the
# DE-NOVO souporcell cluster (singlet assignment), and score cluster separation
# with the silhouette on the SVD space. A negative / near-zero silhouette is the
# design's failure-mode 2 signature (genuinely too few cells / no signal).
#
# NOTE: snakemake-driven script (reads the `snakemake` object). Env:
# envs/embedding.yaml (scipy + scikit-learn + umap-learn + matplotlib + numpy).
# Deterministic: seed 0 for SVD and UMAP. Falls back to the first two SVD
# components if umap-learn cannot be imported.
#
# cellSNP matrix orientation: sites x cells (rows = SNPs, cols = barcodes in
# cellSNP.samples.tsv order) -- verified against a real pileup; we transpose to
# cells x sites for the embedding.
# ===================

import os
import logging as log

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RANDOM_STATE = 0
np.random.seed(RANDOM_STATE)

log.basicConfig(
    filename=snakemake.log[0],
    level=log.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s",
)

f_ad = snakemake.input["ad_mtx"]
f_dp = snakemake.input["dp_mtx"]
f_samples = snakemake.input["samples"]
f_clusters = snakemake.input["clusters"]
f_png = snakemake.output["png"]
f_tsv = snakemake.output["tsv"]

METRIC_COLS = ["cluster", "n_cells", "silhouette"]


def _exists(path):
    return path and os.path.exists(path) and os.path.getsize(path) > 0


def _empty(note):
    """Emit an empty-but-valid metrics TSV and a placeholder PNG."""
    log.warning("embedding: %s", note)
    pd.DataFrame(columns=METRIC_COLS).to_csv(f_tsv, sep="\t", index=False)
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.text(0.5, 0.5, note, ha="center", va="center", wrap=True, fontsize=9)
    ax.axis("off")
    fig.savefig(f_png, dpi=120, bbox_inches="tight")
    plt.close(fig)


if not (_exists(f_ad) and _exists(f_dp)):
    _empty("cellSNP AD/DP matrices missing/empty; no embedding computed")
    raise SystemExit(0)

from scipy.io import mmread

AD = mmread(f_ad).tocsr().astype(np.float64)
DP = mmread(f_dp).tocsr().astype(np.float64)
n_sites, n_cells = AD.shape
log.info("pileup: %d sites x %d cells", n_sites, n_cells)

if n_cells < 3 or n_sites < 2:
    _empty("too few cells/sites for an embedding (%d x %d)" % (n_sites, n_cells))
    raise SystemExit(0)

# sites covered in >= max(20, 3% of cells) cells (DP > 0)
cover = np.asarray((DP > 0).sum(axis=1)).ravel()
min_cover = max(20, int(np.ceil(0.03 * n_cells)))
keep = cover >= min_cover
n_keep = int(keep.sum())
log.info("sites covered in >=%d cells: %d/%d", min_cover, n_keep, n_sites)
if n_keep < 2:
    _empty("only %d sites pass the coverage floor (>=%d cells)" % (n_keep, min_cover))
    raise SystemExit(0)

# cells x sites alt-fraction (0 where DP == 0)
sub_ad = AD[keep, :].toarray()
sub_dp = DP[keep, :].toarray()
with np.errstate(divide="ignore", invalid="ignore"):
    af = np.where(sub_dp > 0, sub_ad / sub_dp, 0.0)
X = af.T  # cells x sites
log.info("alt-fraction matrix: %s", X.shape)

# ---- TruncatedSVD (seed 0) --------------------------------------------------
from sklearn.decomposition import TruncatedSVD

n_comp = int(min(30, X.shape[1] - 1, X.shape[0] - 1))
n_comp = max(2, n_comp)
svd = TruncatedSVD(n_components=n_comp, random_state=RANDOM_STATE)
svd_emb = svd.fit_transform(X)
log.info("SVD embedding: %s (%d comps)", svd_emb.shape, n_comp)

# ---- UMAP (seed 0), fallback to first 2 SVD comps ---------------------------
try:
    import umap
    reducer = umap.UMAP(random_state=RANDOM_STATE)
    emb2 = reducer.fit_transform(svd_emb)
    method = "UMAP"
except Exception as exc:  # noqa: BLE001 - any import/runtime failure -> fallback
    log.warning("umap unavailable (%s); falling back to first 2 SVD components", exc)
    emb2 = svd_emb[:, :2]
    method = "SVD"

# ---- de-novo souporcell cluster labels, aligned to sample order -------------
samples = []
if _exists(f_samples):
    with open(f_samples) as fh:
        samples = [ln.strip() for ln in fh if ln.strip()]
if len(samples) != n_cells:
    log.warning("samples.tsv has %d barcodes but matrix has %d cells; "
                "labels may be unaligned", len(samples), n_cells)

labels = np.array(["NA"] * n_cells, dtype=object)
if _exists(f_clusters) and samples and len(samples) == n_cells:
    cl = pd.read_csv(f_clusters, sep="\t", dtype=str)
    cl = cl.rename(columns={cl.columns[0]: "barcode"})
    status = cl.set_index("barcode")["status"] if "status" in cl.columns else None
    assign = cl.set_index("barcode")["assignment"] if "assignment" in cl.columns else None
    if assign is not None:
        for i, bc in enumerate(samples):
            if bc in assign.index:
                st = status.get(bc) if status is not None else "singlet"
                if st == "singlet":
                    v = assign.get(bc)
                    if v is not None and not pd.isna(v):
                        labels[i] = str(v)
else:
    log.warning("souporcell clusters unavailable/unaligned; points uncoloured")

# ---- silhouette on the SVD space -------------------------------------------
from sklearn.metrics import silhouette_score, silhouette_samples

metric_rows = []
mask = labels != "NA"
Xs = svd_emb[mask]
ls = labels[mask]
uniq = sorted(set(ls.tolist()))
if len(uniq) >= 2 and Xs.shape[0] > len(uniq):
    overall = float(silhouette_score(Xs, ls))
    samp = silhouette_samples(Xs, ls)
    metric_rows.append({"cluster": "overall", "n_cells": int(Xs.shape[0]),
                        "silhouette": round(overall, 4)})
    for c in uniq:
        cm = ls == c
        metric_rows.append({"cluster": str(c), "n_cells": int(cm.sum()),
                            "silhouette": round(float(samp[cm].mean()), 4)})
    log.info("overall silhouette=%.4f over %d clusters", overall, len(uniq))
else:
    log.warning("not enough labelled clusters for silhouette (%d clusters)", len(uniq))

pd.DataFrame(metric_rows, columns=METRIC_COLS).to_csv(f_tsv, sep="\t", index=False)
log.info("wrote %s", f_tsv)

# ---- plot -------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 6))
palette = plt.get_cmap("tab20")
labelled = [c for c in uniq]
color_of = {c: palette(i % 20) for i, c in enumerate(labelled)}
# unassigned/doublet cells first (grey), then coloured clusters on top
na_mask = labels == "NA"
if na_mask.any():
    ax.scatter(emb2[na_mask, 0], emb2[na_mask, 1], s=4, c="#cccccc",
               label="NA", alpha=0.5, linewidths=0)
for c in labelled:
    cm = labels == c
    ax.scatter(emb2[cm, 0], emb2[cm, 1], s=5, color=color_of[c],
               label="cluster %s" % c, alpha=0.8, linewidths=0)
overall_txt = ""
if metric_rows and metric_rows[0]["cluster"] == "overall":
    overall_txt = " (silhouette %.3f)" % metric_rows[0]["silhouette"]
ax.set_title("%s of variant alt-fraction, coloured by de-novo cluster%s"
             % (method, overall_txt))
ax.set_xlabel("%s-1" % method)
ax.set_ylabel("%s-2" % method)
if len(labelled) <= 20:
    ax.legend(fontsize=7, markerscale=2, loc="best", ncol=2)
fig.tight_layout()
fig.savefig(f_png, dpi=120)
plt.close(fig)
log.info("wrote %s", f_png)
