# ===================
# leave_one_out.py
# -------------------
# Cheap, score-level LEAVE-ONE-OUT validator (design doc §E). Reuses the vireosnp
# GPb_diff matrix (full_GPb_diff_scores_cluster_genotype.tsv: rows = known genotype
# donor, cols = snp_cluster_<id>, LOWER = better match) instead of re-running any
# clustering.
#
# For each donor d:
#   * best_cluster   = the cluster d matches best (argmin GPb_diff over d's row);
#   * own_score      = d's GPb_diff at that cluster;
#   * drop d's row, then look at that cluster's NEXT-best donor:
#       next_best_donor / next_best_score = argmin / min over remaining donors;
#   * margin         = next_best_score - own_score;
#   * verdict        = "absorbed" if margin < margin_min (another donor is nearly
#                      as good -> the cluster is NOT uniquely d; removing d it gets
#                      absorbed / not separable) else "distinct" (clear gap).
#
# NOTE: snakemake-driven script (reads the `snakemake` object). Env:
# envs/vireosnp.yaml (pandas). Deterministic. Robust to a missing/empty matrix.
# ===================

import os
import logging as log

import pandas as pd

log.basicConfig(
    filename=snakemake.log[0],
    level=log.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s",
)

f_scores = snakemake.input["scores_full"]
f_out = snakemake.output["tsv"]
margin_min = float(snakemake.params["margin_min"])

OUT_COLS = ["donor", "best_cluster", "own_score", "next_best_donor",
            "next_best_score", "margin", "verdict"]


def _empty(note):
    log.warning("leave_one_out: %s", note)
    pd.DataFrame(columns=OUT_COLS).to_csv(f_out, sep="\t", index=False)


if not (f_scores and os.path.exists(f_scores) and os.path.getsize(f_scores) > 0):
    _empty("scores_full missing/empty")
    raise SystemExit(0)

scores = pd.read_csv(f_scores, sep="\t")
if "known genotype" not in scores.columns:
    _empty("scores_full lacks 'known genotype' column")
    raise SystemExit(0)

scores = scores.set_index("known genotype")
clust_cols = [c for c in scores.columns if c.startswith("snp_cluster_")]
scores = scores[clust_cols].apply(pd.to_numeric, errors="coerce")
if scores.empty or not clust_cols:
    _empty("no snp_cluster_* columns / no donors")
    raise SystemExit(0)

log.info("LOO over %d donors x %d clusters (margin_min=%s)",
         scores.shape[0], scores.shape[1], margin_min)

rows = []
for d in scores.index:
    row = scores.loc[d].dropna()
    if row.empty:
        rows.append({"donor": d, "best_cluster": "", "own_score": float("nan"),
                     "next_best_donor": "", "next_best_score": float("nan"),
                     "margin": float("nan"), "verdict": "distinct"})
        continue
    best_col = row.idxmin()
    own = float(row.min())
    best_cluster = best_col[len("snp_cluster_"):]

    others = scores[best_col].drop(index=d).dropna()
    if others.empty:
        rows.append({"donor": d, "best_cluster": best_cluster,
                     "own_score": round(own, 6), "next_best_donor": "",
                     "next_best_score": float("nan"), "margin": float("inf"),
                     "verdict": "distinct"})
        continue
    nb_donor = others.idxmin()
    nb_score = float(others.min())
    margin = nb_score - own
    verdict = "absorbed" if margin < margin_min else "distinct"
    rows.append({"donor": d, "best_cluster": best_cluster,
                 "own_score": round(own, 6), "next_best_donor": nb_donor,
                 "next_best_score": round(nb_score, 6),
                 "margin": round(margin, 6), "verdict": verdict})

out = pd.DataFrame(rows, columns=OUT_COLS)
out.to_csv(f_out, sep="\t", index=False)
log.info("wrote %s; verdicts=%s", f_out, out["verdict"].value_counts().to_dict())
