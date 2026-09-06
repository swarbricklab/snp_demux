# ===================
# hierarchical_consensus.py
# -------------------
# Sex-aware hierarchical consensus, run ALONGSIDE demuxafy (see
# docs/hierarchical_demux_design.md). The core idea: DE-NOVO souporcell decides
# *which cells group together* (recovery, high trust); GENOTYPES decide *which
# donor a group is* (labelling); SEX breaks low-margin ties. Recovery-confidence
# and label-confidence are reported as SEPARATE columns so a de-novo-recovered
# cell with a weak genotype label is kept, never silently dropped.
#
# NOTE: snakemake-driven script (reads the `snakemake` object), mirroring
# annotate.py / demuxafy_counts_AtLeastHalfSinglet.py. Deterministic (no RNG).
#
# Output columns:
#   barcode, cluster, donor, recovery_confidence, label_confidence, doublet_flag,
#   sex_call, cluster_sex, sex_concordant, denovo_known_agree, status
# status in {singlet, doublet, unassigned, low_cell, low_conf}
# ===================

import os
import math
import logging as log

import numpy as np
import pandas as pd

log.basicConfig(
    filename=snakemake.log[0],
    level=log.INFO,
    format="%(asctime)s:%(levelname)s:%(message)s",
)

# ---- inputs / params -------------------------------------------------------
f_denovo = snakemake.input["denovo_clusters"]
f_known = snakemake.input["known_clusters"]
f_vsnp_annot = snakemake.input["vireosnp_annotated"]
f_scores_full = snakemake.input["scores_full"]
f_vireo = snakemake.input["vireo_donor_ids"]
f_known_gt_vcf = snakemake.input["known_genotypes_vcf"]
f_sex = snakemake.input["sex_joined"]
f_out = snakemake.output["hier"]

margin_min = float(snakemake.params["margin_min"])
sex_weight = float(snakemake.params["sex_weight"])
min_cells_floor = int(snakemake.params["min_cells_floor"])
low_conf_keep = bool(snakemake.params["low_conf_keep"])

# Module-internal sex-calling thresholds (chrX heterozygosity). Not exposed as
# config because they are intrinsic to the biology, not tuning knobs.
CHRX_HET_THRESH = 0.08   # female if chrX het fraction exceeds this
MIN_SEX_SITES = 10       # min called chrX (or chrY) sites to make a call
SEX_CONF_MIN = 0.5       # numeric per-barcode sex confidence considered "high"

log.info("margin_min=%s sex_weight=%s min_cells_floor=%s low_conf_keep=%s",
         margin_min, sex_weight, min_cells_floor, low_conf_keep)


def _read_tsv(path, **kw):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return None
    try:
        return pd.read_csv(path, sep="\t", **kw)
    except Exception as exc:
        log.warning("Could not read %s: %s", path, exc)
        return None


# ---- 1. de-novo souporcell clusters (recovery backbone) --------------------
# souporcell clusters.tsv: barcode, status, assignment, log_prob_singleton,
# log_prob_doublet, cluster0..clusterN. status in {singlet,doublet,unassigned};
# assignment is the cluster id for singlets ("0"), or "0/1" for doublets.
denovo = _read_tsv(f_denovo, dtype=str)
if denovo is None:
    raise RuntimeError("de-novo souporcell clusters.tsv missing/empty: %s" % f_denovo)
denovo = denovo.rename(columns={denovo.columns[0]: "barcode"})
log.info("de-novo barcodes: %d", len(denovo))

sc_status = denovo["status"].str.lower() if "status" in denovo else pd.Series("singlet", index=denovo.index)
# cluster membership: only meaningful for singlets; doublets/unassigned -> NA
cluster = denovo["assignment"].where(sc_status.eq("singlet"))

# recovery_confidence = P(singlet) from souporcell's log-probs (softmax of the
# singleton vs doublet log-likelihoods). High => confidently a clean cluster member.
def _p_singlet(row):
    try:
        lps = float(row["log_prob_singleton"])
        lpd = float(row["log_prob_doublet"])
    except (KeyError, ValueError, TypeError):
        return np.nan
    m = max(lps, lpd)
    return math.exp(lps - m) / (math.exp(lps - m) + math.exp(lpd - m))

if {"log_prob_singleton", "log_prob_doublet"}.issubset(denovo.columns):
    recovery_conf = denovo.apply(_p_singlet, axis=1)
else:
    log.warning("souporcell log-prob columns absent; recovery_confidence set NaN")
    recovery_conf = pd.Series(np.nan, index=denovo.index)

# ---- 2. known-GT souporcell (de-novo vs known agreement QC) -----------------
# Independent cluster ids between the two runs, so we compare PARTITIONS, not ids:
# per de-novo cluster, the majority known-run cluster id; a barcode agrees if its
# known-run cluster matches that majority.
known = _read_tsv(f_known, dtype=str)
denovo_known_agree = pd.Series(pd.NA, index=denovo.index, dtype="object")
if known is not None:
    known = known.rename(columns={known.columns[0]: "barcode"})
    kmap = dict(zip(known["barcode"], known["assignment"])) if "assignment" in known else {}
    kass = denovo["barcode"].map(kmap)
    tmp = pd.DataFrame({"cluster": cluster, "kass": kass, "st": sc_status})
    singlet_mask = tmp["st"].eq("singlet") & tmp["cluster"].notna() & tmp["kass"].notna()
    majority = {}
    for cl, grp in tmp[singlet_mask].groupby("cluster"):
        majority[cl] = grp["kass"].value_counts().idxmax()

    def _agree(r):
        if r["st"] != "singlet" or pd.isna(r["cluster"]) or pd.isna(r["kass"]):
            return pd.NA
        return bool(r["kass"] == majority.get(r["cluster"]))

    denovo_known_agree = tmp.apply(_agree, axis=1)
else:
    log.info("known-GT souporcell clusters absent; denovo_known_agree=NA")

# ---- 3. cluster -> donor label from vireosnp GPb_diff margin ----------------
# scores_full: rows = known genotype (donor), cols = "snp_cluster_<id>". LOWER
# GPb_diff = better match (see annotate.py). best = min; margin = 2nd_min - min.
cluster_label = {}   # cluster id -> best donor
cluster_margin = {}  # cluster id -> margin (best vs 2nd best); inf if single donor
cluster_ranked = {}  # cluster id -> donors sorted best-first (ascending GPb_diff)

scores = _read_tsv(f_scores_full)
if scores is not None and "known genotype" in scores.columns:
    scores = scores.set_index("known genotype")
    for col in scores.columns:
        if not col.startswith("snp_cluster_"):
            continue
        cl = col[len("snp_cluster_"):]
        s = scores[col].dropna().astype(float).sort_values()  # ascending: best first
        if s.empty:
            continue
        cluster_ranked[cl] = list(s.index)
        cluster_label[cl] = s.index[0]
        cluster_margin[cl] = float(s.iloc[1] - s.iloc[0]) if len(s) > 1 else float("inf")
    log.info("vireosnp margins for %d clusters", len(cluster_label))
else:
    log.warning("scores_full missing/empty; falling back to vireosnp clusters_annotated for labels")

# Fallback labels from clusters_annotated.tsv (barcode,status,assignment=donor).
vsnp_annot = _read_tsv(f_vsnp_annot, dtype=str)
annot_label = {}
if vsnp_annot is not None and "assignment" in vsnp_annot.columns:
    bcol = vsnp_annot.columns[0]
    a = vsnp_annot.rename(columns={bcol: "barcode"})
    # map de-novo cluster id -> donor via the annotated assignment (majority)
    j = pd.DataFrame({"barcode": denovo["barcode"], "cluster": cluster}).merge(
        a[["barcode", "assignment"]], on="barcode", how="left")
    for cl, grp in j.dropna(subset=["cluster"]).groupby("cluster"):
        vc = grp["assignment"].dropna().value_counts()
        if not vc.empty:
            annot_label[cl] = vc.idxmax()

# ---- 4. donor -> sex from the panel (chrX heterozygosity / chrY presence) ---
# Principled + self-contained: males are hemizygous on chrX so their chrX GTs are
# ~all homozygous (low het fraction); females carry chrX hets. chrY presence is a
# fallback male signal. No external metadata used.
def donor_sex_from_vcf(path):
    if not path or not os.path.exists(path) or os.path.getsize(path) == 0:
        return {}
    samples = []
    het = {}
    called = {}
    y_called = {}
    with open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.rstrip("\n").split("\t")[9:]
                het = {s: 0 for s in samples}
                called = {s: 0 for s in samples}
                y_called = {s: 0 for s in samples}
                continue
            if not samples:
                continue
            f = line.rstrip("\n").split("\t")
            chrom = f[0].lower().replace("chr", "")
            if chrom not in ("x", "y"):
                continue
            gts = f[9:]
            for s, cell in zip(samples, gts):
                gt = cell.split(":", 1)[0].replace("|", "/")
                alleles = gt.split("/")
                if len(alleles) != 2 or "." in alleles:
                    continue  # missing
                if chrom == "y":
                    y_called[s] += 1
                    continue
                called[s] += 1
                if alleles[0] != alleles[1]:
                    het[s] += 1
    out = {}
    for s in samples:
        if called.get(s, 0) >= MIN_SEX_SITES:
            frac = het[s] / called[s]
            out[s] = "female" if frac > CHRX_HET_THRESH else "male"
        elif y_called.get(s, 0) >= MIN_SEX_SITES:
            out[s] = "male"
        else:
            out[s] = "unknown"
    return out

donor_sex = donor_sex_from_vcf(f_known_gt_vcf)
log.info("donor sexes from panel: %s", donor_sex)

# ---- 5. per-barcode sex + cluster sex --------------------------------------
sex = _read_tsv(f_sex, dtype=str)
if sex is None:
    sex = pd.DataFrame({"barcode": denovo["barcode"], "sex_call": pd.NA, "sex_confidence": pd.NA})
if "sex_call" not in sex.columns:
    sex["sex_call"] = pd.NA
if "sex_confidence" not in sex.columns:
    sex["sex_confidence"] = pd.NA
sex_map = dict(zip(sex["barcode"], sex["sex_call"]))
conf_map = dict(zip(sex["barcode"], sex["sex_confidence"]))

def _high_conf(v):
    if v is None or (isinstance(v, float) and math.isnan(v)) or pd.isna(v):
        return True  # no confidence column -> trust the call
    try:
        return float(v) >= SEX_CONF_MIN
    except (ValueError, TypeError):
        return str(v).lower() in ("high", "confident", "pass", "true")

def _norm_sex(v):
    if v is None or pd.isna(v):
        return None
    v = str(v).lower()
    return v if v in ("male", "female") else None

# cluster_sex: majority of HIGH-confidence, non-no_signal per-barcode calls.
cluster_sex = {}
csx = pd.DataFrame({"barcode": denovo["barcode"], "cluster": cluster})
for cl, grp in csx.dropna(subset=["cluster"]).groupby("cluster"):
    votes = {"male": 0, "female": 0}
    for bc in grp["barcode"]:
        sc = _norm_sex(sex_map.get(bc))
        if sc and _high_conf(conf_map.get(bc)):
            votes[sc] += 1
    if votes["male"] == 0 and votes["female"] == 0:
        cluster_sex[cl] = "unknown"
    elif votes["male"] == votes["female"]:
        cluster_sex[cl] = "unknown"  # tie -> do not constrain
    else:
        cluster_sex[cl] = "male" if votes["male"] > votes["female"] else "female"
log.info("cluster sexes: %s", cluster_sex)

# ---- 3b. SEX-CONSTRAINED re-ranking for low-margin clusters -----------------
# For a low-margin (ambiguous) cluster, exclude candidate donors whose panel sex
# conflicts with the cluster's sex, then re-pick the best. This directly breaks
# the "flat cluster" ties described in the design doc.
cluster_sex_resolved = {}  # cluster -> bool (sex changed/confirmed the label)
for cl, ranked in cluster_ranked.items():
    csex = cluster_sex.get(cl, "unknown")
    cluster_sex_resolved[cl] = False
    if csex == "unknown" or sex_weight <= 0:
        continue
    concordant = [d for d in ranked if donor_sex.get(d, "unknown") in (csex, "unknown")]
    if not concordant:
        continue
    if cluster_margin.get(cl, 0.0) < margin_min and concordant[0] != cluster_label.get(cl):
        # sex re-ranking changed the winner among ambiguous candidates
        cluster_label[cl] = concordant[0]
        cluster_sex_resolved[cl] = True
    elif concordant[0] == cluster_label.get(cl) and donor_sex.get(concordant[0], "unknown") == csex:
        cluster_sex_resolved[cl] = True  # sex confirms the margin winner

# ---- vireo doublet posterior ----------------------------------------------
vireo = _read_tsv(f_vireo, dtype=str)
vireo_pd = {}
if vireo is not None:
    bcol = "cell" if "cell" in vireo.columns else vireo.columns[0]
    if "prob_doublet" in vireo.columns:
        for _, r in vireo.iterrows():
            try:
                vireo_pd[r[bcol]] = float(r["prob_doublet"])
            except (ValueError, TypeError):
                pass
else:
    log.info("vireo donor_ids absent; prob_doublet unavailable")

# ---- cluster cell counts (low_cell flag) -----------------------------------
cluster_n = csx.dropna(subset=["cluster"]).groupby("cluster")["barcode"].count().to_dict()

# ---- 6. assemble per-cell rows ---------------------------------------------
rows = []
for i, bc in enumerate(denovo["barcode"]):
    st_sc = sc_status.iloc[i]
    cl = cluster.iloc[i]
    cl = None if pd.isna(cl) else cl

    # label
    donor = None
    margin = np.nan
    if cl is not None:
        donor = cluster_label.get(cl, annot_label.get(cl))
        margin = cluster_margin.get(cl, np.nan)

    csex = cluster_sex.get(cl, "unknown") if cl is not None else "unknown"
    dsex = donor_sex.get(donor, "unknown") if donor else "unknown"
    bc_sex = _norm_sex(sex_map.get(bc))

    # per-cell sex concordance: cell's own confident sex vs assigned donor sex
    if bc_sex and dsex in ("male", "female") and _high_conf(conf_map.get(bc)):
        sex_conc = bool(bc_sex == dsex)
    else:
        sex_conc = pd.NA

    # label_confidence: margin scaled to [0,1], penalised on donor/cluster sex conflict
    if np.isfinite(margin):
        lc = min(margin / (2.0 * margin_min), 1.0) if margin_min > 0 else (1.0 if margin > 0 else 0.0)
    elif margin == float("inf"):
        lc = 1.0  # single-donor pool -> label is unambiguous
    else:
        lc = 0.0
    if csex in ("male", "female") and dsex in ("male", "female") and csex != dsex:
        lc *= max(0.0, 1.0 - sex_weight)  # donor sex conflicts with cluster sex

    rec = recovery_conf.iloc[i]

    # doublet flag: souporcell doublet OR vireo prob_doublet>0.9 OR sex conflict
    vpd = vireo_pd.get(bc, np.nan)
    doublet_flag = bool(
        st_sc == "doublet"
        or (not np.isnan(vpd) and vpd > 0.9)
        or (sex_conc is False)
    )

    # status
    low_conf_label = cl is not None and np.isfinite(margin) and margin < margin_min \
        and not cluster_sex_resolved.get(cl, False)
    is_low_cell = cl is not None and cluster_n.get(cl, 0) < min_cells_floor

    if st_sc == "unassigned" or cl is None:
        status = "unassigned"
    elif doublet_flag:
        status = "doublet"
        donor = "doublet"
    elif is_low_cell:
        status = "low_cell"
    elif low_conf_label:
        if low_conf_keep:
            status = "low_conf"   # keep recovered cell with its weak label
        else:
            status = "unassigned"
            donor = None
    else:
        status = "singlet"

    rows.append({
        "barcode": bc,
        "cluster": cl,
        "donor": donor,
        "recovery_confidence": round(rec, 4) if not (rec is None or (isinstance(rec, float) and np.isnan(rec))) else np.nan,
        "label_confidence": round(lc, 4),
        "doublet_flag": doublet_flag,
        "sex_call": sex_map.get(bc) if not pd.isna(sex_map.get(bc, pd.NA)) else "",
        "cluster_sex": csex,
        "sex_concordant": sex_conc,
        "denovo_known_agree": denovo_known_agree.iloc[i],
        "status": status,
    })

out = pd.DataFrame(rows, columns=[
    "barcode", "cluster", "donor", "recovery_confidence", "label_confidence",
    "doublet_flag", "sex_call", "cluster_sex", "sex_concordant",
    "denovo_known_agree", "status",
])
out.to_csv(f_out, sep="\t", index=False)
log.info("wrote %s (%d cells)", f_out, len(out))
log.info("status counts: %s", out["status"].value_counts().to_dict())
