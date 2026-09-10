# ===================
# hierarchical_consensus.py
# -------------------
# Sex-aware hierarchical consensus, run ALONGSIDE demuxafy (see
# docs/hierarchical_demux_design.md).
#
# HIERARCHY (INVERTED vs the pilot). The pilot's cluster-level GPb_diff margin
# gate was too conservative and dropped real donors that vireo's per-cell
# assignment clearly recovers (e.g. donor 4713 in Atlas_Pool_7a: vireo found
# ~106 cells the old cluster gate missed). With the WGS-harmonized panels vireo
# per-cell is strong, so:
#
#   1. BACKBONE  = vireo per-cell donor call (labelling).
#   2. RECOVERY  = souporcell DE-NOVO clusters, used ONLY to rescue cells vireo
#                  left "unassigned" (map the cell's de-novo cluster -> donor).
#   3. SEX       = a tie-break / QC ADJUSTMENT (down-weight + doublet signal),
#                  never a silent drop.
#   4. MARGIN    = DEMOTED. The GPb_diff margin only lowers confidence for
#                  already-weak (recovered / low_conf) labels; it can NEVER
#                  override a confident vireo singlet.
#
# recovery_confidence (souporcell P(singlet)) and label_confidence (how sure we
# are of the DONOR) stay SEPARATE columns, so a rescued cell with a weak label
# is kept, never silently dropped.
#
# NOTE: snakemake-driven script (reads the `snakemake` object), mirroring
# annotate.py / demuxafy_counts_AtLeastHalfSinglet.py. Deterministic (no RNG).
#
# Output columns:
#   barcode, cluster, donor, recovery_confidence, label_confidence, doublet_flag,
#   sex_call, cluster_sex, sex_concordant, denovo_known_agree, status
# status in {singlet, recovered, low_conf, doublet, unassigned, low_cell}
#   singlet    - confident vireo per-cell donor call
#   recovered  - vireo-unassigned (or weak) cell rescued via its de-novo cluster
#   low_conf   - vireo assigned a donor but weakly (prob_max < soft threshold)
#                and de-novo did not confirm it
#   doublet    - vireo doublet / high prob_doublet / souporcell doublet / confident
#                sex conflict
#   unassigned - no vireo donor and no de-novo recovery
#   low_cell   - recovered from a de-novo cluster below min_cells_floor
# ===================

import os
import math
import logging as log
from collections import Counter, defaultdict

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

# Module-internal thresholds. Not exposed as config because they are intrinsic
# to the method/biology, not tuning knobs.
CHRX_HET_THRESH = 0.08     # female if chrX het fraction exceeds this
MIN_SEX_SITES = 10         # min called chrX (or chrY) sites to make a call
SEX_CONF_MIN = 0.5         # numeric per-barcode sex confidence considered "high"
VIREO_SINGLET_MIN = 0.9    # soft prob_max threshold for a "confident" vireo call
VIREO_DOUBLET_MAX = 0.9    # prob_doublet above this => doublet
# de-novo tokens that are NOT a real donor label
VIREO_SPECIAL = {"doublet", "unassigned", "", "nan", "na"}

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


def _to_float(v):
    try:
        f = float(v)
        return f if f == f else np.nan  # reject NaN strings
    except (ValueError, TypeError):
        return np.nan


def _is_real_donor(vd):
    """True if a vireo donor_id is an actual donor (not doublet/unassigned/NA)."""
    if vd is None:
        return False
    if isinstance(vd, float) and math.isnan(vd):
        return False
    if pd.isna(vd):
        return False
    return str(vd).strip().lower() not in VIREO_SPECIAL


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

# per-barcode lookups off the de-novo table (used in the union assembly below)
sc_status_map = dict(zip(denovo["barcode"], sc_status))
cluster_map = dict(zip(denovo["barcode"], cluster))          # values may be NaN
recovery_conf_map = dict(zip(denovo["barcode"], recovery_conf))

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
agree_map = dict(zip(denovo["barcode"], denovo_known_agree))

# ---- 3. cluster -> donor label from vireosnp GPb_diff margin ----------------
# scores_full: rows = known genotype (donor), cols = "snp_cluster_<id>". LOWER
# GPb_diff = better match (see annotate.py). best = min; margin = 2nd_min - min.
# In the inverted hierarchy this is now a SECONDARY (fallback) labeller for the
# recovery layer, not the primary gate.
cluster_label = {}   # cluster id -> best donor (vireosnp)
cluster_margin = {}  # cluster id -> margin (best vs 2nd best); inf if single donor
cluster_ranked = {}  # cluster id -> donors sorted best-first (ascending GPb_diff)

scores = _read_tsv(f_scores_full)
scores_donors = set()  # genotyped donors present in the GPb matrix (NOT the by-elimination donor)
if scores is not None and "known genotype" in scores.columns:
    scores = scores.set_index("known genotype")
    scores_donors = set(str(d) for d in scores.index)
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
# conflicts with the cluster's sex, then re-pick the best. This sharpens the
# vireosnp FALLBACK label used by the recovery layer (it no longer gates the
# vireo backbone).
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

# ---- 6. vireo per-cell backbone --------------------------------------------
# donor_ids.tsv columns: cell, donor_id, prob_max, prob_doublet, n_vars,
# best_singlet, best_doublet, doublet_logLikRatio. donor_id is the REAL donor
# (relabelled), or "doublet"/"unassigned". prob_max = labelling confidence,
# prob_doublet = doublet posterior.
vireo = _read_tsv(f_vireo, dtype=str)
vdonor_map = {}   # barcode -> donor_id (real donor / "doublet" / "unassigned")
vpmax_map = {}    # barcode -> prob_max (float)
vpd_map = {}      # barcode -> prob_doublet (float)
vireo_barcodes = []
if vireo is not None:
    bcol = "cell" if "cell" in vireo.columns else vireo.columns[0]
    has_pmax = "prob_max" in vireo.columns
    has_pd = "prob_doublet" in vireo.columns
    has_don = "donor_id" in vireo.columns
    for _, r in vireo.iterrows():
        bc = r[bcol]
        vireo_barcodes.append(bc)
        if has_don:
            vdonor_map[bc] = r["donor_id"]
        if has_pmax:
            vpmax_map[bc] = _to_float(r["prob_max"])
        if has_pd:
            vpd_map[bc] = _to_float(r["prob_doublet"])
    log.info("vireo backbone cells: %d (real=%d)",
             len(vireo_barcodes),
             sum(1 for v in vdonor_map.values() if _is_real_donor(v)))
else:
    log.warning("vireo donor_ids missing/empty; backbone unavailable, relying on de-novo recovery only")

# ---- 6a. de-novo cluster -> donor recovery map -----------------------------
# For cells vireo left "unassigned" we recover the donor from the cell's de-novo
# cluster. Two evidence sources, in order:
#   (a) the MAJORITY real-donor vireo label among that cluster's CONFIDENT members
#       (prob_max >= VIREO_SINGLET_MIN). label_confidence = fraction of the
#       cluster's members that vireo confidently gave that donor (capped < 1).
#   (b) fallback: the vireosnp GPb_diff cluster label (sex-reranked), or the
#       clusters_annotated majority. Given a low, margin-scaled confidence.
cluster_members = defaultdict(list)
for bc, cl in cluster_map.items():
    if cl is not None and not (isinstance(cl, float) and pd.isna(cl)) and not pd.isna(cl):
        cluster_members[cl].append(bc)
cluster_n = {cl: len(m) for cl, m in cluster_members.items()}

recovery_donor = {}       # cluster -> recovered donor
recovery_label_conf = {}  # cluster -> label_confidence for a recovery
for cl, members in cluster_members.items():
    n_members = len(members)
    votes = Counter()
    for bc in members:
        vd = vdonor_map.get(bc)
        vp = vpmax_map.get(bc, np.nan)
        if _is_real_donor(vd) and (vp == vp) and vp >= VIREO_SINGLET_MIN:
            votes[str(vd)] += 1
    if votes:
        # deterministic argmax: highest count, tie-break on donor id string
        donor, cnt = sorted(votes.items(), key=lambda kv: (-kv[1], kv[0]))[0]
        recovery_donor[cl] = donor
        # fraction of the WHOLE cluster confidently given this donor by vireo,
        # capped below 1 so a recovered cell never outranks a native vireo singlet
        recovery_label_conf[cl] = min(cnt / n_members, 0.95) if n_members else 0.0
    elif cl in annot_label and str(annot_label[cl]) not in scores_donors:
        # vireoSNP annotated this cluster as a donor ABSENT from the GPb matrix --
        # the ungenotyped BY-ELIMINATION donor. vireoSNP is the only authoritative
        # source for it; cluster_label (GPb argmin) below can only pick a GENOTYPED
        # donor and would mislabel this leftover cluster as its nearest known donor.
        recovery_donor[cl] = annot_label[cl]
        recovery_label_conf[cl] = 0.4
    elif cl in cluster_label:
        recovery_donor[cl] = cluster_label[cl]
        m = cluster_margin.get(cl, np.nan)
        if m == float("inf"):
            base = 1.0
        elif margin_min > 0 and m == m:
            base = min(m / (2.0 * margin_min), 1.0)
        else:
            base = 0.0
        recovery_label_conf[cl] = min(0.5, 0.5 * base)  # secondary evidence: low
    elif cl in annot_label:
        recovery_donor[cl] = annot_label[cl]
        recovery_label_conf[cl] = 0.3
log.info("de-novo recovery map: %d clusters", len(recovery_donor))

# ---- 6b. assemble per-cell rows --------------------------------------------
# Cell universe = union of de-novo and vireo barcodes. A confident vireo singlet
# with NO souporcell cluster (cluster NA) is still labelled from the backbone.
seen = set()
master_barcodes = []
for bc in denovo["barcode"]:
    if bc not in seen:
        seen.add(bc)
        master_barcodes.append(bc)
for bc in vireo_barcodes:
    if bc not in seen:
        seen.add(bc)
        master_barcodes.append(bc)

rows = []
for bc in master_barcodes:
    st_sc = sc_status_map.get(bc)            # souporcell status (may be None for vireo-only)
    cl = cluster_map.get(bc)
    cl = None if (cl is None or (isinstance(cl, float) and pd.isna(cl)) or pd.isna(cl)) else cl

    vdonor = vdonor_map.get(bc)
    vpmax = vpmax_map.get(bc, np.nan)
    vpd = vpd_map.get(bc, np.nan)

    # ---- pick donor + provisional confidence from backbone, then recovery ----
    donor = None
    label_conf = 0.0
    backbone = "none"
    # 0. BY-ELIMINATION OVERRIDE (takes precedence over the vireo backbone). If
    # souporcell placed this cell in a cluster that vireoSNP annotated as the
    # UNGENOTYPED donor (a donor absent from the GPb matrix, which vireo cannot
    # emit), trust that de-novo/by-elimination label. vireo (genotyped donors
    # only) would otherwise mislabel these cells to their nearest KNOWN donor.
    # Only for souporcell singlets; doublets fall through to the doublet logic.
    byelim_donor = (str(annot_label[cl]) if (cl is not None and cl in annot_label
                    and str(annot_label[cl]) not in scores_donors) else None)
    if byelim_donor is not None and st_sc == "singlet":
        donor = byelim_donor
        label_conf = recovery_label_conf.get(cl, 0.4)
        backbone = "recovered"
    elif _is_real_donor(vdonor):
        # 1. BACKBONE: vireo per-cell donor call
        donor = str(vdonor)
        label_conf = vpmax if (vpmax == vpmax) else 0.0
        backbone = "vireo"
    elif vdonor is not None and str(vdonor).strip().lower() == "doublet":
        backbone = "vireo_doublet"   # resolved as doublet below
    else:
        # vireo "unassigned" OR barcode missing from vireo -> try de-novo recovery
        if cl is not None and cl in recovery_donor:
            donor = str(recovery_donor[cl])
            label_conf = recovery_label_conf.get(cl, 0.0)
            backbone = "recovered"
        else:
            backbone = "none"

    # ---- sex QC (adjustment, not a gate) -------------------------------------
    csex = cluster_sex.get(cl, "unknown") if cl is not None else "unknown"
    dsex = donor_sex.get(donor, "unknown") if donor else "unknown"
    bc_sex = _norm_sex(sex_map.get(bc))
    # per-cell concordance: cell's own confident sex vs assigned donor's panel sex
    if bc_sex and dsex in ("male", "female") and _high_conf(conf_map.get(bc)):
        sex_conc = bool(bc_sex == dsex)
    else:
        sex_conc = pd.NA
    # confident sex conflict -> down-weight label_confidence (never a silent drop)
    if sex_conc is False:
        label_conf *= max(0.0, 1.0 - sex_weight)

    # ---- doublet flag --------------------------------------------------------
    # For a cell resolved via the by-elimination / de-novo recovery override,
    # vireo's own doublet signal is UNRELIABLE: the cell belongs to an ungenotyped
    # donor absent from vireo's reference, so vireo (unable to place it) inflates
    # its doublet posterior and would spuriously drop a real singlet. Trust
    # souporcell's singlet call instead; keep the souporcell-doublet + sex signals.
    vireo_doublet_sig = (
        (vdonor is not None and str(vdonor).strip().lower() == "doublet")
        or ((vpd == vpd) and vpd > VIREO_DOUBLET_MAX)
    )
    doublet_flag = bool(
        (vireo_doublet_sig and backbone != "recovered")
        or (st_sc == "doublet")
        or (sex_conc is False)
    )

    # ---- status --------------------------------------------------------------
    is_low_cell = (cl is not None) and (cluster_n.get(cl, 0) < min_cells_floor)

    if doublet_flag:
        status = "doublet"
        donor = "doublet"
    elif backbone == "vireo":
        if (vpmax == vpmax) and vpmax >= VIREO_SINGLET_MIN:
            # confident vireo singlet: stays singlet regardless of cluster margin
            status = "singlet"
        elif cl is not None and recovery_donor.get(cl) == donor:
            # weak vireo call CONFIRMED by its de-novo cluster -> rescued to confidence
            label_conf = max(label_conf, recovery_label_conf.get(cl, 0.0))
            status = "recovered"
        elif low_conf_keep:
            status = "low_conf"   # keep the weak vireo label
        else:
            status = "unassigned"
            donor = None
    elif backbone == "recovered":
        # vireo-unassigned cell rescued from its de-novo cluster
        status = "low_cell" if is_low_cell else "recovered"
    else:
        status = "unassigned"
        donor = None

    rec = recovery_conf_map.get(bc, np.nan)

    rows.append({
        "barcode": bc,
        "cluster": cl,
        "donor": donor,
        "recovery_confidence": round(rec, 4) if not (rec is None or (isinstance(rec, float) and np.isnan(rec))) else np.nan,
        "label_confidence": round(float(label_conf), 4),
        "doublet_flag": doublet_flag,
        "sex_call": sex_map.get(bc) if not pd.isna(sex_map.get(bc, pd.NA)) else "",
        "cluster_sex": csex,
        "sex_concordant": sex_conc,
        "denovo_known_agree": agree_map.get(bc, pd.NA),
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
