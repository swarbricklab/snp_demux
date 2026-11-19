"""Concordance gate for the by-elimination missing-donor assignment.

When one pool donor has no reference genotype, vireo learns its genotype
(GT_donors.vireo.vcf.gz) and the souporcell annotation hands it the one
leftover cluster, both purely by elimination. Nothing checks that the leftover
cluster really is that donor. This script does.

Per capture:
    concordance          genotype agreement (fraction of shared, called sites
                         with equal alt-allele dosage) between the leftover
                         donor's vireo-learned genotype and the de-novo
                         genotype of the souporcell cluster it was given
    ceiling              median of the same agreement over the pool's KNOWN
                         donors and their matched clusters: what this capture
                         can achieve when the assignment is right
    status               PASS | FAIL | NOT_APPLICABLE (no missing donor)

PASS requires n_sites >= min_sites, concordance >= min_concordance, and
concordance >= ceiling - max_drop (the last only when a ceiling exists).
Measured on real pools: easy 6-plex pools give ~0.93-0.97 (reliable), hard
8-plex pools ~0.45 (chance), with a known-donor ceiling around 0.87.

Rule combine_elimination_gate concatenates the per-capture tables and, for a
missing donor gated in several captures, adds the lowest pairwise agreement
between its vireo-learned genotypes across those captures (a donor whose
learned genotype is not stable across captures is not being learned).
"""
import gzip
import statistics
import sys
from itertools import combinations

COLUMNS = [
    "capture_id",
    "status",
    "missing_donor",
    "snp_cluster",
    "concordance",
    "n_sites",
    "ceiling",
    "n_known",
    "action",
    "reason",
]
COHORT_COLUMNS = COLUMNS + ["cross_capture_min_concordance", "cross_capture_n"]


# ---------------------------------------------------------------- VCF parsing
def _open(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


def _dosage(gt):
    alleles = gt.replace("|", "/").split("/")
    if not alleles or any(a in (".", "") for a in alleles):
        return None
    return sum(a != "0" for a in alleles)


def read_genotypes(path):
    """{sample: {(chrom, pos, ref, alt): dosage}} from a VCF's GT field.

    Returns {} for a file with no #CHROM line (e.g. the plain-text placeholder
    reheader_vireo_vcf writes when vireo learned no genotypes).
    """
    out = None
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            f = line.rstrip("\n").split("\t")
            if line.startswith("#CHROM"):
                names = f[9:]
                out = {n: {} for n in names}
                continue
            if out is None:
                return {}
            key = (f[0], f[1], f[3], f[4])
            fmt = f[8].split(":")
            if "GT" not in fmt:
                continue
            gi = fmt.index("GT")
            for n, cell in zip(names, f[9:]):
                parts = cell.split(":")
                d = _dosage(parts[gi]) if gi < len(parts) else None
                if d is not None:
                    out[n][key] = d
    return out or {}


def concordance(a, b):
    """(fraction of shared called sites with equal dosage, n shared sites)."""
    shared = a.keys() & b.keys()
    if not shared:
        return float("nan"), 0
    agree = sum(a[k] == b[k] for k in shared)
    return agree / len(shared), len(shared)


# --------------------------------------------------------------------- gate
def read_mapping(path):
    """{donor_id: snp_cluster} from final_sample_to_snpcluster_mappings.tsv."""
    out = {}
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        i_s, i_c = header.index("SAMPLE_ID"), header.index("SNP_CLUSTER")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            out[f[i_s]] = f[i_c]
    return out


def evaluate(
    capture_id,
    pool_donors,
    genotyped,
    mapping,
    vireo_gt,
    soup_gt,
    min_concordance,
    max_drop,
    min_sites,
    action,
):
    nan = float("nan")
    row = dict.fromkeys(COLUMNS, "")
    row.update(capture_id=capture_id, action=action, concordance=nan, ceiling=nan)
    row.update(n_sites=0, n_known=0)

    missing = [d for d in pool_donors if d not in genotyped]
    if not missing:
        row.update(status="NOT_APPLICABLE", reason="all donors genotyped")
        return row
    if len(missing) > 1:
        row.update(
            status="FAIL",
            missing_donor=",".join(missing),
            reason="more than one donor lacks a genotype",
        )
        return row
    donor = missing[0]
    row["missing_donor"] = donor

    # Known-donor ceiling: vireo-learned vs souporcell-de-novo for donors whose
    # cluster match was made on genotype, not by elimination.
    known = []
    for d in pool_donors:
        if d == donor or d not in mapping or d not in vireo_gt:
            continue
        c, n = concordance(vireo_gt[d], soup_gt.get(mapping[d], {}))
        if n >= min_sites:
            known.append(c)
    row["n_known"] = len(known)
    if known:
        row["ceiling"] = statistics.median(known)

    cluster = mapping.get(donor)
    row["snp_cluster"] = cluster or ""
    if cluster is None or cluster not in soup_gt:
        row.update(status="FAIL", reason="no souporcell cluster for missing donor")
        return row
    if donor not in vireo_gt:
        row.update(status="FAIL", reason="vireo learned no genotype for missing donor")
        return row

    c, n = concordance(vireo_gt[donor], soup_gt[cluster])
    row.update(concordance=c, n_sites=n)
    reasons = []
    if n < min_sites:
        reasons.append(f"only {n} shared sites (< {min_sites})")
    elif c < min_concordance:
        reasons.append(f"concordance {c:.3f} < {min_concordance}")
    if n >= min_sites and known and c < row["ceiling"] - max_drop:
        reasons.append(
            f"concordance {c:.3f} more than {max_drop} below known-donor ceiling {row['ceiling']:.3f}"
        )
    row["status"] = "FAIL" if reasons else "PASS"
    row["reason"] = "; ".join(reasons)
    return row


def cross_capture(rows, vireo_paths):
    """Lowest pairwise vireo-GT agreement per missing donor across captures."""
    by_donor = {}
    for r in rows:
        if r["status"] in ("PASS", "FAIL") and r["missing_donor"] and "," not in r["missing_donor"]:
            by_donor.setdefault(r["missing_donor"], []).append(r["capture_id"])
    result = {}
    for donor, caps in by_donor.items():
        if len(caps) < 2:
            continue
        gts = {}
        for cap in caps:
            g = read_genotypes(vireo_paths[cap]).get(donor)
            if g:
                gts[cap] = g
        vals = [concordance(gts[a], gts[b])[0] for a, b in combinations(sorted(gts), 2)]
        vals = [v for v in vals if v == v]
        if vals:
            result[donor] = (min(vals), len(gts))
    return result


# ------------------------------------------------------------------------ IO
def _fmt(v):
    if isinstance(v, float):
        return "NA" if v != v else f"{v:.4f}"
    return "NA" if v in ("", None) else str(v)


def write_tsv(rows, path, columns):
    with open(path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for r in rows:
            fh.write("\t".join(_fmt(r.get(c, "")) for c in columns) + "\n")


def read_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        return [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh]


def _shout(msg, log, console):
    """Write to the rule log AND the console: a FAIL must not go unnoticed."""
    banner = "\n".join(["", "!" * 78, msg, "!" * 78, ""])
    print(banner, file=log)
    print(banner, file=console)


def run(sm, log, console):
    if sm.rule == "combine_elimination_gate":
        rows = [r for p in sm.input.gates for r in read_tsv(p)]
        vireo_paths = dict(zip(sm.params.capture_ids, sm.input.vireo))
        xc = cross_capture(rows, vireo_paths)
        for r in rows:
            v = xc.get(r["missing_donor"])
            r["cross_capture_min_concordance"] = v[0] if v else float("nan")
            r["cross_capture_n"] = v[1] if v else ""
        write_tsv(rows, sm.output.tsv, COHORT_COLUMNS)
        return

    p = sm.params
    with open(sm.input.genotyped_donors) as fh:
        genotyped = {line.strip() for line in fh if line.strip()}
    row = evaluate(
        capture_id=sm.wildcards.capture_id,
        pool_donors=p.donors.split(","),
        genotyped=genotyped,
        mapping=read_mapping(sm.input.mapping),
        vireo_gt=read_genotypes(sm.input.vireo_vcf),
        soup_gt=read_genotypes(sm.input.cluster_genotypes),
        min_concordance=float(p.min_concordance),
        max_drop=float(p.max_drop),
        min_sites=int(p.min_sites),
        action=p.action,
    )
    print(row, file=log)
    if row["status"] == "FAIL":
        verb = "WITHHELD (cells set to unassigned)" if p.action == "withhold" else "FLAGGED"
        _shout(
            f"ELIMINATION GATE FAIL [{row['capture_id']}]: by-elimination assignment of "
            f"donor {row['missing_donor']} to souporcell cluster {row['snp_cluster'] or '?'} "
            f"is not supported by its genotype: {row['reason']}. Assignment {verb}.",
            log,
            console,
        )
    write_tsv([row], sm.output.tsv, COLUMNS)


def main(sm):
    stderr = sys.stderr
    with open(sm.log[0], "w") as log:
        sys.stderr = log
        try:
            run(sm, log, console=stderr)
        finally:
            sys.stderr = stderr


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
