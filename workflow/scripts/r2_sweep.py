"""Benchmark the imputation-R2 threshold: sites separating each donor from the rest.

For each threshold in params.imputation.sweep_r2, apply the same rule as
filter_panel.sh (drop imputed sites with R2 < threshold, keep genotyped sites)
to a capture's raw panel and report, per donor:

    n_sites             sites retained (all donors)
    n_imputed           of which imputed
    min_pairwise_sites  retained sites where this donor's called genotype
                        differs from that of its closest pool donor, i.e. the
                        fewest sites separating it from any other pool donor
    closest_donor       that closest pool donor

min_pairwise_sites is what limits how well a donor can be told apart from the
pool. R2 is systematically lower for under-represented ancestries, so an
aggressive cut can thin some donors' separating sites more than others'; the
per-donor spread is what this table is for.

Rule r2_sweep concatenates the per-capture tables.
"""
import gzip
import sys

COLUMNS = [
    "capture_id",
    "min_r2",
    "donor_id",
    "n_sites",
    "n_imputed",
    "min_pairwise_sites",
    "closest_donor",
]


def _open(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


def _info(field):
    out = {}
    for kv in field.split(";"):
        k, _, v = kv.partition("=")
        out[k] = v if _ else True
    return out


def _dosage(gt):
    alleles = gt.replace("|", "/").split("/")
    if any(a in (".", "") for a in alleles):
        return None
    return sum(a != "0" for a in alleles)


def read_panel(path, donors=None):
    """Return (donors, sites); each site is (imputed, r2 or None, [dosage per donor])."""
    declared = set()
    filters = set()
    sites = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID="):
                declared.add(line[len("##INFO=<ID=") :].split(",", 1)[0])
                continue
            if line.startswith("##FILTER=<ID="):
                filters.add(line[len("##FILTER=<ID=") :].split(",", 1)[0])
                continue
            if line.startswith("#CHROM"):
                header = line.rstrip("\n").split("\t")
                names = header[9:]
                keep = [i for i, n in enumerate(names) if donors is None or n in donors]
                donors_out = [names[i] for i in keep]
                continue
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            info = _info(f[7])
            if "IMPUTED" in declared:
                imputed = "IMPUTED" in info
            elif "TYPED" in declared:
                imputed = "TYPED" not in info
            elif "GENOTYPED" in filters:
                # Minimac3 marks genotyped sites in FILTER, not INFO
                imputed = not {"GENOTYPED", "GENOTYPED_ONLY"} & set(f[6].split(";"))
            else:
                imputed = True
            r2 = info.get("R2") if "R2" in declared else None
            r2 = float(r2) if isinstance(r2, str) and r2 not in ("", ".") else None
            gt_idx = f[8].split(":").index("GT")
            calls = [f[9 + i].split(":")[gt_idx] for i in keep]
            sites.append((imputed, r2, [_dosage(c) for c in calls]))
    return donors_out, sites


def retained(site, min_r2):
    """Same rule as filter_panel.sh: drop imputed sites with R2 below threshold."""
    imputed, r2, _ = site
    return not (imputed and r2 is not None and r2 < min_r2)


def sweep(donors, sites, thresholds, capture_id):
    rows = []
    for t in thresholds:
        kept = [s for s in sites if retained(s, t)]
        n_imputed = sum(1 for s in kept if s[0])
        n = len(donors)
        diff = [[0] * n for _ in range(n)]
        for _, _, dos in kept:
            for i in range(n):
                if dos[i] is None:
                    continue
                for j in range(i + 1, n):
                    if dos[j] is not None and dos[j] != dos[i]:
                        diff[i][j] += 1
                        diff[j][i] += 1
        for i, donor in enumerate(donors):
            others = [(diff[i][j], donors[j]) for j in range(n) if j != i]
            sites_min, closest = min(others) if others else (0, "")
            rows.append([capture_id, t, donor, len(kept), n_imputed, sites_min, closest])
    return rows


def write_tsv(rows, path, header=COLUMNS):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(str(x) for x in r) + "\n")


def main(sm):
    sys.stderr = open(sm.log[0], "w")
    if sm.rule == "r2_sweep":
        with open(sm.output.tsv, "w") as out:
            out.write("\t".join(COLUMNS) + "\n")
            for p in sm.input.tables:
                with open(p) as fh:
                    next(fh)
                    out.writelines(fh)
        return
    donors = set(sm.params.donors.split(","))
    names, sites = read_panel(sm.input.vcf, donors)
    rows = sweep(names, sites, [float(t) for t in sm.params.thresholds], sm.wildcards.capture_id)
    write_tsv(rows, sm.output.tsv)


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
