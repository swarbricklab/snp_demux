"""Benchmark the imputation-R2 threshold: donor-discriminating sites per donor.

For each threshold in params.imputation.sweep_r2, apply the same rule as
filter_panel.sh (drop imputed sites with R2 < threshold, keep genotyped sites)
to a capture's raw panel and count, per donor:

    n_sites            sites retained (all donors)
    n_imputed          of which imputed
    n_discriminating   retained sites where this donor has a called genotype
                       that differs from at least one other pool donor's

R2 is systematically lower for under-represented ancestries, so an aggressive
cut can thin some donors' discriminating sites far more than others'; the
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
    "n_discriminating",
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
    sites = []
    with _open(path) as fh:
        for line in fh:
            if line.startswith("##INFO=<ID="):
                declared.add(line[len("##INFO=<ID=") :].split(",", 1)[0])
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
        disc = [0] * len(donors)
        for _, _, dos in kept:
            called = [d for d in dos if d is not None]
            if len(set(called)) < 2:
                continue  # nobody differs from anybody here
            for i, d in enumerate(dos):
                if d is not None and any(o is not None and o != d for o in dos):
                    disc[i] += 1
        for donor, n in zip(donors, disc):
            rows.append([capture_id, t, donor, len(kept), n_imputed, n])
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
