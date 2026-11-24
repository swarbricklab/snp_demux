"""Append vireo's and souporcell's per-barcode posteriors to the consensus table.

Input is the per-barcode, per-method call table (donor space) that
``wrangle_demuxafy_assignments`` turns into ``sample_barcode_assignments.tsv``.
Every input row and column is passed through unchanged; six columns are
appended:

    Vireo_prob_max                vireo donor_ids.tsv ``prob_max``: posterior of
                                  the best singlet donor (joint with the doublet
                                  hypotheses, so prob_max + prob_doublet <= 1)
    Vireo_prob_doublet            vireo ``prob_doublet``: posterior of the most
                                  likely donor PAIR (vireo reports the max over
                                  pairs, not their sum)
    Vireo_best_singlet            vireo ``best_singlet``
    Souporcell_singlet_posterior  troublet ``singlet_posterior``: posterior of the
                                  best cluster, normalised over clusters (i.e.
                                  given that the barcode is a singlet)
    Souporcell_doublet_posterior  troublet ``doublet_posterior``: best donor pair
                                  vs best singlet, under troublet's doublet prior
                                  (0.5 unless --doublet_prior is passed)
    Souporcell_best_singlet       the cluster with the highest of troublet's
                                  per-cluster singlet log-probabilities (the
                                  ``assignment`` of a souporcell singlet), named
                                  after the donor the consensus table gives that
                                  cluster

Posterior values are copied verbatim (no float re-formatting). The status
thresholds both tools apply are on these same numbers (vireo: prob_max < 0.9 ->
unassigned, prob_doublet >= 0.9 -> doublet; troublet: singlet_posterior > 0.9
for a singlet, doublet_posterior >= 0.9 for a doublet).

``*_best_singlet`` holds donor ids of this capture only, which
``wrangle_assignments.R`` then recodes to sample ids along with the
``*_Individual_Assignment`` columns. Anything that is not a donor of the capture
becomes NA: a vireo donor that was never relabelled, a souporcell cluster that
the consensus left unnamed, or one it named inconsistently. Souporcell's cluster
names are read from the table itself (``Souporcell_Cluster`` ->
``Souporcell_Individual_Assignment``), so they follow whichever cluster naming
produced the table; nothing here re-derives them. A named cluster with no
barcode assigned to it in the table cannot be read off and is NA as well.

A barcode missing from vireo's or souporcell's output gets NA in that tool's
columns.

Pure standard library.
"""
import csv
import gzip
import sys
from collections import defaultdict

NA = "NA"
VIREO_COLUMNS = ["Vireo_prob_max", "Vireo_prob_doublet", "Vireo_best_singlet"]
SOUPORCELL_COLUMNS = [
    "Souporcell_singlet_posterior",
    "Souporcell_doublet_posterior",
    "Souporcell_best_singlet",
]
COLUMNS = VIREO_COLUMNS + SOUPORCELL_COLUMNS


def read_table(path, delimiter="\t"):
    """(header, rows) of a (optionally gzipped) delimited file; rows are dicts."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        return list(reader.fieldnames or []), list(reader)


def capture_donors(samples, capture_id):
    """donor_ids declared for a capture in samples.csv."""
    return {r["donor_id"] for r in samples if r["capture_id"] == capture_id}


def souporcell_cluster_names(rows, donors):
    """cluster -> donor, as the consensus table names souporcell's clusters.

    Read from rows whose Souporcell_Individual_Assignment is a donor of the
    capture. A cluster the table gives two different donors (possible with
    Demuxafy's duplicate rows) is ambiguous and left out, as are clusters no
    such row mentions.
    """
    seen = defaultdict(set)
    for r in rows:
        donor = r.get("Souporcell_Individual_Assignment")
        if donor in donors:
            seen[r["Souporcell_Cluster"]].add(donor)
    return {c: next(iter(d)) for c, d in seen.items() if len(d) == 1}


def vireo_posteriors(rows, donors):
    """barcode -> [prob_max, prob_doublet, best_singlet] from donor_ids.tsv."""
    out = {}
    for r in rows:
        if r["cell"] in out:
            continue
        best = r["best_singlet"]
        out[r["cell"]] = [
            r["prob_max"],
            r["prob_doublet"],
            best if best in donors else NA,
        ]
    return out


def souporcell_best_cluster(row, cluster_columns):
    """Index (as a string) of the highest per-cluster log-probability."""
    values = [float(row[c]) for c in cluster_columns]
    return str(max(range(len(values)), key=values.__getitem__))


def souporcell_posteriors(header, rows, names):
    """barcode -> [singlet_posterior, doublet_posterior, best_singlet donor]."""
    for col in ("singlet_posterior", "doublet_posterior"):
        if col not in header:
            raise ValueError(
                f"souporcell clusters.tsv has no {col} column (troublet < 2.5?)"
            )
    cluster_columns = [c for c in header if c.startswith("cluster")]
    if [c[len("cluster"):] for c in cluster_columns] != [
        str(i) for i in range(len(cluster_columns))
    ]:
        raise ValueError(f"unexpected souporcell cluster columns: {cluster_columns}")
    out = {}
    for r in rows:
        if r["barcode"] in out:
            continue
        cluster = souporcell_best_cluster(r, cluster_columns)
        out[r["barcode"]] = [
            r["singlet_posterior"],
            r["doublet_posterior"],
            names.get(cluster, NA),
        ]
    return out


def add_posteriors(rows, vireo, souporcell):
    """Append COLUMNS to each consensus row (in place); returns rows."""
    missing = [NA, NA, NA]
    for r in rows:
        barcode = r["Barcode"]
        r.update(zip(VIREO_COLUMNS, vireo.get(barcode, missing)))
        r.update(zip(SOUPORCELL_COLUMNS, souporcell.get(barcode, missing)))
    return rows


def write_table(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=header, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def run(table, vireo_path, souporcell_path, samples_path, capture_id, out, log):
    header, rows = read_table(table)
    clash = [c for c in COLUMNS if c in header]
    if clash:
        raise ValueError(f"consensus table already has {clash}")
    _, samples = read_table(samples_path, delimiter=",")
    donors = capture_donors(samples, capture_id)

    _, v_rows = read_table(vireo_path)
    vireo = vireo_posteriors(v_rows, donors)
    s_header, s_rows = read_table(souporcell_path)
    names = souporcell_cluster_names(rows, donors)
    souporcell = souporcell_posteriors(s_header, s_rows, names)

    add_posteriors(rows, vireo, souporcell)
    write_table(out, header + COLUMNS, rows)

    print(f"{capture_id}: {len(rows)} barcodes", file=log)
    print(f"souporcell cluster names read from the table: {names}", file=log)
    for tool, calls in (("vireo", vireo), ("souporcell", souporcell)):
        n_missing = sum(r["Barcode"] not in calls for r in rows)
        if n_missing:
            print(f"{n_missing} barcodes have no {tool} row", file=log)


def main(sm):
    with open(sm.log[0], "w") as log:
        run(
            sm.input.table,
            sm.input.vireo,
            sm.input.souporcell,
            sm.input.samples,
            sm.wildcards.capture_id,
            sm.output.table,
            log,
        )


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
