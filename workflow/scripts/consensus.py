"""Per-barcode consensus of the four demultiplexers for one capture.

Replaces Demuxafy's ``Combine_Results.R`` (``--method AtLeastHalfSinglet
--ref Vireo``) and ``add_vireosnp_to_combined_results.R``. This is a direct port
of the code paths that configuration exercises, kept output-equivalent so the
switch can be verified against the golden; behaviour changes come separately.

Steps
  1. Read each method's per-barcode call and normalise it to
     ``<M>_DropletType`` in {singlet, doublet, unassigned}, plus a donor
     (Demuxlet, Vireo) or a cluster id (Freemuxlet, Souporcell).
  2. Map Freemuxlet / Souporcell clusters to donors by agreement with Vireo: a
     cluster takes the Vireo assignment shared by at least ``pct_agreement`` of
     its barcodes. As in Demuxafy, a barcode whose cluster has no such mapping
     (including the pseudo-clusters ``doublet`` and ``unassigned``) becomes a
     ``doublet``.
  3. AtLeastHalfSinglet: a singlet if at least half the methods call a singlet,
     assigned to the donor most methods name (a tie leaves it unassigned);
     otherwise a doublet.
  4. Append souporcell's genotype-annotated call (vireoSNP) as
     ``Souporcell_vireosnp_*``. It does not vote.

Rows follow demuxlet's barcode order (Demuxafy left-joins onto demuxlet); a
barcode missing from another method gets ``unassigned`` for that method.

Where Demuxafy's join emits several rows for one barcode (a cluster that ties
for the mapping, possible when pct_agreement <= 0.5), downstream kept the first;
this port emits that row only.

Pure standard library (runs in any container with python3 >= 3.6).
"""
import csv
import gzip
import re
import sys
from collections import Counter, OrderedDict

METHODS = ("Demuxlet", "Freemuxlet", "Souporcell", "Vireo")
CLUSTER_METHODS = ("Freemuxlet", "Souporcell")
SENTINELS = ("doublet", "unassigned")
POPSCLE_TYPES = {"SNG": "singlet", "DBL": "doublet", "AMB": "unassigned"}
MISSING = None  # a barcode absent from a table / an unmatched join

CONSENSUS = "AtLeastHalfSinglet"
VIREOSNP = "Souporcell_vireosnp"
COLUMNS = [
    "Barcode",
    "Demuxlet_DropletType",
    "Demuxlet_Individual_Assignment",
    "Freemuxlet_Cluster",
    "Freemuxlet_Individual_Assignment",
    "Freemuxlet_DropletType",
    "Souporcell_Cluster",
    "Souporcell_Individual_Assignment",
    "Souporcell_DropletType",
    "Vireo_Individual_Assignment",
    "Vireo_DropletType",
    CONSENSUS + "_DropletType",
    CONSENSUS + "_Individual_Assignment",
    VIREOSNP + "_DropletType",
    VIREOSNP + "_Individual_Assignment",
]


# ---------------------------------------------------------------- reading ----
def read_table(path):
    """Rows of a (optionally gzipped) TSV as dicts of strings."""
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _sentinel_or(dtype, value):
    return dtype if dtype in SENTINELS else value


def read_popscle(rows, id_column):
    """demuxlet .best / freemuxlet .samples.gz -> OrderedDict barcode -> call.

    BEST.GUESS is "a,b,p" for demuxlet and "a,b" for freemuxlet; a singlet's
    donor/cluster is the first field (Demuxafy strips ",b,p" and ",b" with
    regexes; for these formats both reduce to the first field).
    """
    strip = re.compile(r",.+,.\..+") if id_column == "Individual_Assignment" else re.compile(r",.+")
    calls = OrderedDict()
    for r in rows:
        dtype = POPSCLE_TYPES.get(r["DROPLET.TYPE"], r["DROPLET.TYPE"])
        calls[r["BARCODE"]] = {
            "DropletType": dtype,
            id_column: _sentinel_or(dtype, strip.sub("", r["BEST.GUESS"])),
        }
    return calls


def read_souporcell(rows):
    calls = OrderedDict()
    for r in rows:
        calls[r["barcode"]] = {
            "DropletType": r["status"],
            "Cluster": _sentinel_or(r["status"], r["assignment"]),
        }
    return calls


def read_vireo(rows):
    calls = OrderedDict()
    for r in rows:
        donor = r["donor_id"]
        calls[r["cell"]] = {
            "Individual_Assignment": donor,
            "DropletType": donor if donor in SENTINELS else "singlet",
        }
    return calls


# ---------------------------------------------------------------- mapping ----
def cluster_map(barcodes, clusters, reference, pct_agreement):
    """cluster -> reference assignment, by agreement over ``barcodes``.

    For each (cluster, reference) pair: Proportion = pair count / cluster size.
    Pairs with Proportion >= pct_agreement are candidates. Demuxafy keeps every
    candidate (a join on them duplicates barcodes); downstream kept the first row,
    i.e. the candidate with the highest count, ties going to the pair seen *last*
    in barcode order (Demuxafy sorts with rev(order(Count))). That is the one
    returned here. ``MISSING`` is a value like any other (a barcode absent from a
    table), exactly as data.table groups NA.
    """
    totals = Counter()
    pairs = OrderedDict()  # (cluster, ref) -> count, in order of first appearance
    for bc in barcodes:
        key = (clusters.get(bc), reference.get(bc))
        totals[key[0]] += 1
        pairs[key] = pairs.get(key, 0) + 1
    ranked = sorted(
        ((count, i, key) for i, (key, count) in enumerate(pairs.items())),
        reverse=True,
    )
    mapping = {}
    for count, _, (cluster, ref) in ranked:
        if count / totals[cluster] >= pct_agreement and cluster not in mapping:
            mapping[cluster] = ref
    return mapping


def assign_clusters(calls, mapping):
    """Give each barcode of a cluster method the donor its cluster maps to."""
    for call in calls.values():
        donor = mapping.get(call["Cluster"], MISSING)
        call["Individual_Assignment"] = "doublet" if donor is MISSING else donor
        if call["Individual_Assignment"] == "doublet":
            call["DropletType"] = "doublet"


# -------------------------------------------------------------- consensus ----
def plurality(values):
    """The single most frequent donor, ignoring sentinels; None on a tie/none."""
    counts = Counter(v for v in values if v not in SENTINELS)
    if not counts:
        return None
    top = counts.most_common()
    if len(top) > 1 and top[0][1] == top[1][1]:
        return None
    return top[0][0]


def at_least_half_singlet(row):
    """Demuxafy's AtLeastHalfSinglet over METHODS -> (droplet type, assignment)."""
    n_singlet = sum(row[m + "_DropletType"] == "singlet" for m in METHODS)
    if n_singlet < len(METHODS) / 2:
        return "doublet", "doublet"
    donor = plurality(row[m + "_Individual_Assignment"] for m in METHODS)
    return "singlet", "unassigned" if donor is None else donor


# ------------------------------------------------------------------ build ----
def build(demuxlet, freemuxlet, souporcell, vireo, vireosnp, pct_agreement):
    """Inputs are the parsed tables (lists of dicts); returns output rows."""
    calls = {
        "Demuxlet": read_popscle(demuxlet, "Individual_Assignment"),
        "Freemuxlet": read_popscle(freemuxlet, "Cluster"),
        "Souporcell": read_souporcell(souporcell),
        "Vireo": read_vireo(vireo),
    }
    barcodes = list(calls["Demuxlet"])
    vireo_ref = {bc: c["Individual_Assignment"] for bc, c in calls["Vireo"].items()}
    for m in CLUSTER_METHODS:
        clusters = {bc: c["Cluster"] for bc, c in calls[m].items()}
        assign_clusters(calls[m], cluster_map(barcodes, clusters, vireo_ref, pct_agreement))

    snp = {r["barcode"]: r for r in vireosnp}
    out = []
    for bc in barcodes:
        row = {"Barcode": bc}
        for m in METHODS:
            call = calls[m].get(bc, {})
            for field in ("Cluster", "Individual_Assignment", "DropletType"):
                if m in CLUSTER_METHODS or field != "Cluster":
                    value = call.get(field, MISSING)
                    row[m + "_" + field] = "unassigned" if value is MISSING else value
        row[CONSENSUS + "_DropletType"], row[CONSENSUS + "_Individual_Assignment"] = (
            at_least_half_singlet(row)
        )
        s = snp.get(bc)
        row[VIREOSNP + "_DropletType"] = s["status"] if s else "NA"
        row[VIREOSNP + "_Individual_Assignment"] = s["assignment"] if s else "NA"
        out.append(row)
    return out


def write_tsv(rows, path):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main(sm):
    log = open(sm.log[0], "w")
    tables = {k: read_table(sm.input[k]) for k in ("demuxlet", "freemuxlet", "souporcell", "vireo", "vireosnp")}
    rows = build(pct_agreement=float(sm.params.pct_agreement), **tables)
    write_tsv(rows, sm.output.calls)

    only = {r["Barcode"] for r in rows} ^ {r["barcode"] for r in tables["vireosnp"]}
    if only:
        print("WARNING: %d barcode(s) in only one of demuxlet / vireoSNP" % len(only), file=log)
    for m in METHODS + (CONSENSUS,):
        counts = Counter(r[m + "_DropletType"] for r in rows)
        print("%-20s %s" % (m, dict(sorted(counts.items()))), file=log)
    log.close()


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
