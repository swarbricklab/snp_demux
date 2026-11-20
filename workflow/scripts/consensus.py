"""Per-barcode consensus of the four demultiplexers for one capture.

Replaces Demuxafy's ``Combine_Results.R`` and ``add_vireosnp_to_combined_results.R``.

Steps
  1. Read each method's per-barcode call and normalise it to
     ``<M>_DropletType`` in {singlet, doublet, unassigned}, plus a donor
     (Demuxlet, Vireo) or a cluster id (Freemuxlet, Souporcell). A method's own
     ``unassigned`` / ``doublet`` call is kept as it is.
  2. Name Freemuxlet / Souporcell clusters after donors, by ``cluster_naming``:
       genotype (default)
         From the vireoSNP GPb_diff between each cluster's genotype and each
         donor's known genotype: a cluster and a donor that are each other's
         closest match, with the cluster clearly closer to that donor than to
         the next (``min_genotype_margin``). This keeps the cluster methods
         independent of Vireo. See ``genotype_map``.
       vireo_agreement
         Demuxafy's method: a cluster takes the Vireo donor shared by at least
         ``pct_agreement`` of its barcodes. Only real clusters are named (not
         the ``doublet`` / ``unassigned`` calls), only after donors, and only
         when exactly one donor qualifies. See ``cluster_map``.
     Either way, a singlet in an unnamed cluster abstains (``unassigned``): a
     cluster no donor accounts for is as likely to be doublets or noise as a
     donor, so its singlet call is not evidence.
  3. Consensus (``Consensus_*``), by ``method``:
       informative_majority (default)
         Votes are the methods' droplet types; ``unassigned`` abstains.
         Fewer than ``min_informative`` singlet/doublet votes -> unassigned;
         otherwise the majority of singlet vs doublet, a tie -> unassigned.
         A singlet takes the donor most methods name; a tie -> unassigned.
         ``recover_from`` (default souporcell and vireo): a barcode that
         exactly one method calls, as a singlet of a named donor, is
         ``recovered`` with that donor rather than ``unassigned`` when that
         method is listed; ``Consensus_recovered_by`` names it. On bcarc, lone
         souporcell / vireo calls matched an independent call's donor ~85-88%
         of the time, freemuxlet ~72%, demuxlet ~42%. A lone doublet call
         stays ``unassigned``.
       at_least_half_singlet
         Demuxafy's AtLeastHalfSinglet (below).
     Vote counts are written as ``Consensus_n_singlet``, ``Consensus_n_doublet``
     and ``Consensus_n_informative``, so consumers can apply their own rule.
  4. ``AtLeastHalfSinglet_*`` is always written as well: a singlet if at least
     half the methods call a singlet, otherwise a doublet (so ``unassigned``
     counts against singlet). Kept for comparison and for existing consumers.
  5. Append souporcell's genotype-annotated call (vireoSNP) as
     ``Souporcell_vireosnp_*``. It does not vote.

Two Demuxafy behaviours are deliberately not reproduced (see #65): it turned
every barcode whose cluster it could not name, including a method's own
``unassigned`` calls, into a ``doublet``; and AtLeastHalfSinglet has no
``unassigned`` outcome, so barcodes most methods could not call became doublets.
``method: at_least_half_singlet`` restores the second; nothing restores the first.

Rows follow demuxlet's barcode order; a barcode missing from another method gets
``unassigned`` for that method.

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
    "Consensus_DropletType",
    "Consensus_Individual_Assignment",
    "Consensus_n_singlet",
    "Consensus_n_doublet",
    "Consensus_n_informative",
    "Consensus_recovered_by",
]
RULES = ("informative_majority", "at_least_half_singlet")
NAMINGS = ("genotype", "vireo_agreement")


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
    """Name clusters by agreement with ``reference`` over ``barcodes``.

    Returns (mapping cluster -> donor, report rows). A cluster's agreement with a
    donor is the share of the cluster's barcodes (as seen in ``barcodes``) that the
    reference assigns to that donor; the denominator includes barcodes the
    reference calls doublet/unassigned or lacks. A cluster is named only if
    exactly one donor reaches ``pct_agreement`` (two can at <= 0.5).
    """
    totals = Counter()
    pairs = Counter()
    for bc in barcodes:
        cl = clusters.get(bc, MISSING)
        if cl is MISSING or cl in SENTINELS:
            continue
        totals[cl] += 1
        ref = reference.get(bc, MISSING)
        if ref is not MISSING and ref not in SENTINELS:
            pairs[(cl, ref)] += 1
    mapping, report = {}, []
    for cl in sorted(totals, key=lambda x: (-totals[x], x)):
        frac = {ref: n / totals[cl] for (c, ref), n in pairs.items() if c == cl}
        ok = sorted(ref for ref, f in frac.items() if f >= pct_agreement)
        best = max(frac, key=lambda r: (frac[r], r)) if frac else None
        if len(ok) == 1:
            mapping[cl] = ok[0]
        report.append((cl, totals[cl], best, frac.get(best, 0.0),
                       ok[0] if len(ok) == 1 else ("ambiguous" if ok else "unnamed")))
    return mapping, report


def read_scores(rows):
    """GPb_diff matrix (``known genotype`` + ``snp_cluster_<id>`` columns) ->
    {cluster: {donor: score}}; lower is a closer match; NA scores are dropped."""
    scores = OrderedDict()
    for r in rows:
        for col, value in r.items():
            if not col.startswith("snp_cluster_"):
                continue
            try:
                score = float(value)
            except (TypeError, ValueError):
                continue
            if score == score:
                scores.setdefault(col[len("snp_cluster_"):], {})[r["known genotype"]] = score
    return scores


def genotype_map(barcodes, clusters, scores, min_margin, ungenotyped=None):
    """Name clusters by genotype: a cluster and a donor that clearly match each other.

    ``scores`` is {cluster: {donor: GPb_diff}}, lower is closer. A cluster is
    named after its closest donor when
      - the relative margin, (second - best) / second, is at least ``min_margin``,
        so the cluster is clearly closer to that donor than to any other; and
      - the cluster is also that donor's closest cluster (a mutual best match).
    The second condition keeps a doublet cluster, whose genotype is a mixture and
    so lies nearest one of its donors, from taking that donor's name; the cost is
    that the smaller half of a donor split over two clusters stays unnamed.
    Unnamed clusters' singlets abstain (see ``assign_clusters``).

    A pool donor with no known genotype cannot be matched, and its cluster
    tends to sit nearest some genotyped donor instead. Vireo is the only method
    that learns such a donor, so ``ungenotyped`` maps a cluster to an
    ungenotyped donor when Vireo agreement names it so (see ``cluster_map``);
    that name replaces the genotype one. (Naming the one leftover cluster "by
    elimination", as annotate.py does, mislabels a split donor's second cluster.)

    Returns (mapping cluster -> donor, report rows) like ``cluster_map``.
    """
    sizes = Counter()
    for bc in barcodes:
        cl = clusters.get(bc, MISSING)
        if cl is not MISSING and cl not in SENTINELS:
            sizes[cl] += 1
    closest_cluster = {}
    for cl, row in scores.items():
        for donor, score in row.items():
            if donor not in closest_cluster or score < scores[closest_cluster[donor]][donor]:
                closest_cluster[donor] = cl
    mapping, report = {}, []
    for cl in sorted(set(sizes) | set(scores), key=lambda x: (-sizes[x], x)):
        ranked = sorted(scores.get(cl, {}).items(), key=lambda kv: (kv[1], kv[0]))
        if len(ranked) < 2 or ranked[1][1] <= 0:
            report.append([cl, sizes[cl], ranked[0][0] if ranked else None, float("nan"), "unnamed"])
            continue
        (best, b), (_, second) = ranked[0], ranked[1]
        margin = (second - b) / second
        if margin < min_margin:
            named = "unnamed"
        elif closest_cluster[best] != cl:
            named = "unnamed (%s is closer to cluster %s)" % (best, closest_cluster[best])
        else:
            mapping[cl] = named = best
        report.append([cl, sizes[cl], best, margin, named])
    for r in report:
        donor = (ungenotyped or {}).get(r[0])
        if donor is not None:
            mapping[r[0]] = donor
            r[4] = "%s (no genotype; named by vireo agreement)" % donor
    return mapping, [tuple(r) for r in report]


def assign_clusters(calls, mapping):
    """Give each singlet of a cluster method the donor its cluster is named after.

    Doublet / unassigned calls keep that as their assignment; a singlet in an
    unnamed cluster becomes ``unassigned`` (its cluster id is kept).
    """
    for call in calls.values():
        if call["DropletType"] not in SENTINELS:
            donor = mapping.get(call["Cluster"])
            if donor is None:
                call["DropletType"] = "unassigned"
        call["Individual_Assignment"] = (
            call["DropletType"] if call["DropletType"] in SENTINELS else donor
        )


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


def votes(row):
    types = [row[m + "_DropletType"] for m in METHODS]
    return types.count("singlet"), types.count("doublet")


def informative_majority(row, min_informative=2, recover_from=()):
    """Majority of the methods that made a call -> (droplet type, assignment).

    A barcode whose only call is a singlet from a method in ``recover_from``
    (below ``min_informative``) is ``recovered`` with that donor.
    """
    n_singlet, n_doublet = votes(row)
    if n_singlet == 1 and n_doublet == 0 and min_informative > 1:
        lone = next(m for m in METHODS if row[m + "_DropletType"] == "singlet")
        donor = row[lone + "_Individual_Assignment"]
        if lone in recover_from and donor not in SENTINELS:
            return "recovered", donor
    if n_singlet + n_doublet < min_informative or n_singlet == n_doublet:
        return "unassigned", "unassigned"
    if n_doublet > n_singlet:
        return "doublet", "doublet"
    donor = plurality(row[m + "_Individual_Assignment"] for m in METHODS)
    return ("singlet", donor) if donor is not None else ("unassigned", "unassigned")


def at_least_half_singlet(row):
    """Demuxafy's AtLeastHalfSinglet over METHODS -> (droplet type, assignment)."""
    n_singlet = sum(row[m + "_DropletType"] == "singlet" for m in METHODS)
    if n_singlet < len(METHODS) / 2:
        return "doublet", "doublet"
    donor = plurality(row[m + "_Individual_Assignment"] for m in METHODS)
    return "singlet", "unassigned" if donor is None else donor


# ------------------------------------------------------------------ build ----
def build(demuxlet, freemuxlet, souporcell, vireo, vireosnp, pct_agreement,
          method="informative_majority", min_informative=2, log=None,
          cluster_naming="vireo_agreement", scores=None, min_genotype_margin=0.1,
          ungenotyped=(), recover_from=()):
    """Inputs are the parsed tables (lists of dicts); returns output rows.

    ``scores`` maps each cluster method to its parsed GPb_diff matrix (rows of
    dicts); it is required when ``cluster_naming`` is ``genotype``.
    """
    if method not in RULES:
        raise ValueError("consensus method must be one of %s, not %r" % (RULES, method))
    if cluster_naming not in NAMINGS:
        raise ValueError("cluster_naming must be one of %s, not %r" % (NAMINGS, cluster_naming))
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
        if cluster_naming == "genotype":
            agreement = cluster_map(barcodes, clusters, vireo_ref, pct_agreement)[0]
            no_genotype = {cl: d for cl, d in agreement.items() if d in ungenotyped}
            mapping, report = genotype_map(barcodes, clusters, read_scores(scores[m]),
                                           min_genotype_margin, no_genotype)
            how = "by genotype (relative GPb_diff margin >= %s)" % min_genotype_margin
        else:
            mapping, report = cluster_map(barcodes, clusters, vireo_ref, pct_agreement)
            how = "by agreement with vireo (pct_agreement %s)" % pct_agreement
        assign_clusters(calls[m], mapping)
        if log:
            print("%s clusters named %s:" % (m, how), file=log)
            for cl, n, best, value, named in report:
                print("  cluster %-6s n=%-6d best=%s (%.3f) -> %s" % (cl, n, best, value, named), file=log)

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
        if method == "informative_majority":
            dtype, donor = informative_majority(row, min_informative, recover_from)
        else:
            dtype, donor = row[CONSENSUS + "_DropletType"], row[CONSENSUS + "_Individual_Assignment"]
        n_singlet, n_doublet = votes(row)
        row.update({
            "Consensus_DropletType": dtype,
            "Consensus_Individual_Assignment": donor,
            "Consensus_n_singlet": n_singlet,
            "Consensus_n_doublet": n_doublet,
            "Consensus_n_informative": n_singlet + n_doublet,
            "Consensus_recovered_by": next(
                (m for m in METHODS if row[m + "_DropletType"] == "singlet"), "NA"
            ) if dtype == "recovered" else "NA",
        })
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
    scores = {
        "Freemuxlet": read_table(sm.input.freemuxlet_scores),
        "Souporcell": read_table(sm.input.souporcell_scores),
    }
    genotyped = {line.strip() for line in open(sm.input.genotyped_donors) if line.strip()}
    pool = [d for d in str(sm.params.pool_donors).split(",") if d]
    ungenotyped = [d for d in pool if d not in genotyped]
    if ungenotyped:
        print("pool donors without a genotype: %s" % ", ".join(ungenotyped), file=log)
    rows = build(pct_agreement=float(sm.params.pct_agreement), method=sm.params.method,
                 min_informative=int(sm.params.min_informative), log=log,
                 cluster_naming=sm.params.cluster_naming, scores=scores,
                 min_genotype_margin=float(sm.params.min_genotype_margin),
                 ungenotyped=ungenotyped,
                 recover_from=[m.capitalize() for m in sm.params.recover_from], **tables)
    write_tsv(rows, sm.output.calls)

    only = {r["Barcode"] for r in rows} ^ {r["barcode"] for r in tables["vireosnp"]}
    if only:
        print("WARNING: %d barcode(s) in only one of demuxlet / vireoSNP" % len(only), file=log)
    for m in METHODS + (CONSENSUS, "Consensus"):
        counts = Counter(r[m + "_DropletType"] for r in rows)
        print("%-20s %s" % (m, dict(sorted(counts.items()))), file=log)
    by = Counter(r["Consensus_recovered_by"] for r in rows if r["Consensus_DropletType"] == "recovered")
    if by:
        print("%-20s %s" % ("recovered by", dict(sorted(by.items()))), file=log)
    log.close()


if "snakemake" in globals():
    main(snakemake)  # noqa: F821
