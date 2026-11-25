"""Read the calls of any demultiplexing output into one schema; the `calls` command.

Every reader returns one row per barcode with these text columns:

    barcode      the cell barcode as written
    status       singlet, doublet, unassigned or recovered; any other word
                 is kept as written (e.g. "NA" in a file that writes only
                 singlets). hier-v2's "low_cell" becomes unassigned.
    status_raw   the status word as the file wrote it
    donor        the donor of a singlet or recovered cell; "" otherwise
    pair         a doublet's two donors in string order joined by a comma
                 ("d0,d2"); "" when the file names none. vireo writes a best
                 pair for every cell, so its rows carry one whatever the call
    confidence   the file's own score for the call ("" when it has none):
                 vireo prob_max, the consensus label_confidence, the Scorer's margin
    cluster      souporcell's cluster ("3", or "1/4" for a doublet) where the
                 file has one; "" otherwise
    source       the kind the rows were read as

Kinds:

    vireo           donor_ids.tsv. The pair comes from best_doublet, corrected
                    with the header of prob_doublet.tsv.gz (vireoSNP 0.5.7
                    names the wrong pair; see vireo_files). Without that file
                    the written name is used as is.
    souporcell      clusters.tsv (or troublet's output). Clusters become
                    donors through a cluster map; without one a singlet's
                    donor is its cluster number and a doublet's pair its two
                    cluster numbers.
    hier            hier-v2's hierarchical_consensus.tsv.
    assignment      cell_assignment.tsv (barcode, assignment, status). A
                    status other than singlet, doublet, unassigned, recovered
                    or low_cell (snp_demux main writes NA) is split into
                    doublet and unassigned with sample_barcode_assignments.tsv
                    when it is given (Consensus_DropletType if present, else
                    AtLeastHalfSinglet_DropletType).
    demuxafy:TOOL   one tool's columns of sample_barcode_assignments.tsv:
                    TOOL_DropletType and TOOL_Individual_Assignment, e.g.
                    demuxafy:Vireo, demuxafy:Souporcell.
    scorer          the Scorer's per-cell table (vireo's columns plus best_pair
                    and margin).
    generic:barcode=COL,status=COL[,donor=COL][,pair=COL][,confidence=COL]
                    any table; values are taken as written.

Labels (vireo's donor names, Demuxafy's sample ids) become donors through an
optional label map: a table with columns label and donor, or sample_id and
donor_id. A label not in the map is matched against an optional regular
expression whose first group is the donor, and is otherwise kept as written.
"""
import os
import re

import pandas as pd

from demux_tools._io import (UsageError, check_output_path, open_text, read_table, start_clock,
                             write_provenance, write_table)

SCHEMA = ["barcode", "status", "status_raw", "donor", "pair", "confidence", "cluster", "source"]
STATUSES = ("singlet", "doublet", "unassigned", "recovered")
_STATUS_WORDS = {"singlet": "singlet", "doublet": "doublet", "unassigned": "unassigned",
                 "recovered": "recovered", "low_cell": "unassigned"}
KINDS = ("vireo", "souporcell", "hier", "assignment", "demuxafy:TOOL", "scorer", "generic:...")


# ---------------------------------------------------------------------------
# Small readers


def sniff_sep(path):
    """"\\t" when the first line of a text file holds a tab, else ","."""
    with open_text(path) as f:
        first = f.readline()
    return "\t" if "\t" in first else ","


def read_any_table(path):
    """A tab- or comma-separated table (chosen from its first line), every column text."""
    return read_table(path, sep=sniff_sep(path))


def pair_key(a, b):
    """Two donor names in string order joined by a comma."""
    return ",".join(sorted((a, b)))


def read_label_map(path):
    """{label: donor} from a table with columns label and donor, or sample_id and donor_id.

    A label listed twice with two donors raises ValueError.
    """
    t = read_any_table(path)
    for lab, don in (("label", "donor"), ("sample_id", "donor_id")):
        if lab in t.columns and don in t.columns:
            break
    else:
        raise ValueError("%s: needs columns label and donor, or sample_id and donor_id" % path)
    out = {}
    for a, b in zip(t[lab], t[don]):
        if a in out and out[a] != b:
            raise ValueError("%s: label %r maps to %r and %r" % (path, a, out[a], b))
        out[a] = b
    return out


class LabelMapper(object):
    """Turns a label into a donor: the map, then the pattern's first group, then the label itself."""

    def __init__(self, mapping=None, pattern=None):
        self.mapping = dict(mapping or {})
        self.pattern = re.compile(pattern) if pattern else None

    def __call__(self, label):
        if label in self.mapping:
            return self.mapping[label]
        if self.pattern is not None:
            m = self.pattern.match(label)
            if m:
                return m.group(1)
        return label

    def many(self, labels):
        return [self(x) for x in labels]


def make_mapper(label_map=None, label_pattern=None):
    """A LabelMapper from a label-map path (or dict) and a pattern; None and "" mean none."""
    mapping = label_map
    if isinstance(label_map, str):
        mapping = read_label_map(label_map) if label_map else None
    return LabelMapper(mapping, label_pattern or None)


def _status(raw):
    """Normalised status of each raw word (list in, list out)."""
    return [_STATUS_WORDS.get(w.strip().lower(), w) for w in raw]


def _frame(barcode, status_raw, donor, pair, confidence, cluster, source):
    n = len(barcode)
    blank = [""] * n
    df = pd.DataFrame({
        "barcode": list(barcode),
        "status": _status(status_raw),
        "status_raw": list(status_raw),
        "donor": list(donor) if donor is not None else blank,
        "pair": list(pair) if pair is not None else blank,
        "confidence": list(confidence) if confidence is not None else blank,
        "cluster": list(cluster) if cluster is not None else blank,
        "source": [source] * n,
    }, columns=SCHEMA)
    keep_donor = df.status.isin(["singlet", "recovered"])
    df.loc[~keep_donor, "donor"] = ""
    dup = df.barcode.duplicated()
    if dup.any():
        raise ValueError("barcode %s is listed twice" % df.barcode[dup].iloc[0])
    return df


# ---------------------------------------------------------------------------
# One reader per kind


def _pairs_from_names(names):
    out = []
    for w in names:
        parts = w.split(",")
        out.append(pair_key(parts[0], parts[1]) if len(parts) == 2 and all(parts) else "")
    return out


def read_vireo(path, prob_doublet=None, mapper=None):
    """vireo's donor_ids.tsv. prob_doublet: its prob_doublet.tsv.gz, for the true pair."""
    from demux_tools import vireo_files

    t = vireo_files.read_donor_ids(path)
    mapper = mapper or LabelMapper()
    raw = ["doublet" if d == "doublet" else "unassigned" if d == "unassigned" else "singlet" for d in t.donor_id]
    donor = mapper.many(t.donor_id)
    if prob_doublet:
        pairs = vireo_files.true_pairs(prob_doublet, t.best_doublet)
    else:
        pairs = _pairs_from_names(t.best_doublet)
    if mapper.mapping or mapper.pattern is not None:
        pairs = [pair_key(*mapper.many(p.split(","))) if p else "" for p in pairs]
    return _frame(t.cell, raw, donor, pairs, t.prob_max, None, "vireo")


def read_souporcell(path, cluster_map=None):
    """souporcell's clusters.tsv. cluster_map: {cluster: donor}; see the module docstring."""
    t = read_any_table(path)
    for col in ("barcode", "status", "assignment"):
        if col not in t.columns:
            raise ValueError("%s lacks column %s" % (path, col))
    df = _frame(t.barcode, t.status, None, None, None, t.assignment, "souporcell")
    return name_clusters(df, cluster_map)


def name_clusters(calls, cluster_map=None):
    """Set donor (singlets) and pair (doublets "i/j") of souporcell-like calls from their cluster.

    cluster_map: {cluster: donor}. Without a map the cluster numbers stand in
    for donors. A doublet with a cluster the map does not name gets pair "".
    Returns a new DataFrame.
    """
    df = calls.copy()
    donor, pair = [], []
    for st, cl in zip(df.status, df.cluster):
        d, p = "", ""
        if st in ("singlet", "recovered") and cl and "/" not in cl:
            d = cluster_map.get(cl, "") if cluster_map is not None else cl
        elif st == "doublet" and cl.count("/") == 1:
            i, j = cl.split("/")
            if cluster_map is not None:
                di, dj = cluster_map.get(i, ""), cluster_map.get(j, "")
            else:
                di, dj = i, j
            p = pair_key(di, dj) if di and dj else ""
        donor.append(d)
        pair.append(p)
    df["donor"] = donor
    df["pair"] = pair
    return df


def read_hier(path):
    """hier-v2's hierarchical_consensus.tsv (barcode, cluster, donor, ..., status)."""
    t = read_any_table(path)
    for col in ("barcode", "donor", "status"):
        if col not in t.columns:
            raise ValueError("%s lacks column %s" % (path, col))
    conf = t.label_confidence if "label_confidence" in t.columns else None
    cluster = t.cluster if "cluster" in t.columns else None
    return _frame(t.barcode, t.status, t.donor, None, conf, cluster, "hier")


def _sba(path):
    t = read_any_table(path)
    if "Barcode" not in t.columns:
        raise ValueError("%s lacks column Barcode" % path)
    return t


def read_assignment(path, mapper=None, sba=None):
    """cell_assignment.tsv; sba: sample_barcode_assignments.tsv to split a non-standard status."""
    t = read_any_table(path)
    for col in ("barcode", "assignment", "status"):
        if col not in t.columns:
            raise ValueError("%s lacks column %s" % (path, col))
    mapper = mapper or LabelMapper()
    raw = list(t.status)
    if sba:
        s = _sba(sba)
        col = "Consensus_DropletType" if "Consensus_DropletType" in s.columns else "AtLeastHalfSinglet_DropletType"
        if col not in s.columns:
            raise ValueError("%s has neither Consensus_DropletType nor AtLeastHalfSinglet_DropletType" % sba)
        kind = dict(zip(s.Barcode, s[col]))
        raw = [r if r.strip().lower() in _STATUS_WORDS
               else ("doublet" if kind.get(b, "") == "doublet" else "unassigned" if b in kind else r)
               for b, r in zip(t.barcode, raw)]
    df = _frame(t.barcode, raw, mapper.many(t.assignment), None, None, None, "assignment")
    df["status_raw"] = list(t.status)
    return df


def read_demuxafy(path, tool, mapper=None):
    """One tool's calls from sample_barcode_assignments.tsv (TOOL_DropletType, TOOL_Individual_Assignment)."""
    t = _sba(path)
    tcol, acol = tool + "_DropletType", tool + "_Individual_Assignment"
    if tcol not in t.columns or acol not in t.columns:
        tools = sorted(c[:-len("_DropletType")] for c in t.columns if c.endswith("_DropletType"))
        raise ValueError("%s has no columns for %s (tools: %s)" % (path, tool, ", ".join(tools)))
    mapper = mapper or LabelMapper()
    cluster = t[tool + "_Cluster"] if tool + "_Cluster" in t.columns else None
    return _frame(t.Barcode, t[tcol], mapper.many(t[acol]), None, None, cluster, "demuxafy:" + tool)


def read_scorer(path, mapper=None):
    """The Scorer's per-cell table: vireo's 8 columns plus best_pair and margin."""
    t = read_any_table(path)
    bc = "cell" if "cell" in t.columns else "barcode"
    for col in (bc, "donor_id", "best_pair"):
        if col not in t.columns:
            raise ValueError("%s lacks column %s" % (path, col))
    mapper = mapper or LabelMapper()
    raw = ["doublet" if d == "doublet" else "unassigned" if d == "unassigned" else "singlet" for d in t.donor_id]
    pairs = [pair_key(*mapper.many(p.split(","))) if p.count(",") == 1 else "" for p in t.best_pair]
    conf = t.margin if "margin" in t.columns else None
    return _frame(t[bc], raw, mapper.many(t.donor_id), pairs, conf, None, "scorer")


def parse_generic(kind):
    """{"barcode": col, "status": col, ...} from "generic:barcode=..,status=..[,donor=..]"."""
    spec = {}
    for part in kind[len("generic:"):].split(","):
        if "=" not in part:
            raise UsageError("generic kind needs key=column pairs, not %r" % part)
        k, v = part.split("=", 1)
        if k not in ("barcode", "status", "donor", "pair", "confidence", "cluster"):
            raise UsageError("generic kind: unknown key %r" % k)
        spec[k] = v
    for k in ("barcode", "status"):
        if k not in spec:
            raise UsageError("generic kind needs %s=COLUMN" % k)
    return spec


def read_generic(path, spec, mapper=None):
    t = read_any_table(path)
    missing = [c for c in spec.values() if c not in t.columns]
    if missing:
        raise ValueError("%s lacks column(s) %s" % (path, ", ".join(missing)))
    mapper = mapper or LabelMapper()
    get = lambda k: t[spec[k]] if k in spec else None  # noqa: E731
    donor = mapper.many(t[spec["donor"]]) if "donor" in spec else None
    return _frame(t[spec["barcode"]], t[spec["status"]], donor, get("pair"), get("confidence"),
                  get("cluster"), "generic")


def detect_kind(path):
    """The kind of a calls file from its header; raises ValueError when it cannot tell."""
    with open_text(path) as f:
        first = f.readline().rstrip("\r\n")
    cols = set(first.split("\t" if "\t" in first else ","))
    if {"donor_id", "best_doublet", "best_pair"} <= cols:
        return "scorer"
    if {"donor_id", "best_doublet"} <= cols:
        return "vireo"
    if {"barcode", "status", "assignment"} <= cols and ("log_prob_singleton" in cols or "cluster0" in cols):
        return "souporcell"
    if {"barcode", "donor", "status"} <= cols and ("label_confidence" in cols or "recovery_confidence" in cols):
        return "hier"
    if cols == {"barcode", "assignment", "status"}:
        return "assignment"
    if "Barcode" in cols and any(c.endswith("_DropletType") for c in cols):
        raise ValueError("%s is Demuxafy's sample_barcode_assignments.tsv: give kind demuxafy:TOOL" % path)
    raise ValueError("%s: cannot tell the kind of calls from its header" % path)


def read_calls(path, kind=None, label_map=None, label_pattern=None, prob_doublet=None,
               cluster_map=None, sba=None):
    """Read one calls file into the normalised schema (see the module docstring).

    kind: one of KINDS, or None to detect it from the header.
    label_map: path or dict {label: donor}; label_pattern: regular expression
    whose first group is the donor. prob_doublet: vireo's prob_doublet.tsv.gz
    (vireo kind; found beside donor_ids.tsv when not given). cluster_map:
    {cluster: donor} for souporcell. sba: sample_barcode_assignments.tsv for
    the assignment kind. Returns a DataFrame of text in SCHEMA order, rows in
    file order.
    """
    path = os.fspath(path)
    kind = kind or detect_kind(path)
    mapper = make_mapper(label_map, label_pattern)
    if kind == "vireo":
        if prob_doublet is None:
            beside = os.path.join(os.path.dirname(path), "prob_doublet.tsv.gz")
            prob_doublet = beside if os.path.exists(beside) else None
        return read_vireo(path, prob_doublet or None, mapper)
    if kind == "souporcell":
        return read_souporcell(path, cluster_map)
    if kind == "hier":
        df = read_hier(path)
        if mapper.mapping or mapper.pattern is not None:
            df["donor"] = [mapper(d) if d else "" for d in df.donor]
        return df
    if kind == "assignment":
        return read_assignment(path, mapper, sba or None)
    if kind.startswith("demuxafy:"):
        return read_demuxafy(path, kind.split(":", 1)[1], mapper)
    if kind == "scorer":
        return read_scorer(path, mapper)
    if kind.startswith("generic:"):
        return read_generic(path, parse_generic(kind), mapper)
    raise UsageError("unknown kind %r (kinds: %s)" % (kind, ", ".join(KINDS)))


def files_read(path, kind=None, prob_doublet=None, label_map=None, cluster_map=None, sba=None):
    """Every file read_calls reads for these options, for provenance and the output check.

    Includes vireo's prob_doublet.tsv.gz found beside donor_ids.tsv and the
    table or sample_barcode_assignments.tsv behind a cluster map.
    """
    out = [path]
    try:
        kind = kind or detect_kind(path)
    except (ValueError, OSError):
        kind = kind or ""
    if kind == "vireo" and not prob_doublet:
        beside = os.path.join(os.path.dirname(os.fspath(path)), "prob_doublet.tsv.gz")
        prob_doublet = beside if os.path.exists(beside) else None
    out += [prob_doublet, label_map if isinstance(label_map, str) else None, sba]
    if cluster_map and isinstance(cluster_map, str) and not cluster_map.startswith("majority"):
        out.append(cluster_map[4:] if cluster_map.startswith("sba:") else cluster_map)
    return [os.fspath(p) for p in out if p]


# ---------------------------------------------------------------------------
# Naming souporcell clusters


def majority_map(clusters, donors):
    """{cluster: the donor most of its cells carry}, over cells where both are non-empty.

    Ties go to the donor pandas' value_counts lists first.
    """
    d = pd.DataFrame({"cluster": list(clusters), "donor": list(donors)})
    d = d[(d.cluster != "") & (d.donor != "")]
    if d.empty:
        return {}
    return d.groupby("cluster").donor.agg(lambda x: x.value_counts().index[0]).to_dict()


def majority_map_from_calls(clusters_calls, reference_calls):
    """Name souporcell clusters by the donor most of their singlets carry in a reference run.

    Only cells that both call a singlet are used: e.g. the clusters of hier-v2's
    first souporcell run named by its consensus donors.
    """
    m = clusters_calls[["barcode", "status", "cluster"]].merge(
        reference_calls[["barcode", "status", "donor"]], on="barcode", suffixes=("_c", "_r"))
    sing = (m.status_c == "singlet") & (m.status_r == "singlet")
    return majority_map(m.cluster[sing], m.donor[sing])


def cluster_map_from_sba(path, mapper=None, tool="Souporcell"):
    """{cluster: donor} from Demuxafy's own naming in sample_barcode_assignments.tsv.

    Uses TOOL_Cluster and TOOL_Individual_Assignment of the cells whose
    TOOL_DropletType is singlet; a cluster seen with two labels keeps the last.
    """
    t = _sba(path)
    mapper = mapper or LabelMapper()
    sing = t[tool + "_DropletType"] == "singlet"
    return dict(zip(t[tool + "_Cluster"][sing], mapper.many(t[tool + "_Individual_Assignment"][sing])))


def read_cluster_map(path):
    """{cluster: donor} from a table with columns cluster and donor."""
    t = read_any_table(path)
    if "cluster" not in t.columns or "donor" not in t.columns:
        raise ValueError("%s: a cluster map needs columns cluster and donor" % path)
    return dict(zip(t.cluster, t.donor))


def resolve_cluster_map(spec, mapper=None):
    """A cluster map from a spec: "" or None (none), "sba:PATH" or a cluster,donor table.

    "majority" specs need a reference run and are resolved by the caller.
    """
    if not spec:
        return None
    if spec.startswith("sba:"):
        return cluster_map_from_sba(spec[4:], mapper)
    return read_cluster_map(spec)


# ---------------------------------------------------------------------------
# Command: calls


def configure_calls(parser):
    parser.description = ("Read one calls file (vireo, souporcell, hier-v2 consensus, cell_assignment.tsv, "
                          "Demuxafy, the Scorer or any table) and write it in the normalised calls schema: "
                          "barcode, status, status_raw, donor, pair, confidence, cluster, source.")
    parser.add_argument("path", help="the calls file")
    parser.add_argument("--kind", default=None,
                        help="vireo, souporcell, hier, assignment, demuxafy:TOOL, scorer or "
                             "generic:barcode=COL,status=COL[,donor=COL,...] (default: from the header)")
    parser.add_argument("--prob-doublet", default=None, help="vireo's prob_doublet.tsv.gz (default: beside the file)")
    parser.add_argument("--label-map", default=None, help="table label,donor (or sample_id,donor_id)")
    parser.add_argument("--label-pattern", default=None, help="regular expression; group 1 is the donor")
    parser.add_argument("--cluster-map", default=None, help="souporcell: table cluster,donor, or sba:PATH")
    parser.add_argument("--sba", default=None, help="assignment kind: sample_barcode_assignments.tsv")
    parser.add_argument("--out", required=True, help="output table (.gz for gzip)")


def run_calls(args):
    clock = start_clock()
    inputs = files_read(args.path, args.kind, args.prob_doublet, args.label_map, args.cluster_map, args.sba)
    check_output_path(args.out, inputs)
    try:
        mapper = make_mapper(args.label_map, args.label_pattern)
        cmap = resolve_cluster_map(args.cluster_map, mapper)
        calls = read_calls(args.path, args.kind, args.label_map, args.label_pattern, args.prob_doublet,
                           cmap, args.sba)
    except (ValueError, OSError, KeyError, re.error) as e:
        raise UsageError(str(e))
    write_table(calls, args.out)
    counts = calls.status.value_counts().to_dict()
    write_provenance(args.out, args.argv, inputs, clock=clock,
                     extra={"cells": len(calls), "status": counts, "source": calls.source.iloc[0] if len(calls) else ""})
    print("calls: %d cells (%s) in %s" % (len(calls), ", ".join("%s %d" % kv for kv in sorted(counts.items())), args.out))
    return 0

