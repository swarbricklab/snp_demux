"""Compare two runs' calls cell by cell, and check donor labels against per-cell sex calls.

Commands: `compare-runs` and `sex-check`. Both read calls with `calls.read_calls`,
so either side can be any kind of calls file.

compare-runs joins two runs on the barcode (all barcodes of either run) and
writes, per capture:

    PREFIX.cells.tsv.gz     one row per barcode, sorted by barcode within a
                            capture: capture, barcode, status_a, status_b,
                            donor_a, donor_b, pair_a, pair_b, in_a, in_b,
                            same_status, same_donor, same_pair, same_call,
                            transition
    PREFIX.transitions.tsv  capture, status_a, status_b, cells, same_donor,
                            different_donor
    PREFIX.summary.tsv      capture, cells_a, cells_b, cells_both, only_a,
                            only_b, same_call, same_call_pct,
                            singlet_donor_changes, doublets_both,
                            doublet_pairs_compared, doublet_pairs_same
    PREFIX.sex_check.tsv    with sex calls: the sex check of each run
    PREFIX.provenance.json

Definitions, for a barcode in both runs:
- same_status: the two statuses are equal.
- same_donor: the two donor fields are equal ("" for a doublet or an
  unassigned cell, so two doublets have the same donor).
- same_pair: "true" or "false" when both runs name a pair, "" when either
  names none (hier, assignment and Demuxafy calls name no pair; vireo names
  one for every cell).
- same_call: same status and same donor, and for two doublets whose pairs
  both runs name, the same pair.
- transition: "status_a -> status_b"; a barcode missing from a run has status
  "absent" there.
In the transitions table, same_donor and different_donor count the cells
where both runs name a donor; a cell with a donor on one side only is in
neither column. In the summary, doublets_both counts cells both runs call
doublet, doublet_pairs_compared those of them with a pair on both sides, and
doublet_pairs_same those with the same pair.

The sex check counts, per run, status and label sex, the cells whose donor's
sex (from the donor table, or male for --male-donors and female otherwise)
disagrees with the cell's own sex call. Sex calls are compared in lower case.
Only cells whose sex call is male or female are counted as with_sex_call.
"""
import os
import re

import numpy as np
import pandas as pd

from demux_tools._io import (UsageError, check_output_path, start_clock, write_provenance,
                             write_table)
from demux_tools.calls import (files_read, majority_map_from_calls, make_mapper, name_clusters,
                               read_any_table, read_calls, resolve_cluster_map)

CELL_COLUMNS = ["capture", "barcode", "status_a", "status_b", "donor_a", "donor_b", "pair_a", "pair_b",
                "in_a", "in_b", "same_status", "same_donor", "same_pair", "same_call", "transition"]
TRANSITION_COLUMNS = ["capture", "status_a", "status_b", "cells", "same_donor", "different_donor"]
SUMMARY_COLUMNS = ["capture", "cells_a", "cells_b", "cells_both", "only_a", "only_b", "same_call",
                   "same_call_pct", "singlet_donor_changes", "doublets_both", "doublet_pairs_compared",
                   "doublet_pairs_same"]
SEX_COLUMNS = ["capture", "run", "status", "label_sex", "cells", "with_sex_call", "discordant",
               "discordant_share"]
ABSENT = "absent"


def _bool_text(values):
    return np.where(np.asarray(values, dtype=bool), "true", "false")


def compare(a, b, capture=""):
    """Compare calls a and b (normalised schema). Returns (cells, transitions, summary) DataFrames.

    The per-cell rows are sorted by barcode, so every pandas version writes
    them in the same order.
    """
    m = a[["barcode", "status", "donor", "pair"]].merge(
        b[["barcode", "status", "donor", "pair"]], on="barcode", how="outer", suffixes=("_a", "_b"),
        indicator=True, sort=False)
    m = m.sort_values("barcode", kind="mergesort").reset_index(drop=True)
    in_a = m._merge.isin(["both", "left_only"]).to_numpy()
    in_b = m._merge.isin(["both", "right_only"]).to_numpy()
    for col in ("status_a", "status_b", "donor_a", "donor_b", "pair_a", "pair_b"):
        m[col] = m[col].fillna("")
    m.loc[~in_a, "status_a"] = ABSENT
    m.loc[~in_b, "status_b"] = ABSENT
    both = in_a & in_b
    same_status = both & (m.status_a == m.status_b).to_numpy()
    same_donor = both & (m.donor_a == m.donor_b).to_numpy()
    paired = ((m.pair_a != "") & (m.pair_b != "")).to_numpy()
    pair_equal = (m.pair_a == m.pair_b).to_numpy()
    dbl_both = both & (m.status_a == "doublet").to_numpy() & (m.status_b == "doublet").to_numpy()
    same_call = same_status & same_donor & (~(dbl_both & paired) | pair_equal)
    cells = pd.DataFrame({
        "capture": capture, "barcode": m.barcode,
        "status_a": m.status_a, "status_b": m.status_b, "donor_a": m.donor_a, "donor_b": m.donor_b,
        "pair_a": m.pair_a, "pair_b": m.pair_b,
        "in_a": _bool_text(in_a), "in_b": _bool_text(in_b),
        "same_status": _bool_text(same_status), "same_donor": _bool_text(same_donor),
        "same_pair": np.where(paired, _bool_text(pair_equal), ""),
        "same_call": _bool_text(same_call),
        "transition": m.status_a + " -> " + m.status_b,
    }, columns=CELL_COLUMNS)
    named = (m.donor_a != "") & (m.donor_b != "")
    t = pd.DataFrame({"status_a": m.status_a, "status_b": m.status_b,
                      "same": named & (m.donor_a == m.donor_b), "diff": named & (m.donor_a != m.donor_b)})
    g = t.groupby(["status_a", "status_b"], sort=True)
    trans = pd.DataFrame({"cells": g.size(), "same_donor": g["same"].sum().astype(int),
                          "different_donor": g["diff"].sum().astype(int)}).reset_index()
    trans.insert(0, "capture", capture)
    n_both = int(both.sum())
    sing = both & (m.status_a == "singlet").to_numpy() & (m.status_b == "singlet").to_numpy()
    summary = pd.DataFrame([{
        "capture": capture, "cells_a": int(in_a.sum()), "cells_b": int(in_b.sum()), "cells_both": n_both,
        "only_a": int((in_a & ~in_b).sum()), "only_b": int((in_b & ~in_a).sum()),
        "same_call": int(same_call.sum()),
        "same_call_pct": round(100.0 * same_call.sum() / n_both, 4) if n_both else float("nan"),
        "singlet_donor_changes": int((sing & (m.donor_a != m.donor_b).to_numpy()).sum()),
        "doublets_both": int(dbl_both.sum()), "doublet_pairs_compared": int((dbl_both & paired).sum()),
        "doublet_pairs_same": int((dbl_both & paired & pair_equal).sum()),
    }], columns=SUMMARY_COLUMNS)
    return cells, trans[TRANSITION_COLUMNS], summary


# ---------------------------------------------------------------------------
# Sex


def read_sex_calls(path, sex_column="sex_call", barcode_column="barcode"):
    """{barcode: sex call in lower case} from a per-cell table (comma- or tab-separated)."""
    t = read_any_table(path)
    for col in (barcode_column, sex_column):
        if col not in t.columns:
            raise ValueError("%s lacks column %s" % (path, col))
    return dict(zip(t[barcode_column], [x.strip().lower() for x in t[sex_column]]))


def read_donor_sex(path=None, male_donors=None, default_sex=None):
    """A function donor -> "male" | "female" | "unknown".

    path: a donor table with columns donor and sex. male_donors: donor names
    that are male. default_sex: the sex of a donor in neither ("female" when
    only male_donors is given, else "unknown").
    """
    table = {}
    if path:
        t = read_any_table(path)
        if "donor" not in t.columns or "sex" not in t.columns:
            raise ValueError("%s: a donor table needs columns donor and sex" % path)
        table = {d: s.strip().lower() for d, s in zip(t.donor, t.sex)}
    for d in male_donors or ():
        table[d] = "male"
    if default_sex is None:
        default_sex = "female" if (male_donors and not path) else "unknown"

    def sex_of(donor):
        return table.get(donor, default_sex) if donor else ""
    return sex_of


def sex_check(calls, sex_calls, donor_sex, statuses=("singlet", "recovered"), capture="", run=""):
    """Cells whose donor's sex disagrees with their own sex call.

    calls: normalised calls. sex_calls: {barcode: call}. donor_sex: donor ->
    "male" | "female" | "unknown" (e.g. from read_donor_sex). Returns one row
    per status in `statuses` and label sex (male, female) with: cells,
    with_sex_call (call male or female), discordant (the other sex) and
    discordant_share.
    """
    c = calls[calls.status.isin(statuses) & (calls.donor != "")]
    label = np.array([donor_sex(d) for d in c.donor], dtype=object)
    call = np.array([sex_calls.get(b, "").strip().lower() for b in c.barcode], dtype=object)
    has = np.isin(call, ["male", "female"])
    rows = []
    for st in statuses:
        s = (c.status == st).to_numpy()
        for sex in ("male", "female"):
            k = s & (label == sex)
            n = int((k & has).sum())
            bad = int((k & has & (call != sex)).sum())
            rows.append({"capture": capture, "run": run, "status": st, "label_sex": sex, "cells": int(k.sum()),
                         "with_sex_call": n, "discordant": bad,
                         "discordant_share": bad / float(n) if n else float("nan")})
    return pd.DataFrame(rows, columns=SEX_COLUMNS)


# ---------------------------------------------------------------------------
# Command helpers


def _parse_run(text):
    """("kind", "path") from "KIND:PATH" (KIND may itself hold a colon, e.g. demuxafy:Vireo) or "PATH"."""
    if os.path.exists(text):
        return None, text
    for split in range(len(text)):
        if text[split] == ":" and os.path.exists(text[split + 1:]):
            return text[:split], text[split + 1:]
    raise UsageError("no such calls file in %r (use KIND:PATH or PATH)" % text)


_SIDE_KEYS = ("kind", "path", "prob_doublet", "label_map", "label_pattern", "cluster_map", "sba")


def _side(row, side, defaults):
    """The read options of one side of a manifest row, falling back to the command-line defaults."""
    out = {}
    for k in _SIDE_KEYS:
        v = row.get(side + "_" + k, "")
        out[k] = v if v else defaults.get(k, "")
    return out


def _read_side(opts, other=None):
    """Read one side; a cluster map of "majority" names clusters by the other side's singlet donors."""
    mapper = make_mapper(opts["label_map"] or None, opts["label_pattern"] or None)
    spec = opts["cluster_map"]
    cmap = None if spec in ("", "majority") else resolve_cluster_map(spec, mapper)
    calls = read_calls(opts["path"], opts["kind"] or None, opts["label_map"] or None, opts["label_pattern"] or None,
                       opts["prob_doublet"] or None, cmap, opts["sba"] or None)
    if spec == "majority":
        if other is None:
            raise UsageError("cluster map 'majority' needs the other run to be read first")
        calls = name_clusters(calls, majority_map_from_calls(calls, other))
    return calls


def _inputs(opts):
    """Every file one side's options make read_calls read (calls.files_read)."""
    return files_read(opts["path"], opts.get("kind") or None, opts.get("prob_doublet") or None,
                      opts.get("label_map") or None, opts.get("cluster_map") or None, opts.get("sba") or None)


def _add_side_options(parser, side):
    up = side.upper()
    parser.add_argument("--" + side, default=None, metavar="KIND:PATH",
                        help="run %s: calls file, optionally prefixed by its kind" % up)
    parser.add_argument("--%s-prob-doublet" % side, default="", help="vireo prob_doublet.tsv.gz of run %s" % up)
    parser.add_argument("--%s-label-map" % side, default="", help="label map of run %s" % up)
    parser.add_argument("--%s-label-pattern" % side, default="", help="label pattern of run %s" % up)
    parser.add_argument("--%s-cluster-map" % side, default="",
                        help="souporcell cluster map of run %s: table, sba:PATH or majority" % up)
    parser.add_argument("--%s-sba" % side, default="", help="sample_barcode_assignments.tsv of run %s" % up)


def _defaults(args, side):
    d = {}
    for k in _SIDE_KEYS[2:]:
        d[k] = getattr(args, "%s_%s" % (side, k)) or ""
    return d


def _read_manifest(path, required):
    t = read_any_table(path)
    missing = [c for c in required if c not in t.columns]
    if missing:
        raise UsageError("manifest %s lacks column(s) %s" % (path, ", ".join(missing)))
    return t


def _sex_options(parser):
    parser.add_argument("--sex-calls", default=None, help="per-cell sex calls (columns barcode and sex_call)")
    parser.add_argument("--sex-column", default="sex_call", help="sex call column (default sex_call)")
    parser.add_argument("--male-donors", default=None, help="comma-separated male donors; other donors are female")
    parser.add_argument("--donors", default=None, help="donor table with columns donor and sex")
    parser.add_argument("--statuses", default="singlet,recovered", help="statuses checked (default singlet,recovered)")


def _donor_sex(args):
    males = [d for d in (args.male_donors or "").split(",") if d]
    if not males and not args.donors:
        return None
    try:
        return read_donor_sex(args.donors, males)
    except (ValueError, OSError) as e:
        raise UsageError(str(e))


# ---------------------------------------------------------------------------
# Command: compare-runs


def configure_compare_runs(parser):
    parser.description = ("Compare two runs' calls per cell. Writes PREFIX.cells.tsv.gz, PREFIX.transitions.tsv, "
                          "PREFIX.summary.tsv (and PREFIX.sex_check.tsv with sex calls) and PREFIX.provenance.json. "
                          "One capture from --a and --b, or many from --manifest.")
    _add_side_options(parser, "a")
    _add_side_options(parser, "b")
    parser.add_argument("--a-name", default="a", help="name of run A in the sex check (default a)")
    parser.add_argument("--b-name", default="b", help="name of run B in the sex check (default b)")
    parser.add_argument("--capture", default="", help="capture name for the single-capture form")
    parser.add_argument("--manifest", default=None,
                        help="table with columns capture, a_path, b_path and optionally a_kind, b_kind, "
                             "a_/b_prob_doublet, a_/b_label_map, a_/b_label_pattern, a_/b_cluster_map, "
                             "a_/b_sba and sex_calls")
    _sex_options(parser)
    parser.add_argument("--out", required=True, help="output prefix")


def run_compare_runs(args):
    clock = start_clock()
    rows = []
    if args.manifest:
        man = _read_manifest(args.manifest, ["capture", "a_path", "b_path"])
        rows = man.to_dict("records")
    else:
        if not args.a or not args.b:
            raise UsageError("give --a and --b, or --manifest")
        ka, pa = _parse_run(args.a)
        kb, pb = _parse_run(args.b)
        rows = [{"capture": args.capture, "a_kind": ka or "", "a_path": pa, "b_kind": kb or "", "b_path": pb,
                 "sex_calls": args.sex_calls or ""}]
    donor_sex = _donor_sex(args)
    statuses = tuple(s for s in args.statuses.split(",") if s)
    out = args.out
    targets = [out + s for s in (".cells.tsv.gz", ".transitions.tsv", ".summary.tsv", ".sex_check.tsv")]
    inputs = [args.manifest] if args.manifest else []
    cells, trans, summ, sex = [], [], [], []
    for row in rows:
        oa, ob = _side(row, "a", _defaults(args, "a")), _side(row, "b", _defaults(args, "b"))
        sex_path = row.get("sex_calls", "") or args.sex_calls or ""
        read_here = _inputs(oa) + _inputs(ob) + ([sex_path] if sex_path else [])
        inputs += _inputs(oa) + _inputs(ob)
        for t in targets:
            check_output_path(t, read_here)
        try:
            a = _read_side(oa)
            b = _read_side(ob, other=a)
        except (ValueError, OSError, KeyError, re.error) as e:
            raise UsageError("capture %s: %s" % (row.get("capture", ""), e))
        c, t, s = compare(a, b, row.get("capture", ""))
        cells.append(c)
        trans.append(t)
        summ.append(s)
        if sex_path and donor_sex is not None:
            inputs.append(sex_path)
            try:
                calls_sex = read_sex_calls(sex_path, args.sex_column)
            except (ValueError, OSError) as e:
                raise UsageError("capture %s: %s" % (row.get("capture", ""), e))
            for name, x in ((args.a_name, a), (args.b_name, b)):
                sex.append(sex_check(x, calls_sex, donor_sex, statuses, row.get("capture", ""), name))
    d = os.path.dirname(os.path.abspath(out))
    os.makedirs(d, exist_ok=True)
    C = pd.concat(cells, ignore_index=True)
    T = pd.concat(trans, ignore_index=True)
    S = pd.concat(summ, ignore_index=True)
    write_table(C, targets[0])
    write_table(T, targets[1])
    write_table(S, targets[2], float_formats={"same_call_pct": "%.4f"})
    written = targets[:3]
    totals = {"captures": len(S), "cells_both": int(S.cells_both.sum()), "same_call": int(S.same_call.sum()),
              "only_a": int(S.only_a.sum()), "only_b": int(S.only_b.sum()),
              "singlet_donor_changes": int(S.singlet_donor_changes.sum()),
              "doublets_both": int(S.doublets_both.sum()), "doublet_pairs_compared": int(S.doublet_pairs_compared.sum()),
              "doublet_pairs_same": int(S.doublet_pairs_same.sum())}
    if sex:
        X = pd.concat(sex, ignore_index=True)
        write_table(X, targets[3], float_formats={"discordant_share": "%.6g"})
        written.append(targets[3])
    write_provenance(out, args.argv, inputs, clock=clock, extra={"outputs": written, "totals": totals})
    tt = T.groupby(["status_a", "status_b"]).sum(numeric_only=True).reset_index()
    print("compare-runs: %d captures, %d cells in both runs, %d with the same call; only in A %d, only in B %d"
          % (totals["captures"], totals["cells_both"], totals["same_call"], totals["only_a"], totals["only_b"]))
    print(tt.to_string(index=False))
    if sex:
        print(X.groupby(["run", "status", "label_sex"])[["cells", "with_sex_call", "discordant"]].sum()
              .reset_index().to_string(index=False))
    return 0


# ---------------------------------------------------------------------------
# Command: sex-check


def configure_sex_check(parser):
    parser.description = ("Count cells whose donor's sex disagrees with the cell's own sex call, per run, "
                          "status and label sex. One run from --calls, or many from --manifest.")
    parser.add_argument("--calls", default=None, metavar="KIND:PATH", help="calls file, optionally prefixed by its kind")
    parser.add_argument("--prob-doublet", default="", help="vireo prob_doublet.tsv.gz")
    parser.add_argument("--label-map", default="", help="label map (label,donor or sample_id,donor_id)")
    parser.add_argument("--label-pattern", default="", help="regular expression; group 1 is the donor")
    parser.add_argument("--cluster-map", default="", help="souporcell cluster map: table or sba:PATH")
    parser.add_argument("--sba", default="", help="sample_barcode_assignments.tsv (assignment kind)")
    parser.add_argument("--capture", default="", help="capture name for the single-run form")
    parser.add_argument("--run", default="", help="run name for the single-run form")
    parser.add_argument("--manifest", default=None,
                        help="table with columns capture, run, path, sex_calls and optionally kind, prob_doublet, "
                             "label_map, label_pattern, cluster_map, sba")
    _sex_options(parser)
    parser.add_argument("--out", required=True, help="output table")


def run_sex_check(args):
    clock = start_clock()
    donor_sex = _donor_sex(args)
    if donor_sex is None:
        raise UsageError("give --male-donors or --donors")
    defaults = {"prob_doublet": args.prob_doublet, "label_map": args.label_map, "label_pattern": args.label_pattern,
                "cluster_map": args.cluster_map, "sba": args.sba}
    if args.manifest:
        rows = _read_manifest(args.manifest, ["capture", "run", "path", "sex_calls"]).to_dict("records")
    else:
        if not args.calls or not args.sex_calls:
            raise UsageError("give --calls and --sex-calls, or --manifest")
        k, p = _parse_run(args.calls)
        rows = [{"capture": args.capture, "run": args.run, "kind": k or "", "path": p, "sex_calls": args.sex_calls}]
    statuses = tuple(s for s in args.statuses.split(",") if s)
    inputs = [args.manifest] if args.manifest else []
    out = []
    for row in rows:
        opts = {k: (row.get(k, "") or defaults.get(k, "") or "") for k in _SIDE_KEYS}
        if opts["cluster_map"] == "majority":
            raise UsageError("sex-check cannot name clusters by majority; give a cluster map")
        inputs += _inputs(opts) + [row["sex_calls"]]
        check_output_path(args.out, _inputs(opts) + [row["sex_calls"]])
        try:
            calls = _read_side(opts)
            sc = read_sex_calls(row["sex_calls"], args.sex_column)
        except (ValueError, OSError, KeyError, re.error) as e:
            raise UsageError("capture %s: %s" % (row.get("capture", ""), e))
        out.append(sex_check(calls, sc, donor_sex, statuses, row.get("capture", ""), row.get("run", "")))
    X = pd.concat(out, ignore_index=True)
    write_table(X, args.out, float_formats={"discordant_share": "%.6g"})
    tot = X.groupby(["run", "status", "label_sex"])[["cells", "with_sex_call", "discordant"]].sum().reset_index()
    write_provenance(args.out, args.argv, inputs, clock=clock,
                     extra={"rows": len(X), "totals": tot.to_dict("records")})
    print(tot.to_string(index=False))
    return 0

