"""Truth tables for scoring calls: per-cell truth, donor attributes, capture attributes, cell classes.

Nothing here names a cohort. A study passes its own tables:

- Per-cell truth: a manifest with columns capture, path and optionally
  columns (comma-separated; default every column but the barcode). Each path
  is a per-cell table (comma- or tab-separated) with a `barcode` column, such
  as a sex-call table (n_umi, chrY_counts, xist_counts, sex_call) or a cell
  type table (coarse_cell_type). Several rows per capture are joined by
  barcode.
- Donors: columns donor, sex (male, female or empty) and optionally
  immune_only and decoy (true or false). immune_only marks a donor whose
  cells are all immune, such as a blood control; decoy marks a donor that is
  not in the capture.
- Captures: columns capture and optionally pool, capture_set, donors
  (comma-separated: the donors the sample sheet puts in the capture),
  sister_key and genotyped_<run> (comma-separated donors that run had a
  genotype for). An empty pool becomes "donors:" plus the sorted donors, so
  captures without a pool that share their donors share a pool. An empty
  capture_set is "main". An empty sister_key is the pool plus the sorted
  donors.
- Cell classes: `data/cell_classes.tsv` (columns label, class) maps generic
  coarse cell-type labels to immune or non-immune. A label it does not hold is
  "unclear", an empty label "no annotation". Another table can replace it.
"""
import os

import numpy as np
import pandas as pd

from demux_tools._io import read_table
from demux_tools.calls import read_any_table

DEFAULT_CELL_CLASSES = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "cell_classes.tsv")
TRUE_WORDS = ("true", "yes", "1", "y", "t")


def flag(values):
    """Boolean array from text flags (true, yes, 1; anything else is false)."""
    return np.array([str(v).strip().lower() in TRUE_WORDS for v in values], dtype=bool)


def split_list(text):
    """Comma-separated text to a list of non-empty, stripped items."""
    return [x.strip() for x in str(text).split(",") if x.strip()]


def read_cell_classes(path=None):
    """{label: "immune" | "non-immune"} from a table with columns label and class."""
    t = read_table(path or DEFAULT_CELL_CLASSES)
    if "label" not in t.columns or "class" not in t.columns:
        raise ValueError("%s: a cell class table needs columns label and class" % path)
    bad = sorted(set(t["class"]) - {"immune", "non-immune"})
    if bad:
        raise ValueError("%s: classes must be immune or non-immune, not %s" % (path, ", ".join(bad)))
    return dict(zip(t["label"], t["class"]))


def cell_class(labels, classes):
    """Per label: its class, "no annotation" for "", "unclear" for a label not in `classes`."""
    return np.array([classes.get(x, "unclear") if x != "" else "no annotation" for x in labels], dtype=object)


def read_donors(path):
    """Donor table as a DataFrame: donor, sex (lower case), immune_only and decoy (bool)."""
    t = read_any_table(path)
    if "donor" not in t.columns:
        raise ValueError("%s: a donor table needs a donor column" % path)
    if t.donor.duplicated().any():
        raise ValueError("%s: donor %s is listed twice" % (path, t.donor[t.donor.duplicated()].iloc[0]))
    out = pd.DataFrame({"donor": t.donor})
    out["sex"] = [s.strip().lower() for s in t["sex"]] if "sex" in t.columns else ""
    out["immune_only"] = flag(t["immune_only"]) if "immune_only" in t.columns else False
    out["decoy"] = flag(t["decoy"]) if "decoy" in t.columns else False
    return out


def donor_key(donors):
    return "donors:" + ",".join(sorted(donors))


def read_captures(path):
    """Capture table as a DataFrame with capture, pool, capture_set, donors (list), donors_given, sister_key, genotyped_*."""
    t = read_any_table(path)
    if "capture" not in t.columns:
        raise ValueError("%s: a capture table needs a capture column" % path)
    if t.capture.duplicated().any():
        raise ValueError("%s: capture %s is listed twice" % (path, t.capture[t.capture.duplicated()].iloc[0]))
    out = pd.DataFrame({"capture": t.capture})
    donors = [split_list(x) for x in t["donors"]] if "donors" in t.columns else [[] for _ in range(len(t))]
    out["donors"] = donors
    out["donors_given"] = [bool(d) for d in donors]
    pools = t["pool"] if "pool" in t.columns else [""] * len(t)
    out["pool"] = [p if p else donor_key(d) for p, d in zip(pools, donors)]
    sets = t["capture_set"] if "capture_set" in t.columns else [""] * len(t)
    out["capture_set"] = [s if s else "main" for s in sets]
    sisters = t["sister_key"] if "sister_key" in t.columns else [""] * len(t)
    out["sister_key"] = [s if s else p + "|" + ",".join(sorted(d)) for s, p, d in zip(sisters, out.pool, donors)]
    for col in t.columns:
        if col.startswith("genotyped_"):
            out[col] = [split_list(x) for x in t[col]]
    return out


def read_truth_manifest(path):
    """Rows of a per-cell truth manifest: capture, path, columns (list or None)."""
    t = read_any_table(path)
    for col in ("capture", "path"):
        if col not in t.columns:
            raise ValueError("%s: a truth manifest needs columns capture and path" % path)
    cols = t["columns"] if "columns" in t.columns else [""] * len(t)
    return [{"capture": c, "path": p, "columns": split_list(k) or None} for c, p, k in zip(t.capture, t.path, cols)]


def load_cell_truth(rows, capture):
    """Per-cell truth of one capture: one DataFrame keyed by barcode (text), or None when it has none.

    rows: from read_truth_manifest. Files of one capture are joined by barcode
    (outer); a column given by two files raises ValueError.
    """
    out = None
    for r in rows:
        if r["capture"] != capture or not r["path"]:
            continue
        if not os.path.exists(r["path"]):
            raise FileNotFoundError(r["path"])
        t = read_any_table(r["path"])
        if "barcode" not in t.columns:
            raise ValueError("%s lacks a barcode column" % r["path"])
        cols = r["columns"] or [c for c in t.columns if c != "barcode"]
        missing = [c for c in cols if c not in t.columns]
        if missing:
            raise ValueError("%s lacks column(s) %s" % (r["path"], ", ".join(missing)))
        t = t[["barcode"] + cols].drop_duplicates("barcode")
        if out is None:
            out = t
        else:
            clash = [c for c in cols if c in out.columns]
            if clash:
                raise ValueError("capture %s: column(s) %s come from two truth files" % (capture, ", ".join(clash)))
            out = out.merge(t, on="barcode", how="outer")
    return out
