"""Read vireo's output files and recover the true doublet pair; the `vireo-pairs` command.

vireo writes, per capture:

    donor_ids.tsv         cell, donor_id, prob_max, prob_doublet, n_vars,
                          best_singlet, best_doublet, doublet_logLikRatio
    prob_singlet.tsv.gz   cell + one column per donor (vireo's donor order)
    prob_doublet.tsv.gz   cell + one column per donor pair

vireo scores doublets over `itertools.combinations` of its donors. vireoSNP
0.5.7 (and earlier) names the pairs with `itertools.permutations` instead, in
both the prob_doublet.tsv.gz header and `best_doublet`: the name written for
the k-th pair is the k-th permutation. Only pairs that include the first donor
come out right. With donors (a, b, c):

    k   pair scored (combinations)   name written (permutations)
    0   a,b                          a,b
    1   a,c                          a,c
    2   b,c                          b,a

So a cell whose best pair is (b, c) is written "b,a". The header then has
n(n-1) names over n(n-1)/2 value columns. vireoSNP 0.5.9 names pairs with
combinations. 0.5.10 renames donor_ids.tsv's `prob_doublet` to
`prob_max_doublet`; the readers here accept either name.

The true pair is recovered from the header: the permutation names give vireo's
donor order, and the written name's position k in the header gives the pair
combinations(order, 2)[k]. When every pair probability of a cell underflows to
0 (`prob_doublet` written 0.00e+00), vireo's argmax is pair 0 whatever the
data, so no pair can be recovered for that cell.

Pairs are returned as the two donor names in string order joined by a comma,
e.g. "a,c", whatever order vireo used.
"""
import itertools
import os

import numpy as np

from demux_tools._io import (UsageError, check_output_path, open_text, read_table, start_clock,
                             write_provenance, write_table)

DONOR_IDS_COLUMNS = ["cell", "donor_id", "prob_max", "prob_doublet", "n_vars",
                     "best_singlet", "best_doublet", "doublet_logLikRatio"]
PAIRS_COLUMNS = DONOR_IDS_COLUMNS + ["best_pair", "pair_prob_underflow", "vireo_naming"]


def pair_key(a, b):
    """The two donor names in string order, joined by a comma."""
    return ",".join(sorted((a, b)))


def read_donor_ids(path):
    """vireo's donor_ids.tsv as a DataFrame of text (every column str, as written).

    vireoSNP 0.5.10's `prob_max_doublet` column is renamed `prob_doublet`.
    Raises ValueError if any of the 8 columns in DONOR_IDS_COLUMNS is missing.
    Extra columns are kept.
    """
    df = read_table(path)
    if "prob_max_doublet" in df.columns and "prob_doublet" not in df.columns:
        df = df.rename(columns={"prob_max_doublet": "prob_doublet"})
    missing = [c for c in DONOR_IDS_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError("%s lacks vireo column(s) %s" % (path, ", ".join(missing)))
    return df


def read_prob_header(path):
    """Column names after "cell" in prob_singlet.tsv.gz or prob_doublet.tsv.gz (gzip or plain)."""
    with open_text(path) as f:
        names = f.readline().rstrip("\r\n").split("\t")
    if not names or names[0] != "cell":
        raise ValueError("%s: first column is %r, not 'cell'" % (path, names[0] if names else ""))
    return names[1:]


def donor_order(prob_singlet_path):
    """Donors in vireo's order: the prob_singlet.tsv.gz header after "cell"."""
    return read_prob_header(prob_singlet_path)


def _n_from_count(m, per):
    """n such that n(n-1)/per == m, or None."""
    n = int(round((1 + (1 + 4 * per * m) ** 0.5) / 2))
    return n if n >= 2 and n * (n - 1) == per * m else None


def _order_from_header(header, n):
    """vireo's donor order read off a pair header: first name's first donor, then
    the second donor of each of the first n - 1 names (true for both namings)."""
    first = header[0].split(",", 1)
    if len(first) != 2:
        raise ValueError("pair name %r has no comma" % header[0])
    order = [first[0]]
    for name in header[:n - 1]:
        parts = name.split(",", 1)
        if len(parts) != 2:
            raise ValueError("pair name %r has no comma" % name)
        order.append(parts[1])
    return order


def pair_naming(prob_doublet_path, donors=None):
    """How a prob_doublet.tsv.gz header names its pairs.

    donors: vireo's donor order (e.g. `donor_order(prob_singlet)`), optional;
    needed only when a donor name contains a comma.

    Returns a dict:
        kind:   "permutations" (vireoSNP <= 0.5.7) or "combinations" (>= 0.5.9)
        donors: vireo's donor order (list of str)
        header: the pair names as written (list of str)
        combos: the true pair of value column k, as pair_key(a, b), in
                combinations order (list of str, n(n-1)/2 entries)
    Raises ValueError when the header is neither naming of one donor order
    (each candidate is rebuilt in full and compared, so this never guesses).
    """
    header = read_prob_header(prob_doublet_path)
    m = len(header)
    for kind, per, gen in (("permutations", 1, itertools.permutations), ("combinations", 2, itertools.combinations)):
        n = _n_from_count(m, per)
        if n is None:
            continue
        if donors is not None:
            if len(donors) != n:
                continue
            order = list(donors)
        else:
            try:
                order = _order_from_header(header, n)
            except ValueError:
                continue
        if [",".join(p) for p in gen(order, 2)] == header:
            return {"kind": kind, "donors": order, "header": header,
                    "combos": [pair_key(a, b) for a, b in itertools.combinations(order, 2)]}
    raise ValueError("%s: %d pair names are neither vireo's permutation nor its combination naming%s"
                     % (prob_doublet_path, m, "" if donors is None else " of the donors given"))


def true_pairs(prob_doublet_path, written, donors=None):
    """The pair vireo scored for each `best_doublet` name it wrote.

    written: iterable of best_doublet values from the same capture's
    donor_ids.tsv. Returns a list of pair_key strings, one per value; "" for a
    value not in the header (e.g. an empty field). With combination naming the
    result is the written pair in string order. Raises ValueError for a name
    whose header position is beyond the value columns (vireo never writes one).
    Cells whose pair probabilities all underflow get pair 0 here, which means
    nothing; `vireo_pairs` blanks them.
    """
    naming = pair_naming(prob_doublet_path, donors=donors)
    combos = naming["combos"]
    index = {name: k for k, name in enumerate(naming["header"])}
    out = []
    for w in written:
        k = index.get(w)
        if k is None:
            out.append("")
        elif k >= len(combos):
            raise ValueError("%s: best_doublet %r is pair name %d of %d, beyond the %d value columns"
                             % (prob_doublet_path, w, k, len(naming["header"]), len(combos)))
        else:
            out.append(combos[k])
    return out


def read_prob_matrix(path):
    """(cells, values) of prob_singlet.tsv.gz or prob_doublet.tsv.gz.

    cells: list of str; values: float64 array (cells x value columns), parsed
    from the written text ("%.2e"). The header is not used for the width, so a
    permutation-named doublet file reads correctly.
    """
    import pandas as pd

    from demux_tools._io import is_gzip

    path = os.fspath(path)
    df = pd.read_csv(path, sep="\t", header=None, skiprows=1, dtype=str, keep_default_na=False,
                     na_filter=False, compression="gzip" if is_gzip(path) else None, engine="c")
    if df.shape[0] == 0:
        return [], np.zeros((0, 0), dtype=np.float64)
    cells = df[0].tolist()
    values = df.iloc[:, 1:].to_numpy(dtype=str).astype(np.float64)
    return cells, values


def best_pairs_from_matrix(prob_doublet_path, donors=None):
    """Each cell's best pair read from the probability values themselves.

    Returns a DataFrame with columns:
        cell
        matrix_pair  pair_key of the column with the largest written value
                     ("" when every value is 0)
        max_prob     that value (float)
        underflow    True when every pair probability is 0 (no pair is known)
        tied         True when another column has the same written value, so
                     the argmax of the rounded text may differ from vireo's
    vireo's own `best_doublet` comes from unrounded values; prefer `true_pairs`,
    and use this as a cross-check.
    """
    import pandas as pd

    naming = pair_naming(prob_doublet_path, donors=donors)
    cells, values = read_prob_matrix(prob_doublet_path)
    n_pairs = len(naming["combos"])
    if values.shape[0] and values.shape[1] != n_pairs:
        raise ValueError("%s: rows have %d values, the header implies %d pairs"
                         % (prob_doublet_path, values.shape[1], n_pairs))
    if values.shape[0] == 0:
        return pd.DataFrame({"cell": [], "matrix_pair": [], "max_prob": [], "underflow": [], "tied": []})
    k = values.argmax(axis=1)
    top = values[np.arange(len(k)), k]
    underflow = top == 0
    tied = (values == top[:, None]).sum(axis=1) > 1
    combos = np.asarray(naming["combos"], dtype=object)
    pair = np.where(underflow, "", combos[k])
    return pd.DataFrame({"cell": cells, "matrix_pair": pair, "max_prob": top,
                         "underflow": underflow, "tied": tied & ~underflow})


def vireo_pairs(donor_ids_path, prob_doublet_path, prob_singlet_path=None):
    """donor_ids.tsv with the true pair added, plus a summary dict.

    Returns (table, summary). table: vireo's 8 columns as text (as written;
    prob_max_doublet renamed prob_doublet) plus
        best_pair            true pair (pair_key) of best_doublet; "" when
                             every pair probability underflows
        pair_prob_underflow  "true" when prob_doublet is 0 (every pair
                             probability underflowed), else "false"
        vireo_naming         "permutations" or "combinations"
    summary: cells, underflow, naming, donors, renamed (best_pair differs from
    best_doublet as a set), and the matrix cross-check: matrix_agree,
    matrix_tied, matrix_disagree (cells where the argmax of the written values
    is a clear winner but differs from best_pair; must be 0).
    Raises ValueError if the two files list different cells.
    """
    donors = donor_order(prob_singlet_path) if prob_singlet_path else None
    ids = read_donor_ids(donor_ids_path)
    naming = pair_naming(prob_doublet_path, donors=donors)
    pairs = true_pairs(prob_doublet_path, ids["best_doublet"], donors=naming["donors"])
    underflow = ids["prob_doublet"].astype(np.float64).to_numpy() == 0
    best_pair = ["" if u else p for p, u in zip(pairs, underflow)]

    matrix = best_pairs_from_matrix(prob_doublet_path, donors=naming["donors"])
    if matrix["cell"].tolist() != ids["cell"].tolist():
        if sorted(matrix["cell"]) != sorted(ids["cell"]):
            raise ValueError("%s and %s list different cells" % (donor_ids_path, prob_doublet_path))
        matrix = matrix.set_index("cell").loc[ids["cell"]].reset_index()
    if not np.array_equal(matrix["underflow"].to_numpy(dtype=bool), underflow):
        raise ValueError("prob_doublet in %s disagrees with the largest value in %s"
                         % (donor_ids_path, prob_doublet_path))
    clear = ~underflow & ~matrix["tied"].to_numpy(dtype=bool)
    agree = matrix["matrix_pair"].to_numpy(dtype=object) == np.asarray(best_pair, dtype=object)

    table = ids[DONOR_IDS_COLUMNS].copy()
    table["best_pair"] = best_pair
    table["pair_prob_underflow"] = np.where(underflow, "true", "false")
    table["vireo_naming"] = naming["kind"]
    written_sorted = [pair_key(*w.split(",", 1)) if "," in w else w for w in ids["best_doublet"]]
    summary = {
        "cells": int(len(table)),
        "underflow": int(underflow.sum()),
        "naming": naming["kind"],
        "donors": naming["donors"],
        "renamed": int(sum(1 for b, w, u in zip(best_pair, written_sorted, underflow) if not u and b != w)),
        "matrix_agree": int((clear & agree).sum()),
        "matrix_tied": int((~underflow & matrix["tied"].to_numpy(dtype=bool)).sum()),
        "matrix_disagree": int((clear & ~agree).sum()),
    }
    return table, summary


# ---------------------------------------------------------------------------
# Command: vireo-pairs


def configure_vireo_pairs(parser):
    parser.description = (
        "Add the true doublet pair to vireo's donor_ids.tsv. vireoSNP 0.5.7 writes the "
        "wrong pair name in best_doublet for most cells; best_pair is the pair vireo scored. "
        "Writes vireo's 8 columns plus best_pair, pair_prob_underflow and vireo_naming, and "
        "<out>.provenance.json with counts and a cross-check against the probability matrix.")
    src = parser.add_argument_group("input (a vireo folder, or the files)")
    src.add_argument("--vireo-dir", help="vireo output folder (donor_ids.tsv, prob_doublet.tsv.gz, prob_singlet.tsv.gz)")
    src.add_argument("--donor-ids", help="donor_ids.tsv")
    src.add_argument("--prob-doublet", help="prob_doublet.tsv.gz")
    src.add_argument("--prob-singlet", help="prob_singlet.tsv.gz (optional; gives the donor order)")
    parser.add_argument("--out", required=True, help="output table (.tsv or .tsv.gz)")


def run_vireo_pairs(args):
    clock = start_clock()
    donor_ids, prob_doublet, prob_singlet = args.donor_ids, args.prob_doublet, args.prob_singlet
    if args.vireo_dir:
        donor_ids = donor_ids or os.path.join(args.vireo_dir, "donor_ids.tsv")
        prob_doublet = prob_doublet or os.path.join(args.vireo_dir, "prob_doublet.tsv.gz")
        if prob_singlet is None and os.path.exists(os.path.join(args.vireo_dir, "prob_singlet.tsv.gz")):
            prob_singlet = os.path.join(args.vireo_dir, "prob_singlet.tsv.gz")
    if not donor_ids or not prob_doublet:
        raise UsageError("give --vireo-dir, or --donor-ids and --prob-doublet")
    inputs = [p for p in (donor_ids, prob_doublet, prob_singlet) if p]
    for p in inputs:
        if not os.path.isfile(p):
            raise UsageError("input not found: %s" % p)
    check_output_path(args.out, inputs)
    table, summary = vireo_pairs(donor_ids, prob_doublet, prob_singlet)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    if summary["matrix_disagree"]:
        raise ValueError("best_pair disagrees with the probability matrix for %d cells; the naming "
                         "rule does not hold for these files" % summary["matrix_disagree"])
    write_table(table, args.out)
    write_provenance(args.out, args.argv, inputs, clock=clock, extra=summary)
    print("vireo-pairs: %(cells)d cells, %(naming)s naming, %(renamed)d best_doublet names corrected, "
          "%(underflow)d with every pair probability 0, matrix check %(matrix_agree)d agree / "
          "%(matrix_tied)d tied / %(matrix_disagree)d disagree" % summary)
    return 0
