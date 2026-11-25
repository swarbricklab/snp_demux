"""The Scorecard: per-cell tables of several runs and their scores against named truth.

Commands:

    cell-table      one row per cell of every capture, with each run's status,
                    donor and pair beside the per-cell truth
    scorecard       the shares of singlets kept, wrong-donor calls and doublet
                    errors per run and capture set, each against a named truth,
                    with intervals from resampling pools; the composition of
                    each tier by cell class and UMI depth; per-donor recall
    calibrate       the wrong-donor rate by sex against a per-cell score
                    (such as the Scorer's margin), in bins, with Wilson intervals
    doublet-excess  per capture: doublet calls by the tool that made them,
                    against the loading expectation, with ambient estimates

A "run" is one set of calls (calls.read_calls) named in the cell table:
`<run>_status`, `<run>_donor` and `<run>_pair`. A run can name tool runs whose
doublet calls split its doublet tier into groups ("both tools", "vireo only",
...) and give its doublet pairs: `--tools RUN=LABEL:STATUS_RUN[/PAIR_RUN],...`.
The pair of a cell is the pair of the first tool that called it a doublet, or
the first tool's pair when none did.

Statistics. A share is numerator over denominator, both counts of cells. Its
interval is the 2.5th and 97.5th percentile of the share over resamples of
pools: n_boot draws, with replacement, of as many pools as the capture set
has, from numpy's default_rng(seed) started afresh for each capture set. A
resample whose pools hold none of the denominator's cells is left out. No
interval is given when the denominator's cells sit in fewer than min_pools
pools or fewer than 95 per cent of the resamples are left.
"""
import os
import sys
import warnings

import numpy as np
import pandas as pd

from demux_tools import doublets as D
from demux_tools._io import (UsageError, check_output_path, read_table, start_clock, write_provenance,
                             write_table)
from demux_tools.calls import (files_read, majority_map_from_calls, make_mapper, name_clusters, read_any_table,
                               read_calls, read_label_map, resolve_cluster_map)
from demux_tools.truth import (cell_class, load_cell_truth, read_captures, read_cell_classes, read_donors,
                               read_truth_manifest, split_list)

RUN_KEYS = ("kind", "path", "prob_doublet", "label_map", "label_pattern", "cluster_map", "sba")
TIERS = ("singlet", "recovered", "doublet", "unassigned")
NOT_IN_RUN = "not in run"
SCORE_COLUMNS = ["capture_set", "measure", "run", "truth", "value", "low", "high", "numerator", "denominator",
                 "observed_no_chrY", "expected_no_chrY", "calls_in_bins_without_data"]
# Measures that need to know which captures hold a male donor (the capture table's donors).
NEEDS_MALE_PRESENCE = (
    "Wrong donor by sex: male share of singlets labelled a female donor, captures with a male donor",
    "Floor: male share of singlets, captures with no male donor",
    "Ambiguous-sex cells called doublet, captures with a male donor",
    "Ambiguous-sex cells called singlet, captures with a male donor",
    "Ambiguous-sex share of cells with a sex call value, captures with a male donor",
    "Floor: ambiguous-sex share of cells with a sex call value, captures with no male donor",
)


# ---------------------------------------------------------------------------
# Cell table


def read_runs_manifest(path):
    """Rows of a runs manifest: capture, run, path and optionally the read options of calls.read_calls."""
    t = read_any_table(path)
    for col in ("capture", "run", "path"):
        if col not in t.columns:
            raise UsageError("runs manifest %s lacks column %s" % (path, col))
    rows = []
    for r in t.to_dict("records"):
        rows.append(dict({k: r.get(k, "") or "" for k in RUN_KEYS}, capture=r["capture"], run=r["run"]))
    return rows


def _label_maps(rows):
    cache = {}
    for r in rows:
        p = r["label_map"]
        if p and p not in cache:
            cache[p] = read_label_map(p)
    return cache


def read_capture_runs(rows, label_maps=None):
    """{run: calls} of one capture's manifest rows.

    A cluster map "majority:RUN" names souporcell clusters by RUN's singlet
    donors, counting cells both call singlet; RUN must be in the same capture.
    """
    label_maps = label_maps or {}
    out = {}
    later = []
    for r in rows:
        if r["cluster_map"].startswith("majority:"):
            later.append(r)
            continue
        out[r["run"]] = _read_row(r, label_maps)
    for r in later:
        ref = r["cluster_map"].split(":", 1)[1]
        if ref not in out:
            raise UsageError("capture %s: run %s names clusters by run %s, which it does not have"
                             % (r["capture"], r["run"], ref))
        calls = _read_row(dict(r, cluster_map=""), label_maps)
        out[r["run"]] = name_clusters(calls, majority_map_from_calls(calls, out[ref]))
    return out


def _read_row(r, label_maps):
    mapping = label_maps.get(r["label_map"]) if r["label_map"] else None
    mapper = make_mapper(mapping, r["label_pattern"] or None)
    cmap = resolve_cluster_map(r["cluster_map"], mapper)
    return read_calls(r["path"], r["kind"] or None, mapping, r["label_pattern"] or None, r["prob_doublet"] or None,
                      cmap, r["sba"] or None)


def build_cell_table(run_rows, truth_rows=None, captures=None, universe=None, progress=None):
    """One row per cell: capture, pool, capture_set, barcode, truth columns, then each run's status, donor and pair.

    run_rows: read_runs_manifest rows. truth_rows: read_truth_manifest rows.
    captures: read_captures DataFrame (every capture of run_rows must be in
    it); without one a capture is its own pool in set "main". universe: runs
    whose barcodes make the rows (default every run); the other runs are
    joined to them. A run without a cell has "" in its columns.
    """
    runs = list(dict.fromkeys(r["run"] for r in run_rows))
    caps = list(dict.fromkeys(r["capture"] for r in run_rows))
    universe = list(universe or runs)
    cap_info = captures.set_index("capture") if captures is not None else None
    label_maps = _label_maps(run_rows)
    frames = []
    for cap in caps:
        calls = read_capture_runs([r for r in run_rows if r["capture"] == cap], label_maps)
        barcodes = []
        seen = set()
        for run in universe:
            if run in calls:
                for b in calls[run].barcode:
                    if b not in seen:
                        seen.add(b)
                        barcodes.append(b)
        d = pd.DataFrame({"barcode": barcodes})
        if truth_rows:
            t = load_cell_truth(truth_rows, cap)
            if t is not None:
                d = d.merge(t, on="barcode", how="left")
        for run in runs:
            cols = ["%s_%s" % (run, k) for k in ("status", "donor", "pair")]
            if run in calls:
                c = calls[run][["barcode", "status", "donor", "pair"]]
                c.columns = ["barcode"] + cols
                d = d.merge(c, on="barcode", how="left")
            else:
                for col in cols:
                    d[col] = ""
        if cap_info is not None:
            if cap not in cap_info.index:
                raise UsageError("capture %s is not in the capture table" % cap)
            pool, cset = cap_info.at[cap, "pool"], cap_info.at[cap, "capture_set"]
        else:
            pool, cset = cap, "main"
        d.insert(0, "capture_set", cset)
        d.insert(0, "pool", pool)
        d.insert(0, "capture", cap)
        frames.append(d)
        if progress:
            progress("%s %d" % (cap, len(d)))
    A = pd.concat(frames, ignore_index=True, sort=False)
    truth_cols = [c for c in A.columns if c not in ("capture", "pool", "capture_set", "barcode")
                  and not any(c == "%s_%s" % (r, k) for r in runs for k in ("status", "donor", "pair"))]
    order = ["capture", "pool", "capture_set", "barcode"] + truth_cols + [
        "%s_%s" % (r, k) for r in runs for k in ("status", "donor", "pair")]
    return A[order].fillna("")


# ---------------------------------------------------------------------------
# Views of one run


class Columns(object):
    """Names of the per-cell truth columns (all configurable on the command line)."""

    def __init__(self, umis="n_umi", sex="sex_call", chry="chrY_counts", xist="xist_counts",
                 cell_type="coarse_cell_type", ambiguous="ambiguous"):
        self.umis, self.sex, self.chry, self.xist = umis, sex, chry, xist
        self.cell_type, self.ambiguous = cell_type, ambiguous


def parse_tools(specs):
    """{run: [(label, status_run, pair_run), ...]} from "RUN=LABEL:STATUS_RUN[/PAIR_RUN],..." texts."""
    out = {}
    for spec in specs or ():
        if "=" not in spec:
            raise UsageError("--tools needs RUN=LABEL:STATUS_RUN[/PAIR_RUN],..., not %r" % spec)
        run, rest = spec.split("=", 1)
        tools = []
        for part in rest.split(","):
            if ":" not in part:
                raise UsageError("--tools: %r needs LABEL:STATUS_RUN" % part)
            label, runs = part.split(":", 1)
            status_run, _, pair_run = runs.partition("/")
            tools.append((label, status_run, pair_run or status_run))
        if not 1 <= len(tools) <= 2:
            raise UsageError("--tools: a run takes one or two tools")
        out[run] = tools
    return out


def _col(A, name):
    if name not in A.columns:
        raise UsageError("the cell table has no column %s" % name)
    return A[name]


def run_view(A, run, tools=None, cols=None):
    """Per cell for one run: tier, donor, group (doublet tool group), pair and UMI ratio.

    tier: the run's status, "not in run" where it has none. donor: set for
    singlets and recovered cells. group and pair: from `tools` (see the module
    docstring); without tools, group is "" and pair is the run's own.
    """
    cols = cols or Columns()
    tier = _col(A, run + "_status").replace({"": NOT_IN_RUN})
    donor = _col(A, run + "_donor").where(tier.isin(["singlet", "recovered"]), "")
    if tools:
        called = [(_col(A, s + "_status") == "doublet").to_numpy() for _, s, _ in tools]
        group = pd.Series(D.tool_groups(called, [lab for lab, _, _ in tools]), index=A.index)
        pair = pd.Series(D.pick_pair(called, [_col(A, p + "_pair").to_numpy() for _, _, p in tools]), index=A.index)
    else:
        group = pd.Series("", index=A.index)
        pair = _col(A, run + "_pair")
    umis = A["_umis"] if "_umis" in A.columns else pd.Series(np.nan, index=A.index)
    ctype = A[cols.cell_type] if cols.cell_type in A.columns else pd.Series("", index=A.index)
    ratio = D.umi_ratio(umis, (tier == "singlet").to_numpy(), A.capture, ctype)
    return pd.DataFrame({"tier": tier, "donor": donor, "group": group, "pair": pair, "ratio": ratio})


def prepare(A, cols, classes):
    """Numeric truth columns (_umis, _chry, _xist), the sex call (_sex) and the cell class (_class)."""
    A = A.copy()
    for src, dst in ((cols.umis, "_umis"), (cols.chry, "_chry"), (cols.xist, "_xist")):
        A[dst] = pd.to_numeric(A[src], errors="coerce") if src in A.columns else np.nan
    A["_sex"] = A[cols.sex].str.strip().str.lower() if cols.sex in A.columns else ""
    ctype = A[cols.cell_type] if cols.cell_type in A.columns else pd.Series("", index=A.index)
    A["_class"] = cell_class(ctype.tolist(), classes)
    return A


# ---------------------------------------------------------------------------
# Pools and shares


class Pools(object):
    """Pool indices of the cells and one resample draw per capture set."""

    def __init__(self, pool, capture_set, n_boot=1000, seed=1, interval_sets=None, min_pools=5, coverage=0.95,
                 level=0.95):
        self.codes, self.index = np.unique(np.asarray(pool, dtype=object).astype(str), return_inverse=True)
        self.index = np.asarray(self.index).ravel()
        self.capture_set = np.asarray(capture_set, dtype=object)
        self.n_boot, self.min_pools, self.coverage = n_boot, min_pools, coverage
        self.q = (100 * (1 - level) / 2.0, 100 * (1 - (1 - level) / 2.0))
        self.boot = {}
        sets = list(dict.fromkeys(self.capture_set)) if interval_sets is None else list(interval_sets)
        for s in sets:
            pools = np.unique(self.index[self.capture_set == s])
            if len(pools) == 0:
                continue
            rng = np.random.default_rng(seed)
            self.boot[s] = pools[rng.integers(0, len(pools), size=(n_boot, len(pools)))]

    def rate(self, num, den, capture_set):
        """{value, low, high, numerator, denominator} of num over den within one capture set."""
        den = np.asarray(den, dtype=bool) & (self.capture_set == capture_set)
        num = np.asarray(num, dtype=bool) & den
        n = np.bincount(self.index, weights=num.astype(float), minlength=len(self.codes))
        d = np.bincount(self.index, weights=den.astype(float), minlength=len(self.codes))
        low = high = np.nan
        boot = self.boot.get(capture_set)
        if boot is not None:
            db = d[boot].sum(1)
            bt = n[boot].sum(1)[db > 0] / db[db > 0]
            if len(bt) >= self.coverage * self.n_boot and (d > 0).sum() >= self.min_pools:
                low, high = np.percentile(bt, self.q[0]), np.percentile(bt, self.q[1])
        return {"value": n.sum() / max(d.sum(), 1), "low": low, "high": high,
                "numerator": int(n.sum()), "denominator": int(d.sum())}


def _sex_of(donors, default_sex):
    table = dict(zip(donors.donor, donors.sex)) if donors is not None else {}

    def sex_of(d):
        return table.get(d, default_sex) if d else ""
    return sex_of


def male_present(A, captures, sex_of):
    """Per cell: its capture's sample-sheet donors include a male donor; None when that is unknown.

    It is unknown without a capture table, or when the table gives no donors
    for a capture of the cell table.
    """
    if captures is None or "donors_given" not in captures.columns:
        return None
    given = dict(zip(captures.capture, captures.donors_given))
    if not all(given.get(c, False) for c in A.capture.unique()):
        return None
    caps = {c for c, ds in zip(captures.capture, captures.donors) if any(sex_of(d) == "male" for d in ds)}
    return A.capture.isin(caps)


# ---------------------------------------------------------------------------
# Scorecard


def scorecard(A, runs, tools=None, donors=None, captures=None, cols=None, classes=None, default_sex="unknown",
              n_boot=1000, seed=1, interval_sets=None, min_pools=5, level=0.95, umi_bins=None,
              recall_runs=None, notes=None, no_chry_fill="nearest"):
    """Score `runs` of cell table A. Returns (scores, lost_cells, donor_recall) DataFrames.

    See the module docstring and docs/experiment_tools/scorecard.md for every
    measure. donors: read_donors table; captures: read_captures table.
    notes: a list that receives a sentence for each group of measures left
    out and why (e.g. no capture table, so male presence per capture is unknown).
    no_chry_fill: how a UMI bin without male singlets gets its P(no chrY UMI),
    "nearest" or "highest" (doublets.nearest_bin_values).
    """
    notes = notes if notes is not None else []
    cols = cols or Columns()
    classes = classes if classes is not None else read_cell_classes()
    tools = tools or {}
    A = prepare(A, cols, classes)
    sets = list(dict.fromkeys(A.capture_set))
    pools = Pools(A.pool, A.capture_set, n_boot, seed, interval_sets, min_pools, level=level)
    sex_of = _sex_of(donors, default_sex)
    males = [d for d in donors.donor if sex_of(d) == "male"] if donors is not None else []
    immune_only = list(donors.donor[donors.immune_only]) if donors is not None else []
    decoys = list(donors.donor[donors.decoy]) if donors is not None else []
    mp = male_present(A, captures, sex_of)
    mp_known = mp is not None
    mp = mp.to_numpy() if mp_known else np.zeros(len(A), dtype=bool)
    if males and not mp_known:
        notes.append("the capture table gives no donors for every capture, so it is unknown which captures "
                     "hold a male donor; left out: " + "; ".join(NEEDS_MALE_PRESENCE))
    sex = A["_sex"]
    has_sex = sex.isin(["male", "female"]).to_numpy()
    female, male = (sex == "female").to_numpy(), (sex == "male").to_numpy()
    amb = (sex == cols.ambiguous.strip().lower()).to_numpy()
    sex_given = (sex != "").to_numpy()
    nonimmune = (A["_class"] == "non-immune").to_numpy()
    clear = A["_class"].isin(["immune", "non-immune"]).to_numpy()
    rows, lost, views = [], [], {}
    for run in runs:
        v = views[run] = run_view(A, run, tools.get(run), cols)
        tier, donor = v.tier.to_numpy(), v.donor
        ratio = v.ratio.to_numpy()
        lab_sex = np.array([sex_of(d) for d in donor], dtype=object)
        lab_male, lab_female = lab_sex == "male", lab_sex == "female"
        sing, rec, dbl = tier == "singlet", tier == "recovered", tier == "doublet"
        allc = np.ones(len(A), dtype=bool)
        ok_ratio = ~np.isnan(ratio)
        m = [("Singlets kept, share of all cells", "count", sing, allc),
             ("Recovered tier, share of all cells", "count", rec, allc),
             ("Doublet tier, share of all cells", "count", dbl, allc),
             ("Unassigned, share of all cells", "count", tier == "unassigned", allc)]
        if males:
            m.append(("Wrong donor by sex: female share of singlets labelled a male donor", "sex", female,
                      sing & lab_male & has_sex))
            for d in males:
                m.append(("Wrong donor by sex: female share of singlets labelled %s" % d, "sex", female,
                          sing & (donor == d).to_numpy() & has_sex))
            m += [("Wrong donor by sex: female share of recovered cells labelled a male donor", "sex", female,
                   rec & lab_male & has_sex),
                  ("Wrong donor by sex: male share of singlets labelled a female donor, captures with a male donor",
                   "sex", male, sing & lab_female & has_sex & mp),
                  ("Floor: male share of singlets, captures with no male donor", "sex", male, sing & has_sex & ~mp)]
        if immune_only:
            lab_io = donor.isin(immune_only).to_numpy()
            m += [("Wrong donor by cell type: non-immune share of singlets labelled an immune-only donor",
                   "cell type", nonimmune, sing & lab_io & clear),
                  ("Wrong donor by cell type: non-immune share of recovered cells labelled an immune-only donor",
                   "cell type", nonimmune, rec & lab_io & clear)]
        if decoys:
            m.append(("Wrong donor by decoy: share of singlets labelled a decoy donor", "decoy",
                      donor.isin(decoys).to_numpy(), sing))
        if males:
            m += [("Ambiguous-sex cells called doublet, captures with a male donor", "sex", dbl, amb & mp),
                  ("Ambiguous-sex cells called singlet, captures with a male donor", "sex", sing, amb & mp)]
        m += [("UMI ratio of doublet calls: share at 1.5 or more", "UMIs", ok_ratio & (ratio >= 1.5), dbl & ok_ratio),
              ("UMI ratio of singlets: share at 1.5 or more (reference)", "UMIs", ok_ratio & (ratio >= 1.5),
               sing & ok_ratio)]
        has_rec = bool(rec.any())
        if not mp_known:
            m = [x for x in m if x[0] not in NEEDS_MALE_PRESENCE]
        for cset in sets:
            for name, truth, num, den in m:
                if not has_rec and "recovered" in name.lower():
                    continue
                rows.append(dict(capture_set=cset, measure=name, run=run, truth=truth,
                                 **pools.rate(num, den, cset)))
            in_set = (A.capture_set == cset).to_numpy()
            if males:
                keep = in_set & sing & donor.isin(males).to_numpy() & sex_given & A["_umis"].notna().to_numpy()
                p0 = D.p_no_chry(A["_umis"], A["_chry"], keep)
                n, obs, exp, share, filled = D.no_male_share(v.pair, in_set & dbl, sex_given, A["_umis"],
                                                             A["_chry"], p0, males, return_filled=True,
                                                             fill=no_chry_fill)
                if n:
                    rows.append(dict(capture_set=cset, measure="Male-female doublet calls: estimated share with "
                                     "no male cell", run=run, truth="sex and UMIs", value=share, low=np.nan,
                                     high=np.nan, numerator=round((obs - exp) / (1 - exp / n)) if n > exp else np.nan,
                                     denominator=n, observed_no_chrY=obs, expected_no_chrY=round(exp, 1),
                                     calls_in_bins_without_data=filled))
            d_set = in_set & dbl
            if d_set.any():
                rows.append(dict(capture_set=cset, measure="UMI ratio of doublet calls: median", run=run,
                                 truth="UMIs", value=float(np.nanmedian(ratio[d_set])) if ok_ratio[d_set].any()
                                 else np.nan, low=np.nan, high=np.nan, numerator=np.nan,
                                 denominator=int(d_set.sum())))
            keep = in_set & (tier != NOT_IN_RUN)
            lost.append(_composition(A["_class"], tier, keep, run, cset, "cell class"))
            lost.append(_composition(_umi_bins(A["_umis"], umi_bins), tier, keep, run, cset, "UMIs"))
    if males and mp_known:
        for cset in sets:
            for name, num, den in (("Ambiguous-sex share of cells with a sex call value, captures with a male donor",
                                    amb, sex_given & mp),
                                   ("Floor: ambiguous-sex share of cells with a sex call value, captures with no "
                                    "male donor", amb, sex_given & ~mp)):
                rows.append(dict(capture_set=cset, measure=name, run="any", truth="sex",
                                 **pools.rate(num, den, cset)))
    S = pd.DataFrame(rows)
    for c in SCORE_COLUMNS:
        if c not in S.columns:
            S[c] = np.nan
    S = S[S.denominator.fillna(1) > 0][SCORE_COLUMNS]
    order = {name: i for i, name in enumerate(dict.fromkeys(S.measure))}
    S = S.assign(_o=S.measure.map(order), _s=S.capture_set.map({s: i for i, s in enumerate(sets)}))
    S = S.sort_values(["_s", "_o"], kind="mergesort").drop(columns=["_o", "_s"]).reset_index(drop=True)
    L = pd.concat(lost, ignore_index=True, sort=False)
    tiers = [t for t in TIERS if t in L.columns] + sorted(c for c in L.columns if c not in TIERS and c not in
                                                          ("capture_set", "run", "by", "group"))
    L = L[["capture_set", "run", "by", "group"] + tiers]
    recall_runs = list(recall_runs or [r for r in runs if r in views])
    unknown = [r for r in recall_runs if r not in views]
    if unknown:
        raise UsageError("--recall-runs names run(s) not scored: %s" % ", ".join(unknown))
    R = donor_recall(A, views, recall_runs, captures)
    return S, L, R


def _umi_bins(umis, edges=None):
    edges = list(edges or [0, 1000, 3000, 10000, np.inf])
    labels = []
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        if i == 0 and lo == 0:
            labels.append("under {:,.0f}".format(hi))
        elif np.isinf(hi):
            labels.append("{:,.0f} or more".format(lo))
        else:
            labels.append("{:,.0f} to {:,.0f}".format(lo, hi - 1))
    return pd.cut(umis, edges, right=False, labels=labels)


def _composition(group, tier, keep, run, cset, by):
    """Share of each tier's cells in each group (columns sum to 1), over the cells in keep with a group.

    A cell without a group (no UMI count for a UMI bin) is left out.
    """
    g = pd.Series(np.asarray(group, dtype=object), dtype=object)
    keep = np.asarray(keep, dtype=bool) & g.notna().to_numpy()
    g = g[keep]
    t = pd.Series(np.asarray(tier, dtype=object))[keep]
    if g.empty:
        return pd.DataFrame(columns=["capture_set", "run", "by", "group"])
    tab = pd.crosstab(g.astype(str), t, normalize="columns")
    if hasattr(group, "cat"):
        order = [str(c) for c in group.cat.categories if str(c) in tab.index]
        tab = tab.loc[order]
    tab = tab.reset_index()
    tab.columns = ["group"] + list(tab.columns[1:])
    tab.insert(0, "by", by)
    tab.insert(0, "run", run)
    tab.insert(0, "capture_set", cset)
    return tab


def donor_recall(A, views, runs, captures):
    """Per capture set, capture and sample-sheet donor: singlets, recovered cells and share per run, with flags.

    share: the donor's singlets over the capture's singlets. sister_median_share:
    the median share of the same donor in the capture's sister captures (same
    sister_key, same capture set). A run is flagged when its share is under a
    fifth of its sisters' median share or of another run's share, where that
    reference is 0.02 or more. A donor that a run had no genotype for
    (captures column genotyped_<run>) gets "no genotype in the <run> run" in
    place of that run's own flags.
    """
    if captures is None or not runs:
        return pd.DataFrame()
    info = captures.set_index("capture")
    rows = []
    for cap, idx in A.groupby("capture", sort=True).groups.items():
        if cap not in info.index:
            continue
        r0 = {"capture_set": info.at[cap, "capture_set"], "capture": cap, "sister_key": info.at[cap, "sister_key"]}
        for d in info.at[cap, "donors"]:
            row = dict(r0, donor=d)
            for run in runs:
                v = views[run].loc[idx]
                sing = v.tier == "singlet"
                total = int(sing.sum())
                n = int((sing & (v.donor == d)).sum())
                row["singlets_" + run] = n
                if (views[run].tier == "recovered").any():
                    row["recovered_" + run] = int(((v.tier == "recovered") & (v.donor == d)).sum())
                row["share_" + run] = n / float(max(total, 1))
                gcol = "genotyped_" + run
                if gcol in info.columns:
                    row["genotyped_in_" + run] = d in info.at[cap, gcol]
            rows.append(row)
    W = pd.DataFrame(rows)
    if W.empty:
        return W
    for run in runs:
        med = []
        for r in W.itertuples(index=False):
            sis = W[(W.capture_set == r.capture_set) & (W.sister_key == r.sister_key) & (W.capture != r.capture)
                    & (W.donor == r.donor)]
            med.append(sis["share_" + run].median() if len(sis) else np.nan)
        W["sister_median_share_" + run] = med
    flags = []
    for _, r in W.iterrows():
        why = {run: [] for run in runs}
        for run in runs:
            s = r["sister_median_share_" + run]
            if s >= 0.02 and r["share_" + run] < 0.2 * s:
                why[run].append("%s: under a fifth of its sister captures" % run)
            for other in runs:
                if other != run and r["share_" + other] >= 0.02 and r["share_" + run] < 0.2 * r["share_" + other]:
                    why[run].append("%s: under a fifth of the %s run" % (run, other))
        out = []
        tail = []
        for run in runs:
            if ("genotyped_in_" + run) in W.columns and not r["genotyped_in_" + run]:
                tail.append("no genotype in the %s run" % run)
            else:
                out += why[run]
        flags.append("; ".join(out + tail))
    W["flag"] = flags
    return W.drop(columns="sister_key")


# ---------------------------------------------------------------------------
# Calibration


def wilson(k, n, z=1.96):
    """Wilson score interval of k successes in n trials (NaN, NaN when n is 0)."""
    if n == 0:
        return np.nan, np.nan
    p = k / float(n)
    c = (p + z * z / (2 * n)) / (1 + z * z / n)
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / (1 + z * z / n)
    return c - h, c + h


def bin_labels(edges):
    """"0 to 1", "1 to 2.2", ..., "100 or more" for bin edges (the last may be inf)."""
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        out.append("%g or more" % lo if np.isinf(hi) else "%g to %g" % (lo, hi))
    return out


def calibrate(cells, score, label, sex, umis=None, unassigned=None, male_donors=(), edges=None, sex_of=None):
    """Wrong-donor rate by sex in bins of a per-cell score.

    cells: rows to score. score: float array; label: the donor each cell is
    scored as; sex: per-cell sex call; umis: float array (median per bin);
    unassigned: bool array of cells the run did not assign. A cell labelled a
    male donor with a female call is wrong; so is a cell labelled a female
    donor with a male call. Returns one row per bin.
    """
    edges = list(edges or [0, 1, 2.2, 5, 10, 20, 50, 100, np.inf])
    labels = bin_labels(edges)
    sex_of = sex_of or (lambda d: "male" if d in set(male_donors) else "female")
    b = pd.cut(pd.Series(score, dtype=float), edges, right=False, labels=labels)
    lab_sex = np.array([sex_of(d) for d in label], dtype=object)
    sex = np.asarray(sex, dtype=object)
    signal = np.isin(sex, ["male", "female"])
    un = np.zeros(len(sex), dtype=bool) if unassigned is None else np.asarray(unassigned, dtype=bool)
    um = np.full(len(sex), np.nan) if umis is None else np.asarray(umis, dtype=float)
    rows = []
    for lab in labels:
        x = (b == lab).to_numpy()
        m, f = x & (lab_sex == "male"), x & (lab_sex == "female")
        mk, mn = int((m & (sex == "female")).sum()), int((m & signal).sum())
        fk, fn = int((f & (sex == "male")).sum()), int((f & signal).sum())
        ma = m & ~un
        mak, man = int((ma & (sex == "female")).sum()), int((ma & signal).sum())
        lo1, hi1 = wilson(mk, mn)
        lo2, hi2 = wilson(fk, fn)
        rows.append({"bin": lab, "male_label_cells": int(m.sum()), "male_label_with_signal": mn,
                     "female_calls": mk, "female_share": mk / float(mn) if mn else np.nan,
                     "female_share_lo": lo1, "female_share_hi": hi1,
                     "male_label_assigned_with_signal": man, "female_calls_assigned": mak,
                     "female_share_assigned": mak / float(man) if man else np.nan,
                     "female_label_cells": int(f.sum()), "female_label_with_signal": fn, "male_calls": fk,
                     "male_share": fk / float(fn) if fn else np.nan, "male_share_lo": lo2, "male_share_hi": hi2,
                     "unassigned_share": un[x].mean() if x.any() else np.nan,
                     "median_umis": float(np.nanmedian(um[x])) if x.any() and np.isfinite(um[x]).any() else np.nan})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Doublet excess


def doublet_excess(A, runs, tools, reference, cols=None, classes=None, ambient=None):
    """One row per capture: doublets per run and tool group, against the loading expectation.

    Columns per run R: R_doublets, R_doublet (share), R_doublet_<group> and
    R_umi_ratio_<group> per tool group, R_largest_donor_in_pairs and
    R_named_pairs. For the reference run also: expected_total,
    expected_between_donors (its singlets' donor shares),
    expected_between_donors_cells, excess (observed share minus expected
    between donors), R_largest_donor_in_pairs_<group> and, with an ambient
    table, vireo_ambient (median ambient_other of the reference run's
    singlets) and souporcell_ambient.
    """
    cols = cols or Columns()
    classes = classes if classes is not None else read_cell_classes()
    A = prepare(A, cols, classes)
    views = {run: run_view(A, run, tools.get(run), cols) for run in runs}
    amb = ambient or {}
    with warnings.catch_warnings():
        # the median of a group with no UMI ratio is NaN, as it should be
        warnings.filterwarnings("ignore", message="Mean of empty slice")
        return pd.DataFrame(_excess_rows(A, runs, tools, reference, views, amb))


def _excess_rows(A, runs, tools, reference, views, amb):
    rows = []
    for cap, idx in A.groupby("capture", sort=True).groups.items():
        row = {"capture": cap, "pool": A.pool[idx[0]], "capture_set": A.capture_set[idx[0]], "cells": len(idx)}
        for run in runs:
            x = views[run].loc[idx]
            groups = D.group_names([lab for lab, _, _ in tools[run]]) if tools.get(run) else []
            shares = x.donor[x.tier == "singlet"].value_counts(normalize=True)
            dbl = x.tier == "doublet"
            row[run + "_doublets"] = int(dbl.sum())
            row[run + "_doublet"] = dbl.mean()
            for g in groups:
                d = dbl & (x.group == g)
                key = g.replace(" ", "_")
                row["%s_doublet_%s" % (run, key)] = d.mean()
                row["%s_umi_ratio_%s" % (run, key)] = x.ratio[d].median()
            row[run + "_largest_donor_in_pairs"], row[run + "_named_pairs"] = D.involvement(x.pair[dbl], shares)
            if run == reference:
                row["expected_total"] = D.expected_total(len(idx))
                row["expected_between_donors"] = row["expected_total"] * (1 - (shares ** 2).sum())
                row["expected_between_donors_cells"] = row["expected_between_donors"] * len(idx)
                row["excess"] = row[run + "_doublet"] - row["expected_between_donors"]
                for g in groups:
                    row["%s_largest_donor_in_pairs_%s" % (run, g.replace(" ", "_"))] = D.involvement(
                        x.pair[dbl & (x.group == g)], shares)[0]
                if cap in amb:
                    singlets = set(A.barcode[idx][(x.tier == "singlet").to_numpy()])
                    row["vireo_ambient"], row["souporcell_ambient"] = read_ambient(amb[cap], singlets)
        rows.append(row)
    return rows


def read_ambient(paths, singlets):
    """(vireo's median ambient_other over the singlets, souporcell's ambient share) of one capture.

    paths: {"vireo": ambient_summary.csv or "", "souporcell": ambient_rna.txt or ""}.
    """
    import re

    v = s = np.nan
    for key in ("vireo", "souporcell"):
        if paths.get(key) and not os.path.exists(paths[key]):
            raise UsageError("ambient file not found: %s" % paths[key])
    if paths.get("vireo"):
        a = pd.read_csv(paths["vireo"])
        v = a[a.barcode.isin(singlets)].ambient_other.median()
    if paths.get("souporcell"):
        with open(paths["souporcell"]) as f:
            m = re.search(r"estimated as ([0-9.eE+-]+)%", f.read())
        s = float(m.group(1)) / 100 if m else np.nan
    return v, s


def read_ambient_manifest(path):
    """{capture: {"vireo": path, "souporcell": path}} from columns capture, vireo_ambient, souporcell_ambient."""
    t = read_any_table(path)
    if "capture" not in t.columns:
        raise UsageError("ambient manifest %s lacks column capture" % path)
    out = {}
    for r in t.to_dict("records"):
        out[r["capture"]] = {"vireo": r.get("vireo_ambient", "") or "",
                             "souporcell": r.get("souporcell_ambient", "") or ""}
    return out


# ---------------------------------------------------------------------------
# Command-line helpers


def _add_truth_columns(parser):
    parser.add_argument("--umis-column", default="n_umi", help="UMIs per cell (default n_umi)")
    parser.add_argument("--sex-column", default="sex_call", help="sex call per cell (default sex_call)")
    parser.add_argument("--chry-column", default="chrY_counts", help="chrY UMIs per cell (default chrY_counts)")
    parser.add_argument("--cell-type-column", default="coarse_cell_type",
                        help="coarse cell type per cell (default coarse_cell_type)")
    parser.add_argument("--ambiguous-word", default="ambiguous",
                        help="the sex call of a cell with both chrY and XIST signal (default ambiguous)")
    parser.add_argument("--cell-classes", default=None,
                        help="table label,class (immune or non-immune); default tools/demux_tools/data/cell_classes.tsv")


def _columns(args):
    return Columns(args.umis_column, args.sex_column, args.chry_column, "xist_counts", args.cell_type_column,
                   args.ambiguous_word)


def _floats(text):
    return [float(x) for x in text.split(",")]


def _read_cells(path):
    try:
        return read_table(path)
    except (OSError, ValueError) as e:
        raise UsageError(str(e))


# ---------------------------------------------------------------------------
# Command: cell-table


def configure_cell_table(parser):
    parser.description = ("Write one row per cell of every capture: capture, pool, capture_set, barcode, the "
                          "per-cell truth columns, then <run>_status, <run>_donor and <run>_pair for each run.")
    parser.add_argument("--runs", required=True,
                        help="table with columns capture, run, path and optionally kind, prob_doublet, label_map, "
                             "label_pattern, cluster_map (table, sba:PATH or majority:RUN) and sba")
    parser.add_argument("--truth", default=None, help="table with columns capture, path and optionally columns")
    parser.add_argument("--captures", default=None,
                        help="table with column capture and optionally pool, capture_set, donors, sister_key")
    parser.add_argument("--universe", default=None,
                        help="comma-separated runs whose barcodes make the rows (default every run)")
    parser.add_argument("--out", required=True, help="output table (.gz for gzip)")


def run_cell_table(args):
    clock = start_clock()
    try:
        rows = read_runs_manifest(args.runs)
        truth = read_truth_manifest(args.truth) if args.truth else None
        caps = read_captures(args.captures) if args.captures else None
    except (ValueError, OSError) as e:
        raise UsageError(str(e))
    inputs = [args.runs, args.truth, args.captures]
    for r in rows:
        cmap = "" if r["cluster_map"].startswith("majority:") else r["cluster_map"]
        inputs += files_read(r["path"], r["kind"] or None, r["prob_doublet"] or None, r["label_map"] or None,
                             cmap or None, r["sba"] or None)
    for r in truth or []:
        inputs.append(r["path"])
    check_output_path(args.out, inputs)
    try:
        A = build_cell_table(rows, truth, caps, split_list(args.universe) if args.universe else None)
    except (ValueError, OSError, KeyError) as e:
        raise UsageError(str(e))
    write_table(A, args.out)
    runs = list(dict.fromkeys(r["run"] for r in rows))
    summary = {"cells": len(A), "captures": int(A.capture.nunique()),
               "status": {run: A[run + "_status"].replace({"": NOT_IN_RUN}).value_counts().to_dict() for run in runs}}
    write_provenance(args.out, args.argv, [p for p in dict.fromkeys(inputs) if p], clock=clock, extra=summary)
    print("cell-table: %d cells in %d captures, runs %s" % (len(A), summary["captures"], ", ".join(runs)))
    return 0


# ---------------------------------------------------------------------------
# Command: scorecard


def configure_scorecard(parser):
    parser.description = ("Score runs of a cell table against truth. Writes OUT/scorecard.tsv, "
                          "OUT/lost_cells.tsv, OUT/donor_recall.tsv and OUT/provenance.json.")
    parser.add_argument("--cells", required=True, help="cell table from cell-table")
    parser.add_argument("--run", action="append", required=True, help="a run to score; repeatable")
    parser.add_argument("--tools", action="append", default=[],
                        help="RUN=LABEL:STATUS_RUN[/PAIR_RUN][,LABEL:STATUS_RUN[/PAIR_RUN]]: tool runs that split "
                             "RUN's doublets into groups and give their pairs; repeatable")
    parser.add_argument("--donors", default=None, help="table donor, sex, immune_only, decoy")
    parser.add_argument("--default-sex", default="unknown",
                        help="sex of a donor the donor table does not list (default unknown)")
    parser.add_argument("--captures", default=None,
                        help="table capture, pool, capture_set, donors, sister_key, genotyped_<run>")
    parser.add_argument("--n-boot", type=int, default=1000, help="pool resamples (default 1000)")
    parser.add_argument("--seed", type=int, default=1, help="seed of numpy's default_rng (default 1)")
    parser.add_argument("--interval-sets", default=None,
                        help="comma-separated capture sets that get intervals (default every set)")
    parser.add_argument("--min-pools", type=int, default=5, help="pools a share needs for an interval (default 5)")
    parser.add_argument("--umi-bins", default="0,1000,3000,10000,inf", help="UMI bin edges for lost_cells.tsv")
    parser.add_argument("--recall-runs", default=None, help="comma-separated runs for donor_recall.tsv (default all)")
    parser.add_argument("--no-chry-fill", choices=D.FILLS, default="nearest",
                        help="P(no chrY UMI) of a UMI bin without male singlets: the nearest bin's (default) or the "
                             "highest bin's, as the first scorecard did")
    _add_truth_columns(parser)
    parser.add_argument("--out", required=True, help="output folder (created if needed)")


def run_scorecard(args):
    clock = start_clock()
    A = _read_cells(args.cells)
    try:
        donors = read_donors(args.donors) if args.donors else None
        caps = read_captures(args.captures) if args.captures else None
        classes = read_cell_classes(args.cell_classes)
    except (ValueError, OSError) as e:
        raise UsageError(str(e))
    tools = parse_tools(args.tools)
    os.makedirs(args.out, exist_ok=True)
    outs = [os.path.join(args.out, n) for n in ("scorecard.tsv", "lost_cells.tsv", "donor_recall.tsv")]
    inputs = [p for p in (args.cells, args.donors, args.captures, args.cell_classes) if p]
    for o in outs:
        check_output_path(o, inputs)
    notes = []
    S, L, R = scorecard(A, args.run, tools, donors, caps, _columns(args), classes, args.default_sex.lower(),
                        args.n_boot, args.seed, split_list(args.interval_sets) if args.interval_sets else None,
                        args.min_pools, umi_bins=_floats(args.umi_bins),
                        recall_runs=split_list(args.recall_runs) if args.recall_runs else None, notes=notes,
                        no_chry_fill=args.no_chry_fill)
    for n in notes:
        sys.stderr.write("scorecard: %s\n" % n)
    S = S.copy()
    for c in ("numerator", "denominator", "observed_no_chrY", "calls_in_bins_without_data"):
        S[c] = S[c].round().astype("Int64")
    write_table(S, outs[0], float_formats={"value": "%.6g", "low": "%.6g", "high": "%.6g",
                                           "expected_no_chrY": "%.1f"})
    write_table(L, outs[1], float_formats={c: "%.6g" for c in L.columns[4:]})
    fr = {c: "%.6g" for c in R.columns if c.startswith("share_") or c.startswith("sister_median")}
    write_table(R, outs[2], float_formats=fr)
    write_provenance(args.out, args.argv, inputs, clock=clock,
                     extra={"rows": len(S), "runs": list(args.run), "cells": len(A), "outputs": outs, "notes": notes})
    print("scorecard: %d rows for runs %s in %s" % (len(S), ", ".join(args.run), args.out))
    return 0


# ---------------------------------------------------------------------------
# Command: calibrate


def configure_calibrate(parser):
    parser.description = ("Wrong-donor rate by sex in bins of a per-cell score (such as the Scorer's margin), "
                          "with Wilson 95 per cent intervals. Cells come from per-capture tables listed in a "
                          "manifest; sex calls from a truth manifest.")
    parser.add_argument("--cells", required=True, help="table with columns capture and path (per-capture cell tables)")
    parser.add_argument("--truth", required=True, help="truth manifest (capture, path[, columns]) with the sex calls")
    parser.add_argument("--donors", required=True, help="donor table (donor, sex)")
    parser.add_argument("--default-sex", default="unknown", help="sex of a donor not in the donor table")
    parser.add_argument("--captures", default=None, help="capture table, for --sets")
    parser.add_argument("--sets", default=None, help="comma-separated capture sets to keep (default all)")
    parser.add_argument("--barcode-column", default="cell", help="barcode column of the cell tables (default cell)")
    parser.add_argument("--score-column", default="margin", help="score column (default margin)")
    parser.add_argument("--label-column", default="best_singlet", help="donor each cell is scored as")
    parser.add_argument("--call-column", default="donor_id", help="the run's call (default donor_id)")
    parser.add_argument("--exclude", default="doublet", help="comma-separated calls left out (default doublet)")
    parser.add_argument("--unassigned", default="unassigned", help="call of an unassigned cell (default unassigned)")
    parser.add_argument("--umis-column", default="umis", help="UMIs column of the cell tables (default umis)")
    parser.add_argument("--sex-column", default="sex_call", help="sex call column of the truth (default sex_call)")
    parser.add_argument("--bins", default="0,1,2.2,5,10,20,50,100,inf", help="bin edges of the score")
    parser.add_argument("--out", required=True, help="output table")


def run_calibrate(args):
    clock = start_clock()
    try:
        man = read_any_table(args.cells)
        truth = read_truth_manifest(args.truth)
        donors = read_donors(args.donors)
        caps = read_captures(args.captures) if args.captures else None
    except (ValueError, OSError) as e:
        raise UsageError(str(e))
    if "capture" not in man.columns or "path" not in man.columns:
        raise UsageError("%s needs columns capture and path" % args.cells)
    sex_of = _sex_of(donors, args.default_sex.lower())
    keep_sets = split_list(args.sets) if args.sets else None
    cset = dict(zip(caps.capture, caps.capture_set)) if caps is not None else {}
    with_truth = {r["capture"] for r in truth}
    inputs = [args.cells, args.truth, args.donors, args.captures]
    frames, used = [], []
    for cap, path in zip(man.capture, man.path):
        if keep_sets is not None and cset.get(cap, "main") not in keep_sets:
            continue
        if cap not in with_truth:
            continue
        d = read_table(path)
        need = [args.barcode_column, args.score_column, args.label_column, args.call_column]
        missing = [c for c in need if c not in d.columns]
        if missing:
            raise UsageError("%s lacks column(s) %s" % (path, ", ".join(missing)))
        if not any(sex_of(x) == "male" for x in set(d[args.label_column])):
            continue
        try:
            t = load_cell_truth(truth, cap)
        except (ValueError, OSError) as e:
            raise UsageError(str(e))
        if t is not None and args.sex_column not in t.columns:
            raise UsageError("the truth of capture %s has no column %s (columns: %s)"
                             % (cap, args.sex_column, ", ".join(c for c in t.columns if c != "barcode")))
        sex = dict(zip(t.barcode, t[args.sex_column].str.strip().str.lower())) if t is not None else {}
        d = d[~d[args.call_column].isin(split_list(args.exclude))]
        frames.append(pd.DataFrame({
            "capture": cap, "score": pd.to_numeric(d[args.score_column], errors="coerce"),
            "label": d[args.label_column], "sex": [sex.get(b, "") for b in d[args.barcode_column]],
            "unassigned": d[args.call_column] == args.unassigned,
            "umis": pd.to_numeric(d[args.umis_column], errors="coerce") if args.umis_column in d.columns else np.nan}))
        used.append(cap)
        inputs.append(path)
        inputs += [r["path"] for r in truth if r["capture"] == cap and r["path"]]
    if not frames:
        raise UsageError("no capture has cells labelled a male donor and a truth file")
    X = pd.concat(frames, ignore_index=True)
    T = calibrate(X, X.score.to_numpy(), X.label.tolist(), X.sex.to_numpy(), X.umis.to_numpy(),
                  X.unassigned.to_numpy(), edges=_floats(args.bins), sex_of=sex_of)
    check_output_path(args.out, inputs)
    fl = {c: "%.6g" for c in T.columns if T[c].dtype.kind == "f"}
    write_table(T, args.out, float_formats=fl)
    write_provenance(args.out, args.argv, [p for p in inputs if p], clock=clock,
                     extra={"captures": len(used), "cells": len(X)})
    print("calibrate: %d captures, %d cells; %d bins in %s" % (len(used), len(X), len(T), args.out))
    return 0


# ---------------------------------------------------------------------------
# Command: doublet-excess


def configure_doublet_excess(parser):
    parser.description = ("Per capture: doublet calls of each run split by the tool that made them, against "
                          "the loading expectation (0.8 per cent per 1,000 cells times 1 minus the sum of squared "
                          "donor shares), with the largest donor's involvement, UMI ratios and ambient estimates.")
    parser.add_argument("--cells", required=True, help="cell table from cell-table")
    parser.add_argument("--run", action="append", required=True, help="a run; repeatable")
    parser.add_argument("--tools", action="append", default=[], help="as in scorecard")
    parser.add_argument("--reference-run", default=None, help="run for the expectation and excess (default first)")
    parser.add_argument("--ambient", default=None,
                        help="table capture, vireo_ambient (ambient_summary.csv), souporcell_ambient (ambient_rna.txt)")
    _add_truth_columns(parser)
    parser.add_argument("--out", required=True, help="output table")


def run_doublet_excess(args):
    clock = start_clock()
    reference = args.reference_run or args.run[0]
    if reference not in args.run:
        raise UsageError("--reference-run %s is not one of the --run values (%s)" % (reference, ", ".join(args.run)))
    A = _read_cells(args.cells)
    tools = parse_tools(args.tools)
    try:
        amb = read_ambient_manifest(args.ambient) if args.ambient else None
        classes = read_cell_classes(args.cell_classes)
    except (ValueError, OSError) as e:
        raise UsageError(str(e))
    inputs = [args.cells, args.ambient, args.cell_classes]
    if amb:
        for p in amb.values():
            inputs += [x for x in p.values() if x]
    check_output_path(args.out, inputs)
    T = doublet_excess(A, args.run, tools, reference, _columns(args), classes, amb)
    fl = {c: "%.6g" for c in T.columns if T[c].dtype.kind == "f"}
    write_table(T, args.out, float_formats=fl)
    write_provenance(args.out, args.argv, [p for p in inputs if p], clock=clock,
                     extra={"captures": len(T), "runs": list(args.run)})
    print("doublet-excess: %d captures in %s" % (len(T), args.out))
    return 0
