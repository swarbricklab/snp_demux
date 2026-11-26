"""demux_tools.truth, doublets and scorecard: shares, intervals, doublet arithmetic, calibration, the commands."""
import json

import numpy as np
import pandas as pd
import pytest

from demux_tools import doublets as D
from demux_tools import scorecard as S
from demux_tools import truth as T


def _write(path, rows, sep="\t"):
    with open(str(path), "w") as f:
        f.write("\n".join(sep.join(str(x) for x in r) for r in rows) + "\n")
    return str(path)


# ---------------------------------------------------------------------------
# truth


def test_cell_classes_file_is_generic_and_complete():
    classes = T.read_cell_classes()
    assert classes["T cell"] == "immune" and classes["fibroblast"] == "non-immune"
    assert set(classes.values()) == {"immune", "non-immune"}
    assert len(classes) == 33
    out = T.cell_class(["T cell", "epithelial cell", "", "stem cell"], classes)
    assert out.tolist() == ["immune", "non-immune", "no annotation", "unclear"]


def test_captures_pool_falls_back_to_the_donor_set(tmp_path):
    p = _write(tmp_path / "captures.tsv", [["capture", "pool", "capture_set", "donors"],
                                           ["c1", "P1", "", "b,a"], ["c2", "", "stale", "b,a"],
                                           ["c3", "", "", "a,b"]])
    c = T.read_captures(p).set_index("capture")
    assert c.at["c1", "pool"] == "P1" and c.at["c1", "capture_set"] == "main"
    assert c.at["c2", "pool"] == c.at["c3", "pool"] == "donors:a,b"
    assert c.at["c2", "capture_set"] == "stale"
    assert c.at["c1", "sister_key"] == "P1|a,b"


def test_donors_flags_and_truth_join(tmp_path):
    d = T.read_donors(_write(tmp_path / "donors.tsv", [["donor", "sex", "immune_only", "decoy"],
                                                      ["m", "Male", "true", ""], ["f", "female", "", "yes"]]))
    assert d.sex.tolist() == ["male", "female"]
    assert d.immune_only.tolist() == [True, False] and d.decoy.tolist() == [False, True]
    sex = _write(tmp_path / "sex.csv", [["barcode", "n_umi", "sex_call"], ["A", "10", "male"], ["B", "20", "female"]], ",")
    ct = _write(tmp_path / "ct.csv", [["barcode", "coarse_cell_type", "fine"], ["B", "T cell", "x"], ["C", "B cell", "y"]], ",")
    man = _write(tmp_path / "truth.tsv", [["capture", "path", "columns"], ["c1", sex, ""], ["c1", ct, "coarse_cell_type"]])
    t = T.load_cell_truth(T.read_truth_manifest(man), "c1").set_index("barcode")
    assert list(t.columns) == ["n_umi", "sex_call", "coarse_cell_type"]
    assert t.at["B", "coarse_cell_type"] == "T cell" and t.at["A", "sex_call"] == "male"
    assert T.load_cell_truth(T.read_truth_manifest(man), "other") is None


# ---------------------------------------------------------------------------
# doublets


def test_expected_doublets_and_involvement():
    assert D.expected_total(10000) == pytest.approx(0.08)
    assert D.expected_between(10000, [0.5, 0.5]) == pytest.approx(0.04)
    shares = pd.Series([0.5, 0.25, 0.25], index=["a", "b", "c"])
    # 2p(1-p)/(1-sum p^2) = 0.5 / 0.625 = 0.8; observed: a in 16 of 20 pairs
    pairs = ["a,b"] * 16 + ["b,c"] * 4 + ["", "x"]
    ratio, n = D.involvement(pairs, shares)
    assert n == 20 and ratio == pytest.approx(0.8 / 0.8)
    assert np.isnan(D.involvement(pairs[:19], shares)[0])


def test_no_male_share_worked_example():
    # 100 male-female calls at 4,000 UMIs; a real doublet's male half (2,000 UMIs) lacks chrY with p = 0.25.
    # 60 calls hold no male cell (no chrY), 40 do, and 10 of those 40 lack chrY: O = 70, E = 25.
    n = 100
    pairs = pd.Series(["m,f"] * n)
    umis = pd.Series([4000.0] * n)
    chry = pd.Series([0.0] * 70 + [3.0] * 30)
    p0 = pd.Series({8.0: 0.25})          # bin 2000 // 250 = 8
    calls, obs, exp, share = D.no_male_share(pairs, np.ones(n, bool), np.ones(n, bool), umis, chry, p0, ["m"])
    assert (calls, obs) == (100, 70) and exp == pytest.approx(25.0)
    assert share == pytest.approx(0.6)   # (70 - 25) / (100 - 25); (O - E) / n would say 0.45
    # a pair of two male donors is never counted; a bin beyond p0 takes its last value
    calls2, _, exp2, _ = D.no_male_share(pd.Series(["m,m2", "m,f"]), np.ones(2, bool), np.ones(2, bool),
                                         pd.Series([4000.0, 90000.0]), pd.Series([0.0, 0.0]), p0, ["m"], ["m", "m2"])
    assert calls2 == 1 and exp2 == pytest.approx(0.25)


def test_p_no_chry_bins_and_umi_ratio():
    umis = pd.Series([100.0, 300.0, 300.0, 30000.0])
    chry = pd.Series([0.0, 0.0, 2.0, 0.0])
    p = D.p_no_chry(umis, chry, np.array([True, True, True, True]))
    assert p.to_dict() == {0.0: 1.0, 1.0: 0.5, 80.0: 1.0}
    r = D.umi_ratio(pd.Series([100.0, 300.0, 200.0, 50.0]), np.array([True, True, False, True]),
                    pd.Series(["c", "c", "c", "c"]), pd.Series(["T", "T", "T", "B"]))
    assert r.tolist() == pytest.approx([0.5, 1.5, 1.0, 1.0])


def test_tool_groups_and_pair_choice():
    v = np.array([True, True, False, False])
    s = np.array([True, False, True, False])
    assert D.tool_groups([v, s], ["vireo", "souporcell"]).tolist() == [
        "both tools", "vireo only", "souporcell only", "neither vireo nor souporcell"]
    pv = np.array(["a,b", "a,c", "a,d", "a,e"], dtype=object)
    ps = np.array(["x,y", "", "b,c", ""], dtype=object)
    assert D.pick_pair([v, s], [pv, ps]).tolist() == ["a,b", "a,c", "b,c", "a,e"]


# ---------------------------------------------------------------------------
# shares and intervals


def test_rate_and_pool_intervals():
    pools = np.repeat(["p%d" % i for i in range(6)], 10)
    sets = np.array(["main"] * 50 + ["other"] * 10)
    P = S.Pools(pools, sets, n_boot=200, seed=1)
    num = np.zeros(60, bool)
    num[::5] = True
    r = P.rate(num, np.ones(60, bool), "main")
    assert (r["numerator"], r["denominator"]) == (10, 50) and r["value"] == pytest.approx(0.2)
    assert r["low"] <= 0.2 <= r["high"]
    # the same seed gives the same interval; another seed draws other pools
    r2 = S.Pools(pools, sets, n_boot=200, seed=1).rate(num, np.ones(60, bool), "main")
    assert (r["low"], r["high"]) == (r2["low"], r2["high"])
    assert not np.array_equal(S.Pools(pools, sets, n_boot=200, seed=2).boot["main"], P.boot["main"])
    # under five pools: no interval
    few = P.rate(num, np.arange(60) < 40, "main")
    assert np.isnan(few["low"]) and few["denominator"] == 40
    assert np.isnan(P.rate(num, np.ones(60, bool), "other")["low"])   # one pool


def test_wilson_and_bin_labels():
    lo, hi = S.wilson(5, 50)
    assert lo == pytest.approx(0.04348, abs=1e-5) and hi == pytest.approx(0.21360, abs=1e-5)
    assert np.isnan(S.wilson(0, 0)[0])
    assert S.bin_labels([0, 1, 2.2, np.inf]) == ["0 to 1", "1 to 2.2", "2.2 or more"]


def test_calibrate_by_hand():
    score = np.array([0.5, 0.5, 3, 3, 3, 60, 60])
    label = ["m", "m", "m", "f", "f", "m", "f"]
    sex = np.array(["female", "male", "female", "male", "", "male", "female"], dtype=object)
    un = np.array([True, True, False, False, False, False, False])
    out = S.calibrate(None, score, label, sex, np.arange(7.0), un, male_donors=["m"], edges=[0, 1, 5, np.inf])
    r = out.set_index("bin")
    assert r.at["0 to 1", "male_label_with_signal"] == 2 and r.at["0 to 1", "female_calls"] == 1
    assert r.at["0 to 1", "male_label_assigned_with_signal"] == 0
    assert r.at["1 to 5", "female_label_cells"] == 2 and r.at["1 to 5", "female_label_with_signal"] == 1
    assert r.at["1 to 5", "male_share"] == 1.0 and r.at["5 or more", "female_share"] == 0.0
    assert r.at["0 to 1", "unassigned_share"] == 1.0 and r.at["1 to 5", "median_umis"] == 3.0


# ---------------------------------------------------------------------------
# cell table and the commands, on files of three small runs


def _runs(tmp_path):
    """Two captures; runs: cons (hier), vir (vireo), soup (souporcell named by majority of cons)."""
    files = []
    for cap, shift in (("c1", 0), ("c2", 1)):
        d = tmp_path / cap
        d.mkdir()
        bcs = ["%s-%d" % (cap, i) for i in range(8)]
        donors = ["m", "f", "m", "f", "", "", "f", "m"]
        status = ["singlet", "singlet", "singlet", "singlet", "doublet", "unassigned", "recovered", "singlet"]
        hier = [["barcode", "cluster", "donor", "label_confidence", "status"]]
        for b, dn, st in zip(bcs, donors, status):
            hier.append([b, "", dn, "1.0", st])
        _write(d / "hier.tsv", hier)
        soup = [["barcode", "status", "assignment", "log_prob_singleton", "log_prob_doublet"]]
        for i, b in enumerate(bcs):
            cl = {"m": "0", "f": "1"}.get(donors[i], "0")
            st = "doublet" if i == 4 else "singlet"
            soup.append([b, st, "0/1" if i == 4 else cl, "-1", "-2"])
        _write(d / "clusters.tsv", soup)
        files.append((cap, str(d / "hier.tsv"), str(d / "clusters.tsv")))
    rows = [["capture", "run", "kind", "path", "cluster_map"]]
    for cap, h, s in files:
        rows += [[cap, "cons", "hier", h, ""], [cap, "soup", "souporcell", s, "majority:cons"]]
    runs = _write(tmp_path / "runs.tsv", rows)
    truth_rows = [["capture", "path"]]
    for cap, _, _ in files:
        sex = [["barcode", "n_umi", "chrY_counts", "xist_counts", "sex_call", "coarse_cell_type"]]
        calls = ["male", "female", "female", "female", "male", "", "male", "male"]
        for i in range(8):
            sex.append(["%s-%d" % (cap, i), 1000 + 100 * i, 0 if calls[i] != "male" else 2, 1, calls[i],
                        "T cell" if i % 2 else "fibroblast"])
        truth_rows.append([cap, _write(tmp_path / ("%s_truth.csv" % cap), sex, ",")])
    truth = _write(tmp_path / "truth.tsv", truth_rows)
    caps = _write(tmp_path / "captures.tsv", [["capture", "pool", "capture_set", "donors"],
                                             ["c1", "P1", "main", "m,f"], ["c2", "P2", "main", "m,f"]])
    donors = _write(tmp_path / "donors.tsv", [["donor", "sex", "immune_only"], ["m", "male", "true"],
                                             ["f", "female", "false"]])
    return runs, truth, caps, donors


def test_cell_table_names_clusters_by_majority(tmp_path):
    runs, truth, caps, _ = _runs(tmp_path)
    A = S.build_cell_table(S.read_runs_manifest(runs), T.read_truth_manifest(truth), T.read_captures(caps))
    assert len(A) == 16 and A.capture.tolist()[:1] == ["c1"]
    assert list(A.columns[:4]) == ["capture", "pool", "capture_set", "barcode"]
    row = A.set_index("barcode").loc["c1-0"]
    assert (row.cons_status, row.cons_donor, row.soup_donor, row.sex_call) == ("singlet", "m", "m", "male")
    assert A.set_index("barcode").loc["c1-4"].soup_pair == "f,m"


def test_scorecard_by_hand(tmp_path):
    runs, truth, caps, donors = _runs(tmp_path)
    A = S.build_cell_table(S.read_runs_manifest(runs), T.read_truth_manifest(truth), T.read_captures(caps))
    sc, lost, recall = S.scorecard(A, ["cons"], {"cons": [("soup", "soup", "soup")]}, T.read_donors(donors),
                                   T.read_captures(caps), n_boot=50)
    r = sc.set_index(["measure", "run"])
    kept = r.loc[("Singlets kept, share of all cells", "cons")]
    assert (kept.numerator, kept.denominator) == (10, 16)
    fem = r.loc[("Wrong donor by sex: female share of singlets labelled a male donor", "cons")]
    # singlets labelled m: cells 0, 2, 7 per capture; sex calls male, female, male
    assert (fem.numerator, fem.denominator) == (2, 6)
    # no recovered cell is labelled a male donor, so that row has no denominator and is left out
    assert ("Wrong donor by sex: female share of recovered cells labelled a male donor", "cons") not in r.index
    io = r.loc[("Wrong donor by cell type: non-immune share of singlets labelled an immune-only donor", "cons")]
    assert (io.numerator, io.denominator) == (4, 6)   # cells 0 and 2 are fibroblasts
    assert np.isnan(kept.low)   # two pools: no interval
    assert set(lost.by) == {"cell class", "UMIs"}
    assert set(recall.donor) == {"m", "f"} and "flag" in recall.columns


def test_commands_end_to_end(tmp_path, run_tools):
    runs, truth, caps, donors = _runs(tmp_path)
    cells = tmp_path / "cells.tsv.gz"
    res = run_tools("cell-table", "--runs", runs, "--truth", truth, "--captures", caps, "--out", cells)
    assert res.returncode == 0, res.stderr
    out = tmp_path / "score"
    res = run_tools("scorecard", "--cells", cells, "--run", "cons", "--tools", "cons=soup:soup",
                    "--donors", donors, "--captures", caps, "--n-boot", "50", "--out", out)
    assert res.returncode == 0, res.stderr
    sc = pd.read_csv(str(out / "scorecard.tsv"), sep="\t")
    assert "Singlets kept, share of all cells" in set(sc.measure)
    assert json.load(open(str(out / "provenance.json")))["results"]["runs"] == ["cons"]
    res = run_tools("doublet-excess", "--cells", cells, "--run", "cons", "--tools", "cons=soup:soup",
                    "--out", tmp_path / "excess.tsv")
    assert res.returncode == 0, res.stderr
    ex = pd.read_csv(str(tmp_path / "excess.tsv"), sep="\t")
    assert ex.capture.tolist() == ["c1", "c2"] and ex.cons_doublets.tolist() == [1, 1]
    assert ex.expected_total.tolist() == pytest.approx([0.008 * 8 / 1000] * 2)
    scorer = [["cell", "donor_id", "best_singlet", "margin", "umis"]]
    for i in range(8):
        scorer.append(["c1-%d" % i, "doublet" if i == 4 else "m", "m" if i % 2 == 0 else "f", 3 * i, 100])
    man = _write(tmp_path / "scorer_cells.tsv", [["capture", "path"], ["c1", _write(tmp_path / "s1.tsv", scorer)]])
    res = run_tools("calibrate", "--cells", man, "--truth", truth, "--donors", donors, "--bins", "0,5,inf",
                    "--out", tmp_path / "cal.tsv")
    assert res.returncode == 0, res.stderr
    cal = pd.read_csv(str(tmp_path / "cal.tsv"), sep="\t").set_index("bin")
    # male-labelled (even i, not 4): i = 0 (score 0), 2 (6), 6 (18): sex male, female, male
    assert cal.at["0 to 5", "male_label_with_signal"] == 1 and cal.at["5 or more", "female_calls"] == 1


# ---------------------------------------------------------------------------
# Review fixes


def _cells(rows, capture="c1", pool="P1", cset="main", run="r"):
    """A cell table from (barcode, status, donor, pair, sex_call, n_umi, chrY_counts) tuples."""
    t = pd.DataFrame(rows, columns=["barcode", "%s_status" % run, "%s_donor" % run, "%s_pair" % run, "sex_call",
                                    "n_umi", "chrY_counts"])
    t.insert(0, "capture_set", cset)
    t.insert(0, "pool", pool)
    t.insert(0, "capture", capture)
    return t.astype(str)


def _two_capture_table(tmp_path):
    runs, truth, caps, donors = _runs(tmp_path)
    A = S.build_cell_table(S.read_runs_manifest(runs), T.read_truth_manifest(truth), T.read_captures(caps))
    return A, donors


def test_floor_measures_need_male_presence(tmp_path):
    A, donors = _two_capture_table(tmp_path)
    floor = "Floor: male share of singlets, captures with no male donor"
    notes = []
    sc, _, _ = S.scorecard(A, ["cons"], donors=T.read_donors(donors), n_boot=10, notes=notes)
    assert floor not in set(sc.measure) and notes and "unknown which captures" in notes[0]
    assert not any(m in set(sc.measure) for m in S.NEEDS_MALE_PRESENCE)
    no_donors = _write(tmp_path / "caps_nodonors.tsv", [["capture", "pool"], ["c1", "P1"], ["c2", "P2"]])
    sc, _, _ = S.scorecard(A, ["cons"], donors=T.read_donors(donors), captures=T.read_captures(no_donors), n_boot=10)
    assert floor not in set(sc.measure)
    # c2's sheet holds only the female donor: its singlets give the floor
    caps = _write(tmp_path / "caps2.tsv", [["capture", "pool", "donors"], ["c1", "P1", "m,f"], ["c2", "P2", "f"]])
    sc, _, _ = S.scorecard(A, ["cons"], donors=T.read_donors(donors), captures=T.read_captures(caps), n_boot=10)
    r = sc.set_index("measure").loc[floor]
    # c2 singlets with a male or female call: cells 0, 1, 2, 3, 7 (male, female, female, female, male)
    assert (r.numerator, r.denominator) == (2, 5)
    assert "Ambiguous-sex share of cells with a sex call value, captures with a male donor" in set(sc.measure)


def test_scorecard_command_says_why_measures_are_left_out(tmp_path, run_tools):
    runs, truth, caps, donors = _runs(tmp_path)
    cells = tmp_path / "cells.tsv.gz"
    assert run_tools("cell-table", "--runs", runs, "--truth", truth, "--out", cells).returncode == 0
    res = run_tools("scorecard", "--cells", cells, "--run", "cons", "--donors", donors, "--n-boot", "10",
                    "--out", tmp_path / "sc")
    assert res.returncode == 0, res.stderr
    assert "unknown which captures hold a male donor" in res.stderr
    prov = json.load(open(str(tmp_path / "sc" / "provenance.json")))
    assert prov["results"]["notes"]


def test_no_crash_when_every_male_singlet_lacks_chry(tmp_path):
    A, donors = _two_capture_table(tmp_path)
    A["chrY_counts"] = "0"
    sc, _, _ = S.scorecard(A, ["cons"], {"cons": [("soup", "soup", "soup")]}, T.read_donors(donors), n_boot=10)
    r = sc.set_index("measure").loc["Male-female doublet calls: estimated share with no male cell"]
    assert np.isnan(r.value) and pd.isna(r.numerator) and r.denominator == 2


def test_no_crash_when_one_donor_holds_every_singlet():
    rows = [("s%d" % i, "singlet", "m", "", "male", "1000", "1") for i in range(10)]
    rows += [("d%d" % i, "doublet", "", "f,m", "male", "2000", "1") for i in range(20)]
    out = S.doublet_excess(_cells(rows), ["r"], {}, "r")
    assert np.isnan(out.r_largest_donor_in_pairs.iloc[0]) and out.r_named_pairs.iloc[0] == 20
    assert out.expected_between_donors.iloc[0] == 0.0


def test_missing_bins_take_the_nearest_bin():
    p0 = pd.Series([0.1, 0.2, 0.3], index=[4.0, 8.0, 12.0])
    vals, filled = D.nearest_bin_values(np.array([0.0, 6.0, 7.0, 8.0, 10.0, 20.0]), p0)
    assert vals.tolist() == pytest.approx([0.1, 0.1, 0.2, 0.2, 0.2, 0.3]) and filled == 5
    # a call at 1,000 UMIs (half in bin 2) below the lowest observed bin takes bin 4's p, not the highest bin's
    n, obs, exp, share, filled = D.no_male_share(pd.Series(["m,f"]), np.ones(1, bool), np.ones(1, bool),
                                                 pd.Series([1000.0]), pd.Series([0.0]), p0, ["m"],
                                                 return_filled=True)
    assert exp == pytest.approx(0.1) and filled == 1
    vals, filled = D.nearest_bin_values(np.array([0.0, 8.0]), p0, fill="highest")
    assert vals.tolist() == pytest.approx([0.3, 0.2]) and filled == 1


def test_calibrate_refuses_a_missing_sex_column(tmp_path, run_tools):
    runs, truth, caps, donors = _runs(tmp_path)
    scorer = [["cell", "donor_id", "best_singlet", "margin", "umis"]] + [
        ["c1-%d" % i, "m", "m", 3 * i, 100] for i in range(8)]
    man = _write(tmp_path / "sc.tsv", [["capture", "path"], ["c1", _write(tmp_path / "s1.tsv", scorer)]])
    res = run_tools("calibrate", "--cells", man, "--truth", truth, "--donors", donors, "--sex-column", "sexcall",
                    "--out", tmp_path / "cal.tsv")
    assert res.returncode == 2 and "no column sexcall" in res.stderr


def test_doublet_excess_numbers_through_the_command(tmp_path, run_tools):
    runs, truth, caps, donors = _runs(tmp_path)
    cells = tmp_path / "cells.tsv.gz"
    assert run_tools("cell-table", "--runs", runs, "--truth", truth, "--captures", caps, "--out", cells).returncode == 0
    res = run_tools("doublet-excess", "--cells", cells, "--run", "cons", "--tools", "cons=soup:soup",
                    "--out", tmp_path / "ex.tsv")
    assert res.returncode == 0, res.stderr
    ex = pd.read_csv(str(tmp_path / "ex.tsv"), sep="\t").set_index("capture")
    # c1: 8 cells, singlets m, f, m, f, m (shares 0.6 and 0.4); one doublet, called by souporcell too
    total = 0.008 * 8 / 1000
    assert ex.at["c1", "expected_between_donors"] == pytest.approx(total * (1 - 0.6 ** 2 - 0.4 ** 2), rel=1e-5)
    assert ex.at["c1", "cons_doublet"] == pytest.approx(0.125)
    assert ex.at["c1", "cons_doublet_soup"] == pytest.approx(0.125)
    assert ex.at["c1", "excess"] == pytest.approx(0.125 - total * 0.48, rel=1e-5)
    bad = run_tools("doublet-excess", "--cells", cells, "--run", "cons", "--reference-run", "other",
                    "--out", tmp_path / "ex2.tsv")
    assert bad.returncode == 2 and "--reference-run" in bad.stderr


def test_read_ambient(tmp_path):
    v = _write(tmp_path / "ambient_summary.csv", [["barcode", "ambient_other"], ["a", "0.1"], ["b", "0.3"],
                                                  ["c", "0.9"]], ",")
    s = tmp_path / "ambient_rna.txt"
    s.write_text("ambient RNA estimated as 4.5% of reads\n")
    assert S.read_ambient({"vireo": v, "souporcell": str(s)}, {"a", "b"}) == pytest.approx((0.2, 0.045))
    with pytest.raises(S.UsageError):
        S.read_ambient({"vireo": str(tmp_path / "none.csv")}, {"a"})


def test_lost_cells_columns_sum_to_one(tmp_path):
    A, donors = _two_capture_table(tmp_path)
    _, lost, _ = S.scorecard(A, ["cons"], donors=T.read_donors(donors), n_boot=10)
    for _, g in lost.groupby(["capture_set", "run", "by"]):
        sums = g[[c for c in ("singlet", "recovered", "doublet", "unassigned") if c in g.columns]].sum()
        assert sums.to_numpy() == pytest.approx(np.ones(len(sums)))


def test_donor_recall_flags():
    rows = []
    for cap, nx in (("c1", 5), ("c2", 5), ("c3", 0)):
        t = [("%s-x%d" % (cap, i), "singlet", "x", "", "", "", "") for i in range(nx)]
        t += [("%s-y%d" % (cap, i), "singlet", "y", "", "", "", "") for i in range(10 - nx)]
        rows.append(_cells(t, capture=cap, pool="P"))
    A = pd.concat(rows, ignore_index=True)
    caps = pd.DataFrame({"capture": ["c1", "c2", "c3"], "pool": "P", "capture_set": "main",
                         "donors": [["x", "y"]] * 3, "donors_given": True, "sister_key": "k"})
    _, _, rec = S.scorecard(A, ["r"], captures=caps, n_boot=10)
    flags = rec.set_index(["capture", "donor"]).flag
    assert flags[("c3", "x")] == "r: under a fifth of its sister captures"
    assert flags[("c1", "x")] == "" and flags[("c3", "y")] == ""
    caps["genotyped_r"] = [["x", "y"], ["x", "y"], ["y"]]
    _, _, rec = S.scorecard(A, ["r"], captures=caps, n_boot=10)
    assert rec.set_index(["capture", "donor"]).flag[("c3", "x")] == "no genotype in the r run"
    with pytest.raises(S.UsageError):
        S.scorecard(A, ["r"], captures=caps, n_boot=10, recall_runs=["q"])


def test_bad_tables_are_usage_errors(tmp_path, run_tools):
    runs, truth, caps, donors = _runs(tmp_path)
    cells = tmp_path / "cells.tsv.gz"
    assert run_tools("cell-table", "--runs", runs, "--truth", truth, "--out", cells).returncode == 0
    bad_donors = _write(tmp_path / "bad_donors.tsv", [["name", "sex"], ["m", "male"]])
    res = run_tools("scorecard", "--cells", cells, "--run", "cons", "--donors", bad_donors, "--out", tmp_path / "s1")
    assert res.returncode == 2 and "Traceback" not in res.stderr
    res = run_tools("scorecard", "--cells", cells, "--run", "cons", "--recall-runs", "other", "--captures", caps,
                    "--out", tmp_path / "s2")
    assert res.returncode == 2 and "not scored" in res.stderr
    amb = _write(tmp_path / "amb.tsv", [["capture", "vireo_ambient"], ["c1", str(tmp_path / "missing.csv")]])
    res = run_tools("doublet-excess", "--cells", cells, "--run", "cons", "--ambient", amb, "--out", tmp_path / "e")
    assert res.returncode == 2 and "missing.csv" in res.stderr
