"""demux_tools.calls and demux_tools.compare: call readers, run comparison and the sex check."""
import gzip
import itertools
import json

import pandas as pd
import pytest

from demux_tools import calls as C
from demux_tools import compare as K


def _write(path, rows, sep="\t"):
    with open(str(path), "w") as f:
        f.write("\n".join(sep.join(r) for r in rows) + "\n")
    return str(path)


def _vireo(tmp_path, donors=("a", "b", "c")):
    """donor_ids.tsv and prob_doublet.tsv.gz of three cells, with 0.5.7's permutation names.

    Cell 3's best pair is the combination (b, c), which vireo 0.5.7 writes "b,a".
    """
    d = tmp_path / "vireo"
    d.mkdir()
    perms = [",".join(p) for p in itertools.permutations(donors, 2)]
    _write(d / "donor_ids.tsv", [
        ["cell", "donor_id", "prob_max", "prob_doublet", "n_vars", "best_singlet", "best_doublet",
         "doublet_logLikRatio"],
        ["AAA-1", "a", "9.90e-01", "1.00e-02", "40", "a", "a,b", "-5.000"],
        ["CCC-1", "unassigned", "6.00e-01", "1.00e-02", "3", "b", "a,c", "-1.000"],
        ["GGG-1", "doublet", "1.00e-03", "9.90e-01", "80", "b", perms[2], "4.000"],
    ])
    with open(str(d / "prob_doublet.tsv.gz"), "wb") as f:
        f.write(gzip.compress(("\t".join(["cell"] + perms) + "\n").encode()))
    return d


def test_vireo_reader_corrects_the_pair_and_blanks_non_singlet_donors(tmp_path):
    d = _vireo(tmp_path)
    c = C.read_calls(str(d / "donor_ids.tsv"))   # kind and prob_doublet found
    assert list(c.columns) == C.SCHEMA
    assert c.status.tolist() == ["singlet", "unassigned", "doublet"]
    assert c.donor.tolist() == ["a", "", ""]
    assert c.pair.tolist() == ["a,b", "a,c", "b,c"]   # "b,a" at position 2 is combination (b, c)
    assert c.confidence.tolist() == ["9.90e-01", "6.00e-01", "1.00e-03"]
    plain = C.read_calls(str(d / "donor_ids.tsv"), kind="vireo", prob_doublet="")
    assert plain.pair.tolist()[2] == "a,b"   # without the header the written name is used


def test_vireo_reader_maps_labels(tmp_path):
    d = _vireo(tmp_path)
    lm = _write(tmp_path / "map.csv", [["sample_id", "donor_id"], ["a", "D1"], ["b", "D2"], ["c", "D3"]], ",")
    c = C.read_calls(str(d / "donor_ids.tsv"), label_map=lm)
    assert c.donor.tolist() == ["D1", "", ""]
    assert c.pair.tolist() == ["D1,D2", "D1,D3", "D2,D3"]


def _clusters(tmp_path):
    return _write(tmp_path / "clusters.tsv", [
        ["barcode", "status", "assignment", "log_prob_singleton", "log_prob_doublet", "cluster0", "cluster1"],
        ["AAA-1", "singlet", "0", "-1", "-5", "-1", "-9"],
        ["CCC-1", "singlet", "1", "-1", "-5", "-9", "-1"],
        ["GGG-1", "doublet", "1/0", "-9", "-1", "-5", "-5"],
        ["TTT-1", "unassigned", "0/1", "-3", "-3", "-3", "-3"],
    ])


def test_souporcell_clusters_named_by_map_or_kept_as_numbers(tmp_path):
    path = _clusters(tmp_path)
    raw = C.read_calls(path)
    assert raw.source.iloc[0] == "souporcell"
    assert raw.donor.tolist() == ["0", "1", "", ""]
    assert raw.pair.tolist() == ["", "", "0,1", ""]
    named = C.read_calls(path, cluster_map={"0": "b", "1": "a"})
    assert named.donor.tolist() == ["b", "a", "", ""]
    assert named.pair.tolist() == ["", "", "a,b", ""]
    half = C.name_clusters(raw, {"0": "b"})
    assert half.donor.tolist() == ["b", "", "", ""] and half.pair.tolist() == ["", "", "", ""]


def test_majority_map_uses_cells_singlet_in_both(tmp_path):
    soup = C.read_calls(_clusters(tmp_path))
    ref = pd.DataFrame({"barcode": ["AAA-1", "CCC-1", "GGG-1"], "status": ["singlet", "singlet", "singlet"],
                        "donor": ["x", "y", "z"]})
    assert C.majority_map_from_calls(soup, ref) == {"0": "x", "1": "y"}
    assert C.majority_map(["1", "1", "1", "2", ""], ["a", "b", "b", "", "c"]) == {"1": "b"}


def _hier(tmp_path):
    return _write(tmp_path / "hierarchical_consensus.tsv", [
        ["barcode", "cluster", "donor", "recovery_confidence", "label_confidence", "status"],
        ["AAA-1", "0", "b", "1.0", "1.0", "singlet"],
        ["CCC-1", "1", "a", "0.3", "0.5", "recovered"],
        ["GGG-1", "", "", "0.0", "0.0", "doublet"],
        ["TTT-1", "", "", "0.0", "0.0", "low_cell"],
    ])


def test_hier_reader_keeps_recovered_and_folds_low_cell(tmp_path):
    h = C.read_calls(_hier(tmp_path))
    assert h.source.iloc[0] == "hier"
    assert h.status.tolist() == ["singlet", "recovered", "doublet", "unassigned"]
    assert h.status_raw.tolist()[3] == "low_cell"
    assert h.donor.tolist() == ["b", "a", "", ""]


def _sba(tmp_path):
    return _write(tmp_path / "sample_barcode_assignments.tsv", [
        ["Souporcell_Cluster", "Souporcell_Individual_Assignment", "Barcode", "Souporcell_DropletType",
         "Vireo_Individual_Assignment", "Vireo_DropletType", "AtLeastHalfSinglet_DropletType",
         "AtLeastHalfSinglet_Individual_Assignment"],
        ["0", "S_b", "AAA-1", "singlet", "S_b", "singlet", "singlet", "S_b"],
        ["1", "S_a", "CCC-1", "singlet", "unassigned", "unassigned", "unassigned", "unassigned"],
        ["1/0", "doublet", "GGG-1", "doublet", "doublet", "doublet", "doublet", "doublet"],
        ["0/1", "unassigned", "TTT-1", "unassigned", "S_a", "singlet", "unassigned", "unassigned"],
    ])


def test_assignment_and_demuxafy_readers_with_label_map_and_pattern(tmp_path):
    sba = _sba(tmp_path)
    ca = _write(tmp_path / "cell_assignment.tsv", [
        ["barcode", "assignment", "status"], ["AAA-1", "S_b", "singlet"], ["CCC-1", "NA", "NA"],
        ["GGG-1", "NA", "NA"], ["TTT-1", "X_q_9", "singlet"]])
    lm = _write(tmp_path / "labels.tsv", [["label", "donor"], ["S_a", "a"], ["S_b", "b"]])
    c = C.read_calls(ca, label_map=lm, label_pattern=r"X_(.+)_\d+$")
    assert c.source.iloc[0] == "assignment"
    assert c.status.tolist() == ["singlet", "NA", "NA", "singlet"]   # main writes NA; no sba, kept
    assert c.donor.tolist() == ["b", "", "", "q"]
    split = C.read_calls(ca, kind="assignment", label_map=lm, sba=sba)
    assert split.status.tolist() == ["singlet", "unassigned", "doublet", "singlet"]
    assert split.status_raw.tolist() == ["singlet", "NA", "NA", "singlet"]
    v = C.read_calls(sba, kind="demuxafy:Vireo", label_map=lm)
    assert v.status.tolist() == ["singlet", "unassigned", "doublet", "singlet"] and v.donor.tolist() == ["b", "", "", "a"]
    assert C.cluster_map_from_sba(sba, C.make_mapper(lm)) == {"0": "b", "1": "a"}
    with pytest.raises(ValueError):
        C.read_calls(sba)   # needs demuxafy:TOOL
    with pytest.raises(ValueError):
        C.read_calls(sba, kind="demuxafy:Demuxlet")


def test_generic_and_scorer_kinds(tmp_path):
    g = _write(tmp_path / "g.csv", [["bc", "call", "who"], ["A-1", "singlet", "d1"], ["B-1", "doublet", "d2"]], ",")
    c = C.read_calls(g, kind="generic:barcode=bc,status=call,donor=who")
    assert c.barcode.tolist() == ["A-1", "B-1"] and c.donor.tolist() == ["d1", ""]
    s = _write(tmp_path / "scorer.tsv", [
        ["cell", "donor_id", "prob_max", "prob_doublet", "n_vars", "best_singlet", "best_doublet",
         "doublet_logLikRatio", "best_pair", "margin"],
        ["A-1", "d1", "1", "0", "9", "d1", "d2,d0", "-1", "d0,d2", "12.5"]])
    sc = C.read_calls(s)
    assert sc.source.iloc[0] == "scorer" and sc.pair.tolist() == ["d0,d2"] and sc.confidence.tolist() == ["12.5"]


def test_duplicate_barcodes_are_refused(tmp_path):
    g = _write(tmp_path / "g.tsv", [["barcode", "assignment", "status"], ["A-1", "x", "singlet"], ["A-1", "y", "singlet"]])
    with pytest.raises(ValueError):
        C.read_calls(g)


def _calls(rows):
    return pd.DataFrame([dict(zip(["barcode", "status", "donor", "pair"], r)) for r in rows])


def test_compare_transitions_by_hand():
    a = _calls([("1", "singlet", "a", ""), ("2", "singlet", "b", ""), ("3", "doublet", "", "a,b"),
                ("4", "recovered", "c", ""), ("5", "singlet", "a", ""), ("6", "doublet", "", "a,c")])
    b = _calls([("1", "singlet", "a", ""), ("2", "singlet", "c", ""), ("3", "singlet", "a", ""),
                ("4", "singlet", "c", ""), ("6", "doublet", "", "a,b"), ("7", "unassigned", "", "")])
    cells, trans, summ = K.compare(a, b, "cap")
    t = {(r.status_a, r.status_b): (r.cells, r.same_donor, r.different_donor) for r in trans.itertuples()}
    assert t == {("singlet", "singlet"): (2, 1, 1), ("doublet", "singlet"): (1, 0, 0),
                 ("recovered", "singlet"): (1, 1, 0), ("singlet", "absent"): (1, 0, 0),
                 ("doublet", "doublet"): (1, 0, 0), ("absent", "unassigned"): (1, 0, 0)}
    s = summ.iloc[0]
    assert (s.cells_a, s.cells_b, s.cells_both, s.only_a, s.only_b) == (6, 6, 5, 1, 1)
    # same call: cell 1 only (cell 6 is a doublet in both with another pair)
    assert s.same_call == 1 and s.singlet_donor_changes == 1
    assert abs(s.same_call_pct - 20.0) < 1e-9
    row = cells.set_index("barcode").loc["6"]
    assert (row.same_status, row.same_donor, row.same_call) == ("true", "true", "false")
    assert cells.set_index("barcode").loc["7"].transition == "absent -> unassigned"


def test_sex_check_by_hand():
    calls = _calls([("1", "singlet", "m", ""), ("2", "singlet", "m", ""), ("3", "singlet", "m", ""),
                    ("4", "singlet", "f", ""), ("5", "recovered", "m", ""), ("6", "doublet", "", "m,f")])
    sex = {"1": "male", "2": "female", "3": "ambiguous", "4": "male", "5": "female", "6": "female"}
    out = K.sex_check(calls, sex, K.read_donor_sex(male_donors=["m"]), capture="c", run="r")
    got = {(r.status, r.label_sex): (r.cells, r.with_sex_call, r.discordant) for r in out.itertuples()}
    assert got[("singlet", "male")] == (3, 2, 1)
    assert got[("singlet", "female")] == (1, 1, 1)
    assert got[("recovered", "male")] == (1, 1, 1)
    assert got[("recovered", "female")] == (0, 0, 0)


def test_donor_table_sex_and_unknown(tmp_path):
    t = _write(tmp_path / "donors.tsv", [["donor", "sex"], ["m", "Male"], ["f", "female"]])
    sex_of = K.read_donor_sex(t)
    assert (sex_of("m"), sex_of("f"), sex_of("x"), sex_of("")) == ("male", "female", "unknown", "")


def test_compare_runs_command_with_manifest(tmp_path, run_tools):
    d = _vireo(tmp_path)
    hier = _hier(tmp_path)
    sex = _write(tmp_path / "sex_calls.csv", [["barcode", "sex_call"], ["AAA-1", "female"], ["CCC-1", "male"]], ",")
    man = _write(tmp_path / "manifest.tsv", [
        ["capture", "a_kind", "a_path", "b_kind", "b_path", "sex_calls"],
        ["cap1", "hier", hier, "vireo", str(d / "donor_ids.tsv"), sex]])
    out = tmp_path / "out" / "cmp"
    res = run_tools("compare-runs", "--manifest", man, "--male-donors", "b", "--a-name", "consensus",
                    "--b-name", "vireo", "--out", out)
    assert res.returncode == 0, res.stderr
    summ = pd.read_csv(str(out) + ".summary.tsv", sep="\t")
    assert summ.cells_both.tolist() == [3] and summ.only_a.tolist() == [1]
    sex_t = pd.read_csv(str(out) + ".sex_check.tsv", sep="\t")
    r = sex_t[(sex_t.run == "consensus") & (sex_t.status == "singlet") & (sex_t.label_sex == "male")].iloc[0]
    assert (r.cells, r.with_sex_call, r.discordant) == (1, 1, 1)   # AAA-1 labelled b (male) but female
    prov = json.load(open(str(out) + ".provenance.json"))
    assert prov["results"]["totals"]["cells_both"] == 3
    one = run_tools("compare-runs", "--a", "hier:" + hier, "--b", str(d / "donor_ids.tsv"), "--capture", "x",
                    "--out", tmp_path / "one")
    assert one.returncode == 0, one.stderr


def test_calls_and_sex_check_commands(tmp_path, run_tools):
    hier = _hier(tmp_path)
    res = run_tools("calls", hier, "--out", tmp_path / "c.tsv.gz")
    assert res.returncode == 0, res.stderr
    c = pd.read_csv(str(tmp_path / "c.tsv.gz"), sep="\t", dtype=str, keep_default_na=False)
    assert c.status.tolist() == ["singlet", "recovered", "doublet", "unassigned"]
    same = run_tools("calls", hier, "--out", hier)
    assert same.returncode == 2 and "also an input" in same.stderr
    sex = _write(tmp_path / "sex.tsv", [["barcode", "sex_call"], ["AAA-1", "male"], ["CCC-1", "male"]])
    res = run_tools("sex-check", "--calls", hier, "--sex-calls", sex, "--male-donors", "a", "--run", "r",
                    "--out", tmp_path / "sex_out.tsv")
    assert res.returncode == 0, res.stderr
    t = pd.read_csv(str(tmp_path / "sex_out.tsv"), sep="\t")
    got = {(r.status, r.label_sex): (r.with_sex_call, r.discordant) for r in t.itertuples()}
    assert got[("singlet", "female")] == (1, 1) and got[("recovered", "male")] == (1, 0)


# ---------------------------------------------------------------------------
# Review fixes


def test_assignment_with_sba_keeps_recovered(tmp_path):
    ca = _write(tmp_path / "ca.tsv", [["barcode", "assignment", "status"], ["A-1", "S_a", "singlet"],
                                      ["B-1", "S_b", "recovered"], ["C-1", "NA", "doublet"]])
    sba = _write(tmp_path / "sba.tsv", [["Barcode", "Consensus_DropletType", "AtLeastHalfSinglet_DropletType"],
                                        ["A-1", "singlet", "singlet"], ["B-1", "recovered", "unassigned"],
                                        ["C-1", "doublet", "doublet"]])
    c = C.read_calls(ca, kind="assignment", sba=sba)
    assert c.status.tolist() == ["singlet", "recovered", "doublet"]
    assert c.donor.tolist() == ["S_a", "S_b", ""]


def test_pairs_are_compared_only_when_both_runs_name_one():
    a = _calls([("1", "doublet", "", ""), ("2", "doublet", "", "a,b"), ("3", "doublet", "", "a,b"),
                ("4", "singlet", "a", "a,b")])
    b = _calls([("1", "doublet", "", "a,b"), ("2", "doublet", "", "a,b"), ("3", "doublet", "", "a,c"),
                ("4", "singlet", "a", "")])
    cells, _, summ = K.compare(a, b, "cap")
    c = cells.set_index("barcode")
    assert c.same_pair.tolist() == ["", "true", "false", ""]
    assert c.same_call.tolist() == ["true", "true", "false", "true"]
    s = summ.iloc[0]
    assert (s.doublets_both, s.doublet_pairs_compared, s.doublet_pairs_same, s.same_call) == (3, 2, 1, 3)


def test_cells_are_sorted_by_barcode():
    a = _calls([("C", "singlet", "a", ""), ("A", "singlet", "a", "")])
    b = _calls([("B", "singlet", "a", ""), ("A", "singlet", "a", "")])
    cells, _, _ = K.compare(a, b)
    assert cells.barcode.tolist() == ["A", "B", "C"]


def test_sex_calls_are_compared_in_lower_case(tmp_path):
    sex = K.read_sex_calls(_write(tmp_path / "s.csv", [["barcode", "sex_call"], ["1", "Female"], ["2", " MALE "]], ","))
    assert sex == {"1": "female", "2": "male"}
    calls = _calls([("1", "singlet", "m", ""), ("2", "singlet", "m", "")])
    out = K.sex_check(calls, {"1": "Female", "2": "male"}, K.read_donor_sex(male_donors=["m"]))
    r = out[(out.status == "singlet") & (out.label_sex == "male")].iloc[0]
    assert (r.with_sex_call, r.discordant) == (2, 1)


def test_inputs_include_files_found_or_named_indirectly(tmp_path, run_tools):
    d = _vireo(tmp_path)
    ids = str(d / "donor_ids.tsv")
    assert C.files_read(ids) == [ids, str(d / "prob_doublet.tsv.gz")]
    clusters = _clusters(tmp_path)
    cmap = _write(tmp_path / "cmap.tsv", [["cluster", "donor"], ["0", "b"], ["1", "a"]])
    sba = _sba(tmp_path)
    assert C.files_read(clusters, cluster_map=cmap) == [clusters, cmap]
    assert C.files_read(clusters, cluster_map="sba:" + sba) == [clusters, sba]
    res = run_tools("calls", clusters, "--cluster-map", cmap, "--out", cmap)
    assert res.returncode == 2 and "also an input" in res.stderr
    res = run_tools("calls", ids, "--out", tmp_path / "v.tsv")
    assert res.returncode == 0, res.stderr
    prov = json.load(open(str(tmp_path / "v.tsv.provenance.json")))
    assert str(d / "prob_doublet.tsv.gz") in [x["path"] for x in prov["inputs"]]


def test_bad_tables_are_usage_errors(tmp_path, run_tools):
    hier = _hier(tmp_path)
    donors = _write(tmp_path / "donors.tsv", [["donor", "gender"], ["a", "male"]])
    sex = _write(tmp_path / "sex.tsv", [["barcode", "sex"], ["AAA-1", "male"]])
    res = run_tools("compare-runs", "--a", hier, "--b", hier, "--sex-calls", sex, "--donors", donors,
                    "--out", tmp_path / "o")
    assert res.returncode == 2 and "Traceback" not in res.stderr
    res = run_tools("compare-runs", "--a", hier, "--b", hier, "--sex-calls", sex, "--male-donors", "a",
                    "--out", tmp_path / "o2")
    assert res.returncode == 2 and "sex_call" in res.stderr and "Traceback" not in res.stderr
    res = run_tools("calls", hier, "--label-map", donors, "--out", tmp_path / "c.tsv")
    assert res.returncode == 2 and "Traceback" not in res.stderr
