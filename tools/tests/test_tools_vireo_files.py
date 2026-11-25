"""demux_tools.vireo_files: pair naming, true-pair recovery and the vireo-pairs command."""
import gzip
import itertools
import json
import os

import numpy as np
import pytest

from demux_tools import vireo_files as vf


def _gz(path, lines):
    with open(str(path), "wb") as f:
        f.write(gzip.compress(("\n".join(lines) + "\n").encode()))
    return str(path)


def _vireo_dir(root, donors, cells, naming="permutations", ids_column="prob_doublet"):
    """Write donor_ids.tsv, prob_singlet.tsv.gz and prob_doublet.tsv.gz as vireo would.

    cells: list of (barcode, singlet probabilities, pair probabilities in
    combinations order). best_doublet is the `naming` name at the argmax.
    """
    root.mkdir(parents=True, exist_ok=True)
    gen = itertools.permutations if naming == "permutations" else itertools.combinations
    names = [",".join(p) for p in gen(donors, 2)]
    rows = ["\t".join(["cell", "donor_id", "prob_max", ids_column, "n_vars", "best_singlet",
                       "best_doublet", "doublet_logLikRatio"])]
    singlet = ["\t".join(["cell"] + list(donors))]
    doublet = ["\t".join(["cell"] + names)]
    for barcode, ps, pd_ in cells:
        ps, pd_ = np.asarray(ps, dtype=float), np.asarray(pd_, dtype=float)
        rows.append("\t".join([barcode, donors[int(ps.argmax())], "%.2e" % ps.max(), "%.2e" % pd_.max(), "50",
                               donors[int(ps.argmax())], names[int(pd_.argmax())], "%.3f" % -3.0]))
        singlet.append("\t".join([barcode] + ["%.2e" % x for x in ps]))
        doublet.append("\t".join([barcode] + ["%.2e" % x for x in pd_]))
    with open(str(root / "donor_ids.tsv"), "w") as f:
        f.write("\n".join(rows) + "\n")
    _gz(root / "prob_singlet.tsv.gz", singlet)
    _gz(root / "prob_doublet.tsv.gz", doublet)
    return root


def _header_file(tmp_path, names):
    return _gz(tmp_path / "prob_doublet.tsv.gz", ["\t".join(["cell"] + names)])


def test_name_written_at_combination_k_is_the_kth_permutation(tmp_path):
    donors = ["a", "b", "c"]
    perms = [",".join(p) for p in itertools.permutations(donors, 2)]
    combos = [",".join(p) for p in itertools.combinations(donors, 2)]
    path = _header_file(tmp_path, perms)
    naming = vf.pair_naming(path)
    assert naming["kind"] == "permutations"
    assert naming["donors"] == donors
    assert naming["combos"] == combos
    # vireo writes perms[k] for the k-th combination
    assert vf.true_pairs(path, [perms[k] for k in range(len(combos))]) == combos
    assert vf.true_pairs(path, ["b,a", "a,c", "", "zz"]) == ["b,c", "a,c", "", ""]


def test_pairs_come_back_in_string_order_for_any_donor_order(tmp_path):
    donors = ["z", "a", "m"]   # not sorted: vireo keeps the VCF's order
    path = _header_file(tmp_path, [",".join(p) for p in itertools.permutations(donors, 2)])
    assert vf.pair_naming(path)["combos"] == ["a,z", "m,z", "a,m"]
    assert vf.true_pairs(path, ["a,z"]) == ["a,m"]   # name 2 is combination 2: (a, m)


def test_combination_headers_pass_through(tmp_path):
    donors = ["d1", "d0", "d2", "d3"]
    names = [",".join(p) for p in itertools.combinations(donors, 2)]
    path = _header_file(tmp_path, names)   # 6 names: also n(n-1) for n = 3
    naming = vf.pair_naming(path)
    assert naming["kind"] == "combinations" and naming["donors"] == donors
    assert vf.true_pairs(path, names) == [vf.pair_key(*n.split(",")) for n in names]


def test_six_permutation_names_are_three_donors(tmp_path):
    path = _header_file(tmp_path, [",".join(p) for p in itertools.permutations(["x", "y", "w"], 2)])
    assert vf.pair_naming(path)["kind"] == "permutations"


def test_two_donors(tmp_path):
    assert vf.pair_naming(_header_file(tmp_path, ["p,q", "q,p"]))["combos"] == ["p,q"]
    os.remove(str(tmp_path / "prob_doublet.tsv.gz"))
    assert vf.pair_naming(_header_file(tmp_path, ["q,p"]))["kind"] == "combinations"


def test_donor_names_with_commas_need_the_donor_list(tmp_path):
    donors = ["a,1", "b", "c"]
    path = _header_file(tmp_path, [",".join(p) for p in itertools.permutations(donors, 2)])
    with pytest.raises(ValueError):
        vf.pair_naming(path)
    assert vf.pair_naming(path, donors=donors)["combos"] == ["a,1,b", "a,1,c", "b,c"]


def test_bad_headers_raise(tmp_path):
    with pytest.raises(ValueError):
        vf.pair_naming(_header_file(tmp_path, ["a,b", "b,c", "c,a", "x,y", "y,z"]))
    path = _header_file(tmp_path, ["a,b", "a,c", "b,a", "b,c", "c,a", "c,b"])
    with pytest.raises(ValueError, match="beyond"):
        vf.true_pairs(path, ["c,b"])   # position 5, but there are 3 value columns


def test_read_donor_ids_accepts_prob_max_doublet(tmp_path):
    d = _vireo_dir(tmp_path / "v", ["a", "b", "c"], [("AAAC-1", [0.9, 0.05, 0.05], [0.1, 0.0, 0.2])],
                   ids_column="prob_max_doublet")
    df = vf.read_donor_ids(str(d / "donor_ids.tsv"))
    assert list(df.columns) == vf.DONOR_IDS_COLUMNS
    assert df["prob_doublet"].tolist() == ["2.00e-01"]


def _cells():
    # singlet probabilities over (a, b, c); pair probabilities over (a,b), (a,c), (b,c)
    return [
        ("AAAC-1", [0.98, 0.01, 0.01], [0.001, 0.002, 0.95]),   # best pair b,c: written "b,a"
        ("AAAG-1", [0.2, 0.7, 0.1], [0.6, 0.1, 0.05]),          # best pair a,b: right name
        ("AAAT-1", [0.5, 0.25, 0.25], [0.0, 0.0, 0.0]),         # every pair probability 0
        ("AACA-1", [0.4, 0.3, 0.3], [0.3, 0.3, 0.1]),           # tie in the written values
    ]


def test_vireo_pairs_table(tmp_path):
    d = _vireo_dir(tmp_path / "v", ["a", "b", "c"], _cells())
    table, summary = vf.vireo_pairs(str(d / "donor_ids.tsv"), str(d / "prob_doublet.tsv.gz"),
                                    str(d / "prob_singlet.tsv.gz"))
    assert list(table.columns) == vf.PAIRS_COLUMNS
    assert table["best_doublet"].tolist() == ["b,a", "a,b", "a,b", "a,b"]
    assert table["best_pair"].tolist() == ["b,c", "a,b", "", "a,b"]
    assert table["pair_prob_underflow"].tolist() == ["false", "false", "true", "false"]
    assert set(table["vireo_naming"]) == {"permutations"}
    assert summary == {"cells": 4, "underflow": 1, "naming": "permutations", "donors": ["a", "b", "c"],
                       "renamed": 1, "matrix_agree": 2, "matrix_tied": 1, "matrix_disagree": 0}


def test_best_pairs_from_matrix(tmp_path):
    d = _vireo_dir(tmp_path / "v", ["a", "b", "c"], _cells())
    m = vf.best_pairs_from_matrix(str(d / "prob_doublet.tsv.gz"))
    assert m["matrix_pair"].tolist() == ["b,c", "a,b", "", "a,b"]
    assert m["underflow"].tolist() == [False, False, True, False]
    assert m["tied"].tolist() == [False, False, False, True]
    assert m["max_prob"].tolist() == [0.95, 0.6, 0.0, 0.3]


def test_vireo_pairs_command(tmp_path, run_tools):
    d = _vireo_dir(tmp_path / "v", ["a", "b", "c"], _cells())
    before = sorted(os.listdir(str(d)))
    out = tmp_path / "out" / "pairs.tsv.gz"   # the folder is created
    res = run_tools("vireo-pairs", "--vireo-dir", d, "--out", out)
    assert res.returncode == 0, res.stderr
    assert "1 best_doublet names corrected" in res.stdout
    assert sorted(os.listdir(str(d))) == before   # nothing written beside the inputs
    from demux_tools._io import read_table

    t = read_table(str(out))
    assert t["best_pair"].tolist() == ["b,c", "a,b", "", "a,b"]
    rec = json.loads(open(str(out) + ".provenance.json").read())
    assert rec["results"]["matrix_disagree"] == 0
    assert len(rec["inputs"]) == 3
    # the files can be named one by one; refusing to overwrite an input
    res = run_tools("vireo-pairs", "--donor-ids", d / "donor_ids.tsv", "--prob-doublet",
                    d / "prob_doublet.tsv.gz", "--out", d / "donor_ids.tsv")
    assert res.returncode == 2 and "also an input" in res.stderr
    res = run_tools("vireo-pairs", "--donor-ids", d / "donor_ids.tsv", "--out", tmp_path / "x.tsv")
    assert res.returncode == 2


def test_vireo_pairs_refuses_a_wrong_naming_rule(tmp_path):
    # best_doublet names a pair that the matrix contradicts: the rule would not hold
    d = _vireo_dir(tmp_path / "v", ["a", "b", "c"], [("AAAC-1", [0.98, 0.01, 0.01], [0.001, 0.002, 0.95])])
    text = open(str(d / "donor_ids.tsv")).read().replace("b,a", "a,c")
    open(str(d / "donor_ids.tsv"), "w").write(text)
    table, summary = vf.vireo_pairs(str(d / "donor_ids.tsv"), str(d / "prob_doublet.tsv.gz"))
    assert summary["matrix_disagree"] == 1


@pytest.mark.vireo
def test_true_pairs_on_files_written_by_vireo(tmp_path):
    """vireoSNP's own writer: 0.5.7 names pairs by permutations, 0.5.9 by combinations.
    Either way best_pair is the argmax of the pair probabilities."""
    from vireoSNP.utils.io_utils import write_donor_id

    rng = np.random.RandomState(7)
    donors = ["d3", "d1", "d0", "d2", "d4"]   # vireo's order need not be sorted
    n_cells, n_pairs = 300, 10
    ll_s = rng.normal(size=(n_cells, 5)) * 5
    ll_d = rng.normal(size=(n_cells, n_pairs)) * 5
    ll_d[:20] -= 900   # pair probabilities that underflow to 0
    id_prob = np.exp(ll_s - ll_s.max(1, keepdims=True))
    id_prob /= id_prob.sum(1, keepdims=True)
    doublet_prob = np.exp(ll_d - np.maximum(ll_s.max(1), ll_d.max(1))[:, None] - 1)
    doublet_prob[:20] = 0.0
    res = {"ID_prob": id_prob, "doublet_prob": doublet_prob, "LB_doublet": 0.0, "theta_shapes": None,
           "doublet_LLR": rng.normal(size=n_cells), "ambient_Psi": None}
    out = tmp_path / "vireo"
    out.mkdir()
    cells = ["C%03d-1" % k for k in range(n_cells)]
    cwd = os.getcwd()
    os.chdir(str(out))   # vireo gzips through a shell command
    try:
        write_donor_id(str(out), donors, cells, np.full(n_cells, 50), res)
    finally:
        os.chdir(cwd)
    table, summary = vf.vireo_pairs(str(out / "donor_ids.tsv"), str(out / "prob_doublet.tsv.gz"),
                                    str(out / "prob_singlet.tsv.gz"))
    combos = [vf.pair_key(a, b) for a, b in itertools.combinations(donors, 2)]
    expected = [combos[k] for k in doublet_prob.argmax(1)]
    keep = doublet_prob.max(1) > 0
    assert summary["underflow"] == 20
    assert [p for p, k in zip(table["best_pair"], keep) if k] == [p for p, k in zip(expected, keep) if k]
    assert set(table.loc[~keep, "best_pair"]) == {""}
    assert summary["matrix_disagree"] == 0
    import vireoSNP

    old = tuple(int(x) for x in vireoSNP.__version__.split(".")[:3]) < (0, 5, 9)
    assert summary["naming"] == ("permutations" if old else "combinations")
