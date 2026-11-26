"""demux_tools.vireo_port: the ported vireoSNP functions give vireoSNP's numbers."""
import gzip
import itertools
import os

import numpy as np
import pytest

from demux_tools import vireo_port as vp


def _random_ids(rng, n, pool):
    return ["chr%d_%d_A_G" % (rng.randint(1, 3), rng.randint(0, pool)) for _ in range(n)]


def test_match_pairs_the_kth_copies():
    ref = ["b", "a", "b", "c", "b"]
    new = ["b", "x", "b", "a"]
    out = vp.match(ref, new)
    # sorted ref: a(1) b(0) b(2) b(4) c(3); sorted new: a(3) b(0) b(2) x(1)
    assert out.tolist() == [0, 3, 2, -1, -1]
    assert vp.match(ref, new, uniq_ref_only=False).tolist() == [0, 3, 0, -1, 0]
    assert vp.match([], new).tolist() == [] and vp.match(ref, []).tolist() == [-1] * 5


@pytest.mark.vireo
@pytest.mark.parametrize("seed", [0, 1, 2])
def test_match_equals_vireo(seed):
    from vireoSNP.utils.vireo_base import match

    rng = np.random.RandomState(seed)
    ref = _random_ids(rng, 400, 150)            # many repeats, beyond insertion-sort sizes
    new = _random_ids(rng, 120, 150)
    for uniq in (True, False):
        theirs = np.array([-1 if x is None else int(x) for x in match(ref, new, uniq_ref_only=uniq)])
        assert vp.match(ref, new, uniq_ref_only=uniq).tolist() == theirs.tolist()


@pytest.mark.vireo
def test_match_variants_equals_match_donor_vcf():
    from vireoSNP.utils import io_utils
    from vireoSNP.utils.io_utils import match_donor_VCF

    rng = np.random.RandomState(3)
    cases = [(["1_%d_A_G" % i for i in rng.randint(0, 50, 80)], ["1_%d_A_G" % i for i in range(0, 50, 3)]),
             (["1_%d_A_G" % i for i in range(10, 20)], ["chr1_%d_A_G" % i for i in range(5, 30)])]
    if hasattr(io_utils, "match_SNPs"):     # 0.5.9; 0.5.7's "chr" retries stop with a TypeError on these
        cases += [(["1_%d_A_G" % i for i in range(20)], ["chr1_%d_A_G" % i for i in range(5, 30)]),
                  (["chr1_%d_A_G" % i for i in range(20)], ["1_%d_A_G" % i for i in range(5, 30)])]
    for cells, donors in cases:
        cell_dat = {"variants": list(cells), "AD": np.arange(len(cells))[:, None], "DP": np.arange(len(cells))[:, None],
                    "FixedINFO": {}}
        donor_vcf = {"variants": list(donors), "FixedINFO": {}, "GenoINFO": {"GT": list(range(len(donors)))}}
        cell_dat, donor_vcf = match_donor_VCF(cell_dat, donor_vcf)
        idx_cell, idx_donor, _ = vp.match_variants(cells, donors)
        assert idx_cell.tolist() == cell_dat["AD"].ravel().tolist()
        assert idx_donor.tolist() == donor_vcf["GenoINFO"]["GT"]


def test_genotype_probabilities_is_parse_gt_code_per_value():
    codes = np.array([["0/1", "./."], ["1|1", "0"], ["1", "."]], dtype=object)
    probs, missing = vp.genotype_probabilities(codes, "GT")
    for i in range(3):
        for j in range(2):
            assert np.array_equal(probs[i, j], vp.parse_gt_code(codes[i, j], "GT"))
    assert missing.tolist() == [[False, True], [False, False], [False, True]]
    gp = np.array([["0.1,0.2,0.7", "1e-3,0.5,0.499"], ["0.3,0.3,0.3", "."]], dtype=object)
    probs, _ = vp.genotype_probabilities(gp, "GP")
    for i in range(2):
        for j in range(2):
            assert np.array_equal(probs[i, j], vp.parse_gt_code(gp[i, j], "GP"))


@pytest.mark.vireo
def test_add_doublet_gt_equals_vireo_in_any_chunks():
    from vireoSNP.utils.vireo_doublet import add_doublet_GT

    rng = np.random.RandomState(5)
    G = vp.normalize(rng.rand(40, 5, 3))
    theirs = add_doublet_GT(G)[:, 5:, :]
    assert np.array_equal(vp.add_doublet_gt(G), theirs)
    for chunk in (1, 3, 7):
        parts = [vp.add_doublet_gt(G, np.arange(s, min(10, s + chunk))) for s in range(0, 10, chunk)]
        assert np.array_equal(np.concatenate(parts, axis=1), theirs)


@pytest.mark.scipy
def test_posteriors_sum_to_one_and_use_vireos_prior():
    rng = np.random.RandomState(6)
    ll_s = -rng.rand(50, 4) * 30
    ll_d = -rng.rand(50, 6) * 30
    p_s, p_d, llr = vp.posteriors(ll_s, ll_d, 0.2)
    np.testing.assert_allclose(p_s.sum(1) + p_d.sum(1), 1.0, atol=1e-12)
    np.testing.assert_allclose(llr, ll_d.max(1) - ll_s.max(1))
    # equal log-likelihoods: the prior alone, (1 - rate) / 4 per donor and rate / 6 per pair
    p_s, p_d, _ = vp.posteriors(np.zeros((1, 4)), np.zeros((1, 6)), 0.2)
    np.testing.assert_allclose(p_s, 0.8 / 4)
    np.testing.assert_allclose(p_d, 0.2 / 6)


def test_write_vireo_files_names_pairs_by_version(tmp_path):
    donors = ["a", "b", "c"]
    ps = np.array([[0.05, 0.05, 0.9], [0.3, 0.3, 0.3]])
    pd_ = np.array([[0.0, 0.0, 0.0], [0.0, 0.0, 0.1]])
    llr = np.array([-3.0, 1.25])
    for naming, col, best in (("permutations", "prob_doublet", "b,a"), ("combinations", "prob_max_doublet", "b,c")):
        out = vp.write_vireo_files(str(tmp_path / naming), donors, ["x-1", "y-1"], np.array([20, 5]), ps, pd_, llr,
                                   naming=naming, prob_doublet_column=col)
        lines = open(out["donor_ids"]).read().splitlines()
        assert lines[0].split("\t")[3] == col
        assert lines[1] == "x-1\tc\t9.00e-01\t0.00e+00\t20\tc\ta,b\t-3.000"
        assert lines[2] == "y-1\tunassigned\t3.00e-01\t1.00e-01\t5\ta\t%s\t1.250" % best
        header = gzip.open(out["prob_doublet"], "rt").readline().rstrip("\n").split("\t")
        gen = itertools.permutations if naming == "permutations" else itertools.combinations
        assert header[1:] == [",".join(p) for p in gen(donors, 2)]
        assert open(out["summary"]).read() == "Var1\tFreq\nc\t1\nunassigned\t1\n"
        assert os.path.getsize(out["prob_singlet"]) > 0
