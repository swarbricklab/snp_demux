"""demux_tools.scorer: vireo's fixed-mode arithmetic, the commands, and identity with vireo itself."""
import gzip
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from demux_tools import scorer as S
from demux_tools import synth
from demux_tools.cellsnp import Pileup

pytestmark = pytest.mark.scipy

HERE = os.path.dirname(os.path.abspath(__file__))
THETA = ([[0.01, 0.5, 0.99]], [[2000.0, 2000.0, 2000.0]])


def _given_fit(capture, theta=THETA):
    return S.Fit(theta[0], theta[1], doublet_prior=S.default_doublet_prior(capture.n_cells))


def _load(cap, **kw):
    return S.load_capture("t", cellsnp=cap["cellsnp"], vcf=cap["vcf"], **kw)


def _vireo_version():
    import vireoSNP

    return vireoSNP.__version__


def _run_vireo(cellsnp, vcf, out_dir, tmp_path, tag="GT"):
    """vireoSNP's own command line in fixed mode, one process, no plots."""
    code = ("import sys; from vireoSNP.vireo import main; sys.argv = ['vireo'] + sys.argv[1:]; main()")
    env = dict(os.environ, MPLCONFIGDIR=str(tmp_path / "mpl"))
    res = subprocess.run([sys.executable, "-c", code, "--cellData", cellsnp, "--donorFile", vcf, "--outDir",
                          out_dir, "--genoTag", tag, "--nproc", "1", "--noPlot"], stdout=subprocess.PIPE,
                         stderr=subprocess.PIPE, universal_newlines=True, env=env)
    assert res.returncode == 0, res.stdout + res.stderr
    return out_dir


# ---------------------------------------------------------------------------
# Arithmetic


def test_loglik_singlet_is_the_sum_of_site_terms(synthetic_capture):
    cap = synthetic_capture(n_donors=4, n_sites=120, n_cells=12, blank_rate=0.2, tags=("GT",), seed=1)
    capture = _load(cap)
    fit = _given_fit(capture, ([[0.03, 0.45, 0.97]], [[500.0, 80.0, 300.0]]))
    ll = S.loglik_singlet(capture.AD, capture.DP, capture.GPb, fit.theta_mean, fit.theta_sum)
    for j in (0, 5, 12):
        t = S.cell_site_terms(capture, fit, capture.barcodes[j])
        sums = np.array([t["term__" + d].sum() for d in capture.donors])
        np.testing.assert_allclose(sums, ll[j], rtol=0, atol=1e-10)
        assert (t["alt"] + t["ref"] > 0).all()


def test_doublet_chunks_give_identical_values(synthetic_capture):
    cap = synthetic_capture(n_donors=6, n_sites=80, n_cells=20, blank_rate=0.1, tags=("GP",), gp_noise=0.2, seed=2)
    capture = _load(cap, geno_tag="GP")
    args = (capture.AD, capture.DP, capture.GPb, THETA[0], THETA[1])
    whole = S.loglik_doublet(*args, pair_chunk=1000)
    assert whole.shape == (21, 15)
    for chunk in (1, 2, 4):
        assert np.array_equal(S.loglik_doublet(*args, pair_chunk=chunk), whole)


def test_posteriors_and_doublet_prior(synthetic_capture):
    capture = _load(synthetic_capture(n_donors=3, n_sites=50, n_cells=9, seed=3))
    fit = _given_fit(capture)
    assert fit.doublet_prior == pytest.approx(10 / 100000)
    assert S.default_doublet_prior(20000) == 0.2 and S.default_doublet_prior(200000) == 0.5
    sc = S.score(capture, fit)
    np.testing.assert_allclose(sc.prob_singlet.sum(1) + sc.prob_doublet.sum(1), 1.0, atol=1e-12)
    assert sc.source == "scorer" and sc.check == {}


def test_recovers_truth_and_flags_a_50_50_doublet():
    sc_ = synth.make_capture(n_donors=4, n_sites=300, n_cells=200, seed=7, doublets=((0, 1, 0.5), (2, 3, 0.5)))
    import tempfile

    with tempfile.TemporaryDirectory() as d:
        paths = synth.write_capture(sc_, d)
        capture = S.load_capture("t", cellsnp=paths["cellsnp"], vcf=paths["vcf"])
        res = S.score(capture, _given_fit(capture))
    t = res.table
    truth = sc_.truth
    singlet = (truth["kind"] == "singlet").to_numpy() & (t["n_vars"].to_numpy() >= 10)
    assert singlet.sum() > 150
    assert (t["best_singlet"].to_numpy()[singlet] == truth["donor"].to_numpy()[singlet]).all()
    assert (t["donor_id"].to_numpy()[singlet] == truth["donor"].to_numpy()[singlet]).all()
    doublets = (truth["kind"] == "doublet").to_numpy()
    assert (t["doublet_margin"].to_numpy()[doublets] < 0).all()
    assert t["best_pair"].tolist()[-2:] == ["d0,d1", "d2,d3"]
    assert (t["margin"].to_numpy() >= 0).all()


def _capture_with_an_empty_cell(tmp_path, n_donors):
    sc_ = synth.make_capture(n_donors=n_donors, n_sites=40, n_cells=6, seed=8, doublets=())
    p = sc_.pileup
    layers = {}
    for name, (r, c, v) in p.layers.items():
        keep = c != 0
        layers[name] = (r[keep], c[keep], v[keep])
    sc_ = sc_._replace(pileup=Pileup(p.chrom, p.pos, p.ref, p.alt, p.barcodes, layers=layers))
    return synth.write_capture(sc_, str(tmp_path / ("d%d" % n_donors)))


@pytest.mark.parametrize("n_donors", [1, 2])
def test_one_and_two_donors_and_a_cell_without_reads(tmp_path, n_donors):
    paths = _capture_with_an_empty_cell(tmp_path, n_donors)
    capture = S.load_capture("t", cellsnp=paths["cellsnp"], vcf=paths["vcf"])
    res = S.score(capture, _given_fit(capture))
    t = res.table
    assert t.loc[0, "n_vars"] == 0 and t.loc[0, "donor_id"] == "unassigned" and np.isnan(t.loc[0, "fit_per_umi"])
    assert t.loc[0, "ll_best"] == 0.0
    assert res.ll_doublet.shape[1] == n_donors * (n_donors - 1) // 2
    if n_donors == 1:
        assert (t["prob_max"] == 1.0).all() and (t["prob_doublet"] == 0.0).all()
        assert set(t["best_pair"]) == {""} and t["ll_second"].isna().all()
    else:
        assert set(t["best_pair"]) == {"d0,d1"}
    out = S.write_cells(res, str(tmp_path / "cells.tsv.gz"))
    back = pd.read_csv(out, sep="\t", dtype=str, keep_default_na=False)
    assert list(back.columns) == S.CELL_COLUMNS and back.loc[0, "fit_per_umi"] == ""


def test_genotype_text():
    gp = np.array([[1.0, 0, 0], [0, 0, 1.0], [1 / 3] * 3, [0.1, 0.8, 0.1]])
    assert S.genotype_text(gp, [False, False, True, False]) == ["0/0", "1/1", "./.", "0.1,0.8,0.1"]


# ---------------------------------------------------------------------------
# Loading and the cache


def test_cache_hit_reads_no_matrix_and_mtime_changes_the_key(synthetic_capture, tmp_path, monkeypatch):
    cap = synthetic_capture(n_donors=3, n_sites=60, n_cells=10, blank_rate=0.1, seed=9)
    first = _load(cap, cache_dir=str(tmp_path / "cache"), info=("AF",))
    assert first.load_info["cache"] == "miss"

    def refuse(*a, **k):
        raise AssertionError("read_mtx called on a cache hit")

    monkeypatch.setattr(S, "read_mtx", refuse)
    second = _load(cap, cache_dir=str(tmp_path / "cache"), info=("AF",))
    assert second.load_info["cache"] == "hit"
    for name in ("AD", "DP"):
        a, b = getattr(first, name), getattr(second, name)
        assert np.array_equal(a.indptr, b.indptr) and np.array_equal(a.indices, b.indices)
        assert np.array_equal(a.data, b.data) and a.data.dtype == b.data.dtype == np.int64
    assert np.array_equal(first.GPb, second.GPb) and np.array_equal(first.missing, second.missing)
    assert first.variants == second.variants and first.barcodes == second.barcodes
    assert np.array_equal(first.info["AF"], second.info["AF"], equal_nan=True)
    assert first.counts == second.counts
    st = os.stat(cap["vcf"])
    os.utime(cap["vcf"], ns=(st.st_atime_ns, st.st_mtime_ns + 10 ** 9))
    with pytest.raises(AssertionError, match="read_mtx called"):
        _load(cap, cache_dir=str(tmp_path / "cache"), info=("AF",))


def test_sites_outside_the_panel_and_barcode_subset(synthetic_capture, tmp_path):
    cap = synthetic_capture(n_donors=3, n_sites=50, n_cells=10, vcf_rows=list(range(0, 50, 2)), seed=10)
    capture = _load(cap)
    assert capture.n_sites == 25 and capture.counts["sites_dropped"]["not_in_panel"] == 25
    assert capture.pileup_rows.tolist() == list(range(0, 50, 2))
    keep = capture.barcodes[3:6]
    (tmp_path / "bc.txt").write_text("\n".join(reversed(keep)) + "\n")
    sub = _load(cap, barcodes=str(tmp_path / "bc.txt"))
    assert sub.barcodes == keep
    assert np.array_equal(sub.DP.toarray(), capture.DP.toarray()[:, 3:6])


def test_no_matching_site_and_unknown_donor(synthetic_capture, tmp_path):
    cap = synthetic_capture(n_donors=2, n_sites=20, n_cells=4, seed=11)
    other = synthetic_capture("other", n_donors=2, n_sites=20, n_cells=4, seed=12)
    with pytest.raises(ValueError, match="no pileup site matches"):
        S.load_capture("t", cellsnp=cap["cellsnp"], vcf=other["vcf"])
    with pytest.raises(S.UsageError, match="not in the VCF"):
        _load(cap, donors=["d0", "zz"])


def test_donor_orders(synthetic_capture, tmp_path):
    cap = synthetic_capture(n_donors=3, n_sites=30, n_cells=6, seed=13, donor_names=["m", "z", "a"])
    assert _load(cap).donors == ["m", "z", "a"]
    assert _load(cap, donor_order="sorted").donors == ["a", "m", "z"]
    assert _load(cap, donors=["z", "a"], donor_order="list").donors == ["z", "a"]
    assert _load(cap, donors=["a", "z"]).donors == ["z", "a"]          # the VCF's order, as vireo
    header = tmp_path / "prob_singlet.tsv.gz"
    header.write_bytes(gzip.compress(b"cell\ta\tz\tm\n"))
    loaded = _load(cap, donor_order="prob-header:" + str(header))
    assert loaded.donors == ["a", "z", "m"]
    order = [["m", "z", "a"].index(d) for d in loaded.donors]
    assert np.array_equal(loaded.GPb, _load(cap).GPb[:, order, :])


def _vcf_with_repeats(sc_, path, n_repeats=60, n_unpiled=3000, seed=1):
    """The capture's donor VCF plus a second copy of `n_repeats` piled records with every genotype changed and
    `n_unpiled` records at sites with no reads, sorted by position (the review's experiment A)."""
    tmp = path + ".plain.vcf"
    synth.write_vcf(sc_, tmp)
    lines = open(tmp).read().splitlines()
    os.remove(tmp)
    header = [x for x in lines if x.startswith("#")]
    body = [x.split("\t") for x in lines if not x.startswith("#")]
    flip = {"0/0": "0/1", "0/1": "1/1", "1/1": "0/0", "./.": "./."}
    rng = np.random.RandomState(seed)
    extra = []
    for k in sorted(rng.choice(len(body), n_repeats, replace=False)):
        r = list(body[k])
        r[9:] = [flip[x] for x in r[9:]]
        extra.append(r)
    for i in range(n_unpiled):
        extra.append(["chr3", str(10000 + 10 * i), ".", "A", "G", ".", "PASS", "AF=0.5", "GT"] + ["0/1"] * len(sc_.donors))
    records = body + extra
    order = sorted(range(len(records)), key=lambda k: (records[k][0], int(records[k][1]), k))
    with gzip.open(path, "wt") as f:
        f.write("\n".join(header) + "\n")
        for k in order:
            f.write("\t".join(records[k]) + "\n")
    return path


def test_repeated_donor_records_follow_vireos_match_on_every_record(tmp_path):
    from demux_tools import vcf as vcf_mod
    from demux_tools.vireo_port import match

    sc_ = synth.make_capture(n_donors=4, n_sites=400, n_cells=40, seed=7, depth=(1, 3))
    paths = synth.write_capture(sc_, str(tmp_path / "cap"))
    vcf = _vcf_with_repeats(sc_, str(tmp_path / "dup.vcf.gz"))
    capture = S.load_capture("t", cellsnp=paths["cellsnp"], vcf=vcf)
    assert capture.counts["vcf_ids_repeated_at_piled_sites"] == 60 and capture.counts["matched_on"] == "every SNP record"
    # vireo: match the pileup ids against every SNP record id of the file, then read those records
    ids, rows = vcf_mod.snp_record_ids(vcf)
    pos = match(capture.variants, ids)
    every = vcf_mod.read_genotypes(vcf)
    by_row = {int(r): k for k, r in enumerate(every.rows)}
    expected = every.gp[[by_row[int(rows[p])] for p in pos]]
    assert np.array_equal(capture.GPb, expected)


def test_matching_on_every_record_equals_the_fast_path_without_repeats(tmp_path):
    sc_ = synth.make_capture(n_donors=3, n_sites=300, n_cells=20, seed=8, blank_rate=0.1)
    paths = synth.write_capture(sc_, str(tmp_path / "cap"), vcf_rows=[i for i in range(300) if i % 5])
    vcf = _vcf_with_repeats(sc_, str(tmp_path / "unpiled.vcf.gz"), n_repeats=0)
    fast = S.load_capture("t", cellsnp=paths["cellsnp"], vcf=vcf)
    full = S.load_capture("t", cellsnp=paths["cellsnp"], vcf=vcf, match="all")
    assert fast.counts["matched_on"] == "records at piled sites" and full.counts["matched_on"] == "every SNP record"
    assert fast.variants == full.variants and np.array_equal(fast.pileup_rows, full.pileup_rows)
    assert np.array_equal(fast.vcf_rows, full.vcf_rows)
    assert np.array_equal(fast.GPb, full.GPb) and np.array_equal(fast.missing, full.missing)
    for name in ("AD", "DP"):
        assert (getattr(fast, name) != getattr(full, name)).nnz == 0


@pytest.mark.vireo
def test_repeated_donor_records_match_vireos_loader(tmp_path):
    sc_ = synth.make_capture(n_donors=4, n_sites=400, n_cells=60, seed=7, doublets=((0, 1, 0.5),), depth=(1, 3))
    paths = synth.write_capture(sc_, str(tmp_path / "cap"))
    vcf = _vcf_with_repeats(sc_, str(tmp_path / "dup.vcf.gz"))
    native = S.load_capture("t", cellsnp=paths["cellsnp"], vcf=vcf)
    vireo = S.load_capture("t", cellsnp=paths["cellsnp"], vcf=vcf, loader="vireo")
    assert native.variants == list(vireo.variants)
    assert np.array_equal(native.GPb, vireo.GPb) and np.array_equal(native.missing, vireo.missing)


def test_cache_write_failure_warns_and_a_corrupt_entry_is_a_miss(synthetic_capture, tmp_path):
    from demux_tools import cache

    cap = synthetic_capture(n_donors=3, n_sites=40, n_cells=6, seed=17)
    locked = tmp_path / "locked"
    locked.mkdir()
    os.chmod(str(locked), 0o500)
    try:
        if os.access(str(locked), os.W_OK):
            pytest.skip("running as a user who can write a read-only folder")
        res = S.load_capture("t", cellsnp=cap["cellsnp"], vcf=cap["vcf"], cache_dir=str(locked / "sub"))
        assert res.load_info["cache"] == "miss, not written" and res.n_sites == 40
    finally:
        os.chmod(str(locked), 0o700)
    first = S.load_capture("t", cellsnp=cap["cellsnp"], vcf=cap["vcf"], cache_dir=str(tmp_path / "c"))
    entry = cache.entry_path(str(tmp_path / "c"), first.load_info["cache_key"])
    with open(entry, "r+b") as f:
        f.truncate(100)
    again = S.load_capture("t", cellsnp=cap["cellsnp"], vcf=cap["vcf"], cache_dir=str(tmp_path / "c"))
    assert again.load_info["cache"] == "miss" and np.array_equal(again.GPb, first.GPb)


def test_long_pair_names_are_cut_as_vireo_cuts_them(tmp_path):
    names = ["donor_%s_%d" % ("x" * 50, k) for k in range(3)]
    sc_ = synth.make_capture(n_donors=3, n_sites=60, n_cells=9, seed=18, donor_names=names)
    paths = synth.write_capture(sc_, str(tmp_path / "cap"))
    capture = S.load_capture("t", cellsnp=paths["cellsnp"], vcf=paths["vcf"])
    res = S.score(capture, _given_fit(capture))
    assert set(len(x) for x in res.table["best_doublet"]) == {100}
    assert all(len(x) > 100 for x in res.table["best_pair"])


# ---------------------------------------------------------------------------
# Commands without vireoSNP


def test_site_terms_with_a_held_fit_and_the_n_donor_guard(synthetic_capture, run_tools, tmp_path):
    cap = synthetic_capture(n_donors=3, n_sites=60, n_cells=8, seed=14)
    fit_json = tmp_path / "held.fit.json"
    fit_json.write_text(json.dumps({"theta_mean": THETA[0][0], "theta_sum": THETA[1][0]}))
    capture = _load(cap)
    ll = S.loglik_singlet(capture.AD, capture.DP, capture.GPb, THETA[0], THETA[1])
    res = run_tools("site-terms", "--cellsnp", cap["cellsnp"], "--vcf", cap["vcf"], "--cell", capture.barcodes[2],
                    "--cell", capture.barcodes[4], "--fit", fit_json, "--out", tmp_path / "terms.tsv")
    assert res.returncode == 0, res.stderr
    t = pd.read_csv(str(tmp_path / "terms.tsv"), sep="\t", dtype={"cell": str})
    for j in (2, 4):
        rows = t[t["cell"] == capture.barcodes[j]]
        sums = [rows["term__" + d].sum() for d in capture.donors]
        np.testing.assert_allclose(sums, ll[j], atol=1e-4)
    res = run_tools("score", "--cellsnp", cap["cellsnp"], "--vcf", cap["vcf"], "--n-donor", "4",
                    "--out", tmp_path / "s")
    assert res.returncode == 2 and "differs from the 3 donors" in res.stderr


def test_validate_vireo_counts_differences(tmp_path, run_tools):
    ref = tmp_path / "vireo"
    ref.mkdir()
    header = "cell\tdonor_id\tprob_max\tprob_doublet\tn_vars\tbest_singlet\tbest_doublet\tdoublet_logLikRatio\n"
    rows = ["a-1\td0\t1.00e+00\t1.00e-05\t30\td0\td0,d1\t-5.000", "b-1\td1\t9.95e-01\t2.00e-03\t25\td1\td0,d1\t-2.500"]
    (ref / "donor_ids.tsv").write_text(header + "\n".join(rows) + "\n")
    ours = tmp_path / "ours.cells.tsv.gz"
    changed = [rows[0], rows[1].replace("9.95e-01", "9.96e-01")]
    ours.write_bytes(gzip.compress((header + "\n".join(changed) + "\n").encode()))
    res = run_tools("validate-vireo", "--scored", tmp_path / "ours", "--vireo-dir", ref, "--out", tmp_path / "v.json")
    assert res.returncode == 1 and "1 of 2 cells identical" in res.stdout
    v = json.load(open(str(tmp_path / "v.json")))
    assert v["differ_prob_max"] == 1 and v["maxabs_prob_max"] == pytest.approx(0.001)
    res = run_tools("validate-vireo", "--scored", ours, "--vireo-dir", ref, "--tolerance", "0.002")
    assert res.returncode == 0, res.stdout


def test_vireo_out_is_an_output_and_names_must_be_file_names(synthetic_capture, run_tools, tmp_path):
    cap = synthetic_capture(n_donors=3, n_sites=40, n_cells=6, seed=19)
    out = tmp_path / "vireo_out"
    out.mkdir()
    (out / "donor_ids.tsv").write_text("cell\n")
    res = run_tools("score", "--cellsnp", cap["cellsnp"], "--vcf", cap["vcf"], "--vireo-out", out,
                    "--out", tmp_path / "s")
    assert res.returncode == 2 and "already holds a donor_ids.tsv" in res.stderr
    (tmp_path / "m.tsv").write_text("capture\tcellsnp\tvcf\n../x\t%s\t%s\n" % (cap["cellsnp"], cap["vcf"]))
    res = run_tools("score-many", "--manifest", tmp_path / "m.tsv", "--out-dir", tmp_path / "o")
    assert res.returncode == 2 and "cannot be used as a file name" in res.stderr


# ---------------------------------------------------------------------------
# Against vireoSNP


@pytest.mark.vireo
@pytest.mark.parametrize("tag,noise,blank", [("GT", 0.0, 0.15), ("GP", 0.2, 0.1), ("PL", 0.1, 0.0)])
def test_native_loader_equals_vireos(synthetic_capture, tag, noise, blank):
    cap = synthetic_capture(n_donors=5, n_sites=150, n_cells=25, blank_rate=blank, tags=("GT", "GP", "PL"),
                            gp_noise=noise, vcf_rows=[i for i in range(150) if i % 7], seed=15)
    a = _load(cap, geno_tag=tag)
    b = _load(cap, geno_tag=tag, loader="vireo")
    assert a.variants == b.variants and a.barcodes == b.barcodes and a.donors == b.donors
    for name in ("AD", "DP"):
        x, y = getattr(a, name), getattr(b, name)
        assert x.format == y.format == "csc" and x.data.dtype == y.data.dtype
        assert np.array_equal(x.indptr, y.indptr) and np.array_equal(x.indices, y.indices)
        assert np.array_equal(x.data, y.data)
    assert np.array_equal(a.GPb, b.GPb) and np.array_equal(a.missing, b.missing)
    assert np.array_equal(a.pileup_rows, b.pileup_rows)


@pytest.mark.vireo
@pytest.mark.parametrize("tag", ["GT", "GP"])
def test_score_reproduces_vireos_command_line(synthetic_capture, run_tools, tmp_path, tag):
    cap = synthetic_capture(n_donors=4, n_sites=250, n_cells=60, blank_rate=0.1, tags=("GT", "GP"), gp_noise=0.1,
                            doublets=((0, 1, 0.5), (2, 3, 0.5), (1, 2, 0.8)), seed=16)
    ref = _run_vireo(cap["cellsnp"], cap["vcf"], str(tmp_path / "vireo"), tmp_path, tag=tag)
    version = _vireo_version()
    res = run_tools("score", "--cellsnp", cap["cellsnp"], "--vcf", cap["vcf"], "--geno-tag", tag,
                    "--vireo-format", version, "--vireo-out", tmp_path / "ours", "--out", tmp_path / "s")
    assert res.returncode == 0, res.stderr
    res = run_tools("validate-vireo", "--scored", tmp_path / "s", "--vireo-dir", ref, "--files", tmp_path / "ours")
    assert res.returncode == 0, res.stdout + res.stderr
    rec = json.load(open(str(tmp_path / "s.fit.json")))
    assert rec["values_from"] == "vireo"
    assert rec["check"] == {"maxdiff_prob_singlet": 0.0, "maxdiff_prob_doublet": 0.0, "maxdiff_doublet_llr": 0.0}
    assert rec["provenance"]["argv"][1] == "score" and len(rec["n_iter"]) == 2


@pytest.mark.vireo
def test_long_donor_names_reproduce_vireos_text(run_tools, tmp_path):
    names = ["donor_%s_%d" % ("x" * 50, k) for k in range(3)]
    cap = synth.build_on_disk(str(tmp_path / "cap"), n_donors=3, n_sites=200, n_cells=30, seed=22, donor_names=names,
                              doublets=((0, 1, 0.5), (1, 2, 0.5)))
    ref = _run_vireo(cap["cellsnp"], cap["vcf"], str(tmp_path / "vireo"), tmp_path)
    res = run_tools("score", "--cellsnp", cap["cellsnp"], "--vcf", cap["vcf"], "--vireo-format", _vireo_version(),
                    "--vireo-out", tmp_path / "ours", "--out", tmp_path / "s")
    assert res.returncode == 0, res.stderr
    res = run_tools("validate-vireo", "--scored", tmp_path / "s", "--vireo-dir", ref, "--files", tmp_path / "ours")
    assert res.returncode == 0, res.stdout + res.stderr


@pytest.mark.vireo
def test_score_many_writes_summary_and_validation(synthetic_capture, run_tools, tmp_path):
    caps = [synthetic_capture("c%d" % k, n_donors=3, n_sites=100 + 40 * k, n_cells=20, seed=20 + k) for k in range(2)]
    ref = _run_vireo(caps[0]["cellsnp"], caps[0]["vcf"], str(tmp_path / "vireo0"), tmp_path)
    lines = ["capture\tcellsnp\tvcf\tvireo_dir",
             "c0\t%s\t%s\t%s" % (caps[0]["cellsnp"], caps[0]["vcf"], ref),
             "c1\t%s\t%s\t" % (caps[1]["cellsnp"], caps[1]["vcf"])]
    (tmp_path / "m.tsv").write_text("\n".join(lines) + "\n")
    res = run_tools("score-many", "--manifest", tmp_path / "m.tsv", "--out-dir", tmp_path / "out", "--jobs", "2",
                    "--vireo-format", _vireo_version(), "--cache-dir", tmp_path / "cache")
    assert res.returncode == 0, res.stdout + res.stderr
    summary = pd.read_csv(str(tmp_path / "out" / "summary.tsv"), sep="\t")
    assert summary["capture"].tolist() == ["c1", "c0"]          # largest DP matrix first
    assert (summary["status"] == "ok").all() and (summary["cache"] == "miss").all()
    validation = pd.read_csv(str(tmp_path / "out" / "validation.tsv"), sep="\t")
    assert validation["capture"].tolist() == ["c0"] and bool(validation["ok"].iloc[0])
    assert os.path.exists(str(tmp_path / "out" / "c1" / "base.cells.tsv.gz"))
    lines.append("c2\t%s\t%s\t" % (str(tmp_path / "absent"), caps[1]["vcf"]))
    (tmp_path / "m2.tsv").write_text("\n".join(lines) + "\n")
    res = run_tools("score-many", "--manifest", tmp_path / "m2.tsv", "--out-dir", tmp_path / "out2",
                    "--vireo-format", _vireo_version(), "--cache-dir", tmp_path / "cache")
    assert res.returncode == 1
    summary = pd.read_csv(str(tmp_path / "out2" / "summary.tsv"), sep="\t")
    assert summary.set_index("capture").loc["c2", "status"].startswith("error")
    assert (summary.set_index("capture").loc[["c0", "c1"], "cache"] == "hit").all()


# The cells table of this capture at the defaults, written in the vireo image
# (vireoSNP 0.5.7) when the Scorer was added. Later changes to the Scorer must
# leave it byte for byte unchanged (compared decompressed).
DEFAULT_TABLE = os.path.join(HERE, "data", "scorer_default.cells.tsv")


def default_capture(tmp_path):
    return synth.build_on_disk(str(tmp_path / "default"), n_donors=4, n_sites=300, n_cells=60, blank_rate=0.1,
                               doublets=((0, 1, 0.5), (2, 3, 0.5)), seed=21)


@pytest.mark.vireo
def test_default_output_is_unchanged(run_tools, tmp_path):
    cap = default_capture(tmp_path)
    res = run_tools("score", "--cellsnp", cap["cellsnp"], "--vcf", cap["vcf"], "--out", tmp_path / "s")
    assert res.returncode == 0, res.stderr
    written = gzip.open(str(tmp_path / "s.cells.tsv.gz"), "rb").read()
    assert written == open(DEFAULT_TABLE, "rb").read()


@pytest.mark.parametrize("blocked,needed", [(("scipy", "vireoSNP"), "scipy"), (("vireoSNP",), "vireoSNP")])
def test_a_missing_library_is_named(synthetic_capture, run_tools, tmp_path, blocked, needed):
    cap = synthetic_capture(n_donors=2, n_sites=20, n_cells=4, seed=20)
    res = run_tools("score", "--cellsnp", cap["cellsnp"], "--vcf", cap["vcf"], "--out", tmp_path / "s",
                    blocked=blocked)
    assert res.returncode == 2 and needed in res.stderr and "Traceback" not in res.stderr
