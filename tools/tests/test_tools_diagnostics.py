"""demux_tools.diagnostics, the genotype overrides and the variants: arithmetic checked on small cases."""
import json
import os

import numpy as np
import pandas as pd
import pytest

from demux_tools import diagnostics as dg
from demux_tools import scorer as S
from demux_tools import sites, synth

pytestmark = pytest.mark.scipy

PLUG_IN = ([[0.01, 0.5, 0.99]], [[1e8, 1e8, 1e8]])     # concentrations so large that E[log theta] = log theta


def _tiny(ad, dp, gp, missing=None, donors=None):
    """A Capture from dense (sites x cells) counts and (sites, donors, 3) genotype probabilities."""
    from scipy import sparse

    ad, dp, gp = np.asarray(ad), np.asarray(dp), np.asarray(gp, dtype=float)
    n_sites, n_cells = dp.shape
    donors = donors or ["d%d" % k for k in range(gp.shape[1])]
    missing = np.zeros(gp.shape[:2], dtype=bool) if missing is None else np.asarray(missing, dtype=bool)
    return S.Capture("tiny", ["c%d-1" % k for k in range(n_cells)], ["chr1_%d_A_G" % (100 + k) for k in range(n_sites)],
                     sparse.csc_matrix(ad.astype(np.int64)), sparse.csc_matrix(dp.astype(np.int64)), donors, gp,
                     missing, np.arange(n_sites), np.arange(n_sites), {}, {}, [], {})


def _load(cap, **kw):
    return S.load_capture("t", cellsnp=cap["cellsnp"], vcf=cap["vcf"], **kw)


def _scored(cap, theta=([[0.03, 0.5, 0.97]], [[400.0, 300.0, 400.0]]), **kw):
    capture = _load(cap, **kw)
    fit = S.Fit(theta[0], theta[1], doublet_prior=S.default_doublet_prior(capture.n_cells))
    return capture, fit, S.score(capture, fit)


def test_one_read_blank_costs_mean_field_against_exact():
    third = [1 / 3] * 3
    capture = _tiny([[1]], [[1]], [[third, [0, 0, 1.0]]], missing=[[True, False]])
    fit = S.Fit(*PLUG_IN, doublet_prior=1e-5)
    scores = S.score(capture, fit)
    assert -scores.ll_singlet[0, 0] == pytest.approx(-(np.log(0.01) + np.log(0.5) + np.log(0.99)) / 3, abs=1e-6)
    assert -scores.ll_singlet[0, 0] == pytest.approx(1.77, abs=0.005)
    table, ll = dg.marginal_singlet(capture, fit, scores.ll_singlet, np.array([0]))
    assert -ll[0, 0] == pytest.approx(np.log(2), abs=1e-6)                 # log(1 / mean(0.01, 0.5, 0.99))
    assert ll[0, 1] == scores.ll_singlet[0, 1]                              # hard genotype: unchanged
    assert table.loc[0, "blank_umis"] == 1


def test_hwe_fill_at_five_per_cent():
    capture = _tiny([[0]], [[1]], [[[1 / 3] * 3]], missing=[[True]])
    filled, counts = S.fill_missing(capture, how="hwe", af=np.array([0.05]))
    np.testing.assert_allclose(filled.GPb[0, 0], [0.9025, 0.095, 0.0025])
    assert counts == {"blanks": 1, "filled": 1, "left_uniform": 0} and filled.missing[0, 0]
    fit = S.Fit(*PLUG_IN, doublet_prior=0.0)
    cost = -S.loglik_singlet(filled.AD, filled.DP, filled.GPb, fit.theta_mean, fit.theta_sum)[0, 0]
    assert cost == pytest.approx(0.09, abs=0.005)
    dropped, counts = S.fill_missing(capture, how="drop-site")
    assert dropped.n_sites == 0 and counts["sites_dropped"] == 1


def test_harden_keeps_blanks():
    gp = [[[0.1, 0.8, 0.1], [1 / 3] * 3]]
    hard, counts = S.harden(_tiny([[0]], [[1]], gp, missing=[[False, True]]))
    np.testing.assert_array_equal(hard.GPb[0, 0], [0, 1, 0])
    np.testing.assert_array_equal(hard.GPb[0, 1], [1 / 3] * 3)
    assert counts == {"hardened": 1}


def test_group_margins_add_up_to_the_margin(synthetic_capture):
    cap = synthetic_capture(n_donors=4, n_sites=150, n_cells=30, blank_rate=0.2, seed=31)
    capture, fit, scores = _scored(cap)
    best, second = dg.best_and_second(scores)
    labels = sites.missing_pattern(capture.missing, capture.donors)
    t = dg.group_margins(capture, fit, best, second, labels, "missing")
    margin_cols = [c for c in t.columns if c.startswith("margin__")]
    np.testing.assert_allclose(t[margin_cols].sum(axis=1), scores.table["margin"], atol=1e-9)
    umi_cols = [c for c in t.columns if c.startswith("umis__")]
    np.testing.assert_array_equal(t[umi_cols].sum(axis=1), scores.table["umis"])
    site_cols = [c for c in t.columns if c.startswith("nsites__")]
    np.testing.assert_array_equal(t[site_cols].sum(axis=1), scores.table["n_vars"])


def test_exact_marginal_is_mean_field_without_blanks_and_never_below(synthetic_capture):
    hard = synthetic_capture("hard", n_donors=3, n_sites=80, n_cells=12, seed=32)
    capture, fit, scores = _scored(hard)
    _, ll = dg.marginal_singlet(capture, fit, scores.ll_singlet, dg.best_and_second(scores)[0])
    assert np.array_equal(ll, scores.ll_singlet)
    soft = synthetic_capture("soft", n_donors=3, n_sites=80, n_cells=12, blank_rate=0.3, tags=("GP",),
                             gp_noise=0.3, seed=32)
    capture, fit, scores = _scored(soft, geno_tag="GP")
    table, ll = dg.marginal_singlet(capture, fit, scores.ll_singlet, dg.best_and_second(scores)[0])
    assert (ll >= scores.ll_singlet - 1e-9).all() and (ll > scores.ll_singlet + 1e-6).any()
    np.testing.assert_allclose(table["margin_marginal"], table["ll_best_marginal"] - table["ll_second_marginal"])


def test_pairwise_without_blanks_is_the_best_singlet(synthetic_capture):
    cap = synthetic_capture(n_donors=4, n_sites=100, n_cells=20, seed=33)
    capture, fit, scores = _scored(cap)
    t = dg.pairwise_shared(capture, fit)
    assert t["pairwise_best"].tolist() == scores.table["best_singlet"].tolist()
    assert (t["pairwise_wins"] == 3).all()
    np.testing.assert_allclose(t["pairwise_margin"], scores.table["margin"], atol=1e-9)
    np.testing.assert_array_equal(t["pairwise_shared_umis"], scores.table["umis"])


def test_pairwise_leaves_a_blank_out_of_the_pairs_it_affects(synthetic_capture):
    cap = synthetic_capture(n_donors=3, n_sites=120, n_cells=15, blank_rate=0.25, seed=34)
    capture, fit, scores = _scored(cap)
    t = dg.pairwise_shared(capture, fit)
    for c in range(capture.n_cells):
        b = capture.donors.index(t.loc[c, "pairwise_best"])
        margins = []
        for e in range(3):
            if e == b:
                continue
            rows = np.flatnonzero(~capture.missing[:, b] & ~capture.missing[:, e])
            ll = S.loglik_singlet(capture.AD[rows], capture.DP[rows], capture.GPb[rows], fit.theta_mean,
                                  fit.theta_sum)
            margins.append(ll[c, b] - ll[c, e])
        assert t.loc[c, "pairwise_margin"] == pytest.approx(min(margins), abs=1e-9)


def test_top_site_deltas_add_up_to_the_margin(synthetic_capture):
    cap = synthetic_capture(n_donors=3, n_sites=60, n_cells=10, blank_rate=0.1, seed=35)
    capture, fit, scores = _scored(cap)
    best, second = dg.best_and_second(scores)
    gene = np.array(["g%d" % (k // 10) for k in range(capture.n_sites)])
    long, per_cell = dg.top_sites(capture, fit, best, second, n=1000, gene=gene)
    sums = long.groupby("cell", sort=False)["delta"].sum()
    np.testing.assert_allclose(sums.loc[capture.barcodes].to_numpy(), scores.table["margin"], atol=1e-9)
    last = long.groupby("cell", sort=False)["cum_share"].last()
    np.testing.assert_allclose(last.loc[capture.barcodes].to_numpy(), 1.0, atol=1e-9)
    assert (long.groupby("cell")["rank"].min() == 1).all()
    first = long[long["rank"] == 1].set_index("cell")["delta"]
    np.testing.assert_allclose(per_cell["top1_site_share"],
                               first.loc[capture.barcodes].to_numpy() / scores.table["margin"].to_numpy())
    assert set(per_cell["top_gene"]) <= set(gene)
    np.testing.assert_allclose(per_cell["top10_genes_share"], 1.0)       # six genes: all of them
    rows = long[long["cell"] == capture.barcodes[0]].copy()
    rows["gene"] = gene[[capture.variants.index(v) for v in rows["variant"]]]
    by_gene = rows.groupby("gene")["delta"].sum()
    assert per_cell.loc[0, "top_gene"] == by_gene.idxmax()
    assert per_cell.loc[0, "top_gene_share"] == pytest.approx(by_gene.max() / scores.table.loc[0, "margin"])


def test_cheap_columns_and_expected_fit(synthetic_capture):
    cap = synthetic_capture(n_donors=3, n_sites=400, n_cells=30, depth=(3, 8), doublets=(), seed=36)
    capture, fit, scores = _scored(cap, theta=([[0.01, 0.5, 0.99]], [[1e6, 1e6, 1e6]]))
    best, second = dg.best_and_second(scores)
    t = dg.cheap_columns(capture, fit, best, second)
    np.testing.assert_array_equal(t["alt_umis"], np.asarray(capture.AD.sum(axis=0)).ravel())
    assert (t["n_sites_diff"] <= scores.table["n_vars"]).all() and (t["umis_diff"] <= scores.table["umis"]).all()
    gap = scores.table["fit_per_umi"] - t["fit_expected_per_umi"]
    assert abs(gap.mean()) < 0.05


def test_decoy_donors_are_appended(synthetic_capture, tmp_path):
    cap = synthetic_capture(n_donors=3, n_sites=200, n_cells=30, seed=37)
    other = synth.make_capture(n_donors=5, n_sites=200, n_cells=1, seed=37, donor_names=["x%d" % k for k in range(5)])
    decoys = synth.write_vcf(other, str(tmp_path / "decoys.vcf.gz"), rows=list(range(150)))
    capture = _load(cap)
    more, counts = S.add_donors(capture, decoys, donors=["x1", "x3"])
    assert more.donors == capture.donors + ["x1", "x3"] and more.GPb.shape == (200, 5, 3)
    assert counts["sites_in_extra_vcf"] == 150 and more.missing[150:, 3:].all() and not more.missing[:150, 3:].any()
    np.testing.assert_array_equal(more.GPb[:, :3], capture.GPb)
    fit = S.Fit([[0.01, 0.5, 0.99]], [[1000.0, 1000.0, 1000.0]], doublet_prior=0.0003)
    calls = S.score(more, fit).table["donor_id"]
    assert not calls.isin(["x1", "x3"]).any()
    with pytest.raises(S.UsageError, match="already in the capture"):
        S.add_donors(capture, cap["vcf"], donors=["d0"])


def test_theta_specs():
    assert S.parse_theta("refit") == ("refit", None) and S.parse_theta("hold:a.json") == ("hold", "a.json")
    assert S.parse_theta("fixed:0.01,0.5,0.99:10,20,30") == ("fixed", ([0.01, 0.5, 0.99], [10.0, 20.0, 30.0]))
    with pytest.raises(S.UsageError):
        S.parse_theta("fixed:0.1,0.5")


def test_variants_read_the_capture_once(synthetic_capture, tmp_path, monkeypatch):
    cap = synthetic_capture(n_donors=3, n_sites=120, n_cells=20, blank_rate=0.2, seed=38)
    calls = []
    real = S.load_capture

    def counting(*a, **k):
        calls.append(k.get("info"))
        return real(*a, **k)

    monkeypatch.setattr(S, "load_capture", counting)
    variants = [{"variant": "fixed", "theta": "fixed:0.02,0.5,0.98:500,500,500"},
                {"variant": "held", "theta": "base", "group_by": "@missing;@info:AF:0,0.5,1"},
                {"variant": "hwe", "theta": "base", "fill_blanks": "hwe:AF", "marginal": "true",
                 "top_sites": "3", "extra_columns": "true"},
                {"variant": "thin", "theta": "base", "thin_bp": "1500", "pairwise_shared": "true"}]
    out = tmp_path / "v"
    out.mkdir()
    records = S.run_variants("t", {"cellsnp": cap["cellsnp"], "vcf": cap["vcf"]}, str(out), variants,
                             cache_dir=str(tmp_path / "cache"))
    assert len(calls) == 1 and calls[0] == ["AF"]
    assert [r["theta_mean"] for r in records[1:]] == [records[0]["theta_mean"]] * 3
    held = pd.read_csv(str(out / "held.cells.tsv.gz"), sep="\t")
    assert "margin__missing__none" in held.columns and "margin__AF__[0,0.5)" in held.columns
    base = pd.read_csv(str(out / "fixed.cells.tsv.gz"), sep="\t", dtype=str, keep_default_na=False)
    assert list(base.columns) == S.CELL_COLUMNS
    assert pd.read_csv(str(out / "held.cells.tsv.gz"), sep="\t", dtype=str,
                       keep_default_na=False)[S.CELL_COLUMNS].equals(base)
    hwe = pd.read_csv(str(out / "hwe.cells.tsv.gz"), sep="\t")
    assert {"ll_best_marginal", "blank_umis", "top1_site_share", "alt_umis"} <= set(hwe.columns)
    assert os.path.exists(str(out / "hwe.sites.tsv.gz"))
    rec = json.load(open(str(out / "hwe.fit.json")))
    assert rec["variant"]["fill_blanks"]["filled"] > 0 and rec["variant"]["options"]["theta"] == "base"
    thin = json.load(open(str(out / "thin.fit.json")))
    assert thin["n_sites_used"] < 120
    assert "pairwise_best" in pd.read_csv(str(out / "thin.cells.tsv.gz"), sep="\t").columns


def test_default_score_has_the_d2_columns_and_record(synthetic_capture, run_tools, tmp_path):
    cap = synthetic_capture(n_donors=3, n_sites=60, n_cells=10, seed=39)
    fit_json = tmp_path / "f.json"
    fit_json.write_text(json.dumps({"theta_mean": PLUG_IN[0][0], "theta_sum": [500.0] * 3}))
    res = run_tools("score", "--cellsnp", cap["cellsnp"], "--vcf", cap["vcf"], "--theta", "hold:%s" % fit_json,
                    "--out", tmp_path / "s")
    assert res.returncode == 0, res.stderr
    t = pd.read_csv(str(tmp_path / "s.cells.tsv.gz"), sep="\t")
    assert list(t.columns) == S.CELL_COLUMNS
    rec = json.load(open(str(tmp_path / "s.fit.json")))
    assert rec["backend"] == "held" and rec["variant"]["options"]["theta"].startswith("hold:")
    assert sorted(rec["seconds"]) == ["diagnostics", "fit", "load", "overrides", "score", "write"]


def test_pairwise_with_more_than_64_donors():
    # blank patterns used to be bits of an int64, which loses donors 64 and above (the review's experiment C)
    from scipy import sparse

    rng = np.random.RandomState(0)
    n_sites, n_cells, n = 300, 12, 70
    gt = rng.randint(0, 3, size=(n_sites, n))
    gp = np.eye(3)[gt].astype(float)
    missing = np.zeros((n_sites, n), dtype=bool)
    missing[:150, 66] = True
    gp[missing] = 1 / 3
    dp = rng.randint(0, 3, size=(n_sites, n_cells))
    cell_donor = np.arange(n_cells) % 3 + 64
    ad = rng.binomial(dp, np.array([0.01, 0.5, 0.99])[gt[:, cell_donor]])
    cap = S.Capture("c", ["c%d" % k for k in range(n_cells)], ["chr1_%d_A_G" % k for k in range(n_sites)],
                    sparse.csc_matrix(ad), sparse.csc_matrix(dp), ["d%d" % k for k in range(n)], gp, missing,
                    np.arange(n_sites), np.arange(n_sites), {}, {}, [], {})
    fit = S.Fit([[0.01, 0.5, 0.99]], [[1000.0] * 3], doublet_prior=1e-4)
    t = dg.pairwise_shared(cap, fit)
    for c in range(n_cells):
        b = cap.donors.index(t.loc[c, "pairwise_best"])
        margins = []
        for e in range(n):
            if e == b:
                continue
            rows = np.flatnonzero(~missing[:, b] & ~missing[:, e])
            ll = S.loglik_singlet(cap.AD[rows], cap.DP[rows], cap.GPb[rows], fit.theta_mean, fit.theta_sum)
            margins.append(ll[c, b] - ll[c, e])
        assert t.loc[c, "pairwise_margin"] == pytest.approx(min(margins), abs=1e-9)


def test_top_sites_at_a_zero_margin_and_sites_without_a_gene():
    # two donors with the same genotypes everywhere: every margin is 0
    gp = [[[0, 1.0, 0], [0, 1.0, 0]]] * 4
    capture = _tiny([[1, 0], [0, 1], [1, 1], [0, 0]], [[1, 1], [1, 1], [2, 1], [1, 0]], gp)
    fit = S.Fit(*PLUG_IN, doublet_prior=1e-5)
    scores = S.score(capture, fit)
    best, second = dg.best_and_second(scores)
    long, per_cell = dg.top_sites(capture, fit, best, second, n=5, gene=np.array(["", "", "g1", ""]))
    assert long["cum_share"].isna().all() and per_cell["top1_site_share"].isna().all()
    assert per_cell["top_gene"].tolist() == ["g1", "g1"]
    # different genotypes: the gene shares leave the sites without a gene out
    gp = [[[1.0, 0, 0], [0, 0, 1.0]]] * 4
    capture = _tiny([[0], [0], [0], [0]], [[1], [1], [2], [1]], gp)
    scores = S.score(capture, fit)
    best, second = dg.best_and_second(scores)
    _, per_cell = dg.top_sites(capture, fit, best, second, n=5, gene=np.array(["", "", "g1", ""]))
    assert per_cell.loc[0, "top_gene"] == "g1"
    assert per_cell.loc[0, "top_gene_share"] == pytest.approx(2 / 5)


def test_theta_fixed_values_are_checked_and_names_are_file_names(tmp_path):
    for bad in ("fixed:0,0.5,0.99", "fixed:0.01,0.5,1", "fixed:0.01,0.5,0.99:10,0,10"):
        with pytest.raises(S.UsageError, match="strictly between 0 and 1"):
            S.parse_theta(bad)
    (tmp_path / "v.tsv").write_text("variant\ttheta\nbase\trefit\n../x\tbase\n")
    with pytest.raises(S.UsageError, match="cannot be used as a file name"):
        S.read_variants(str(tmp_path / "v.tsv"))

