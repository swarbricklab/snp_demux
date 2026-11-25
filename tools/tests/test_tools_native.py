"""The native fit against vireoSNP's, held allele ratios, and score-many with variants."""
import json
import os

import numpy as np
import pandas as pd
import pytest

from demux_tools import scorer as S

pytestmark = [pytest.mark.scipy, pytest.mark.vireo]


def _load(cap, **kw):
    return S.load_capture("t", cellsnp=cap["cellsnp"], vcf=cap["vcf"], **kw)


@pytest.mark.parametrize("seed,tag,blank,noise", [(41, "GT", 0.0, 0.0), (42, "GT", 0.2, 0.0), (43, "GP", 0.1, 0.2),
                                                  (44, "GT", 0.3, 0.0)])
def test_native_fit_matches_vireo(synthetic_capture, seed, tag, blank, noise):
    cap = synthetic_capture(n_donors=5, n_sites=250, n_cells=80, blank_rate=blank, tags=("GT", "GP"), gp_noise=noise,
                            doublets=((0, 1, 0.5), (2, 4, 0.7), (3, 1, 0.5)), seed=seed)
    capture = _load(cap, geno_tag=tag)
    vireo = S.fit(capture, backend="vireo")
    native = S.fit(capture, backend="native")
    np.testing.assert_allclose(native.theta_mean, vireo.theta_mean, rtol=1e-9)
    np.testing.assert_allclose(native.theta_sum, vireo.theta_sum, rtol=1e-9)
    assert native.n_iter == vireo.n_iter
    assert native.elbo == pytest.approx(vireo.elbo, rel=1e-9)
    a, b = S.score(capture, vireo), S.score(capture, native)
    for c in ("donor_id", "best_singlet", "best_doublet", "best_pair", "second_donor", "n_vars"):
        assert a.table[c].tolist() == b.table[c].tolist(), c
    np.testing.assert_allclose(b.prob_singlet, a.prob_singlet, atol=1e-9)
    np.testing.assert_allclose(b.prob_doublet, a.prob_doublet, atol=1e-9)


def test_holding_vireos_ratios_gives_its_likelihoods(synthetic_capture):
    cap = synthetic_capture(n_donors=4, n_sites=200, n_cells=50, blank_rate=0.1, seed=45)
    capture = _load(cap)
    vireo = S.fit(capture)
    held = S.fit(capture, theta=(vireo.theta_mean, vireo.theta_sum))
    a, b = S.score(capture, vireo), S.score(capture, held)
    assert np.array_equal(a.ll_singlet, b.ll_singlet) and np.array_equal(a.ll_doublet, b.ll_doublet)
    assert a.source == "vireo" and b.source == "scorer"
    assert a.table["donor_id"].tolist() == b.table["donor_id"].tolist()
    assert a.check["maxdiff_prob_singlet"] == 0.0


def test_score_many_variants(synthetic_capture, run_tools, tmp_path):
    import sys

    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from test_tools_scorer import _run_vireo, _vireo_version

    caps = [synthetic_capture("c%d" % k, n_donors=3, n_sites=150, n_cells=30, blank_rate=0.2, seed=46 + k)
            for k in range(2)]
    ref = _run_vireo(caps[0]["cellsnp"], caps[0]["vcf"], str(tmp_path / "vireo0"), tmp_path)
    (tmp_path / "m.tsv").write_text("capture\tcellsnp\tvcf\tvireo_dir\n" + "".join(
        "c%d\t%s\t%s\t%s\n" % (k, c["cellsnp"], c["vcf"], ref if k == 0 else "") for k, c in enumerate(caps)))
    (tmp_path / "v.tsv").write_text(
        "variant\ttheta\tgroup_by\tfill_blanks\tbackend\n"
        "base\trefit\t\t\t\n"
        "missing\tbase\t@missing\t\t\n"
        "hwe\tbase\t@missing\thwe:AF\t\n"
        "native\trefit\t\t\tnative\n")
    res = run_tools("score-many", "--manifest", tmp_path / "m.tsv", "--variants", tmp_path / "v.tsv", "--out-dir",
                    tmp_path / "out", "--vireo-format", _vireo_version(), "--cache-dir", tmp_path / "cache")
    assert res.returncode == 0, res.stdout + res.stderr
    summary = pd.read_csv(str(tmp_path / "out" / "summary.tsv"), sep="\t")
    assert summary["variant"].tolist() == ["base", "missing", "hwe", "native"] * 2
    assert (summary["status"] == "ok").all()
    assert set(summary["backend"]) == {"vireo", "held", "native"}
    validation = pd.read_csv(str(tmp_path / "out" / "validation.tsv"), sep="\t")
    assert validation["capture"].tolist() == ["c0"] and bool(validation["ok"].iloc[0])
    for k in range(2):
        base = pd.read_csv(str(tmp_path / "out" / ("c%d" % k) / "base.cells.tsv.gz"), sep="\t", dtype=str,
                           keep_default_na=False)
        native = pd.read_csv(str(tmp_path / "out" / ("c%d" % k) / "native.cells.tsv.gz"), sep="\t", dtype=str,
                             keep_default_na=False)
        assert base["donor_id"].tolist() == native["donor_id"].tolist()
        rec = json.load(open(str(tmp_path / "out" / ("c%d" % k) / "hwe.fit.json")))
        assert rec["variant"]["fill_blanks"]["filled"] > 0
