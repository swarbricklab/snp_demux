"""demux_tools.synth: synthetic captures are deterministic and follow vireo's read model."""
import gzip
import hashlib
import os

import numpy as np
import pytest

from demux_tools import synth
from demux_tools.cellsnp import read_cellsnp


def _content_digest(paths):
    """sha1 over the decompressed text of every written file (gzip bytes may vary by zlib build)."""
    h = hashlib.sha1()
    files = [os.path.join(paths["cellsnp"], f) for f in sorted(os.listdir(paths["cellsnp"]))]
    for p in files + [paths["vcf"], paths["truth"]]:
        raw = open(p, "rb").read()
        h.update(gzip.decompress(raw) if raw[:2] == b"\x1f\x8b" else raw)
    return h.hexdigest()


def test_same_seed_same_bytes(synthetic_capture):
    a = synthetic_capture("a", n_donors=4, n_sites=50, n_cells=30, blank_rate=0.1, oth_rate=0.05,
                          tags=("GT", "GP", "PL"), gp_noise=0.1, seed=11)
    b = synthetic_capture("b", n_donors=4, n_sites=50, n_cells=30, blank_rate=0.1, oth_rate=0.05,
                          tags=("GT", "GP", "PL"), gp_noise=0.1, seed=11)
    for key in ("vcf", "truth"):
        assert open(a[key], "rb").read() == open(b[key], "rb").read()
    for name in os.listdir(a["cellsnp"]):
        assert open(os.path.join(a["cellsnp"], name), "rb").read() == open(os.path.join(b["cellsnp"], name), "rb").read()
    c = synthetic_capture("c", n_donors=4, n_sites=50, n_cells=30, seed=12)
    assert _content_digest(c) != _content_digest(a)


# The default capture's content, the same with numpy 1.21 to 2.x and pandas 1.3 to 2.2
# (RandomState's stream is fixed). A change here changes every test built on it.
DEFAULT_DIGEST = "49b3ca1628a6c4bae936fc418808cddab339bc1b"


def test_default_capture_is_the_same_everywhere(synthetic_capture):
    assert _content_digest(synthetic_capture()) == DEFAULT_DIGEST


def test_shapes_truth_and_polymorphic_sites():
    sc = synth.make_capture(n_donors=3, n_sites=30, n_cells=9, doublets=((2, 0, 0.7), (0, 1, 0.5)), seed=2)
    assert sc.donors == ["d0", "d1", "d2"]
    assert sc.pileup.shape == (30, 11)
    assert sc.genotypes.shape == (30, 3) and sc.genotypes.dtype == np.int8
    assert not (sc.genotypes == sc.genotypes[:, :1]).all(axis=1).any()
    t = sc.truth
    assert t["barcode"].tolist() == list(sc.pileup.barcodes)
    assert t["kind"].tolist() == ["singlet"] * 9 + ["doublet"] * 2
    assert t["donor"].tolist()[:4] == ["d0", "d1", "d2", "d0"]
    assert t["pair"].tolist()[9:] == ["d0,d2", "d0,d1"]
    assert t["donor_a"].tolist()[9:] == ["d2", "d0"] and t["share_a"].tolist()[9:] == ["0.7", "0.5"]
    assert len(set(sc.pileup.barcodes)) == 11
    assert all(b.endswith("-1") and len(b) == 18 for b in sc.pileup.barcodes)
    chrom = list(sc.pileup.chrom)
    assert chrom == sorted(chrom) and set(chrom) == {"chr1", "chr2"}
    for c in ("chr1", "chr2"):
        pos = sc.pileup.pos[np.asarray(chrom) == c]
        assert (np.diff(pos) > 0).all()


def test_counts_follow_the_genotypes():
    theta = (0.01, 0.5, 0.99)
    sc = synth.make_capture(n_donors=2, n_sites=400, n_cells=40, depth=(5, 5), coverage=1.0, theta=theta,
                            doublets=(), seed=4)
    ad, dp = sc.pileup.to_dense("AD"), sc.pileup.to_dense("DP")
    assert (dp == 5).all()
    for c in range(40):
        g = sc.genotypes[:, c % 2]
        for k in range(3):
            frac = ad[g == k, c].sum() / float(dp[g == k, c].sum())
            assert abs(frac - theta[k]) < 0.15


def test_blank_and_vcf_text(synthetic_capture):
    built = synthetic_capture(n_donors=3, n_sites=200, n_cells=3, blank_rate=0.25, tags=("GT", "GP", "PL"),
                              gp_noise=0.3, seed=8)
    sc = built["capture"]
    assert 0.15 < sc.blank.mean() < 0.35
    lines = gzip.decompress(open(built["vcf"], "rb").read()).decode().splitlines()
    header = [l for l in lines if l.startswith("#CHROM")][0].split("\t")
    assert header[9:] == ["d0", "d1", "d2"]
    records = [l.split("\t") for l in lines if not l.startswith("#")]
    assert len(records) == 200
    for s, rec in enumerate(records):
        assert rec[8] == "GT:GP:PL"
        for d, field in enumerate(rec[9:]):
            if sc.blank[s, d]:
                assert field == "./.:.:."
            else:
                gt, gp, pl = field.split(":")
                g = int(sc.genotypes[s, d])
                assert gt == synth.GT_TEXT[g]
                probs = [float(x) for x in gp.split(",")]
                assert int(np.argmax(probs)) == g and abs(sum(probs) - 1) < 1e-5
                assert pl.split(",")[g] == "0"


def test_vcf_rows_and_reading_back(synthetic_capture):
    built = synthetic_capture(n_sites=10, vcf_rows=[0, 2, 4])
    lines = [l for l in gzip.decompress(open(built["vcf"], "rb").read()).decode().splitlines()
             if not l.startswith("#")]
    p = built["capture"].pileup
    assert [l.split("\t")[1] for l in lines] == [str(p.pos[k]) for k in (0, 2, 4)]
    q = read_cellsnp(built["cellsnp"])
    assert q.variants == p.variants


def test_bad_arguments():
    with pytest.raises(ValueError):
        synth.make_capture(doublets=((0, 0, 0.5),))
    with pytest.raises(ValueError):
        synth.make_capture(n_donors=2, donor_names=["x", "x"])
    with pytest.raises(ValueError):
        synth.make_capture(depth=(0, 3))
    with pytest.raises(ValueError):
        synth.write_vcf(synth.make_capture(), "unused.vcf", tags=("DS",))


def test_default_donor_names_sort_like_their_index():
    names = synth.default_donor_names(12)
    assert names[:2] == ["d00", "d01"] and names == sorted(names)
    assert synth.default_donor_names(1) == ["d0"]


def test_synth_command(run_tools, tmp_path):
    out = tmp_path / "synthetic"
    res = run_tools("synth", "--out", out, "--donors", "4", "--sites", "30", "--cells", "8",
                    "--doublet", "0:3:0.6", "--doublet", "1:2:0.5", "--tags", "GT,GP", "--blank-rate", "0.1")
    assert res.returncode == 0, res.stderr
    assert sorted(os.listdir(str(out))) == ["cellsnp", "donors.vcf.gz", "provenance.json", "truth.tsv"]
    sc = synth.make_capture(n_donors=4, n_sites=30, n_cells=8, doublets=((0, 3, 0.6), (1, 2, 0.5)),
                            blank_rate=0.1)
    p = read_cellsnp(str(out / "cellsnp"), layers=("AD", "DP"))
    assert p.variants == sc.pileup.variants
    assert np.array_equal(p.to_dense("AD"), sc.pileup.to_dense("AD"))
    res = run_tools("synth", "--out", tmp_path / "bad", "--doublet", "0:9:0.5")
    assert res.returncode == 2 and "bad doublet" in res.stderr
