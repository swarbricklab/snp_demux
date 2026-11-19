import gzip
import math
import random

import pytest

from conftest import load_script

eg = load_script("elimination_gate")

N_SITES = 400
KNOWN = ["D1", "D2", "D3"]
POOL = KNOWN + ["DX"]


def _truth(seed=0):
    rng = random.Random(seed)
    return {d: [rng.choice([0, 1, 2]) for _ in range(N_SITES)] for d in POOL}


def _noisy(dosages, err, seed):
    """Copy of dosages with a fraction `err` of sites replaced by a different call."""
    rng = random.Random(seed)
    return [rng.choice([x for x in (0, 1, 2) if x != d]) if rng.random() < err else d for d in dosages]


def _gt(d):
    return {0: "0/0", 1: "0/1", 2: "1/1", None: "./."}[d]


def _write_vcf(path, samples, dosages, gz=False):
    lines = ["##fileformat=VCFv4.2", "\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + samples)]
    for i in range(N_SITES):
        calls = [_gt(dosages[s][i]) + ":5" for s in samples]
        lines.append("\t".join(["chr1", str(100 + i), ".", "A", "G", ".", "PASS", ".", "GT:DP"] + calls))
    text = "\n".join(lines) + "\n"
    if gz:
        with gzip.open(path, "wt") as fh:
            fh.write(text)
    else:
        path.write_text(text)


def _scenario(tmp_path, dx_cluster_err, known_err=0.08):
    """vireo learned all four donors; souporcell clusters 0-3 are noisy copies.

    DX's cluster carries `dx_cluster_err` disagreement with DX's learned
    genotype; known donors' clusters carry `known_err` (the ceiling ~ 1 - known_err).
    """
    t = _truth()
    vireo = {d: t[d] for d in POOL}
    soup = {str(i): _noisy(t[d], known_err if d != "DX" else dx_cluster_err, seed=i + 1) for i, d in enumerate(POOL)}
    _write_vcf(tmp_path / "vireo.vcf.gz", POOL, vireo, gz=True)
    _write_vcf(tmp_path / "soup.vcf", list(soup), soup)
    (tmp_path / "map.tsv").write_text(
        "SAMPLE_ID\tSNP_CLUSTER\tGPb_diff_score\nD1\t0\t0.1\nD2\t1\t0.1\nD3\t2\t0.1\nDX\t3\t\n"
    )
    return (
        eg.read_mapping(tmp_path / "map.tsv"),
        eg.read_genotypes(tmp_path / "vireo.vcf.gz"),
        eg.read_genotypes(tmp_path / "soup.vcf"),
    )


def _eval(mapping, vireo, soup, genotyped=KNOWN, **kw):
    args = dict(min_concordance=0.8, max_drop=0.1, min_sites=100, action="flag")
    args.update(kw)
    return eg.evaluate("cap", POOL, set(genotyped), mapping, vireo, soup, **args)


def test_concordance_ignores_missing_calls():
    c, n = eg.concordance({1: 0, 2: 1, 3: 2}, {2: 1, 3: 0, 4: 2})
    assert (c, n) == (0.5, 2)
    c, n = eg.concordance({}, {1: 0})
    assert math.isnan(c) and n == 0


def test_dosage_parsing():
    assert [eg._dosage(g) for g in ("0/0", "0|1", "1/0", "1|1", "./.", ".", "0/2")] == [0, 1, 1, 2, None, None, 1]


def test_reliable_pool_passes(tmp_path):
    row = _eval(*_scenario(tmp_path, dx_cluster_err=0.05))
    assert row["status"] == "PASS", row
    assert row["missing_donor"] == "DX" and row["snp_cluster"] == "3"
    assert row["concordance"] == pytest.approx(0.95, abs=0.03)
    assert row["ceiling"] == pytest.approx(0.92, abs=0.03)
    assert row["n_known"] == 3 and row["n_sites"] == N_SITES


def test_noise_level_agreement_fails(tmp_path):
    # Leftover cluster is essentially unrelated to the learned genotype (~0.45).
    row = _eval(*_scenario(tmp_path, dx_cluster_err=0.55))
    assert row["status"] == "FAIL"
    assert "concordance" in row["reason"]


def test_below_ceiling_margin_fails_even_above_floor(tmp_path):
    # 0.82 clears the 0.8 floor but sits > 0.1 below a ~0.98 ceiling.
    row = _eval(*_scenario(tmp_path, dx_cluster_err=0.18, known_err=0.02))
    assert row["status"] == "FAIL"
    assert "ceiling" in row["reason"]


def test_too_few_sites_fails(tmp_path):
    row = _eval(*_scenario(tmp_path, dx_cluster_err=0.05), min_sites=N_SITES + 1)
    assert row["status"] == "FAIL"
    assert "shared sites" in row["reason"]
    assert row["n_known"] == 0  # known donors below min_sites don't set a ceiling


def test_all_genotyped_is_not_applicable(tmp_path):
    row = _eval(*_scenario(tmp_path, 0.05), genotyped=POOL)
    assert row["status"] == "NOT_APPLICABLE"


def test_vireo_placeholder_fails_cleanly(tmp_path):
    mapping, _, soup = _scenario(tmp_path, 0.05)
    (tmp_path / "ph.vcf.gz").write_text("No donor vcf file produced by vireo\n")
    vireo = eg.read_genotypes(tmp_path / "ph.vcf.gz")
    assert vireo == {}
    row = _eval(mapping, vireo, soup)
    assert row["status"] == "FAIL" and "vireo learned no genotype" in row["reason"]


def test_two_missing_donors_fail(tmp_path):
    row = _eval(*_scenario(tmp_path, 0.05), genotyped=["D1", "D2"])
    assert row["status"] == "FAIL" and row["missing_donor"] == "D3,DX"


def test_cross_capture_stability(tmp_path):
    t = _truth()
    paths = {}
    for cap, err in (("c1", 0.0), ("c2", 0.1), ("c3", 0.5)):
        d = tmp_path / cap
        d.mkdir()
        _write_vcf(d / "v.vcf.gz", ["DX"], {"DX": _noisy(t["DX"], err, seed=7)}, gz=True)
        paths[cap] = d / "v.vcf.gz"
    rows = [{"capture_id": c, "status": "PASS", "missing_donor": "DX"} for c in paths]
    rows.append({"capture_id": "c4", "status": "NOT_APPLICABLE", "missing_donor": ""})
    xc = eg.cross_capture(rows, paths)
    lo, n = xc["DX"]
    assert n == 3
    assert lo == pytest.approx(0.5, abs=0.07)  # c3 is the unstable one


def test_tsv_roundtrip(tmp_path):
    row = _eval(*_scenario(tmp_path, 0.05))
    eg.write_tsv([row], tmp_path / "g.tsv", eg.COLUMNS)
    back = eg.read_tsv(tmp_path / "g.tsv")
    assert back[0]["status"] == "PASS"
    assert back[0]["reason"] == "NA"
    assert list(back[0]) == eg.COLUMNS


def test_snakemake_entrypoints(tmp_path, capsys):
    from types import SimpleNamespace as NS

    _scenario(tmp_path, dx_cluster_err=0.55)
    (tmp_path / "genotyped.txt").write_text("\n".join(KNOWN) + "\n")
    sm = NS(
        rule="elimination_gate",
        log=[str(tmp_path / "gate.log")],
        wildcards=NS(capture_id="007"),
        input=NS(
            genotyped_donors=str(tmp_path / "genotyped.txt"),
            mapping=str(tmp_path / "map.tsv"),
            vireo_vcf=str(tmp_path / "vireo.vcf.gz"),
            cluster_genotypes=str(tmp_path / "soup.vcf"),
        ),
        output=NS(tsv=str(tmp_path / "gate.tsv")),
        params=NS(donors=",".join(POOL), min_concordance=0.8, max_drop=0.1, min_sites=100, action="withhold"),
    )
    eg.main(sm)
    row = eg.read_tsv(tmp_path / "gate.tsv")[0]
    assert row["capture_id"] == "007" and row["status"] == "FAIL" and row["action"] == "withhold"
    assert "ELIMINATION GATE FAIL [007]" in (tmp_path / "gate.log").read_text()
    assert "WITHHELD" in capsys.readouterr().err  # shouted on the console too

    cohort = NS(
        rule="combine_elimination_gate",
        log=[str(tmp_path / "cohort.log")],
        input=NS(gates=[str(tmp_path / "gate.tsv")], vireo=[str(tmp_path / "vireo.vcf.gz")]),
        params=NS(capture_ids=["007"]),
        output=NS(tsv=str(tmp_path / "cohort.tsv")),
    )
    eg.main(cohort)
    back = eg.read_tsv(tmp_path / "cohort.tsv")
    assert list(back[0]) == eg.COHORT_COLUMNS
    assert back[0]["cross_capture_min_concordance"] == "NA"  # DX gated in one capture only
