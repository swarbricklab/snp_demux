import shutil
import subprocess
from pathlib import Path

import pytest

from conftest import SCRIPTS, load_script

r2 = load_script("r2_sweep")

HEADER = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=1000000>
##INFO=<ID=R2,Number=1,Type=Float,Description="r2">
{flags}##INFO=<ID=MAF,Number=1,Type=Float,Description="maf">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tD1\tD2\tD3
"""
FLAGS = (
    '##INFO=<ID=TYPED,Number=0,Type=Flag,Description="typed">\n'
    '##INFO=<ID=IMPUTED,Number=0,Type=Flag,Description="imputed">\n'
)
# id: what it exercises
#  a imputed, high R2       -> kept
#  b imputed, low R2        -> dropped at 0.8
#  c typed, R2=1, low MAF   -> kept (dropped only by a MAF cut)
#  d typed, low R2          -> kept: genotyped sites are never R2-filtered
#  e imputed, no R2 value   -> kept
#  f imputed, R2=0.6, all donors identical -> kept at 0.5, never discriminating
RECORDS = """chr1\t100\ta\tA\tG\t.\tPASS\tR2=0.95;IMPUTED;MAF=0.2\tGT\t0|0\t0|1\t1|1
chr1\t200\tb\tA\tG\t.\tPASS\tR2=0.5;IMPUTED;MAF=0.2\tGT\t0|0\t0|1\t0|0
chr1\t300\tc\tA\tG\t.\tPASS\tR2=1;TYPED;MAF=0.01\tGT\t0|1\t0|1\t0|0
chr1\t400\td\tA\tG\t.\tPASS\tR2=0.3;TYPED;MAF=0.3\tGT\t0|1\t1|1\t0|0
chr1\t500\te\tA\tG\t.\tPASS\tIMPUTED;MAF=0.3\tGT\t./.\t0|1\t0|0
chr1\t600\tf\tA\tG\t.\tPASS\tR2=0.6;IMPUTED;MAF=0.3\tGT\t0|1\t1|0\t0|1
"""
# The same records the Minimac3 way: genotyped sites are marked in FILTER
# (GENOTYPED, or GENOTYPED_ONLY for sites not in the reference), not in INFO.
FLAGS_MINIMAC3 = (
    '##FILTER=<ID=GENOTYPED,Description="genotyped">\n'
    '##FILTER=<ID=GENOTYPED_ONLY,Description="genotyped only">\n'
)
RECORDS_MINIMAC3 = (
    RECORDS.replace(";IMPUTED", "")
    .replace("\tIMPUTED;", "\t")
    .replace("PASS\tR2=1;TYPED", "GENOTYPED\tR2=1")
    .replace("PASS\tR2=0.3;TYPED", "GENOTYPED_ONLY\tR2=0.3")
)


def _write(tmp_path, flags=FLAGS, records=RECORDS):
    p = tmp_path / "panel.vcf"
    p.write_text(HEADER.format(flags=flags) + records)
    return p


def _ids(path):
    out = []
    for line in Path(path).read_text().splitlines():
        if not line.startswith("#"):
            out.append(line.split("\t")[2])
    return out


def test_sweep_counts(tmp_path):
    donors, sites = r2.read_panel(_write(tmp_path))
    assert donors == ["D1", "D2", "D3"]
    rows = {(r[1], r[2]): r for r in r2.sweep(donors, sites, [0.5, 0.8], "cap")}
    # 0.8 keeps a, c, d, e (b and f are imputed with R2 < 0.8). Pairwise sites
    # differing: D1-D2 at a, d (2); D1-D3 at a, c, d (3; D1 missing at e);
    # D2-D3 at a, c, d, e (4).
    assert rows[(0.8, "D1")][3:] == [4, 2, 2, "D2"]
    assert rows[(0.8, "D2")][3:] == [4, 2, 2, "D1"]
    assert rows[(0.8, "D3")][3:] == [4, 2, 3, "D1"]
    # 0.5 additionally keeps b (R2=0.5, not < 0.5), where D2 differs from D1 and
    # D3, and f (all het: no one differs). D1-D2 3, D1-D3 3, D2-D3 5; a tie goes
    # to the first donor in panel order.
    assert rows[(0.5, "D1")][3:] == [6, 4, 3, "D2"]
    assert rows[(0.5, "D2")][3:] == [6, 4, 3, "D1"]
    assert rows[(0.5, "D3")][3:] == [6, 4, 3, "D1"]


def test_sweep_varies_by_donor():
    """A donor close to another has fewer separating sites than a distant one."""
    donors = ["near1", "near2", "far"]
    # near1 and near2 differ at one site; far differs from both at all three
    sites = [(True, 0.9, [0, 0, 2]), (True, 0.9, [0, 1, 2]), (True, 0.9, [1, 1, 0])]
    rows = {r[2]: r[5:] for r in r2.sweep(donors, sites, [0.8], "cap")}
    assert rows == {"near1": [1, "near2"], "near2": [1, "near1"], "far": [3, "near1"]}


def test_sweep_restricts_to_pool_donors(tmp_path):
    donors, _ = r2.read_panel(_write(tmp_path), donors={"D1", "D3"})
    assert donors == ["D1", "D3"]


def test_minimac3_filter_marks_genotyped_sites(tmp_path):
    _, sites = r2.read_panel(_write(tmp_path, FLAGS_MINIMAC3, RECORDS_MINIMAC3))
    assert [s[0] for s in sites] == [True, True, False, False, True, True]
    kept = [s for s in sites if r2.retained(s, 0.8)]
    assert len(kept) == 4  # a, c, d, e: genotyped c and d are kept


def test_without_flags_every_site_is_imputed(tmp_path):
    _, sites = r2.read_panel(_write(tmp_path, flags=""))
    assert all(s[0] for s in sites)
    kept = [s for s in sites if r2.retained(s, 0.8)]
    assert len(kept) == 3  # a, c, e (d is now R2-filtered too)


bcftools = pytest.mark.skipif(
    not (shutil.which("bcftools") and shutil.which("bgzip")),
    reason="bcftools/bgzip not installed",
)


def _filter(tmp_path, min_r2, min_maf, flags=FLAGS, records=RECORDS):
    vcf = _write(tmp_path, flags, records)
    subprocess.run(["bgzip", "-f", str(vcf)], check=True)
    out = tmp_path / "out.vcf.gz"
    res = subprocess.run(
        ["bash", str(SCRIPTS / "filter_panel.sh"), f"{vcf}.gz", str(out), str(min_r2), str(min_maf)],
        capture_output=True,
        text=True,
    )
    return res, out


def _ids_gz(path):
    out = subprocess.run(["bcftools", "view", "-H", str(path)], capture_output=True, text=True, check=True)
    return [line.split("\t")[2] for line in out.stdout.splitlines()]


@bcftools
@pytest.mark.parametrize(
    "min_r2,min_maf,flags,records,expected",
    [
        (0.8, 0, FLAGS, RECORDS, ["a", "c", "d", "e"]),
        (0.5, 0, FLAGS, RECORDS, ["a", "b", "c", "d", "e", "f"]),
        (0.8, 0.05, FLAGS, RECORDS, ["a", "d", "e"]),  # MAF cut reaches genotyped c
        (0, 0, FLAGS, RECORDS, ["a", "b", "c", "d", "e", "f"]),
        (0.8, 0, "", RECORDS, ["a", "c", "e"]),  # no flags: all sites R2-filtered
        (0.8, 0, FLAGS_MINIMAC3, RECORDS_MINIMAC3, ["a", "c", "d", "e"]),
    ],
)
def test_filter_panel_sh(tmp_path, min_r2, min_maf, flags, records, expected):
    res, out = _filter(tmp_path, min_r2, min_maf, flags, records)
    assert res.returncode == 0, res.stderr
    assert _ids_gz(out) == expected
    assert Path(f"{out}.csi").exists()


@bcftools
def test_filter_matches_sweep_rule(tmp_path):
    """filter_panel.sh and r2_sweep.retained() must agree on what survives."""
    for name, flags, records in (("m4", FLAGS, RECORDS), ("m3", FLAGS_MINIMAC3, RECORDS_MINIMAC3), ("none", "", RECORDS)):
        for t in (0.3, 0.5, 0.8, 0.9):
            d = tmp_path / f"{name}_{t}"
            d.mkdir()
            res, out = _filter(d, t, 0, flags, records)
            assert res.returncode == 0, res.stderr
            ids = [line.split("\t")[2] for line in records.splitlines()]
            _, sites = r2.read_panel(_write(d, flags, records))
            expected = [i for i, s in zip(ids, sites) if r2.retained(s, t)]
            assert _ids_gz(out) == expected, (name, t)


@bcftools
def test_maf_without_header_fails(tmp_path):
    vcf = tmp_path / "p.vcf"
    vcf.write_text(HEADER.format(flags=FLAGS).replace('##INFO=<ID=MAF,Number=1,Type=Float,Description="maf">\n', "") + RECORDS.replace(";MAF=", ";XMAF="))
    subprocess.run(["bgzip", "-f", str(vcf)], check=True)
    res = subprocess.run(
        ["bash", str(SCRIPTS / "filter_panel.sh"), f"{vcf}.gz", str(tmp_path / "o.vcf.gz"), "0.8", "0.05"],
        capture_output=True, text=True,
    )
    assert res.returncode != 0
    assert "declares no INFO/MAF" in res.stderr
