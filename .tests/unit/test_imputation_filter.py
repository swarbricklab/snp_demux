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


def _write(tmp_path, flags=FLAGS):
    p = tmp_path / "panel.vcf"
    p.write_text(HEADER.format(flags=flags) + RECORDS)
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
    # 0.8 keeps a, c, d, e (b and f are imputed with R2 < 0.8)
    assert rows[(0.8, "D1")][3:] == [4, 2, 3]  # D1 differs at a, c, d; missing at e
    assert rows[(0.8, "D2")][3:] == [4, 2, 4]
    assert rows[(0.8, "D3")][3:] == [4, 2, 4]
    # 0.5 additionally keeps b (R2=0.5, not < 0.5) and f (all het: no one differs)
    assert rows[(0.5, "D1")][3:] == [6, 4, 4]
    assert rows[(0.5, "D2")][3:] == [6, 4, 5]


def test_sweep_restricts_to_pool_donors(tmp_path):
    donors, _ = r2.read_panel(_write(tmp_path), donors={"D1", "D3"})
    assert donors == ["D1", "D3"]


def test_without_flags_every_site_is_imputed(tmp_path):
    _, sites = r2.read_panel(_write(tmp_path, flags=""))
    assert all(s[0] for s in sites)
    kept = [s for s in sites if r2.retained(s, 0.8)]
    assert len(kept) == 3  # a, c, e (d is now R2-filtered too)


bcftools = pytest.mark.skipif(
    not (shutil.which("bcftools") and shutil.which("bgzip")),
    reason="bcftools/bgzip not installed",
)


def _filter(tmp_path, min_r2, min_maf, flags=FLAGS):
    vcf = _write(tmp_path, flags)
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
    "min_r2,min_maf,flags,expected",
    [
        (0.8, 0, FLAGS, ["a", "c", "d", "e"]),
        (0.5, 0, FLAGS, ["a", "b", "c", "d", "e", "f"]),
        (0.8, 0.05, FLAGS, ["a", "d", "e"]),
        (0, 0, FLAGS, ["a", "b", "c", "d", "e", "f"]),
        (0.8, 0, "", ["a", "c", "e"]),  # no TYPED/IMPUTED: all sites R2-filtered
    ],
)
def test_filter_panel_sh(tmp_path, min_r2, min_maf, flags, expected):
    res, out = _filter(tmp_path, min_r2, min_maf, flags)
    assert res.returncode == 0, res.stderr
    assert _ids_gz(out) == expected
    assert Path(f"{out}.csi").exists()


@bcftools
def test_filter_matches_sweep_rule(tmp_path):
    """filter_panel.sh and r2_sweep.retained() must agree on what survives."""
    for t in (0.3, 0.5, 0.8, 0.9):
        d = tmp_path / str(t)
        d.mkdir()
        res, out = _filter(d, t, 0)
        assert res.returncode == 0, res.stderr
        ids = [line.split("\t")[2] for line in RECORDS.splitlines()]
        _, sites = r2.read_panel(_write(d))
        expected = [i for i, s in zip(ids, sites) if r2.retained(s, t)]
        assert _ids_gz(out) == expected


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
