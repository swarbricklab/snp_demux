"""demux_tools.sites: site lists, masks, thinning, GTF regions, site groups and the sites command."""
import gzip

import numpy as np
import pandas as pd
import pytest

from demux_tools import sites


VARIANTS = ["chr1_100_A_G", "chr1_150_C_T", "chr1_400_G_A", "chr1_1000_T_C", "chr2_100_A_C", "chr2_5000_G_T"]


def test_site_lists_by_id_by_position_and_from_a_vcf(tmp_path):
    (tmp_path / "ids.tsv").write_text("variant\textra\nchr1_150_C_T\tx\nchr2_100_A_G\ty\n")
    (tmp_path / "pos.tsv").write_text("chrom\tpos\nchr1\t150\nchr2\t100\n")
    (tmp_path / "full.tsv").write_text("chrom\tpos\tref\talt\nchr1\t150\tC\tT\nchr2\t100\tA\tG\n")
    vcf = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t400\t.\tG\tA\t.\t.\t.\n"
    (tmp_path / "s.vcf.gz").write_bytes(gzip.compress(vcf.encode()))
    assert sites.in_list(VARIANTS, str(tmp_path / "ids.tsv")).tolist() == [0, 1, 0, 0, 0, 0]
    assert sites.in_list(VARIANTS, str(tmp_path / "pos.tsv")).tolist() == [0, 1, 0, 0, 1, 0]
    assert sites.in_list(VARIANTS, str(tmp_path / "full.tsv")).tolist() == [0, 1, 0, 0, 0, 0]
    assert sites.in_list(VARIANTS, str(tmp_path / "s.vcf.gz")).tolist() == [0, 0, 1, 0, 0, 0]
    m = sites.mask(VARIANTS, keep=str(tmp_path / "pos.tsv"), exclude=str(tmp_path / "ids.tsv"))
    assert m.tolist() == [0, 0, 0, 0, 1, 0]
    (tmp_path / "bad.tsv").write_text("x\n1\n")
    with pytest.raises(sites.UsageError):
        sites.read_site_list(str(tmp_path / "bad.tsv"))


def test_thin_first_and_deepest():
    # chr1 sites at 100, 150, 400, 1000; window 100 bp
    assert sites.thin(VARIANTS, 100).tolist() == [1, 0, 1, 1, 1, 1]
    depth = np.array([1, 9, 1, 1, 1, 1])
    assert sites.thin(VARIANTS, 100, keep="deepest", depth=depth).tolist() == [0, 1, 1, 1, 1, 1]
    assert sites.thin(VARIANTS, 300, keep="deepest", depth=depth).tolist() == [0, 1, 0, 1, 1, 1]
    assert sites.thin(VARIANTS, 250).tolist() == [1, 0, 1, 1, 1, 1]
    assert sites.thin(VARIANTS, 300).tolist() == [1, 0, 0, 1, 1, 1]
    with pytest.raises(sites.UsageError):
        sites.thin(VARIANTS, 10, keep="deepest")


GTF = "\n".join([
    "#!genome-build test",
    'chr1\tt\tgene\t90\t500\t.\t+\t.\tgene_id "G1"; gene_name "ONE";',
    'chr1\tt\texon\t90\t120\t.\t+\t.\tgene_id "G1"; gene_name "ONE"; exon_number "1";',
    'chr1\tt\texon\t380\t400\t.\t+\t.\tgene_id "G1"; gene_name "ONE"; exon_number "2";',
    'chr1\tt\tgene\t140\t160\t.\t-\t.\tgene_id "G2"; gene_name "TWO";',
    'chr1\tt\ttranscript\t140\t160\t.\t-\t.\tgene_id "G2"; gene_name "TWO";',
    'chr2\tt\tgene\t4000\t6000\t.\t+\t.\tgene_id "G3";',
]) + "\n"


def test_gtf_regions(tmp_path):
    (tmp_path / "g.gtf").write_text(GTF)
    t = sites.annotate_regions(VARIANTS, str(tmp_path / "g.gtf"))
    assert t["region"].tolist() == ["exonic", "intronic", "exonic", "intergenic", "intergenic", "intronic"]
    assert t["gene_id"].tolist() == ["G1", "G1", "G1", "", "", "G3"]
    assert t["gene_name"].tolist() == ["ONE", "ONE", "ONE", "", "", "G3"]      # no gene_name: the id
    assert t["n_genes"].tolist() == [1, 2, 1, 0, 0, 1]


def test_info_bins_missing_patterns_and_groups():
    assert sites.info_bins([0.1, 0.3, 0.8, 1.0, np.nan, 2], ["0", "0.3", "0.8", "1"]).tolist() == \
        ["[0,0.3)", "[0.3,0.8)", "[0.8,1]", "[0.8,1]", "NA", "NA"]
    missing = np.array([[0, 0, 0], [1, 0, 0], [1, 0, 1], [1, 0, 0], [0, 1, 0]], dtype=bool)
    assert sites.missing_pattern(missing, ["a", "b", "c"]).tolist() == ["none", "a", "a+c", "a", "b"]
    assert sites.missing_pattern(missing, ["a", "b", "c"], top=1).tolist() == ["none", "a", "other", "a", "other"]
    name, labels = sites.groups_from("@chrom", VARIANTS)
    assert name == "chrom" and labels.tolist()[-1] == "chr2"
    table = pd.DataFrame({"variant": VARIANTS[:2], "region": ["exonic", "intronic"]})
    name, labels = sites.groups_from("region", VARIANTS, site_table=table)
    assert labels.tolist() == ["exonic", "intronic", "NA", "NA", "NA", "NA"]
    name, labels = sites.groups_from("@info:R2:0,0.5,1", VARIANTS, info={"R2": np.array([0.2, 0.9] + [np.nan] * 4)})
    assert name == "R2" and labels.tolist()[:3] == ["[0,0.5)", "[0.5,1]", "NA"]
    with pytest.raises(sites.UsageError):
        sites.groups_from("@info:R2:0,1", VARIANTS, info={})


def test_sites_command(synthetic_capture, run_tools, tmp_path):
    cap = synthetic_capture(n_donors=3, n_sites=40, n_cells=6, seed=1)
    base = cap["cellsnp"] + "/cellSNP.base.vcf.gz"
    gtf = tmp_path / "g.gtf"
    gtf.write_text('chr1\tt\tgene\t10000\t20000\t.\t+\t.\tgene_id "G1"; gene_name "ONE";\n'
                   'chr1\tt\texon\t10000\t12000\t.\t+\t.\tgene_id "G1"; gene_name "ONE";\n')
    res = run_tools("sites", "annotate", "--cellsnp", cap["cellsnp"], "--gtf", gtf, "--panel", cap["vcf"],
                    "--info", "AF", "--out", tmp_path / "a.tsv.gz")
    assert res.returncode == 0, res.stderr
    a = pd.read_csv(str(tmp_path / "a.tsv.gz"), sep="\t", dtype=str, keep_default_na=False)
    assert list(a.columns) == ["variant", "chrom", "pos", "ref", "alt", "region", "gene_id", "gene_name", "n_genes",
                               "AF", "in_panel"]
    assert len(a) == 40 and set(a["in_panel"]) == {"true"} and "exonic" in set(a["region"])
    res = run_tools("sites", "thin", "--base", base, "--window-bp", "1500", "--out", tmp_path / "t.tsv")
    assert res.returncode == 0, res.stderr
    t = pd.read_csv(str(tmp_path / "t.tsv"), sep="\t")
    assert 0 < len(t) < 40
    res = run_tools("sites", "thin", "--cellsnp", cap["cellsnp"], "--window-bp", "1500", "--keep", "deepest",
                    "--out", tmp_path / "td.tsv")
    assert res.returncode == 0, res.stderr
    res = run_tools("sites", "from-panel", "--vcf", cap["vcf"], "--min", "AF=0.5", "--snps-only",
                    "--out", tmp_path / "p.tsv")
    assert res.returncode == 0, res.stderr
    p = pd.read_csv(str(tmp_path / "p.tsv"), sep="\t")
    assert (p["AF"] >= 0.5).all() and 0 < len(p) < 40


def test_sites_commands_take_a_vcf_without_samples(run_tools, tmp_path):
    # the review's experiment D: FORMAT and samples absent
    text = ("##fileformat=VCFv4.2\n##INFO=<ID=AF,Number=A,Type=Float,Description=\"x\">\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t.\tA\tG\t.\tPASS\tAF=0.2\nchr1\t200\t.\tC\tT\t.\tPASS\tAF=0.4\n"
            "chr1\t300\t.\tAC\tT\t.\tPASS\tAF=0.1\n")
    (tmp_path / "s.vcf").write_text(text)
    res = run_tools("sites", "from-panel", "--vcf", tmp_path / "s.vcf", "--snps-only", "--min", "AF=0.3",
                    "--out", tmp_path / "p.tsv")
    assert res.returncode == 0, res.stderr
    assert pd.read_csv(str(tmp_path / "p.tsv"), sep="\t")["variant"].tolist() == ["chr1_200_C_T"]
    res = run_tools("sites", "annotate", "--base", tmp_path / "s.vcf", "--panel", tmp_path / "s.vcf", "--info", "AF",
                    "--out", tmp_path / "a.tsv")
    assert res.returncode == 0, res.stderr
    a = pd.read_csv(str(tmp_path / "a.tsv"), sep="\t", dtype=str, keep_default_na=False)
    assert a["in_panel"].tolist() == ["true", "true", "false"] and a["AF"].tolist()[:2] == ["0.2", "0.4"]


def test_gtf_with_other_chromosome_names_warns(tmp_path, capsys):
    (tmp_path / "g.gtf").write_text('1\tt\tgene\t90\t500\t.\t+\t.\tgene_id "G1";\n')
    t = sites.annotate_regions(VARIANTS, str(tmp_path / "g.gtf"))
    assert set(t["region"]) == {"intergenic"}
    assert "no site's chromosome" in capsys.readouterr().err

