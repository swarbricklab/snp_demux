"""demux_tools.vcf: donor genotypes read as vireo reads them, and vireo's site matching."""
import gzip

import numpy as np
import pytest

from demux_tools import vcf
from demux_tools.vireo_port import parse_gt_code

THIRD = [1 / 3, 1 / 3, 1 / 3]


def _vcf(tmp_path, records, samples=("a", "b"), name="d.vcf.gz", extra_header=()):
    lines = ["##fileformat=VCFv4.2"] + list(extra_header)
    lines.append("\t".join(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] + list(samples)))
    for r in records:
        lines.append("\t".join(r))
    path = tmp_path / name
    data = ("\n".join(lines) + "\n").encode()
    path.write_bytes(gzip.compress(data) if name.endswith(".gz") else data)
    return str(path)


def _rec(chrom, pos, ref, alt, fmt, *fields, info="."):
    return [chrom, str(pos), ".", ref, alt, ".", "PASS", info, fmt] + list(fields)


def test_gt_codes_follow_vireo(tmp_path):
    recs = [_rec("chr1", 10, "A", "G", "GT", "0|1", "1/1"),
            _rec("chr1", 20, "C", "T", "GT", "./.", "."),
            _rec("chr1", 30, "C", "T", "GT", "0", "1"),          # haploid: 0 and 2 ALT alleles
            _rec("chr1", 40, "C", "T", "GT", ".|.", "1|0")]
    g = vcf.read_genotypes(_vcf(tmp_path, recs))
    assert g.variants == ["chr1_10_A_G", "chr1_20_C_T", "chr1_30_C_T", "chr1_40_C_T"]
    assert g.samples == ["a", "b"]
    np.testing.assert_array_equal(g.gp[0], [[0, 1, 0], [0, 0, 1]])
    np.testing.assert_array_equal(g.gp[1], [THIRD, THIRD])
    np.testing.assert_array_equal(g.gp[2], [[1, 0, 0], [0, 0, 1]])
    np.testing.assert_array_equal(g.gp[3], [THIRD, [0, 1, 0]])
    np.testing.assert_array_equal(g.missing, [[False, False], [True, True], [False, False], [True, False]])
    assert g.counts == {"records": 4, "snps": 4, "tagged": 4, "kept": 4}
    np.testing.assert_array_equal(g.rows, [0, 1, 2, 3])


def test_gp_and_pl(tmp_path):
    recs = [_rec("1", 5, "A", "C", "GP:PL", "0.1,0.2,0.7:0,10,100", ".:."),
            _rec("1", 6, "A", "C", "GP:PL", "0.2,0.2,0.2:30,0,30", "0,0,1:255,255,0")]
    path = _vcf(tmp_path, recs)
    gp = vcf.read_genotypes(path, tag="GP")
    np.testing.assert_allclose(gp.gp[0, 0], [0.1, 0.2, 0.7])
    np.testing.assert_allclose(gp.gp[1, 0], [1 / 3] * 3)          # divided by its sum
    np.testing.assert_array_equal(gp.missing, [[False, True], [False, False]])
    pl = vcf.read_genotypes(path, tag="PL")
    expected = 10 ** (-0.1 * np.array([0.0, 10, 100]) - 0.025)
    np.testing.assert_allclose(pl.gp[0, 0], expected / expected.sum())
    for i, code in ((0, "0,10,100"), (1, "30,0,30")):
        np.testing.assert_array_equal(pl.gp[i, 0], parse_gt_code(code, "PL"))


def test_indels_and_multiallelic_records_are_dropped_but_alt_dot_kept(tmp_path):
    recs = [_rec("1", 5, "AT", "A", "GT", "0/1", "0/0"),
            _rec("1", 6, "A", "C,G", "GT", "0/1", "0/0"),
            _rec("1", 7, "A", ".", "GT", "0/0", "0/0"),     # vireo keeps it: both fields one character
            _rec("1", 8, "A", "G", "GT", "0/1", "1/1")]
    g = vcf.read_genotypes(_vcf(tmp_path, recs))
    assert g.variants == ["1_7_A_.", "1_8_A_G"]
    assert g.counts["records"] == 4 and g.counts["snps"] == 2
    np.testing.assert_array_equal(g.rows, [2, 3])


def test_multi_key_format_and_a_record_without_the_tag(tmp_path):
    recs = [_rec("1", 5, "A", "G", "GT:DS:GP", "0|1:1.0:0.05,0.9,0.05", "1|1:2:0,0.01,0.99"),
            _rec("1", 6, "A", "G", "GT:DS", "0|0:0", "0|1:1")]
    path = _vcf(tmp_path, recs)
    g = vcf.read_genotypes(path, tag="GT")
    np.testing.assert_array_equal(g.gp[:, :, :].argmax(axis=2), [[1, 2], [0, 1]])
    p = vcf.read_genotypes(path, tag="GP")
    np.testing.assert_allclose(p.gp[0, 0], [0.05, 0.9, 0.05])
    np.testing.assert_array_equal(p.gp[1], [THIRD, THIRD])      # no GP in FORMAT: blank for everyone
    assert p.counts["tagged"] == 1 and p.missing[1].all()


def test_donor_selection_and_order(tmp_path):
    path = _vcf(tmp_path, [_rec("1", 5, "A", "G", "GT", "0/0", "0/1", "1/1")], samples=("x", "y", "z"))
    g = vcf.read_genotypes(path, samples=["z", "x"])
    assert g.samples == ["z", "x"] and g.header_samples == ["x", "y", "z"]
    np.testing.assert_array_equal(g.gp[0].argmax(axis=1), [2, 0])
    with pytest.raises(ValueError, match="no sample w"):
        vcf.read_genotypes(path, samples=["w"])


def test_too_few_records_with_the_tag_is_vireos_error(tmp_path):
    recs = [_rec("1", p, "A", "G", "GT", "0/1", "0/0") for p in range(1, 11)]
    path = _vcf(tmp_path, recs)
    g = vcf.read_genotypes(path, tag="GP")
    assert g.counts["tagged"] == 0
    with pytest.raises(ValueError, match="fewer than 10 per cent"):
        vcf.check_tag(g, path)
    vcf.check_tag(vcf.read_genotypes(path, tag="GT"))


@pytest.mark.parametrize("fmt,field,tag", [("GT", "0/.", "GT"), ("GP", "0.5,0.5", "GP"), ("GT", "2/2", "GT"),
                                           ("GT:GP", "0/1", "GP")])
def test_fields_vireo_cannot_read_raise(tmp_path, fmt, field, tag):
    path = _vcf(tmp_path, [_rec("1", 5, "A", "G", fmt, field, field)])
    with pytest.raises(ValueError):
        vcf.read_genotypes(path, tag=tag)


def test_keep_ids_parses_only_those_records_and_counts_all(tmp_path):
    recs = [_rec("1", p, "A", "G", "GT", "0/1", "0/0") for p in range(1, 6)]
    g = vcf.read_genotypes(_vcf(tmp_path, recs), keep_ids=["1_2_A_G", "1_4_A_G", "9_9_A_G"], chunksize=2)
    assert g.variants == ["1_2_A_G", "1_4_A_G"]
    np.testing.assert_array_equal(g.rows, [1, 3])
    assert g.counts == {"records": 5, "snps": 5, "tagged": 5, "kept": 2}


def test_info_fields(tmp_path):
    recs = [_rec("1", 5, "A", "G", "GT", "0/1", "0/0", info="AF=0.25;R2=0.9;IMPUTED"),
            _rec("1", 6, "A", "G", "GT", "0/1", "0/0", info="MAF=0.1;AF=."),
            _rec("1", 7, "A", "G", "GT", "0/1", "0/0", info=".")]
    g = vcf.read_genotypes(_vcf(tmp_path, recs, name="plain.vcf"), info=("AF", "R2", "IMPUTED"))
    np.testing.assert_array_equal(g.info["AF"], [0.25, np.nan, np.nan])
    np.testing.assert_array_equal(g.info["R2"], [0.9, np.nan, np.nan])
    np.testing.assert_array_equal(g.info["IMPUTED"], [1.0, np.nan, np.nan])


def test_match_uses_one_copy_of_a_repeated_pileup_row():
    # vireo matches the copy that numpy's sort of the whole list puts first; with so few ids the sort is
    # stable, so here that is the first copy
    cells = ["c_1_A_G", "c_2_A_G", "c_1_A_G", "c_3_A_G"]
    donors = ["c_3_A_G", "c_1_A_G", "c_9_A_G"]
    idx_cell, idx_donor, how = vcf.match_sites(cells, donors)
    np.testing.assert_array_equal(idx_cell, [0, 3])
    np.testing.assert_array_equal(idx_donor, [1, 0])
    assert how == "as is"


def test_match_chr_fallbacks_both_ways_and_no_match():
    idx_cell, idx_donor, how = vcf.match_sites(["1_5_A_G", "2_6_C_T"], ["chr2_6_C_T"])
    assert how == "chr added to pileup" and idx_cell.tolist() == [1] and idx_donor.tolist() == [0]
    idx_cell, idx_donor, how = vcf.match_sites(["chr1_5_A_G", "chr2_6_C_T"], ["1_5_A_G"])
    assert how == "chr added to donors" and idx_cell.tolist() == [0] and idx_donor.tolist() == [0]
    idx_cell, _, how = vcf.match_sites(["1_5_A_G"], ["1_5_A_C"])
    assert how == "none" and len(idx_cell) == 0


def _bgzf(data):
    """BGZF bytes for `data`, as bgzip writes them: gzip members of at most 65,280 bytes of input, each with a
    BC extra field holding its size, then an empty end-of-file member."""
    import struct
    import zlib

    out = b""
    for chunk in [data[i:i + 65280] for i in range(0, len(data), 65280)] + [b""]:
        comp = zlib.compressobj(6, zlib.DEFLATED, -15)
        body = comp.compress(chunk) + comp.flush()
        header = struct.pack("<BBBBIBBH", 31, 139, 8, 4, 0, 0, 255, 6) + b"BC" + struct.pack("<HH", 2, len(body) + 25)
        out += header + body + struct.pack("<II", zlib.crc32(chunk) & 0xffffffff, len(chunk))
    return out


def test_bgzf_input_is_read_like_gzip(tmp_path):
    recs = [_rec("chr%d" % (1 + k % 2), 100 + k, "A", "G", "GT", "0/1", "1/1") for k in range(4000)]
    plain = _vcf(tmp_path, recs, name="plain.vcf")
    data = open(plain, "rb").read()
    (tmp_path / "b.vcf.gz").write_bytes(_bgzf(data))
    assert len(data) > 65280 * 2                          # several BGZF blocks
    a, b = vcf.read_genotypes(plain), vcf.read_genotypes(str(tmp_path / "b.vcf.gz"), chunksize=1000)
    assert a.variants == b.variants and np.array_equal(a.gp, b.gp)
    ids, rows = vcf.snp_record_ids(str(tmp_path / "b.vcf.gz"), chunksize=999)
    assert ids == a.variants and rows.tolist() == list(range(4000))


def test_a_vcf_without_samples_gives_sites_and_info(tmp_path):
    text = ("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
            "chr1\t100\t.\tA\tG\t.\tPASS\tAF=0.2\nchr1\t200\t.\tC\tT\t.\tPASS\tAF=0.4\n"
            "chr1\t300\t.\tAC\tT\t.\tPASS\tAF=0.1\n")
    (tmp_path / "s.vcf").write_text(text)
    g = vcf.read_genotypes(str(tmp_path / "s.vcf"), info=("AF",))
    assert g.variants == ["chr1_100_A_G", "chr1_200_C_T"] and g.samples == [] and g.gp.shape == (2, 0, 3)
    np.testing.assert_array_equal(g.info["AF"], [0.2, 0.4])
    assert g.counts == {"records": 3, "snps": 2, "tagged": 0, "kept": 2}


def test_keep_rows_selects_records_by_line(tmp_path):
    recs = [_rec("1", p, "A", "G", "GT", "0/1", "0/0") for p in range(1, 8)]
    recs.insert(3, _rec("1", 50, "AT", "G", "GT", "0/1", "0/0"))           # line 3: not a SNP
    path = _vcf(tmp_path, recs)
    ids, rows = vcf.snp_record_ids(path)
    assert rows.tolist() == [0, 1, 2, 4, 5, 6, 7] and ids[3] == "1_4_A_G"
    g = vcf.read_genotypes(path, keep_rows=[7, 1, 3], chunksize=3)
    assert g.variants == ["1_2_A_G", "1_7_A_G"] and g.rows.tolist() == [1, 7]


@pytest.mark.vireo
def test_reader_equals_vireos_parser(synthetic_capture):
    from vireoSNP.utils.vcf_utils import load_VCF, parse_donor_GPb

    cap = synthetic_capture(n_donors=5, n_sites=200, n_cells=10, blank_rate=0.2, tags=("GT", "GP", "PL"),
                            gp_noise=0.13, seed=4)
    for tag in ("GT", "GP", "PL"):
        ours = vcf.read_genotypes(cap["vcf"], tag=tag, chunksize=37)
        theirs = load_VCF(cap["vcf"], biallelic_only=True, sparse=False, format_list=[tag])
        assert ours.variants == theirs["variants"] and ours.samples == theirs["samples"]
        expected = parse_donor_GPb(theirs["GenoINFO"][tag], tag)
        assert np.array_equal(ours.gp, expected), tag
