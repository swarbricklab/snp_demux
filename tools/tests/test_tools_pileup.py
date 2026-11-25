"""demux_tools.pileup: subset, compare and join pileup folders of both engines."""
import gzip
import json
import os
import shutil

import numpy as np
import pytest

from demux_tools import pileup as P
from demux_tools.cellsnp import Pileup, write_cellsnp
from demux_tools.synth import make_capture


def _souporcell(pile, out_dir, panel_extra="\tGT\t0/1\t1/1"):
    """Write `pile` as a souporcell engine folder: variants.vcf, barcodes.tsv, ref.mtx and alt.mtx.

    Like vartrix: "real" banner, a "% written by sprs" line, space-separated, sorted by row then column,
    and an explicit zero in alt.mtx where a cell has only reference UMIs.
    """
    os.makedirs(out_dir, exist_ok=True)
    dp, ad = pile.to_dense("DP"), pile.to_dense("AD")
    with open(os.path.join(out_dir, "variants.vcf"), "w") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\td0\td1\n")
        for c, p, r, a in zip(pile.chrom, pile.pos, pile.ref, pile.alt):
            f.write("%s\t%d\t.\t%s\t%s\t.\tPASS\tAF=0.5%s\n" % (c, p, r, a, panel_extra))
    with open(os.path.join(out_dir, "barcodes.tsv"), "w") as f:
        f.write("".join(b + "\n" for b in pile.barcodes))
    rows, cols = np.nonzero(dp)
    for name, m in (("ref.mtx", dp - ad), ("alt.mtx", ad)):
        with open(os.path.join(out_dir, name), "w") as f:
            f.write("%%MatrixMarket matrix coordinate real general\n% written by sprs\n")
            f.write("%d %d %d\n" % (pile.n_variants, pile.n_barcodes, len(rows)))
            for r, c in zip(rows, cols):
                f.write("%d %d %d\n" % (r + 1, c + 1, m[r, c]))
    return out_dir


def _capture(tmp_path, name="cap", seed=0, n_sites=30, n_cells=12, prefix="", sites_from=None, shared=0):
    """A synthetic capture in both engines' folders; with sites_from, its first `shared` sites are that pileup's."""
    sc = make_capture(n_donors=3, n_sites=n_sites, n_cells=n_cells, seed=seed, oth_rate=0.2, coverage=0.5)
    pile = sc.pileup
    if prefix or sites_from is not None:
        chrom, pos, ref, alt = pile.chrom.copy(), pile.pos.copy(), pile.ref.copy(), pile.alt.copy()
        if sites_from is not None:
            for arr, src in ((chrom, sites_from.chrom), (pos, sites_from.pos), (ref, sites_from.ref),
                             (alt, sites_from.alt)):
                arr[:shared] = src[:shared]
        pile = Pileup(chrom, pos, ref, alt, [prefix + b for b in pile.barcodes], layers=pile.layers)
    cs = str(tmp_path / name / "cellsnp")
    write_cellsnp(pile, cs)
    sp = _souporcell(pile, str(tmp_path / name / "souporcell"))
    return pile, cs, sp


def _read(path):
    op = gzip.open if open(path, "rb").read(2) == b"\x1f\x8b" else open
    with op(path, "rb") as f:
        return f.read()


def test_occurrence_keys():
    assert P.occurrence_keys([]) == []
    assert P.occurrence_keys(["a", "b", "a"]) == ["a#0", "b#0", "a#1"]


def test_bgzf_is_gzip_with_the_eof_block(tmp_path):
    data = ("x" * 70000 + "\n").encode()
    b = P.bgzf_bytes(data)
    assert gzip.decompress(b) == data and b.endswith(P.BGZF_EOF) and b == P.bgzf_bytes(data)
    assert b[12:16] == b"BC\x02\x00"


def test_site_subset_equals_a_fresh_pileup_of_those_sites(tmp_path):
    pile, cs, sp = _capture(tmp_path)
    keep = np.arange(pile.n_variants) % 3 != 0
    fresh = str(tmp_path / "fresh")
    write_cellsnp(pile.subset(rows=keep), fresh)
    sites = tmp_path / "sites.tsv"
    sites.write_text("CHROM\tPOS\tREF\tALT\n" + "".join(
        "%s\t%d\t%s\t%s\n" % t for t, k in zip(zip(pile.chrom, pile.pos, pile.ref, pile.alt), keep) if k))
    out = str(tmp_path / "sub")
    counts = P.subset(cs, out, sites=P.read_site_list(str(sites)))
    assert counts["rows_out"] == int(keep.sum()) and counts["barcodes_out"] == pile.n_barcodes
    table, verdict, notes = P.compare(out, fresh)
    assert verdict == "same sites and counts", notes
    for f in ("cellSNP.samples.tsv", "cellSNP.tag.AD.mtx", "cellSNP.tag.DP.mtx", "cellSNP.tag.OTH.mtx"):
        assert _read(os.path.join(out, f)) == _read(os.path.join(fresh, f)), f
    assert _read(os.path.join(out, "cellSNP.base.vcf.gz")) == _read(os.path.join(fresh, "cellSNP.base.vcf.gz"))


def test_barcode_subset_with_min_count_equals_a_fresh_pileup(tmp_path):
    pile, cs, sp = _capture(tmp_path)
    cols = np.arange(pile.n_barcodes) % 2 == 0
    sub = pile.subset(cols=cols)
    total = sub.to_dense("DP").sum(1) + sub.to_dense("OTH").sum(1)
    fresh_pile = sub.subset(rows=total >= 4)
    fresh = str(tmp_path / "fresh")
    write_cellsnp(fresh_pile, fresh)
    bl = tmp_path / "bc.tsv.gz"
    bl.write_bytes(gzip.compress("".join(b + "\n" for b in pile.barcodes[cols][::-1]).encode()))
    out = str(tmp_path / "sub")
    counts = P.subset(cs, out, barcodes=P.read_barcode_list(str(bl)), min_count=4)
    assert counts["info"] == "recompute" and counts["rows_out"] == fresh_pile.n_variants
    assert P.compare(out, fresh)[1] == "same sites and counts"
    for f in ("cellSNP.samples.tsv", "cellSNP.tag.DP.mtx", "cellSNP.tag.OTH.mtx"):
        assert _read(os.path.join(out, f)) == _read(os.path.join(fresh, f)), f
    # the INFO sums are those of the kept barcodes, as a fresh pileup writes them
    assert _read(os.path.join(out, "cellSNP.base.vcf.gz")) == _read(os.path.join(fresh, "cellSNP.base.vcf.gz"))
    # columns_of_b does the same cut inside compare
    assert P.compare(cs, fresh, columns_of_b=4)[1] == "same sites and counts"


def test_souporcell_subset_keeps_records_and_rows_aligned(tmp_path):
    pile, cs, sp = _capture(tmp_path)
    keep_keys = set("%s:%d:%s:%s" % t for t in list(zip(pile.chrom, pile.pos, pile.ref, pile.alt))[5:12])
    out = str(tmp_path / "sub")
    P.subset(sp, out, sites=keep_keys, barcodes=list(pile.barcodes[:6]))
    p = P.PileupFolder(out)
    assert p.engine == "souporcell" and len(p.sites) == 7 and len(p.barcodes) == 6
    assert p.sites.keys == ["%s:%d:%s:%s" % t for t in list(zip(pile.chrom, pile.pos, pile.ref, pile.alt))[5:12]]
    assert p.sites.records[0].endswith("GT\t0/1\t1/1")     # genotype columns kept as read
    alt = open(os.path.join(out, "alt.mtx")).read().split("\n")
    assert alt[0].endswith("real general") and alt[1] == "% written by sprs" and " " in alt[2]
    # same counts as souporcell files written from the subset in memory, explicit zeros included
    fresh = _souporcell(pile.subset(rows=np.arange(5, 12), cols=np.arange(6)), str(tmp_path / "fresh"))
    for f in ("alt.mtx", "ref.mtx", "barcodes.tsv", "variants.vcf"):
        assert _read(os.path.join(out, f)) == _read(os.path.join(fresh, f)), f
    with pytest.raises(P.UsageError):
        P.subset(sp, str(tmp_path / "x"), min_count=2)


def test_unchanged_files_are_copied_and_repeated_sites_kept(tmp_path):
    pile, cs, _ = _capture(tmp_path)
    rep = Pileup(np.r_[pile.chrom, pile.chrom[:2]], np.r_[pile.pos, pile.pos[:2]], np.r_[pile.ref, pile.ref[:2]],
                 np.r_[pile.alt, pile.alt[:2]], pile.barcodes,
                 layers={k: (v[0], v[1], v[2]) for k, v in pile.layers.items()})
    d = str(tmp_path / "rep")
    write_cellsnp(rep, d)
    out = str(tmp_path / "same")
    P.subset(d, out)
    for f in os.listdir(d):
        assert open(os.path.join(d, f), "rb").read() == open(os.path.join(out, f), "rb").read(), f
    p = P.PileupFolder(out)
    assert len(p.sites) == pile.n_variants + 2
    key0 = p.sites.keys[0]
    out2 = str(tmp_path / "one")
    P.subset(d, out2, sites={key0})
    assert P.PileupFolder(out2).sites.keys == [key0, key0]   # both copies of the site
    _, verdict, notes = P.compare(d, out)
    assert verdict == "same sites and counts" and any("repeated site rows: A 2, B 2" in n for n in notes)


def test_compare_reports_differences(tmp_path):
    pile, cs, sp = _capture(tmp_path)
    changed = Pileup(pile.chrom, pile.pos, pile.ref, pile.alt, pile.barcodes, layers=dict(pile.layers))
    r, c, v = pile.layers["DP"]
    changed.set_layer("DP", r, c, v + (np.arange(len(v)) == 0))
    d = str(tmp_path / "changed")
    write_cellsnp(changed, d)
    table, verdict, notes = P.compare(cs, d)
    t = table.set_index("layer")
    assert verdict == "differ" and t.at["DP", "equal"] == t.at["DP", "in_both"] - 1
    assert t.at["AD", "equal"] == t.at["AD", "entries_a"]
    fewer = str(tmp_path / "fewer")
    write_cellsnp(pile.subset(rows=np.arange(1, pile.n_variants)), fewer)
    table, verdict, notes = P.compare(cs, fewer)
    assert table.set_index("layer").at["sites", "only_a"] == 1 and "sites only in A 1, only in B 0" in notes
    # at the sites of the smaller pileup the two agree
    sites = P.read_site_list(os.path.join(fewer, "cellSNP.base.vcf.gz"))
    assert P.compare(cs, fewer, sites=sites)[1] == "same sites and counts"
    with pytest.raises(P.UsageError):
        P.compare(cs, sp)


def test_explicit_zeros_equal_missing_entries_unless_kept(tmp_path):
    pile, _, sp = _capture(tmp_path)
    d = str(tmp_path / "nozero")
    os.makedirs(d)
    for f in ("variants.vcf", "barcodes.tsv", "ref.mtx"):
        open(os.path.join(d, f), "wb").write(open(os.path.join(sp, f), "rb").read())
    lines = open(os.path.join(sp, "alt.mtx")).read().split("\n")
    body = [l for l in lines[3:] if l and not l.endswith(" 0")]
    n = lines[2].split()
    open(os.path.join(d, "alt.mtx"), "w").write("\n".join(lines[:2] + ["%s %s %d" % (n[0], n[1], len(body))] + body) + "\n")
    assert P.compare(sp, d)[1] == "same sites and counts"
    assert P.compare(sp, d, keep_zeros=True)[1] == "differ"


@pytest.mark.parametrize("engine", ["cellsnp", "souporcell"])
@pytest.mark.parametrize("how", ["union", "intersection"])
def test_join_then_split_gives_the_originals(tmp_path, engine, how):
    a_pile, a_cs, a_sp = _capture(tmp_path, "a", seed=1, n_sites=25)
    b_pile, b_cs, b_sp = _capture(tmp_path, "b", seed=2, n_sites=20, sites_from=a_pile, shared=15)
    a, b = (a_cs, b_cs) if engine == "cellsnp" else (a_sp, b_sp)
    out = str(tmp_path / "joined")
    counts = P.join([(a, "A_"), (b, "B_")], out, how)
    j = P.PileupFolder(out)
    assert len(j.barcodes) == a_pile.n_barcodes + b_pile.n_barcodes and j.barcodes[0].startswith("A_")
    for src, prefix, pile in ((a, "A_", a_pile), (b, "B_", b_pile)):
        bl = [prefix + x for x in pile.barcodes]
        back = str(tmp_path / ("back_" + prefix))
        sites = P.read_site_list(os.path.join(src, P.ENGINES[engine]["sites"]))
        P.subset(out, back, sites=sites, barcodes=bl, strip_prefix=prefix)
        if how == "union":
            assert P.compare(src, back)[1] == "same sites and counts"
        else:
            common = P.read_site_list(os.path.join(out, P.ENGINES[engine]["sites"]))
            assert P.compare(src, back, sites=common)[1] == "same sites and counts"
    ka, kb = set(P.PileupFolder(a).sites.keys), set(P.PileupFolder(b).sites.keys)
    assert counts["rows"] == (len(ka | kb) if how == "union" else len(ka & kb))
    assert len(ka & kb) >= 15
    if engine == "souporcell":
        assert all(len(r.split("\t")) == 8 for r in j.sites.records)
    with pytest.raises(P.UsageError):
        P.join([(a, ""), (a, "")], str(tmp_path / "clash"), how)   # the same barcodes twice


def test_join_keeps_each_inputs_order_at_one_position(tmp_path):
    # a panel can list two alleles of one position in either order (here G>C before G>A)
    sc = make_capture(n_donors=2, n_sites=4, n_cells=6, seed=3, coverage=0.9, doublets=())
    a = Pileup(["chr1"] * 4, [100, 200, 200, 300], ["G"] * 4, ["T", "C", "A", "G"], sc.pileup.barcodes,
               layers=sc.pileup.layers)
    b = Pileup(["chr1"] * 3, [50, 200, 400], ["G"] * 3, ["T", "A", "G"], ["b%d" % i for i in range(6)],
               layers={k: (r[r < 3], c[r < 3], v[r < 3]) for k, (r, c, v) in sc.pileup.layers.items()})
    da, db = str(tmp_path / "a"), str(tmp_path / "b")
    write_cellsnp(a, da)
    write_cellsnp(b, db)
    out = str(tmp_path / "j")
    P.join([(da, "A_"), (db, "B_")], out, "union")
    keys = P.PileupFolder(out).sites.keys
    assert keys == ["chr1:50:G:T", "chr1:100:G:T", "chr1:200:G:C", "chr1:200:G:A", "chr1:300:G:G", "chr1:400:G:G"]
    back = str(tmp_path / "back")
    P.subset(out, back, sites=P.read_site_list(os.path.join(da, "cellSNP.base.vcf.gz")),
             barcodes=["A_" + x for x in a.barcodes], strip_prefix="A_")
    for f in ("cellSNP.samples.tsv", "cellSNP.tag.AD.mtx", "cellSNP.tag.DP.mtx", "cellSNP.tag.OTH.mtx",
              "cellSNP.base.vcf.gz"):
        assert _read(os.path.join(back, f)) == _read(os.path.join(da, f)), f


def test_commands(tmp_path, run_tools):
    pile, cs, sp = _capture(tmp_path)
    bl = tmp_path / "bc.txt"
    bl.write_text("".join(b + "\n" for b in pile.barcodes[:5]))
    res = run_tools("pileup-subset", cs, tmp_path / "sub", "--barcodes", bl, "--min-count", "2")
    assert res.returncode == 0, res.stderr
    prov = json.load(open(str(tmp_path / "sub" / "provenance.json")))
    assert prov["results"]["barcodes_out"] == 5
    res = run_tools("pileup-compare", cs, tmp_path / "sub", "--columns-of-b", "2", "--check", "--out", tmp_path / "c.tsv")
    assert res.returncode == 0, res.stdout + res.stderr
    assert "same sites and counts" in res.stdout
    res = run_tools("pileup-compare", cs, tmp_path / "sub", "--check", "--out", tmp_path / "c2.tsv")
    assert res.returncode == 1 and "differ" in res.stdout
    res = run_tools("pileup-join", "--in", cs, "--prefix", "x_", "--in", cs, "--prefix", "y_", tmp_path / "j")
    assert res.returncode == 0, res.stderr
    assert len(P.PileupFolder(str(tmp_path / "j")).barcodes) == 2 * pile.n_barcodes
    res = run_tools("pileup-subset", cs, cs)
    assert res.returncode == 2 and "also an input" in res.stderr


# ---------------------------------------------------------------------------
# Review fixes

BANNER = "%%MatrixMarket matrix coordinate integer general\n%\n"
VCF_HEAD = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"


def _hand_cellsnp(d):
    """Three sites, three barcodes, written by hand as cellsnp-lite writes them."""
    os.makedirs(d)
    with open(os.path.join(d, "cellSNP.base.vcf.gz"), "wb") as f:
        f.write(gzip.compress((VCF_HEAD + "chr1\t100\t.\tA\tG\t.\tPASS\tAD=1;DP=3;OTH=0\n"
                               "chr1\t200\t.\tC\tT\t.\tPASS\tAD=2;DP=2;OTH=1\n"
                               "chr2\t50\t.\tG\tA\t.\tPASS\tAD=0;DP=4;OTH=0\n").encode()))
    files = {"cellSNP.samples.tsv": "AAA-1\nCCC-1\nGGG-1\n",
             "cellSNP.tag.AD.mtx": BANNER + "3\t3\t2\n1\t2\t1\n2\t3\t2\n",
             "cellSNP.tag.DP.mtx": BANNER + "3\t3\t5\n1\t1\t1\n1\t2\t2\n2\t3\t2\n3\t1\t3\n3\t3\t1\n",
             "cellSNP.tag.OTH.mtx": BANNER + "3\t3\t1\n2\t1\t1\n"}
    for name, text in files.items():
        with open(os.path.join(d, name), "w") as f:
            f.write(text)
    return d


def test_subset_matches_hand_written_files(tmp_path):
    src = _hand_cellsnp(str(tmp_path / "in"))
    out = str(tmp_path / "out")
    P.subset(src, out, sites={"chr1:100:A:G", "chr2:50:G:A"}, barcodes=["GGG-1", "AAA-1"])
    assert open(os.path.join(out, "cellSNP.samples.tsv")).read() == "AAA-1\nGGG-1\n"
    assert open(os.path.join(out, "cellSNP.tag.DP.mtx")).read() == BANNER + "2\t2\t3\n1\t1\t1\n2\t1\t3\n2\t2\t1\n"
    assert open(os.path.join(out, "cellSNP.tag.AD.mtx")).read() == BANNER + "2\t2\t0\n"
    assert open(os.path.join(out, "cellSNP.tag.OTH.mtx")).read() == BANNER + "2\t2\t0\n"
    assert _read(os.path.join(out, "cellSNP.base.vcf.gz")).decode() == (
        VCF_HEAD + "chr1\t100\t.\tA\tG\t.\tPASS\tAD=0;DP=1;OTH=0\nchr2\t50\t.\tG\tA\t.\tPASS\tAD=0;DP=4;OTH=0\n")
    sp = str(tmp_path / "sp")
    os.makedirs(sp)
    real = "%%MatrixMarket matrix coordinate real general\n% written by sprs\n"
    head = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tx\n"
    recs = ["chr1\t100\t.\tA\tG\t.\tPASS\t.\tGT\t0/1", "chr1\t100\t.\tA\tT\t.\tPASS\t.\tGT\t0/0",
            "chr1\t300\t.\tC\tG\t.\tPASS\t.\tGT\t1/1"]
    for name, text in (("variants.vcf", head + "\n".join(recs) + "\n"), ("barcodes.tsv", "B1\nB2\n"),
                       ("alt.mtx", real + "3 2 4\n1 1 0\n1 2 2\n2 1 1\n3 2 0\n"),
                       ("ref.mtx", real + "3 2 4\n1 1 3\n1 2 0\n2 1 1\n3 2 5\n")):
        with open(os.path.join(sp, name), "w") as f:
            f.write(text)
    out2 = str(tmp_path / "out2")
    P.subset(sp, out2, sites={"chr1:100:A:T", "chr1:300:C:G"}, barcodes=["B2"])
    assert open(os.path.join(out2, "variants.vcf")).read() == head + recs[1] + "\n" + recs[2] + "\n"
    assert open(os.path.join(out2, "barcodes.tsv")).read() == "B2\n"
    assert open(os.path.join(out2, "alt.mtx")).read() == real + "2 1 1\n2 1 0\n"
    assert open(os.path.join(out2, "ref.mtx")).read() == real + "2 1 1\n2 1 5\n"


def test_a_different_barcode_list_differs(tmp_path, run_tools):
    src = _hand_cellsnp(str(tmp_path / "a"))
    b = str(tmp_path / "b")
    shutil.copytree(src, b)
    with open(os.path.join(b, "cellSNP.samples.tsv"), "a") as f:
        f.write("TTT-1\n")                 # one more barcode, with no entries
    for name in ("cellSNP.tag.AD.mtx", "cellSNP.tag.DP.mtx", "cellSNP.tag.OTH.mtx"):
        p = os.path.join(b, name)
        text = open(p).read().replace("3\t3\t", "3\t4\t", 1)
        open(p, "w").write(text)
    table, verdict, notes = P.compare(src, b)
    assert verdict == "differ" and "barcodes only in A 0, only in B 1" in notes
    res = run_tools("pileup-compare", src, b, "--check", "--out", tmp_path / "c.tsv")
    assert res.returncode == 1 and "differ" in res.stdout


def test_links_in_the_output_folder_are_replaced_not_written_through(tmp_path, run_tools):
    src = _hand_cellsnp(str(tmp_path / "in"))
    out = tmp_path / "out"
    out.mkdir()
    targets = tmp_path / "cache"
    targets.mkdir()
    (targets / "w").write_text("cache object\n")
    (targets / "r").write_text("read-only cache object\n")
    os.chmod(str(targets / "r"), 0o444)
    os.symlink(str(targets / "w"), str(out / "cellSNP.samples.tsv"))     # copied unchanged in this run
    os.symlink(str(targets / "r"), str(out / "cellSNP.tag.AD.mtx"))
    res = run_tools("pileup-subset", src, out)
    assert res.returncode == 2 and "--force" in res.stderr           # a folder with files is refused
    res = run_tools("pileup-subset", src, out, "--force")
    assert res.returncode == 0, res.stderr
    for name in ("cellSNP.samples.tsv", "cellSNP.tag.AD.mtx"):
        assert not os.path.islink(str(out / name))
        assert open(str(out / name)).read() == open(os.path.join(src, name)).read()
    assert (targets / "w").read_text() == "cache object\n"
    assert (targets / "r").read_text() == "read-only cache object\n"


def test_force_removes_files_the_run_would_leave_stale(tmp_path):
    src = _hand_cellsnp(str(tmp_path / "in"))
    out = str(tmp_path / "out")
    P.subset(src, out)
    no_oth = str(tmp_path / "no_oth")
    shutil.copytree(src, no_oth)
    os.remove(os.path.join(no_oth, "cellSNP.tag.OTH.mtx"))
    with open(os.path.join(out, "notes.txt"), "w") as f:
        f.write("not a pileup file\n")
    with pytest.raises(P.UsageError):
        P.subset(no_oth, out)
    P.subset(no_oth, out, force=True)
    assert not os.path.exists(os.path.join(out, "cellSNP.tag.OTH.mtx"))
    assert os.path.exists(os.path.join(out, "notes.txt"))


def test_subset_reports_requested_and_found_and_refuses_an_empty_result(tmp_path, run_tools):
    src = _hand_cellsnp(str(tmp_path / "in"))
    counts = P.subset(src, str(tmp_path / "o1"), sites={"chr1:100:A:G", "1:100:A:G"},
                      barcodes=["AAA-1", "AAA"])
    assert (counts["sites_requested"], counts["sites_found"]) == (2, 1)
    assert (counts["barcodes_requested"], counts["barcodes_found"]) == (2, 1)
    with pytest.raises(P.UsageError):
        P.subset(src, str(tmp_path / "o2"), sites={"1:100:A:G"})
    assert not os.path.exists(str(tmp_path / "o2"))
    counts = P.subset(src, str(tmp_path / "o3"), sites={"1:100:A:G"}, allow_empty=True)
    assert counts["rows_out"] == 0
    sites = tmp_path / "sites.tsv"
    sites.write_text("1\t100\tA\tG\n")
    res = run_tools("pileup-subset", src, tmp_path / "o4", "--sites", sites)
    assert res.returncode == 2 and "sites requested 1, found 0" in res.stderr and "--allow-empty" in res.stderr
    assert "None" not in res.stderr
    res = run_tools("pileup-subset", src, tmp_path / "o5", "--barcodes", tmp_path / "missing.txt")
    assert res.returncode == 2 and "Traceback" not in res.stderr


def test_join_reports_rows_and_refuses_a_missing_layer(tmp_path, run_tools):
    a = _hand_cellsnp(str(tmp_path / "a"))
    b = str(tmp_path / "b")
    os.makedirs(b)
    with open(os.path.join(b, "cellSNP.base.vcf.gz"), "wb") as f:
        f.write(gzip.compress((VCF_HEAD + "chr1\t100\t.\tA\tG\t.\tPASS\tAD=1;DP=1\n"
                               "chr3\t10\t.\tT\tC\t.\tPASS\tAD=0;DP=1\n").encode()))
    for name, text in (("cellSNP.samples.tsv", "AAA-1\n"), ("cellSNP.tag.AD.mtx", BANNER + "2\t1\t1\n1\t1\t1\n"),
                       ("cellSNP.tag.DP.mtx", BANNER + "2\t1\t2\n1\t1\t1\n2\t1\t1\n")):
        with open(os.path.join(b, name), "w") as f:
            f.write(text)
    with pytest.raises(P.UsageError):
        P.join([(a, "A_"), (b, "B_")], str(tmp_path / "j0"))           # b lacks OTH
    counts = P.join([(a, "A_"), (b, "B_")], str(tmp_path / "j"), layers=["AD", "DP"])
    assert counts["input_rows"] == [3, 2] and counts["rows_in_every_input"] == 1
    assert counts["rows_only_in_input"] == [2, 1] and counts["rows"] == 4
    assert not os.path.exists(str(tmp_path / "j" / "cellSNP.tag.OTH.mtx"))
    res = run_tools("pileup-join", "--in", a, "--prefix", "A_", "--in", b, "--prefix", "B_", "--layers", "AD,DP",
                    tmp_path / "j2")
    assert res.returncode == 0, res.stderr
    assert "rows in every input: 1" in res.stdout and "2 of them in no other input" in res.stdout
    prov = json.load(open(str(tmp_path / "j2" / "provenance.json")))
    assert prov["results"]["rows_only_in_input"] == [2, 1]


def test_columns_of_b_needs_dp_and_oth(tmp_path, run_tools):
    a = _hand_cellsnp(str(tmp_path / "a"))
    res = run_tools("pileup-compare", a, a, "--columns-of-b", "2", "--layers", "AD", "--out", tmp_path / "c.tsv")
    assert res.returncode == 2 and "Traceback" not in res.stderr and "DP" in res.stderr
