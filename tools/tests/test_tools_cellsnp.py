"""demux_tools.cellsnp: MatrixMarket reading with a row filter, the cellsnp folder, the writer."""
import gzip
import os

import numpy as np
import pytest

from demux_tools import cellsnp
from demux_tools.cellsnp import Pileup


def _write(path, text):
    with open(str(path), "w") as f:
        f.write(text)
    return str(path)


MTX = """%%MatrixMarket matrix coordinate integer general
%
5\t3\t7
1\t1\t2
1\t3\t1
2\t2\t4
3\t1\t1
4\t3\t9
5\t1\t3
5\t2\t1
"""


def _full(path, **kw):
    shape, i, j, v = cellsnp.read_mtx(path, **kw)
    return shape, list(zip(i.tolist(), j.tolist(), v.tolist()))


def test_read_mtx_is_zero_based_and_typed(tmp_path):
    p = _write(tmp_path / "m.mtx", MTX)
    shape, i, j, v = cellsnp.read_mtx(p)
    assert shape == (5, 3)
    assert i.dtype == j.dtype == v.dtype == np.int32
    assert list(zip(i, j, v))[:2] == [(0, 0, 2), (0, 2, 1)]


@pytest.mark.parametrize("chunksize", [1, 2, 3, 100])
def test_read_mtx_keep_rows_equals_full_read_then_filter(tmp_path, chunksize):
    p = _write(tmp_path / "m.mtx", MTX)
    _, full = _full(p)
    keep = [4, 0, 2, 2]   # a set: order and repeats ignored
    shape, kept = _full(p, keep_rows=keep, chunksize=chunksize)
    assert shape == (5, 3)
    assert kept == [e for e in full if e[0] in set(keep)]
    mask = np.zeros(5, dtype=bool)
    mask[[0, 2, 4]] = True
    assert _full(p, keep_rows=mask, chunksize=chunksize)[1] == kept
    assert _full(p, keep_rows=[], chunksize=chunksize)[1] == []


def test_read_mtx_on_a_synthetic_matrix_with_many_chunks(tmp_path):
    rng = np.random.RandomState(1)
    dense = rng.binomial(3, 0.2, size=(200, 30))
    r, c = np.nonzero(dense)
    path = cellsnp.write_mtx(str(tmp_path / "big.mtx.gz"), dense.shape, r, c, dense[r, c])
    keep = rng.random_sample(200) < 0.3
    shape, i, j, v = cellsnp.read_mtx(path, keep_rows=keep, chunksize=37)
    got = np.zeros_like(dense)
    got[i, j] = v
    assert shape == dense.shape
    assert np.array_equal(got, dense * keep[:, None])


def test_read_mtx_rejects_bad_files(tmp_path):
    short = _write(tmp_path / "short.mtx", MTX.replace("5\t3\t7", "5\t3\t8"))
    with pytest.raises(ValueError, match="8 entries"):
        cellsnp.read_mtx(short)
    outside = _write(tmp_path / "outside.mtx", MTX.replace("4\t3\t9", "6\t3\t9"))
    with pytest.raises(ValueError, match="outside"):
        cellsnp.read_mtx(outside)
    real = _write(tmp_path / "real.mtx", MTX.replace("integer", "real"))
    with pytest.raises(ValueError, match="expected"):
        cellsnp.read_mtx(real)
    with pytest.raises(ValueError):
        cellsnp.read_mtx(_write(tmp_path / "m.mtx", MTX), keep_rows=[5])


def test_read_mtx_empty_matrix(tmp_path):
    p = _write(tmp_path / "e.mtx", "%%MatrixMarket matrix coordinate integer general\n%\n4\t2\t0\n")
    shape, i, j, v = cellsnp.read_mtx(p, keep_rows=[1])
    assert shape == (4, 2) and len(i) == len(j) == len(v) == 0


def _tiny_pileup():
    # Row 1 repeats row 0's site, as cellsnp-lite does for a panel with a repeated record.
    layers = {
        "AD": ([0, 1, 2, 0], [1, 1, 0, 0], [1, 1, 2, 3]),
        "DP": ([2, 0, 1, 0, 3], [0, 1, 1, 0, 1], [2, 2, 2, 5, 1]),
        "OTH": ([3], [1], [1]),
    }
    return Pileup(["chr1", "chr1", "chr2", "chr2"], [100, 100, 50, 60], ["A", "A", "C", "G"],
                  ["G", "G", "T", "T"], ["AAAC-1", "TTTG-1"], layers=layers)


def test_write_cellsnp_writes_cellsnp_text(tmp_path):
    paths = cellsnp.write_cellsnp(_tiny_pileup(), str(tmp_path / "p"))
    assert sorted(os.listdir(str(tmp_path / "p"))) == sorted(cellsnp.FILES.values())
    with open(paths["DP"]) as f:
        assert f.read() == ("%%MatrixMarket matrix coordinate integer general\n%\n4\t2\t5\n"
                            "1\t1\t5\n1\t2\t2\n2\t2\t2\n3\t1\t2\n4\t2\t1\n")
    with open(paths["samples"]) as f:
        assert f.read() == "AAAC-1\nTTTG-1\n"
    base = gzip.decompress(open(paths["base"], "rb").read()).decode()
    assert base.splitlines()[2:4] == ["chr1\t100\t.\tA\tG\t.\tPASS\tAD=4;DP=7;OTH=0",
                                      "chr1\t100\t.\tA\tG\t.\tPASS\tAD=1;DP=2;OTH=0"]


def test_round_trip_keeps_repeated_sites_and_counts(tmp_path):
    p = _tiny_pileup()
    cellsnp.write_cellsnp(p, str(tmp_path / "p"))
    q = cellsnp.read_cellsnp(str(tmp_path / "p"), layers=("AD", "DP", "OTH"))
    assert q.variants == ["chr1_100_A_G", "chr1_100_A_G", "chr2_50_C_T", "chr2_60_G_T"]
    assert list(q.barcodes) == ["AAAC-1", "TTTG-1"]
    for name in ("AD", "DP", "OTH"):
        assert np.array_equal(q.to_dense(name), p.to_dense(name))
    assert q.rows.tolist() == [0, 1, 2, 3]


def test_read_cellsnp_keep_rows_renumbers_and_maps_back(tmp_path):
    p = _tiny_pileup()
    cellsnp.write_cellsnp(p, str(tmp_path / "p"))
    q = cellsnp.read_cellsnp(str(tmp_path / "p"), keep_rows=[3, 1])
    assert q.rows.tolist() == [1, 3]
    assert q.variants == ["chr1_100_A_G", "chr2_60_G_T"]
    assert np.array_equal(q.to_dense("DP"), p.to_dense("DP")[[1, 3]])
    assert "OTH" not in q.layers


def test_explicit_files_override_the_folder(tmp_path):
    p = _tiny_pileup()
    paths = cellsnp.write_cellsnp(p, str(tmp_path / "p"))
    other = tmp_path / "object_named_by_md5"
    other.write_bytes(open(paths["AD"], "rb").read())
    resolved = cellsnp.resolve_cellsnp(base=paths["base"], samples=paths["samples"], ad=str(other), dp=paths["DP"])
    assert resolved["AD"] == str(other) and resolved["OTH"] is None
    q = cellsnp.read_cellsnp(resolved)
    assert np.array_equal(q.to_dense("AD"), p.to_dense("AD"))
    with pytest.raises(FileNotFoundError):
        cellsnp.resolve_cellsnp(cellsnp=str(tmp_path / "missing"))
    with pytest.raises(FileNotFoundError):
        cellsnp.read_cellsnp(resolved, layers=("OTH",))


def test_reader_checks_sizes_against_base_and_barcodes(tmp_path):
    paths = cellsnp.write_cellsnp(_tiny_pileup(), str(tmp_path / "p"))
    with open(paths["samples"], "a") as f:
        f.write("EXTRA-1\n")
    with pytest.raises(ValueError, match="barcode list 3"):
        cellsnp.read_cellsnp(str(tmp_path / "p"))


def test_writer_sorts_drops_zeros_and_refuses_repeats(tmp_path):
    path = cellsnp.write_mtx(str(tmp_path / "m.mtx"), (3, 3), [2, 0, 1, 0], [0, 2, 1, 1], [1, 5, 0, 2])
    with open(path) as f:
        assert f.read().splitlines()[2:] == ["3\t3\t3", "1\t2\t2", "1\t3\t5", "3\t1\t1"]
    with pytest.raises(ValueError, match="repeats"):
        cellsnp.write_mtx(str(tmp_path / "d.mtx"), (2, 2), [0, 0], [1, 1], [1, 2])
    with pytest.raises(KeyError):
        cellsnp.write_cellsnp(_tiny_pileup(), str(tmp_path / "q"), layers=("AD", "XX"))


def test_subset_renumbers_rows_and_columns():
    p = _tiny_pileup()
    q = p.subset(rows=[3, 0], cols=np.array([False, True]))
    assert q.shape == (2, 1)
    assert q.rows.tolist() == [3, 0]
    assert np.array_equal(q.to_dense("DP"), p.to_dense("DP")[[3, 0]][:, [1]])
    assert q.subset(rows=[1]).rows.tolist() == [0]
    with pytest.raises(ValueError):
        p.subset(rows=[0, 0])


def test_variant_id_is_vireos():
    assert cellsnp.variant_id("chr1", 10, "A", "G") == "chr1_10_A_G"
    assert cellsnp.variant_id(["1", "X"], np.array([5, 6]), ["C", "T"], ["A", "G"]) == ["1_5_C_A", "X_6_T_G"]


def test_synthetic_capture_round_trips(synthetic_capture):
    built = synthetic_capture(n_donors=4, n_sites=60, n_cells=25, oth_rate=0.1, seed=3)
    p = built["capture"].pileup
    q = cellsnp.read_cellsnp(built["cellsnp"], layers=("AD", "DP", "OTH"))
    assert q.variants == p.variants and list(q.barcodes) == list(p.barcodes)
    for name in ("AD", "DP", "OTH"):
        assert np.array_equal(q.to_dense(name), p.to_dense(name))


@pytest.mark.scipy
def test_to_csc_matches_dense():
    p = _tiny_pileup()
    assert np.array_equal(p.to_csc("DP").toarray(), p.to_dense("DP"))
    assert p.to_csc("DP").format == "csc"


@pytest.mark.vireo
def test_reader_matches_vireos_read_cellsnp(synthetic_capture):
    from vireoSNP.utils.io_utils import read_cellSNP

    built = synthetic_capture(n_donors=3, n_sites=40, n_cells=12, seed=5)
    theirs = read_cellSNP(built["cellsnp"], layers=["AD", "DP", "OTH"])
    ours = cellsnp.read_cellsnp(built["cellsnp"], layers=("AD", "DP", "OTH"))
    assert list(theirs["variants"]) == ours.variants
    assert list(theirs["samples"]) == list(ours.barcodes)
    for name in ("AD", "DP", "OTH"):
        assert np.array_equal(theirs[name].toarray(), ours.to_dense(name))
