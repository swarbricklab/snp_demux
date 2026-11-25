"""demux_tools._io: tables, deterministic gzip, atomic writes, provenance."""
import gzip
import json
import os

import numpy as np
import pandas as pd
import pytest

from demux_tools import _io


def _frame():
    return pd.DataFrame({
        "cell": ["AAAC-1", "AAAG-1", "AAAT-1"],
        "donor": ["0012", "NA", ""],
        "prob": [0.123456, float("nan"), 1e-300],
        "n": [1, 20, 300],
    })


def test_gzip_round_trip_gives_identical_bytes(tmp_path):
    a, b = tmp_path / "a.tsv.gz", tmp_path / "b.tsv.gz"
    _io.write_table(_frame(), str(a), float_formats={"prob": "%.2e"})
    _io.write_table(_frame(), str(b), float_formats={"prob": "%.2e"})
    raw = a.read_bytes()
    assert raw == b.read_bytes()
    assert raw[:2] == b"\x1f\x8b"
    assert raw[4:8] == b"\x00\x00\x00\x00"   # mtime 0
    assert raw[3] & 0x08 == 0               # no file name stored
    back = _io.read_table(str(a))
    assert back.to_dict("list") == {
        "cell": ["AAAC-1", "AAAG-1", "AAAT-1"],
        "donor": ["0012", "NA", ""],
        "prob": ["1.23e-01", "", "1.00e-300"],
        "n": ["1", "20", "300"],
    }
    # A second write and read gives the same text again.
    _io.write_table(back, str(b))
    assert gzip.decompress(b.read_bytes()) == gzip.decompress(raw)


def test_read_table_detects_gzip_by_magic_not_name(tmp_path):
    p = tmp_path / "no_extension"
    p.write_bytes(gzip.compress(b"a\tb\n007\tNA\n"))
    assert _io.read_table(str(p)).to_dict("list") == {"a": ["007"], "b": ["NA"]}
    q = tmp_path / "plain.tsv.gz"   # the name says gzip, the bytes are plain
    q.write_text("a\tb\n1\t2\n")
    assert _io.read_table(str(q)).to_dict("list") == {"a": ["1"], "b": ["2"]}


def test_open_text_reads_gzip_and_plain(tmp_path):
    with _io.open_text(str(tmp_path / "x.txt.gz"), "wt") as f:
        f.write("one\ntwo\n")
    assert _io.read_lines(str(tmp_path / "x.txt.gz")) == ["one", "two"]
    with _io.open_text(str(tmp_path / "x.txt"), "wt") as f:
        f.write("one\n")
    with _io.open_text(str(tmp_path / "x.txt")) as f:
        assert f.read() == "one\n"
    with pytest.raises(ValueError):
        _io.open_text(str(tmp_path / "x.txt"), "a")


def test_failed_write_leaves_no_file(tmp_path):
    target = tmp_path / "out.tsv.gz"
    with pytest.raises(RuntimeError):
        with _io.open_text(str(target), "wt") as f:
            f.write("partial\n")
            raise RuntimeError("stop")
    assert os.listdir(str(tmp_path)) == []
    target.write_text("old\n")
    with pytest.raises(RuntimeError):
        with _io.open_text(str(target), "wt") as f:
            f.write("new\n")
            raise RuntimeError("stop")
    assert target.read_text() == "old\n"
    assert os.listdir(str(tmp_path)) == ["out.tsv.gz"]


def test_format_column_blanks_missing_values():
    assert _io.format_column([0.5, None, float("nan"), np.float32(2)], "%.2e") == ["5.00e-01", "", "", "2.00e+00"]


def test_fingerprint_follows_size_and_mtime(tmp_path):
    p = tmp_path / "f.txt"
    p.write_text("abc")
    first = _io.fingerprint([str(p)])
    assert _io.fingerprint([str(p)]) == first
    st = os.stat(str(p))
    os.utime(str(p), ns=(st.st_atime_ns, st.st_mtime_ns + 10 ** 9))
    assert _io.fingerprint([str(p)]) != first
    with pytest.raises(FileNotFoundError):
        _io.fingerprint([str(tmp_path / "missing")])


def test_check_output_path_refuses_an_input(tmp_path):
    p = tmp_path / "in.tsv"
    p.write_text("x\n")
    link = tmp_path / "link.tsv"
    os.symlink(str(p), str(link))
    with pytest.raises(_io.UsageError):
        _io.check_output_path(str(link), [str(p)])
    _io.check_output_path(str(tmp_path / "out.tsv"), [str(p), None])


def test_provenance_sidecar(tmp_path):
    inp = tmp_path / "in.tsv"
    inp.write_text("x\n")
    out = tmp_path / "out.tsv"
    out.write_text("y\n")
    clock = _io.start_clock()
    path = _io.write_provenance(str(out), ["demux_tools", "demo"], [str(inp)], clock=clock, extra={"n": np.int64(3)})
    assert path == str(out) + ".provenance.json"
    rec = json.loads(open(path).read())
    for key in ("tool", "version", "git_sha", "tools_path", "argv", "cwd", "versions", "inputs",
                "inputs_fingerprint", "started", "seconds_wall", "seconds_cpu", "peak_rss_mb", "results"):
        assert key in rec
    assert rec["results"] == {"n": 3}
    assert rec["inputs"][0]["size"] == 2
    assert rec["inputs_fingerprint"] == _io.fingerprint([str(inp)])
    assert set(rec["versions"]) == {"python", "numpy", "pandas", "scipy", "vireoSNP"}
    assert _io.sidecar_path(str(tmp_path)) == os.path.join(str(tmp_path), "provenance.json")


def _git_layout(root, head, refs=None, packed=None):
    gitdir = root / ".git"
    (gitdir / "refs" / "heads").mkdir(parents=True)
    (gitdir / "HEAD").write_text(head + "\n")
    for ref, sha in (refs or {}).items():
        (gitdir / ref).parent.mkdir(parents=True, exist_ok=True)
        (gitdir / ref).write_text(sha + "\n")
    if packed:
        (gitdir / "packed-refs").write_text("# pack-refs with: peeled\n" + "".join(
            "%s %s\n" % (sha, ref) for ref, sha in packed.items()))
    return gitdir


def test_git_sha_reads_loose_packed_detached_and_worktree(tmp_path):
    sha1, sha2, sha3 = "a" * 40, "b" * 40, "c" * 40
    loose = tmp_path / "loose"
    _git_layout(loose, "ref: refs/heads/main", refs={"refs/heads/main": sha1})
    (loose / "tools" / "demux_tools").mkdir(parents=True)
    assert _io.git_sha(str(loose / "tools" / "demux_tools" / "_io.py")) == sha1

    packed = tmp_path / "packed"
    gitdir = _git_layout(packed, "ref: refs/heads/topic", packed={"refs/heads/topic": sha2})
    assert _io.git_sha(str(packed / "file")) == sha2

    detached = tmp_path / "detached"
    _git_layout(detached, sha3)
    assert _io.git_sha(str(detached / "file")) == sha3

    worktree = tmp_path / "worktree"
    worktree.mkdir()
    wt_gitdir = gitdir / "worktrees" / "wt"
    wt_gitdir.mkdir(parents=True)
    (wt_gitdir / "HEAD").write_text("ref: refs/heads/topic\n")
    (wt_gitdir / "commondir").write_text("../..\n")
    (worktree / ".git").write_text("gitdir: %s\n" % wt_gitdir)
    assert _io.git_sha(str(worktree / "file")) == sha2

    unresolved = tmp_path / "unresolved"
    _git_layout(unresolved, "ref: refs/heads/none")
    assert _io.git_sha(str(unresolved / "file")) is None


def test_library_versions_never_import_scipy():
    import sys

    loaded = "scipy" in sys.modules
    versions = _io.library_versions()
    assert ("scipy" in sys.modules) == loaded
    assert versions["numpy"] == np.__version__
