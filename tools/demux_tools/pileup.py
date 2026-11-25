"""Compare, subset and join pileup folders of either engine; commands pileup-compare, pileup-subset, pileup-join.

Two engines' folders are read and written:

    cellsnp      cellSNP.base.vcf.gz (one record per row), cellSNP.samples.tsv
                 (one barcode per column), cellSNP.tag.AD.mtx, cellSNP.tag.DP.mtx,
                 cellSNP.tag.OTH.mtx
    souporcell   variants.vcf (one record per row; souporcell's
                 --known_genotypes matches records to rows by order),
                 barcodes.tsv, ref.mtx, alt.mtx

A site is CHROM, POS, REF and ALT of its record. A site written twice keeps
both rows; the k-th copy of a site matches the k-th copy in another pileup.

Files are written in the layout and entry order the pileup tools write:

- Matrices keep their header lines as read (banner and comment lines) and the
  separator of their size line and entries. Entries keep their file order
  with rows and columns renumbered, so a matrix sorted by row then column
  stays sorted; explicit zero entries (vartrix writes them) are kept.
- VCF header lines and records are copied as read. A cellsnp record's INFO
  (AD, DP and OTH summed over the barcodes) is recomputed when the barcodes
  change. A gzip input is written as BGZF, which htslib reads and indexes.
- The barcode list keeps the input's order and is written one per line.
- A file whose content does not change is copied byte for byte.

What a subset can and cannot reproduce:

- A site's counts do not depend on which other sites were piled, and a
  barcode's counts do not depend on which other barcodes were piled, so a
  subset of a pileup has the counts a fresh pileup of the same BAM, barcode
  list, tool version and parameters would have at those sites and barcodes.
- A site that was not piled cannot be added, and a site the tool dropped
  (cellsnp's --minCOUNT over the barcodes piled) cannot be brought back.
- cellsnp drops a site whose UMIs over the piled barcodes fall under
  --minCOUNT. After a barcode subset, `--min-count` applies the same rule to
  the kept barcodes (DP + OTH summed over them). cellsnp's --minMAF is not
  re-applied (OTH does not say which other allele was seen).
- souporcell's rows come from a scan of the BAM. When sc_pileup skips
  souporcell's remapping (skip_remap: true, its default) the scan reads the
  whole Cell Ranger BAM and does not depend on the barcodes, so its rows are
  never filtered by count. With the remapping it depends on them.
- Output folders: a folder that holds files is refused unless force; every
  file is written under a temporary name and renamed into place, so a link
  in the folder is replaced, never written through.
"""
import gzip
import os
import shutil
import struct
import zlib

import numpy as np
import pandas as pd

from demux_tools._io import (UsageError, atomic_path, check_output_path, is_gzip, open_text, read_lines,
                             start_clock, write_provenance, write_table)

ENGINES = {
    "cellsnp": {"sites": "cellSNP.base.vcf.gz", "barcodes": "cellSNP.samples.tsv",
                "layers": (("AD", "cellSNP.tag.AD.mtx"), ("DP", "cellSNP.tag.DP.mtx"),
                           ("OTH", "cellSNP.tag.OTH.mtx"))},
    "souporcell": {"sites": "variants.vcf", "barcodes": "barcodes.tsv",
                   "layers": (("alt", "alt.mtx"), ("ref", "ref.mtx"))},
}
BGZF_BLOCK = 0xff00
BGZF_EOF = bytes.fromhex("1f8b08040000000000ff0600424302001b0003000000000000000000")
COMPARE_COLUMNS = ["layer", "entries_a", "entries_b", "in_both", "equal", "only_a", "only_b"]


def detect_engine(folder):
    """"cellsnp" or "souporcell" from the files in a folder; UsageError when neither fits."""
    for name, spec in ENGINES.items():
        if os.path.exists(os.path.join(folder, spec["barcodes"])) and os.path.exists(
                os.path.join(folder, spec["sites"])):
            return name
    raise UsageError("%s holds neither a cellsnp nor a souporcell pileup" % folder)


# ---------------------------------------------------------------------------
# BGZF


def bgzf_bytes(data, level=6):
    """BGZF blocks of `data` (bytes) with the end-of-file block; deterministic for one zlib."""
    out = []
    for start in range(0, len(data), BGZF_BLOCK):
        chunk = data[start:start + BGZF_BLOCK]
        c = zlib.compressobj(level, zlib.DEFLATED, -15)
        cdata = c.compress(chunk) + c.flush()
        bsize = len(cdata) + 25
        out.append(b"\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC\x02\x00" + struct.pack("<H", bsize)
                   + cdata + struct.pack("<II", zlib.crc32(chunk) & 0xffffffff, len(chunk)))
    out.append(BGZF_EOF)
    return b"".join(out)


def _write_bytes(path, data, gz):
    with atomic_path(path) as tmp:
        with open(tmp, "wb") as f:
            f.write(bgzf_bytes(data) if gz else data)
    return path


def _copy(src, dst):
    """Copy a file into place through a temporary name; a link at dst is replaced, never written through."""
    with atomic_path(dst) as tmp:
        shutil.copyfile(src, tmp)
    return dst


# Every file the pileup commands write into an output folder.
OUTPUT_FILES = sorted(set([spec["sites"] for spec in ENGINES.values()] + [spec["barcodes"] for spec in ENGINES.values()]
                          + [fn for spec in ENGINES.values() for _, fn in spec["layers"]] + ["provenance.json"]))


def prepare_out_dir(out_dir, writing, force=False):
    """Create out_dir. A folder with files in it is refused unless force.

    With force, the files these commands write that this run will not write
    again (e.g. a matrix of a layer the new pileup lacks, or provenance.json)
    are removed first, so none is left stale. A link is removed, never followed.
    """
    if os.path.isdir(out_dir) and os.listdir(out_dir):
        if not force:
            raise UsageError("output folder %s is not empty; give --force to replace the pileup files in it" % out_dir)
        for name in OUTPUT_FILES:
            path = os.path.join(out_dir, name)
            if name not in writing and os.path.lexists(path):
                os.unlink(path)
    try:
        os.makedirs(out_dir, exist_ok=True)
    except OSError as e:
        raise UsageError("cannot create output folder %s: %s" % (out_dir, e))


# ---------------------------------------------------------------------------
# Reading


class Matrix(object):
    """A MatrixMarket coordinate matrix of integers: header lines, separators, shape, 0-based entries."""

    def __init__(self, header, size_sep, sep, shape, i, j, v):
        self.header, self.size_sep, self.sep, self.shape = header, size_sep, sep, shape
        self.i, self.j, self.v = i, j, v


def read_matrix(path):
    """Read a coordinate matrix (integer, or real with integral values as vartrix writes them)."""
    header = []
    with open_text(path) as f:
        banner = f.readline()
        if not banner.lower().startswith("%%matrixmarket matrix coordinate"):
            raise ValueError("%s: not a MatrixMarket coordinate matrix" % path)
        header.append(banner)
        size = None
        for line in f:
            if line.startswith("%"):
                header.append(line)
                continue
            size = line
            break
        first = f.readline()
    if size is None:
        raise ValueError("%s: no size line" % path)
    n_rows, n_cols, nnz = (int(x) for x in size.split())
    size_sep = "\t" if "\t" in size else " "
    sep = ("\t" if "\t" in first else " ") if first.strip() else size_sep
    if nnz == 0:
        z = np.zeros(0, dtype=np.int64)
        return Matrix(header, size_sep, sep, (n_rows, n_cols), z, z.copy(), z.copy())
    df = pd.read_csv(path, sep=r"\s+", header=None, skiprows=len(header) + 1, names=["i", "j", "v"],
                     dtype={"i": np.int64, "j": np.int64, "v": str}, engine="c",
                     compression="gzip" if is_gzip(path) else None)
    if len(df) != nnz:
        raise ValueError("%s: size line says %d entries, file has %d" % (path, nnz, len(df)))
    vf = pd.to_numeric(df.v).to_numpy()
    v = vf.astype(np.int64)
    if not np.array_equal(v, vf):
        raise ValueError("%s: a value is not a whole number" % path)
    i, j = df.i.to_numpy() - 1, df.j.to_numpy() - 1
    if len(i) and (i.min() < 0 or i.max() >= n_rows or j.min() < 0 or j.max() >= n_cols):
        raise ValueError("%s: an entry lies outside the %d x %d matrix" % (path, n_rows, n_cols))
    return Matrix(header, size_sep, sep, (n_rows, n_cols), i, j, v)


def write_matrix(path, m, i, j, v, shape):
    """Write entries (0-based, in the order given) with m's header lines and separators."""
    with atomic_path(path) as tmp:
        with open(tmp, "w", newline="\n") as f:
            f.write("".join(m.header))
            f.write(m.size_sep.join(str(x) for x in (shape[0], shape[1], len(v))) + "\n")
            if len(v):
                pd.DataFrame({"i": i + 1, "j": j + 1, "v": v}).to_csv(f, sep=m.sep, header=False, index=False)
    return path


class Sites(object):
    """A VCF read as text: header lines, record lines, site keys and whether it was gzip."""

    def __init__(self, header, records, gz):
        self.header, self.records, self.gz = header, records, gz
        fields = [r.split("\t", 5) for r in records]
        for k, f in enumerate(fields):
            if len(f) < 5:
                raise ValueError("record %d has fewer than five columns" % (k + 1))
        self.chrom = [f[0] for f in fields]
        self.pos = [f[1] for f in fields]
        self.keys = ["%s:%s:%s:%s" % (f[0], f[1], f[3], f[4]) for f in fields]

    def __len__(self):
        return len(self.records)


def read_vcf_text(path):
    """Sites of a VCF (gzip, BGZF or plain), every line kept as written."""
    gz = is_gzip(path)
    opener = gzip.open if gz else open
    with opener(path, "rt") as f:
        text = f.read()
    lines = text.split("\n")
    if lines and lines[-1] == "":
        lines.pop()
    n = 0
    while n < len(lines) and lines[n].startswith("#"):
        n += 1
    return Sites(lines[:n], lines[n:], gz)


def occurrence_keys(keys):
    """Each site key with its copy number: "key#0" for the first copy, "key#1" for a repeat."""
    keys = list(keys)
    if not keys:
        return []
    s = pd.Series(keys, dtype=object)
    return (s + "#" + s.groupby(s, sort=False).cumcount().astype(str)).tolist()


class PileupFolder(object):
    """A pileup folder in memory: engine, sites, barcodes and matrices by layer."""

    def __init__(self, folder, layers=None):
        self.folder = folder
        self.engine = detect_engine(folder)
        spec = ENGINES[self.engine]
        self.sites = read_vcf_text(os.path.join(folder, spec["sites"]))
        self.barcodes = [b for b in read_lines(os.path.join(folder, spec["barcodes"]))]
        self.files = {}
        self.matrices = {}
        for name, fname in spec["layers"]:
            p = os.path.join(folder, fname)
            if not os.path.exists(p) or (layers is not None and name not in layers):
                continue
            self.files[name] = p
            m = read_matrix(p)
            if m.shape != (len(self.sites), len(self.barcodes)):
                raise ValueError("%s is %d x %d; %d records and %d barcodes" % (
                    p, m.shape[0], m.shape[1], len(self.sites), len(self.barcodes)))
            self.matrices[name] = m

    def inputs(self):
        spec = ENGINES[self.engine]
        return [os.path.join(self.folder, spec["sites"]), os.path.join(self.folder, spec["barcodes"])] + list(
            self.files.values())


def read_site_list(path):
    """Site keys from a VCF (records' CHROM, POS, REF, ALT) or a table of CHROM POS REF ALT (header optional)."""
    with open_text(path) as f:
        first = f.readline()
    if first.startswith("#"):
        from demux_tools.cellsnp import read_base

        chrom, pos, ref, alt = read_base(path)
        return set("%s:%d:%s:%s" % t for t in zip(chrom, pos, ref, alt))
    out = set()
    for k, line in enumerate(read_lines(path)):
        f = line.split("\t") if "\t" in line else line.split()
        if len(f) < 4:
            continue
        if k == 0 and not f[1].isdigit():
            continue
        out.add("%s:%s:%s:%s" % (f[0], f[1], f[2], f[3]))
    return out


def read_barcode_list(path):
    """Barcodes, one per line (gzip by magic bytes); empty lines skipped."""
    return [b.strip() for b in read_lines(path) if b.strip()]


# ---------------------------------------------------------------------------
# Subset


def _renumber(keep):
    new = np.full(len(keep), -1, dtype=np.int64)
    new[keep] = np.arange(int(keep.sum()), dtype=np.int64)
    return new


def cellsnp_info(p, row_keep, col_keep):
    """"AD=..;DP=..;OTH=.." per kept row, summed over the kept columns, as cellsnp-lite writes INFO."""
    sums = {}
    for name in ("AD", "DP", "OTH"):
        s = np.zeros(len(row_keep), dtype=np.int64)
        if name in p.matrices:
            m = p.matrices[name]
            k = col_keep[m.j]
            np.add.at(s, m.i[k], m.v[k])
        sums[name] = s[row_keep]
    names = [n for n in ("AD", "DP", "OTH") if n in p.matrices]
    return [";".join("%s=%d" % (n, sums[n][r]) for n in names) for r in range(int(row_keep.sum()))]


def subset(in_dir, out_dir, sites=None, barcodes=None, min_count=None, strip_prefix=None, info="auto",
           allow_empty=False, force=False):
    """Write the rows of `sites` (keys, or None for all) and columns of `barcodes` (or None) of a pileup.

    min_count: cellsnp only; drop rows whose DP + OTH over the kept barcodes is
    under it. strip_prefix: removed from every kept barcode (each must start
    with it). info: cellsnp INFO, "recompute" (sums over the kept barcodes),
    "keep" (as read) or "auto" (recompute when the barcodes change).
    A result with no row or no barcode raises UsageError unless allow_empty;
    force: see prepare_out_dir. Returns a dict of counts, with the sites and
    barcodes requested and found in the pileup.
    """
    p = PileupFolder(in_dir)
    spec = ENGINES[p.engine]
    n_rows, n_cols = len(p.sites), len(p.barcodes)
    row_keep = np.ones(n_rows, dtype=bool)
    if sites is not None:
        row_keep = np.array([k in sites for k in p.sites.keys], dtype=bool)
    col_keep = np.ones(n_cols, dtype=bool)
    found = {"sites_requested": None, "sites_found": None, "barcodes_requested": None, "barcodes_found": None}
    if sites is not None:
        keyset = set(p.sites.keys)
        found["sites_requested"] = len(sites)
        found["sites_found"] = sum(1 for k in sites if k in keyset)
    if barcodes is not None:
        wanted = set(barcodes)
        col_keep = np.array([b in wanted for b in p.barcodes], dtype=bool)
        found["barcodes_requested"] = len(wanted)
        found["barcodes_found"] = len(wanted & set(p.barcodes))
    if min_count is not None:
        if p.engine != "cellsnp":
            raise UsageError("--min-count applies to cellsnp pileups only")
        total = np.zeros(n_rows, dtype=np.int64)
        for name in ("DP", "OTH"):
            if name not in p.matrices:
                raise UsageError("--min-count needs the %s matrix" % name)
            m = p.matrices[name]
            k = col_keep[m.j]
            np.add.at(total, m.i[k], m.v[k])
        row_keep &= total >= min_count
    kept_barcodes = [b for b, k in zip(p.barcodes, col_keep) if k]
    if strip_prefix:
        bad = [b for b in kept_barcodes if not b.startswith(strip_prefix)]
        if bad:
            raise UsageError("barcode %s does not start with %s" % (bad[0], strip_prefix))
        kept_barcodes = [b[len(strip_prefix):] for b in kept_barcodes]
    if (not row_keep.any() or not col_keep.any()) and not allow_empty:
        asked = ["%s requested %d, found %d" % (w, found[w + "_requested"], found[w + "_found"])
                 for w in ("sites", "barcodes") if found[w + "_requested"] is not None]
        raise UsageError("the subset would hold %d sites and %d barcodes (%s); nothing written; check the names, "
                         "or give --allow-empty" % (int(row_keep.sum()), int(col_keep.sum()), "; ".join(asked)))
    for f in [spec["sites"], spec["barcodes"]] + [fn for _, fn in spec["layers"]]:
        check_output_path(os.path.join(out_dir, f), p.inputs())
    prepare_out_dir(out_dir, [spec["sites"], spec["barcodes"]] + [os.path.basename(f) for f in p.files.values()],
                    force)
    rows_same = bool(row_keep.all())
    cols_same = bool(col_keep.all()) and not strip_prefix
    if info == "auto":
        info = "keep" if bool(col_keep.all()) else "recompute"
    # sites file
    src = os.path.join(in_dir, spec["sites"])
    dst = os.path.join(out_dir, spec["sites"])
    if p.engine == "cellsnp" and info == "recompute":
        new_info = cellsnp_info(p, row_keep, col_keep)
        recs = []
        for rec, inf in zip([r for r, k in zip(p.sites.records, row_keep) if k], new_info):
            f = rec.split("\t")
            f[7] = inf
            recs.append("\t".join(f))
    else:
        recs = [r for r, k in zip(p.sites.records, row_keep) if k]
    text = "\n".join(p.sites.header + recs) + "\n"
    if rows_same and recs == p.sites.records:
        _copy(src, dst)
    else:
        _write_bytes(dst, text.encode("utf-8"), p.sites.gz)
    # barcodes
    src = os.path.join(in_dir, spec["barcodes"])
    dst = os.path.join(out_dir, spec["barcodes"])
    if cols_same:
        _copy(src, dst)
    else:
        with atomic_path(dst) as tmp:
            with open(tmp, "w", newline="\n") as f:
                f.write("".join(b + "\n" for b in kept_barcodes))
    # matrices
    rnew, cnew = _renumber(row_keep), _renumber(col_keep)
    entries = {}
    for name, fname in spec["layers"]:
        if name not in p.matrices:
            continue
        m = p.matrices[name]
        dst = os.path.join(out_dir, fname)
        if rows_same and bool(col_keep.all()):
            _copy(p.files[name], dst)
            entries[name] = len(m.v)
            continue
        k = row_keep[m.i] & col_keep[m.j]
        write_matrix(dst, m, rnew[m.i[k]], cnew[m.j[k]], m.v[k], (int(row_keep.sum()), int(col_keep.sum())))
        entries[name] = int(k.sum())
    out = {"engine": p.engine, "rows_in": n_rows, "rows_out": int(row_keep.sum()), "barcodes_in": n_cols,
           "barcodes_out": int(col_keep.sum()), "entries_out": entries, "info": info if p.engine == "cellsnp" else ""}
    out.update(found)
    return out


# ---------------------------------------------------------------------------
# Compare


def _entries(p, name, site_code, bc_code, n_bc, row_keep, col_keep, keep_zeros):
    m = p.matrices[name]
    k = row_keep[m.i] & col_keep[m.j]
    if not keep_zeros:
        k &= m.v != 0
    key = site_code[m.i[k]] * n_bc + bc_code[m.j[k]]
    return key, m.v[k], int(((m.v == 0) & row_keep[m.i] & col_keep[m.j]).sum())


def compare(a_dir, b_dir, sites=None, columns_of_b=None, keep_zeros=False, layers=None):
    """Compare two pileups as sets of (site, barcode, count). Returns (table, verdict, notes).

    The verdict is "same sites and counts" when the two hold the same sites,
    the same barcodes and equal entries in every matrix, and "differ" otherwise.

    sites: keys to compare at (both pileups cut to them). columns_of_b:
    cellsnp only; cut A to B's barcodes and to sites with DP + OTH over them
    of at least this count. keep_zeros: count explicit zero entries as
    entries (by default a zero and a missing entry are the same count).
    """
    a, b = PileupFolder(a_dir, layers), PileupFolder(b_dir, layers)
    if a.engine != b.engine:
        raise UsageError("A is a %s pileup and B a %s pileup" % (a.engine, b.engine))
    notes = []
    ka, kb = occurrence_keys(a.sites.keys), occurrence_keys(b.sites.keys)
    row_a = np.ones(len(ka), dtype=bool)
    row_b = np.ones(len(kb), dtype=bool)
    if sites is not None:
        row_a = np.array([k in sites for k in a.sites.keys], dtype=bool)
        row_b = np.array([k in sites for k in b.sites.keys], dtype=bool)
    col_a = np.ones(len(a.barcodes), dtype=bool)
    col_b = np.ones(len(b.barcodes), dtype=bool)
    if columns_of_b is not None:
        if a.engine != "cellsnp":
            raise UsageError("--columns-of-b applies to cellsnp pileups only")
        lacking = [n for n in ("DP", "OTH") if n not in a.matrices]
        if lacking:
            raise UsageError("--columns-of-b needs A's DP and OTH matrices; --layers leaves out %s"
                             % " and ".join(lacking))
        want = set(b.barcodes)
        col_a = np.array([x in want for x in a.barcodes], dtype=bool)
        total = np.zeros(len(ka), dtype=np.int64)
        for name in ("DP", "OTH"):
            m = a.matrices[name]
            k = col_a[m.j]
            np.add.at(total, m.i[k], m.v[k])
        row_a &= total >= columns_of_b
        notes.append("A cut to B's %d barcodes and to %d sites with %d UMIs or more"
                     % (int(col_a.sum()), int(row_a.sum()), columns_of_b))
    sa = [k for k, r in zip(ka, row_a) if r]
    sb = [k for k, r in zip(kb, row_b) if r]
    rep_a = sum(1 for k in sa if not k.endswith("#0"))
    rep_b = sum(1 for k in sb if not k.endswith("#0"))
    universe = pd.Index(list(dict.fromkeys(sa + sb)))
    code_a = np.full(len(ka), -1, dtype=np.int64)
    code_a[row_a] = universe.get_indexer(sa)
    code_b = np.full(len(kb), -1, dtype=np.int64)
    code_b[row_b] = universe.get_indexer(sb)
    bc_a = [x for x, c in zip(a.barcodes, col_a) if c]
    bcs = pd.Index(list(dict.fromkeys(bc_a + b.barcodes)))
    bcode_a = np.full(len(a.barcodes), -1, dtype=np.int64)
    bcode_a[col_a] = bcs.get_indexer(bc_a)
    bcode_b = bcs.get_indexer(b.barcodes).astype(np.int64)
    set_a, set_b = set(sa), set(sb)
    both_sites = len(set_a & set_b)
    rows = [{"layer": "sites", "entries_a": len(sa), "entries_b": len(sb), "in_both": both_sites, "equal": both_sites,
             "only_a": len(set_a - set_b), "only_b": len(set_b - set_a)}]
    ba, bb = set(bc_a), set(b.barcodes)
    rows.append({"layer": "barcodes", "entries_a": len(ba), "entries_b": len(bb), "in_both": len(ba & bb),
                 "equal": len(ba & bb), "only_a": len(ba - bb), "only_b": len(bb - ba)})
    if rep_a or rep_b:
        notes.append("repeated site rows: A %d, B %d" % (rep_a, rep_b))
    same = set_a == set_b
    if not same:
        notes.append("sites only in A %d, only in B %d" % (len(set_a - set_b), len(set_b - set_a)))
    names = [n for n in a.matrices if n in b.matrices]
    for missing in sorted(set(a.matrices) ^ set(b.matrices)):
        notes.append("layer %s is in one pileup only" % missing)
        same = False
    for name in names:
        key_a, v_a, z_a = _entries(a, name, code_a, bcode_a, len(bcs), row_a, col_a, keep_zeros)
        key_b, v_b, z_b = _entries(b, name, code_b, bcode_b, len(bcs), row_b, col_b, keep_zeros)
        for key, who in ((key_a, "A"), (key_b, "B")):
            if len(np.unique(key)) != len(key):
                raise ValueError("%s's %s matrix lists one site and barcode twice" % (who, name))
        common, ia, ib = np.intersect1d(key_a, key_b, assume_unique=True, return_indices=True)
        equal = int((v_a[ia] == v_b[ib]).sum())
        row = {"layer": name, "entries_a": len(key_a), "entries_b": len(key_b), "in_both": len(common),
               "equal": equal, "only_a": len(key_a) - len(common), "only_b": len(key_b) - len(common)}
        rows.append(row)
        if not (row["entries_a"] == row["entries_b"] == row["in_both"] == equal):
            same = False
            notes.append("%s entries A %d, B %d, in both %d, equal %d" % (name, len(key_a), len(key_b),
                                                                          len(common), equal))
        if (z_a or z_b) and not keep_zeros:
            notes.append("%s explicit zero entries left out: A %d, B %d" % (name, z_a, z_b))
    if ba != bb:
        same = False
        notes.append("barcodes only in A %d, only in B %d" % (len(ba - bb), len(bb - ba)))
    verdict = "same sites and counts" if same else "differ"
    return pd.DataFrame(rows, columns=COMPARE_COLUMNS), verdict, notes


# ---------------------------------------------------------------------------
# Join


def _chrom_rank(chroms):
    order = {}
    for c in chroms:
        if c not in order:
            order[c] = len(order)
    return order


def join(inputs, out_dir, sites="union", layers=None, force=False):
    """Join pileups of one engine into one folder; barcodes get each input's prefix.

    inputs: [(folder, prefix)]. sites: "union" (every site of any input,
    ordered by chromosome, in order of first appearance, and position; at
    one position the records keep the order of the input that lists them
    first, so each input's own order survives) or "intersection" (sites in
    every input, in the first input's order). With "union", the cells of one
    input have no counts at a site only another input piled: such a site reads
    as uncovered for them, not as unknown. A cellsnp record's INFO is summed
    over all barcodes. A souporcell variants.vcf keeps its first eight columns
    only: the inputs' genotype columns belong to different donors.
    layers: matrices to join (default every layer any input has; an input
    lacking one raises UsageError). force: see prepare_out_dir.
    Returns a dict of counts, with each input's rows, the rows every input
    has and the rows only one input has.
    """
    ps = [PileupFolder(d) for d, _ in inputs]
    engine = ps[0].engine
    if any(p.engine != engine for p in ps):
        raise UsageError("all inputs must be pileups of one engine")
    spec = ENGINES[engine]
    known = [name for name, _ in spec["layers"]]
    requested = list(layers) if layers else [n for n in known if any(n in p.matrices for p in ps)]
    bad = [n for n in requested if n not in known]
    if bad:
        raise UsageError("%s pileups have no layer %s (layers: %s)" % (engine, ", ".join(bad), ", ".join(known)))
    for (d, _), p in zip(inputs, ps):
        lacking = [n for n in requested if n not in p.matrices]
        if lacking:
            raise UsageError("input %s lacks layer %s; give --layers to join only the layers every input has"
                             % (d, ", ".join(lacking)))
    barcodes = []
    for p, (_, prefix) in zip(ps, inputs):
        barcodes += [prefix + x for x in p.barcodes]
    if len(set(barcodes)) != len(barcodes):
        raise UsageError("two inputs share a barcode after their prefixes; give distinct --prefix values")
    okeys = [occurrence_keys(p.sites.keys) for p in ps]
    if sites == "intersection":
        common = set(okeys[0])
        for k in okeys[1:]:
            common &= set(k)
        site_list = [k for k in okeys[0] if k in common]
    elif sites == "union":
        rank = _chrom_rank([c for p in ps for c in p.sites.chrom])
        seen = {}
        for p, keys in zip(ps, okeys):
            for k, c, pos in zip(keys, p.sites.chrom, p.sites.pos):
                if k not in seen:
                    seen[k] = (rank[c], int(pos), len(seen))
        site_list = sorted(seen, key=lambda k: seen[k])
    else:
        raise UsageError("--sites must be union or intersection")
    index = pd.Index(site_list)
    wanted = set(site_list)
    records = {}
    for p, keys in zip(ps, okeys):
        for k, rec in zip(keys, p.sites.records):
            if k not in records and k in wanted:
                records[k] = rec
    out_rows = len(site_list)
    sets = [set(k) for k in okeys]
    every = set.intersection(*sets)
    only = []
    for n, k in enumerate(sets):
        others = set().union(*(x for m, x in enumerate(sets) if m != n))
        only.append(len(k - others))
    for f in [spec["sites"], spec["barcodes"]] + [fn for _, fn in spec["layers"]]:
        for p in ps:
            check_output_path(os.path.join(out_dir, f), p.inputs())
    layer_files = dict(spec["layers"])
    prepare_out_dir(out_dir, [spec["sites"], spec["barcodes"]] + [layer_files[n] for n in requested], force)
    # matrices: rows remapped, columns offset, sorted by row then column
    joined = {}
    maps = [index.get_indexer(k).astype(np.int64) for k in okeys]
    for name, fname in spec["layers"]:
        if name not in requested:
            continue
        parts_i, parts_j, parts_v = [], [], []
        off = 0
        for p, rmap in zip(ps, maps):
            m = p.matrices[name]
            r = rmap[m.i]
            k = r >= 0
            parts_i.append(r[k])
            parts_j.append(m.j[k] + off)
            parts_v.append(m.v[k])
            off += len(p.barcodes)
        i, j, v = np.concatenate(parts_i), np.concatenate(parts_j), np.concatenate(parts_v)
        order = np.lexsort((j, i))
        joined[name] = (i[order], j[order], v[order])
        write_matrix(os.path.join(out_dir, fname), ps[0].matrices[name], i[order], j[order], v[order],
                     (out_rows, len(barcodes)))
    # sites file
    recs = [records[k] for k in site_list]
    header = list(ps[0].sites.header)
    if engine == "cellsnp":
        sums = {}
        for name in ("AD", "DP", "OTH"):
            if name in joined:
                s = np.zeros(out_rows, dtype=np.int64)
                np.add.at(s, joined[name][0], joined[name][2])
                sums[name] = s
        names = [n for n in ("AD", "DP", "OTH") if n in sums]
        new = []
        for r, rec in enumerate(recs):
            f = rec.split("\t")
            f[7] = ";".join("%s=%d" % (n, sums[n][r]) for n in names)
            new.append("\t".join(f))
        recs = new
    else:
        recs = ["\t".join(r.split("\t")[:8]) for r in recs]
        header = [("\t".join(h.split("\t")[:8]) if h.startswith("#CHROM") else h) for h in header]
    _write_bytes(os.path.join(out_dir, spec["sites"]), ("\n".join(header + recs) + "\n").encode("utf-8"),
                 ps[0].sites.gz)
    with atomic_path(os.path.join(out_dir, spec["barcodes"])) as tmp:
        with open(tmp, "w", newline="\n") as f:
            f.write("".join(x + "\n" for x in barcodes))
    return {"engine": engine, "rows": out_rows, "barcodes": len(barcodes), "inputs": len(ps),
            "entries": {k: int(len(v[2])) for k, v in joined.items()}, "layers": requested,
            "input_rows": [len(k) for k in okeys], "rows_in_every_input": len(every), "rows_only_in_input": only}


# ---------------------------------------------------------------------------
# Commands


def configure_pileup_compare(parser):
    parser.description = ("Compare two pileup folders of one engine (cellsnp or souporcell) as sets of (site, "
                          "barcode, count), ignoring row and column order. Writes a table of sites, barcodes and "
                          "each matrix's entries, and prints a verdict: 'same sites and counts' or the differences.")
    parser.add_argument("a", help="pileup folder A")
    parser.add_argument("b", help="pileup folder B")
    parser.add_argument("--sites", default=None, help="compare only at these sites (VCF or CHROM POS REF ALT table)")
    parser.add_argument("--columns-of-b", type=int, default=None, metavar="MINCOUNT",
                        help="cellsnp: cut A to B's barcodes and to sites with DP + OTH of at least MINCOUNT over them")
    parser.add_argument("--layers", default=None, help="comma-separated matrices to compare (default all)")
    parser.add_argument("--keep-zeros", action="store_true",
                        help="count explicit zero entries as entries (default: a zero equals a missing entry)")
    parser.add_argument("--check", action="store_true",
                        help="exit with status 1 unless the verdict is 'same sites and counts' (same sites, "
                             "barcodes and entries)")
    parser.add_argument("--out", required=True, help="output table")


def run_pileup_compare(args):
    clock = start_clock()
    layers = [x for x in args.layers.split(",") if x] if args.layers else None
    try:
        sites = read_site_list(args.sites) if args.sites else None
        table, verdict, notes = compare(args.a, args.b, sites, args.columns_of_b, args.keep_zeros, layers)
    except (ValueError, OSError) as e:
        raise UsageError(str(e))
    inputs = []
    for d in (args.a, args.b):
        spec = ENGINES[detect_engine(d)]
        inputs += [os.path.join(d, spec["sites"]), os.path.join(d, spec["barcodes"])]
        inputs += [os.path.join(d, fn) for _, fn in spec["layers"] if os.path.exists(os.path.join(d, fn))]
    if args.sites:
        inputs.append(args.sites)
    check_output_path(args.out, inputs)
    write_table(table, args.out)
    line = verdict if not notes else verdict + "; " + "; ".join(notes)
    write_provenance(args.out, args.argv, inputs, clock=clock, extra={"verdict": verdict, "notes": notes})
    print(table.to_string(index=False))
    print(line)
    return 1 if (args.check and verdict != "same sites and counts") else 0


def configure_pileup_subset(parser):
    parser.description = ("Cut a pileup folder (cellsnp or souporcell) to a site list, a barcode list or both, "
                          "and write it in the layout the pileup tool writes (OUT/provenance.json beside). A site "
                          "that was not piled cannot be added.")
    parser.add_argument("input", help="pileup folder")
    parser.add_argument("out", help="output folder (created if needed)")
    parser.add_argument("--sites", default=None, help="sites to keep: a VCF, or a CHROM POS REF ALT table")
    parser.add_argument("--barcodes", default=None, help="barcodes to keep, one per line (gzip allowed)")
    parser.add_argument("--min-count", type=int, default=None,
                        help="cellsnp: drop sites whose DP + OTH over the kept barcodes is under this "
                             "(cellsnp-lite's --minCOUNT)")
    parser.add_argument("--strip-prefix", default=None, help="remove this prefix from every kept barcode")
    parser.add_argument("--info", choices=("auto", "recompute", "keep"), default="auto",
                        help="cellsnp INFO sums: recompute over the kept barcodes, keep as read, or auto "
                             "(recompute when the barcodes change; default)")
    parser.add_argument("--allow-empty", action="store_true", help="write the subset even when it holds no site "
                        "or no barcode (refused by default)")
    parser.add_argument("--force", action="store_true",
                        help="write into a folder that holds files, removing pileup files this run does not write")


def run_pileup_subset(args):
    clock = start_clock()
    try:
        sites = read_site_list(args.sites) if args.sites else None
        barcodes = read_barcode_list(args.barcodes) if args.barcodes else None
        counts = subset(args.input, args.out, sites, barcodes, args.min_count, args.strip_prefix, args.info,
                        args.allow_empty, args.force)
    except (ValueError, OSError) as e:
        raise UsageError(str(e))
    spec = ENGINES[counts["engine"]]
    inputs = [os.path.join(args.input, spec["sites"]), os.path.join(args.input, spec["barcodes"])]
    inputs += [os.path.join(args.input, fn) for _, fn in spec["layers"]
               if os.path.exists(os.path.join(args.input, fn))]
    inputs += [p for p in (args.sites, args.barcodes) if p]
    write_provenance(args.out, args.argv, inputs, clock=clock, extra=counts)
    print("pileup-subset: %s pileup, %d of %d sites and %d of %d barcodes kept in %s" % (
        counts["engine"], counts["rows_out"], counts["rows_in"], counts["barcodes_out"], counts["barcodes_in"],
        args.out))
    for what in ("sites", "barcodes"):
        if counts[what + "_requested"] is not None:
            print("  %s requested %d, found in the pileup %d" % (what, counts[what + "_requested"],
                                                               counts[what + "_found"]))
    return 0


def configure_pileup_join(parser):
    parser.description = ("Join pileup folders of one engine into one in-silico pool: barcodes take each "
                          "input's prefix, sites are the union or the intersection. OUT/provenance.json beside.")
    parser.add_argument("--in", dest="inputs", action="append", required=True, help="a pileup folder; repeatable")
    parser.add_argument("--prefix", action="append", required=True, help="barcode prefix of each --in, in order")
    parser.add_argument("--sites", choices=("union", "intersection"), default="union",
                        help="sites of the joined pileup (default union: a capture's cells then have no counts at "
                             "sites only another capture piled)")
    parser.add_argument("--layers", default=None,
                        help="comma-separated matrices to join (default every layer; an input lacking one stops)")
    parser.add_argument("--force", action="store_true",
                        help="write into a folder that holds files, removing pileup files this run does not write")
    parser.add_argument("out", help="output folder (created if needed)")


def run_pileup_join(args):
    clock = start_clock()
    if len(args.inputs) != len(args.prefix):
        raise UsageError("give one --prefix per --in")
    if len(args.inputs) < 2:
        raise UsageError("give two or more --in folders")
    layers = [x for x in args.layers.split(",") if x] if args.layers else None
    try:
        counts = join(list(zip(args.inputs, args.prefix)), args.out, args.sites, layers, args.force)
    except (ValueError, OSError) as e:
        raise UsageError(str(e))
    inputs = []
    for d in args.inputs:
        spec = ENGINES[detect_engine(d)]
        inputs += [os.path.join(d, spec["sites"]), os.path.join(d, spec["barcodes"])]
        inputs += [os.path.join(d, fn) for _, fn in spec["layers"] if os.path.exists(os.path.join(d, fn))]
    write_provenance(args.out, args.argv, inputs, clock=clock, extra=counts)
    print("pileup-join: %d %s pileups, %d sites (%s), %d barcodes in %s" % (
        counts["inputs"], counts["engine"], counts["rows"], args.sites, counts["barcodes"], args.out))
    for d, n, o in zip(args.inputs, counts["input_rows"], counts["rows_only_in_input"]):
        print("  %s: %d rows, %d of them in no other input" % (d, n, o))
    print("  rows in every input: %d" % counts["rows_in_every_input"])
    if args.sites == "union" and any(counts["rows_only_in_input"]):
        print("  at a row only one input has, the other inputs' cells have no counts")
    return 0
