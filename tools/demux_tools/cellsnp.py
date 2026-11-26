"""Read and write cellsnp-lite pileups with numpy and pandas (no scipy needed).

A cellsnp-lite folder holds:

    cellSNP.base.vcf.gz   one VCF record per variant row (CHROM, POS, REF, ALT)
    cellSNP.samples.tsv   one barcode per line, one per matrix column
    cellSNP.tag.AD.mtx    alternative-allele UMIs     \
    cellSNP.tag.DP.mtx    REF + ALT UMIs               > MatrixMarket coordinate integer,
    cellSNP.tag.OTH.mtx   UMIs with any other allele  /  rows = variants, columns = barcodes

Conventions in this module:
- Row and column indices in memory are 0-based int32; the files are 1-based.
- A site that appears twice in the base VCF keeps both rows (nothing is merged).
- `keep_rows` filters rows while the matrix is read, so a large matrix is never
  held whole; kept rows stay in file order.
- Matrices are written as cellsnp-lite writes them: sorted by row, then column,
  tab-separated, zero entries left out.
"""
import os

import numpy as np

from demux_tools._io import is_gzip, open_text, read_lines

FILES = {
    "base": "cellSNP.base.vcf.gz",
    "samples": "cellSNP.samples.tsv",
    "AD": "cellSNP.tag.AD.mtx",
    "DP": "cellSNP.tag.DP.mtx",
    "OTH": "cellSNP.tag.OTH.mtx",
}
LAYERS = ("AD", "DP", "OTH")
MTX_BANNER = "%%MatrixMarket matrix coordinate integer general"
_INT32_MAX = np.iinfo(np.int32).max


def variant_id(chrom, pos, ref, alt):
    """vireo's variant id "CHROM_POS_REF_ALT" (vireoSNP's load_VCF).

    Scalars give one str; equal-length sequences give a list of str.
    """
    if isinstance(chrom, str):
        return "%s_%s_%s_%s" % (chrom, pos, ref, alt)
    return ["%s_%s_%s_%s" % t for t in zip(chrom, pos, ref, alt)]


class Pileup(object):
    """A cellsnp-lite pileup in memory.

    Attributes:
        chrom, ref, alt: numpy object arrays of str, one per variant row.
        pos: numpy int64 array, one per variant row.
        barcodes: numpy object array of str, one per column.
        layers: {name: (rows, cols, vals)}, int32 numpy arrays, 0-based
            coordinates (COO). Readers keep file order; any order is valid.
        rows: numpy int64 array, each variant's 0-based row in the source it
            was read or subset from (arange for a full read).
    Derived: n_variants, n_barcodes, shape, variants (vireo ids, list of str).
    """

    def __init__(self, chrom, pos, ref, alt, barcodes, layers=None, rows=None):
        self.chrom = np.asarray(chrom, dtype=object)
        self.pos = np.asarray(pos, dtype=np.int64)
        self.ref = np.asarray(ref, dtype=object)
        self.alt = np.asarray(alt, dtype=object)
        self.barcodes = np.asarray(barcodes, dtype=object)
        n = len(self.chrom)
        if not (len(self.pos) == len(self.ref) == len(self.alt) == n):
            raise ValueError("chrom, pos, ref and alt differ in length")
        self.rows = np.arange(n, dtype=np.int64) if rows is None else np.asarray(rows, dtype=np.int64)
        if len(self.rows) != n:
            raise ValueError("rows has %d entries for %d variants" % (len(self.rows), n))
        self.layers = {}
        for name, coo in (layers or {}).items():
            self.set_layer(name, *coo)
        self._variants = None

    @property
    def n_variants(self):
        return len(self.chrom)

    @property
    def n_barcodes(self):
        return len(self.barcodes)

    @property
    def shape(self):
        return (self.n_variants, self.n_barcodes)

    @property
    def variants(self):
        if self._variants is None:
            self._variants = variant_id(self.chrom, self.pos, self.ref, self.alt)
        return self._variants

    def set_layer(self, name, rows, cols, vals):
        """Store a layer from 0-based COO arrays (converted to int32, bounds checked)."""
        r = np.asarray(rows, dtype=np.int64)
        c = np.asarray(cols, dtype=np.int64)
        v = np.asarray(vals, dtype=np.int64)
        if not (len(r) == len(c) == len(v)):
            raise ValueError("layer %s: rows, cols and vals differ in length" % name)
        if len(r) and (r.min() < 0 or r.max() >= self.n_variants or c.min() < 0 or c.max() >= self.n_barcodes):
            raise ValueError("layer %s: an index is outside the %d x %d pileup" % ((name,) + self.shape))
        if len(v) and (v.max() > _INT32_MAX or v.min() < -_INT32_MAX):
            raise ValueError("layer %s: a value does not fit in int32" % name)
        self.layers[name] = (r.astype(np.int32), c.astype(np.int32), v.astype(np.int32))

    def to_dense(self, name, dtype=np.int64):
        """Layer `name` as a dense (n_variants, n_barcodes) array; repeated entries are summed."""
        r, c, v = self.layers[name]
        out = np.zeros(self.shape, dtype=dtype)
        np.add.at(out, (r, c), v)
        return out

    def to_csc(self, name, dtype=np.int64):
        """Layer `name` as a scipy.sparse.csc_matrix (scipy is imported here, not at module load)."""
        from scipy import sparse

        r, c, v = self.layers[name]
        return sparse.csc_matrix((v.astype(dtype), (r, c)), shape=self.shape)

    def subset(self, rows=None, cols=None):
        """A new Pileup with the given rows and columns, renumbered.

        rows, cols: None (keep all), a boolean mask, or 0-based indices without
        repeats; indices set the new order. Entries outside the kept rows or
        columns are dropped; kept entries stay in their current order. The new
        `rows` attribute maps back to the source of this pileup.
        """
        row_idx = _as_index(rows, self.n_variants, "rows")
        col_idx = _as_index(cols, self.n_barcodes, "cols")
        row_new = _renumber(row_idx, self.n_variants)
        col_new = _renumber(col_idx, self.n_barcodes)
        layers = {}
        for name, (r, c, v) in self.layers.items():
            nr, nc = row_new[r], col_new[c]
            keep = (nr >= 0) & (nc >= 0)
            layers[name] = (nr[keep], nc[keep], v[keep])
        return Pileup(self.chrom[row_idx], self.pos[row_idx], self.ref[row_idx], self.alt[row_idx],
                      self.barcodes[col_idx], layers=layers, rows=self.rows[row_idx])


def _as_index(sel, n, what):
    if sel is None:
        return np.arange(n, dtype=np.int64)
    a = np.asarray(sel)
    if a.dtype == bool:
        if len(a) != n:
            raise ValueError("%s mask has %d entries for %d" % (what, len(a), n))
        return np.flatnonzero(a).astype(np.int64)
    a = a.astype(np.int64)
    if len(a) and (a.min() < 0 or a.max() >= n):
        raise ValueError("%s index outside 0..%d" % (what, n - 1))
    if len(np.unique(a)) != len(a):
        raise ValueError("%s index repeats an entry" % what)
    return a


def _renumber(idx, n):
    new = np.full(n, -1, dtype=np.int64)
    new[idx] = np.arange(len(idx), dtype=np.int64)
    return new


def _row_mask(keep_rows, n_rows):
    """None, or a boolean mask of length n_rows from a mask or 0-based indices (repeats allowed)."""
    if keep_rows is None:
        return None
    a = np.asarray(keep_rows)
    if a.dtype == bool:
        if len(a) != n_rows:
            raise ValueError("keep_rows mask has %d entries for %d rows" % (len(a), n_rows))
        return a
    a = a.astype(np.int64)
    if len(a) and (a.min() < 0 or a.max() >= n_rows):
        raise ValueError("keep_rows index outside 0..%d" % (n_rows - 1))
    mask = np.zeros(n_rows, dtype=bool)
    mask[a] = True
    return mask


# ---------------------------------------------------------------------------
# Paths


def resolve_cellsnp(cellsnp=None, base=None, samples=None, ad=None, dp=None, oth=None):
    """Paths of one pileup's files: {"base", "samples", "AD", "DP", "OTH"} -> path or None.

    `cellsnp` is a cellsnp-lite folder; explicit file arguments override its
    files, so DVC objects (files named by their md5) can be read in place with
    no folder of links. "base" and "samples" must exist; a layer file absent
    from both is None. Raises FileNotFoundError.
    """
    out = {}
    explicit = {"base": base, "samples": samples, "AD": ad, "DP": dp, "OTH": oth}
    for key, name in FILES.items():
        path = explicit[key]
        if path is None and cellsnp is not None:
            candidate = os.path.join(os.fspath(cellsnp), name)
            path = candidate if os.path.exists(candidate) else None
        if path is not None:
            path = os.fspath(path)
            if not os.path.exists(path):
                raise FileNotFoundError(path)
        out[key] = path
    for key in ("base", "samples"):
        if out[key] is None:
            raise FileNotFoundError("cellsnp %s file not given and not found in %r" % (key, cellsnp))
    return out


# ---------------------------------------------------------------------------
# Reading


def mtx_header(path):
    """(header_lines, n_rows, n_cols, nnz) of a MatrixMarket coordinate integer file.

    header_lines counts the banner, comment lines and the size line. Raises
    ValueError for any other MatrixMarket kind.
    """
    with open_text(path) as f:
        banner = f.readline()
        words = banner.lower().split()
        if words[:5] != ["%%matrixmarket", "matrix", "coordinate", "integer", "general"]:
            raise ValueError("%s: expected '%s', found %r" % (path, MTX_BANNER, banner.strip()))
        n = 1
        for line in f:
            n += 1
            s = line.strip()
            if not s or s.startswith("%"):
                continue
            parts = s.split()
            if len(parts) != 3:
                raise ValueError("%s: bad size line %r" % (path, s))
            return n, int(parts[0]), int(parts[1]), int(parts[2])
    raise ValueError("%s: no size line" % path)


def read_mtx(path, keep_rows=None, chunksize=5000000):
    """Read a MatrixMarket coordinate integer matrix with pandas' C reader.

    keep_rows: None, a boolean mask over rows, or 0-based row indices (a set:
    order and repeats ignored). Only entries in kept rows are returned; row
    numbers are NOT renumbered.

    Returns (shape, i, j, v): shape = (n_rows, n_cols) from the size line;
    i, j = 0-based int32 row and column of each kept entry, v = int32 values,
    in file order. Raises ValueError if the entry count differs from the size
    line or an index is out of range. Gzip input is accepted (magic bytes).
    """
    import pandas as pd

    path = os.fspath(path)
    n_header, n_rows, n_cols, nnz = mtx_header(path)
    mask = _row_mask(keep_rows, n_rows)
    empty = np.zeros(0, dtype=np.int32)
    if nnz == 0:
        return (n_rows, n_cols), empty, empty.copy(), empty.copy()
    reader = pd.read_csv(path, sep=r"\s+", header=None, skiprows=n_header, names=["i", "j", "v"],
                         dtype=np.int64, chunksize=chunksize, engine="c",
                         compression="gzip" if is_gzip(path) else None)
    parts_i, parts_j, parts_v = [], [], []
    total = 0
    try:
        for chunk in reader:
            a = chunk.to_numpy(dtype=np.int64)
            total += a.shape[0]
            if a.shape[0] == 0:
                continue
            i = a[:, 0] - 1
            j = a[:, 1] - 1
            if i.min() < 0 or i.max() >= n_rows or j.min() < 0 or j.max() >= n_cols:
                raise ValueError("%s: an entry lies outside the %d x %d matrix" % (path, n_rows, n_cols))
            v = a[:, 2]
            if v.max() > _INT32_MAX or v.min() < -_INT32_MAX:
                raise ValueError("%s: a value does not fit in int32" % path)
            if mask is not None:
                m = mask[i]
                i, j, v = i[m], j[m], v[m]
            parts_i.append(i.astype(np.int32))
            parts_j.append(j.astype(np.int32))
            parts_v.append(v.astype(np.int32))
    finally:
        reader.close()
    if total != nnz:
        raise ValueError("%s: size line says %d entries, file has %d" % (path, nnz, total))
    if not parts_i:
        return (n_rows, n_cols), empty, empty.copy(), empty.copy()
    return (n_rows, n_cols), np.concatenate(parts_i), np.concatenate(parts_j), np.concatenate(parts_v)


def read_base(path):
    """(chrom, pos, ref, alt) of every record of a VCF, in file order, repeats kept.

    chrom, ref, alt: numpy object arrays of str; pos: int64. Reads only the
    first five columns' CHROM, POS, REF and ALT with pandas' C reader.
    """
    import pandas as pd

    path = os.fspath(path)
    n_header = 0
    with open_text(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            n_header += 1
        else:
            empty = np.zeros(0, dtype=object)
            return empty, np.zeros(0, dtype=np.int64), empty.copy(), empty.copy()
    df = pd.read_csv(path, sep="\t", header=None, skiprows=n_header, usecols=[0, 1, 3, 4],
                     dtype={0: str, 1: np.int64, 3: str, 4: str}, keep_default_na=False, na_filter=False,
                     compression="gzip" if is_gzip(path) else None, engine="c")
    return (df[0].to_numpy(dtype=object), df[1].to_numpy(dtype=np.int64),
            df[3].to_numpy(dtype=object), df[4].to_numpy(dtype=object))


def read_samples(path):
    """Barcodes of cellSNP.samples.tsv, one per line, as a list of str."""
    return [line.strip() for line in read_lines(path)]


def read_cellsnp(paths, layers=("AD", "DP"), keep_rows=None, chunksize=5000000):
    """Read a cellsnp-lite pileup into a Pileup.

    paths: a cellsnp-lite folder, or the dict of `resolve_cellsnp`.
    layers: matrices to read (each must exist).
    keep_rows: None, a boolean mask over the base VCF's records, or their
        0-based indices (a set). The Pileup then holds the kept rows only, in
        file order, renumbered; `pileup.rows` gives each one's row in the files.
    Raises ValueError when a matrix's size line disagrees with the base VCF or
    the barcode list.
    """
    if not isinstance(paths, dict):
        paths = resolve_cellsnp(cellsnp=paths)
    chrom, pos, ref, alt = read_base(paths["base"])
    barcodes = read_samples(paths["samples"])
    n_rows, n_cols = len(chrom), len(barcodes)
    mask = _row_mask(keep_rows, n_rows)
    kept = np.arange(n_rows, dtype=np.int64) if mask is None else np.flatnonzero(mask).astype(np.int64)
    new_row = _renumber(kept, n_rows)
    coo = {}
    for name in layers:
        if paths.get(name) is None:
            raise FileNotFoundError("cellsnp layer %s not found" % name)
        shape, i, j, v = read_mtx(paths[name], keep_rows=mask, chunksize=chunksize)
        if shape != (n_rows, n_cols):
            raise ValueError("%s is %d x %d; the base VCF has %d records and the barcode list %d"
                             % (paths[name], shape[0], shape[1], n_rows, n_cols))
        coo[name] = (new_row[i], j, v)
    return Pileup(chrom[kept], pos[kept], ref[kept], alt[kept], barcodes, layers=coo, rows=kept)


# ---------------------------------------------------------------------------
# Writing


def _sorted_entries(name, r, c, v):
    keep = v != 0
    r, c, v = r[keep].astype(np.int64), c[keep].astype(np.int64), v[keep].astype(np.int64)
    order = np.lexsort((c, r))
    r, c, v = r[order], c[order], v[order]
    if len(r) > 1:
        same = (r[1:] == r[:-1]) & (c[1:] == c[:-1])
        if same.any():
            k = int(np.flatnonzero(same)[0])
            raise ValueError("layer %s repeats entry (row %d, column %d)" % (name, r[k], c[k]))
    return r, c, v


def write_mtx(path, shape, rows, cols, vals, name="matrix"):
    """Write 0-based COO entries as cellsnp-lite does: banner, "%" line, size line,
    then 1-based "row<TAB>col<TAB>value" sorted by row then column; zeros left
    out; a repeated (row, col) raises ValueError. Atomic; gzip if ".gz"."""
    import pandas as pd

    r, c, v = _sorted_entries(name, np.asarray(rows), np.asarray(cols), np.asarray(vals))
    with open_text(path, "wt") as f:
        f.write("%s\n%%\n%d\t%d\t%d\n" % (MTX_BANNER, shape[0], shape[1], len(v)))
        if len(v):
            pd.DataFrame({"i": r + 1, "j": c + 1, "v": v}).to_csv(f, sep="\t", header=False, index=False)
    return os.fspath(path)


def write_cellsnp(pileup, out_dir, layers=None):
    """Write `pileup` as a cellsnp-lite folder; returns {key: path} like `resolve_cellsnp`.

    layers: names to write (default: every layer the pileup holds; a name it
    lacks raises KeyError). The base VCF carries INFO "AD=..;DP=..;OTH=.." sums
    of the written layers, as cellsnp-lite writes it, and is plain gzip with
    mtime 0 (not BGZF, so it cannot be indexed). Every file is written
    atomically; the folder is created if needed.
    """
    import pandas as pd

    out_dir = os.fspath(out_dir)
    os.makedirs(out_dir, exist_ok=True)
    names = [k for k in LAYERS if k in pileup.layers] + sorted(k for k in pileup.layers if k not in LAYERS)
    if layers is not None:
        missing = [k for k in layers if k not in pileup.layers]
        if missing:
            raise KeyError("pileup has no layer %s" % ", ".join(missing))
        names = list(layers)
    paths = {"base": os.path.join(out_dir, FILES["base"]), "samples": os.path.join(out_dir, FILES["samples"])}
    with open_text(paths["samples"], "wt") as f:
        for b in pileup.barcodes:
            f.write("%s\n" % b)
    info_parts = []
    for name in [k for k in LAYERS if k in names]:
        r, _, v = pileup.layers[name]
        sums = np.zeros(pileup.n_variants, dtype=np.int64)
        np.add.at(sums, r.astype(np.int64), v.astype(np.int64))
        info_parts.append([("%s=%d" % (name, s)) for s in sums])
    info = [";".join(t) for t in zip(*info_parts)] if info_parts else ["."] * pileup.n_variants
    base = pd.DataFrame({"CHROM": pileup.chrom, "POS": pileup.pos, "ID": ".", "REF": pileup.ref,
                         "ALT": pileup.alt, "QUAL": ".", "FILTER": "PASS", "INFO": info})
    with open_text(paths["base"], "wt") as f:
        f.write("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        if len(base):
            base.to_csv(f, sep="\t", header=False, index=False)
    for name in names:
        file_name = FILES.get(name, "cellSNP.tag.%s.mtx" % name)
        paths[name] = write_mtx(os.path.join(out_dir, file_name), pileup.shape, *pileup.layers[name], name=name)
    for name in LAYERS:
        paths.setdefault(name, None)
    return paths
