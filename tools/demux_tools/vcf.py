"""Donor genotypes from a VCF, read as vireo reads them, and matched to pileup sites.

vireo (vireoSNP 0.5.7, `load_VCF(biallelic_only=True)` and `parse_donor_GPb`)
reads a donor VCF as follows, and `read_genotypes` gives the same numbers:

- A record is kept when REF and ALT are each at most one character
  (`len(REF) > 1 or len(ALT) > 1` drops it), so indels and multi-allelic
  records are dropped and an ALT of "." is kept.
- The record's id is "CHROM_POS_REF_ALT", from the text of those columns.
- The genotype is the FORMAT field named by the tag (GT, GP or PL); a record
  whose FORMAT lacks the tag counts as missing for every donor.
- ".", "./." and ".|." are a third each for 0, 1 and 2 ALT alleles.
- GT adds its first and last characters, so "0|1" and "1/0" are 1 and a
  haploid "1" is 2. GP is the three numbers. PL is turned into probabilities
  with vireo's formula. Each is divided by its sum.
- vireo stops when the tag is present in fewer than 10 per cent of the kept
  records; `check_tag` raises the same error.

Differences from vireo:
- The file is read with pandas' C reader in chunks, and only the records
  asked for (`keep_ids` or `keep_rows`) are parsed further, so a panel of
  millions of records is never held whole. This changes no value vireo would
  compute for those records.
- gzip is detected by its first bytes, so a file need not end in ".gz".
- A field vireo cannot read (e.g. "0/." or a GP with two numbers) raises
  ValueError naming the value; vireo stops with an exception there. A record
  that is not parsed is not checked, so a bad field at an unpiled record
  stops vireo but not the Scorer.
- A VCF without sample columns can be read for its sites and INFO fields.
"""
import csv
import os
import re

import numpy as np

from demux_tools._io import is_gzip, open_text
from demux_tools.cellsnp import variant_id


class Genotypes(object):
    """Donor genotypes for the kept biallelic SNP records of a VCF.

    Attributes:
        variants: list of ids "CHROM_POS_REF_ALT", one per kept record.
        chrom, ref, alt: object arrays of str; pos: object array of the POS text.
        samples: donor names, in the order of the second axis of `gp`.
        gp: float64 (n_records, n_donors, 3), vireo's genotype probabilities.
        missing: bool (n_records, n_donors), True for ".", "./." and ".|.".
        info: {field: float64 array}, INFO values (NaN where absent or ".").
        rows: int64, each kept record's 0-based index among the VCF's records.
        header_samples: every sample name in the header, in file order.
        counts: {"records", "snps", "tagged", "kept"}: data lines; records
            vireo keeps (biallelic SNPs); those with the tag in FORMAT; those
            parsed (in keep_ids).
        tag: the FORMAT field read.
    """

    def __init__(self, variants, chrom, pos, ref, alt, samples, gp, missing, info, rows, header_samples,
                 counts, tag):
        self.variants = list(variants)
        self.chrom = np.asarray(chrom, dtype=object)
        self.pos = np.asarray(pos, dtype=object)
        self.ref = np.asarray(ref, dtype=object)
        self.alt = np.asarray(alt, dtype=object)
        self.samples = list(samples)
        self.gp = gp
        self.missing = missing
        self.info = info
        self.rows = np.asarray(rows, dtype=np.int64)
        self.header_samples = list(header_samples)
        self.counts = counts
        self.tag = tag

    def __len__(self):
        return len(self.variants)


def read_header_columns(path):
    """(n_header_lines, columns) of a VCF: lines starting with "#", and the names on the #CHROM line."""
    n = 0
    columns = None
    with open_text(path) as f:
        for line in f:
            if not line.startswith("#"):
                break
            n += 1
            if line.startswith("#CHROM"):
                columns = line.rstrip().split("\t")
    if columns is None:
        raise ValueError("%s has no #CHROM header line" % path)
    return n, columns


def read_header(path):
    """(n_header_lines, header_samples) of a VCF: lines starting with "#", and the names after FORMAT."""
    n, columns = read_header_columns(path)
    return n, columns[9:]


def snp_record_ids(path, chunksize=500000):
    """(ids, rows) of every record vireo keeps (REF and ALT at most one character), in file order.

    ids: list of "CHROM_POS_REF_ALT"; rows: int64, each record's 0-based index
    among the file's data lines. Reads only CHROM, POS, REF and ALT.
    """
    import pandas as pd

    path = os.fspath(path)
    n_header, _ = read_header_columns(path)
    ids, rows = [], []
    offset = 0
    reader = pd.read_csv(path, sep="\t", header=None, skiprows=n_header, usecols=[0, 1, 3, 4], dtype=str,
                         keep_default_na=False, na_filter=False, quoting=csv.QUOTE_NONE,
                         compression="gzip" if is_gzip(path) else None, engine="c", chunksize=chunksize)
    try:
        for chunk in reader:
            snp = ((chunk[3].str.len() <= 1) & (chunk[4].str.len() <= 1)).to_numpy()
            sub = chunk[snp]
            ids.extend((sub[0] + "_" + sub[1] + "_" + sub[3] + "_" + sub[4]).tolist())
            rows.append(offset + np.flatnonzero(snp))
            offset += len(chunk)
    finally:
        reader.close()
    return ids, (np.concatenate(rows).astype(np.int64) if rows else np.zeros(0, dtype=np.int64))


def _subfield(values, k):
    """The k-th ':'-separated part of each string in a pandas Series (None where absent)."""
    if k == 0:
        return values.str.partition(":")[0]
    parts = values.str.split(":", n=k + 1, expand=True)
    if parts.shape[1] <= k:
        return values.map(lambda x: None)
    return parts[k]


def _info_values(info_text, field):
    """Float values of INFO `field` (NaN where absent, "." or not a number); a flag gives 1.0."""
    import pandas as pd

    pattern = r"(?:^|;)(%s)(?:=([^;]*))?(?:;|$)" % re.escape(field)
    found = info_text.str.extract(pattern)
    present = found[0].notna().to_numpy()
    values = pd.to_numeric(found[1], errors="coerce").to_numpy(dtype=np.float64)
    flag = present & found[1].isna().to_numpy()
    values[flag] = 1.0
    return values


def read_genotypes(path, samples=None, tag="GT", keep_ids=None, info=(), chunksize=200000, keep_rows=None):
    """Read donor genotypes from a VCF (plain or gzip) as vireo would.

    samples: donor names to read, in the order wanted (default every header
        sample in file order); a name not in the header raises ValueError.
    tag: "GT", "GP" or "PL".
    keep_ids: None (keep every biallelic SNP record) or a collection of ids;
        only records whose id is in it are parsed and returned.
    keep_rows: None, or a collection of 0-based data-line indices (as in
        `Genotypes.rows`); only those records are parsed. Applied with keep_ids.
    Counts still cover the whole file.
    info: INFO fields to return as float arrays.
    Returns a Genotypes. Records keep file order. A VCF without sample
    columns (or without FORMAT) gives sites and INFO with no donors.
    """
    import pandas as pd

    from demux_tools.vireo_port import genotype_probabilities

    if tag not in ("GT", "GP", "PL"):
        raise ValueError("tag must be GT, GP or PL, not %r" % (tag,))
    path = os.fspath(path)
    n_header, columns = read_header_columns(path)
    header_samples = columns[9:]
    has_format = len(columns) > 8
    if samples is None:
        samples = list(header_samples)
    samples = list(samples)
    missing_names = [s for s in samples if s not in header_samples]
    if missing_names:
        raise ValueError("%s has no sample %s" % (path, ", ".join(missing_names)))
    if len(set(samples)) != len(samples):
        raise ValueError("a donor is named twice: %s" % ",".join(samples))
    sample_cols = [9 + header_samples.index(s) for s in samples]
    usecols = sorted(set([0, 1, 3, 4, 7] + ([8] if has_format else []) + sample_cols))
    keep = None
    if keep_ids is not None:
        keep = pd.Index(pd.unique(np.asarray(list(keep_ids), dtype=object)))
        keep.get_indexer(keep[:1])   # builds the hash table once
    keep_row_mask = None
    if keep_rows is not None:
        keep_rows = np.asarray(list(keep_rows), dtype=np.int64)
        keep_row_mask = np.zeros(int(keep_rows.max()) + 1 if len(keep_rows) else 0, dtype=bool)
        keep_row_mask[keep_rows] = True

    counts = {"records": 0, "snps": 0, "tagged": 0, "kept": 0}
    parts = {"chrom": [], "pos": [], "ref": [], "alt": [], "info": [], "codes": [], "rows": []}
    reader = pd.read_csv(path, sep="\t", header=None, skiprows=n_header, usecols=usecols, dtype=str,
                         keep_default_na=False, na_filter=False, quoting=csv.QUOTE_NONE,
                         compression="gzip" if is_gzip(path) else None, engine="c", chunksize=chunksize)
    offset = 0
    try:
        for chunk in reader:
            n = len(chunk)
            counts["records"] += n
            ref = chunk[3]
            alt = chunk[4]
            snp = ((ref.str.len() <= 1) & (alt.str.len() <= 1)).to_numpy()
            counts["snps"] += int(snp.sum())
            fmt = chunk[8] if has_format else pd.Series([""] * n, index=chunk.index)
            fmt_tag_index = {}
            for f in pd.unique(fmt[snp]) if snp.any() else []:
                keys = f.split(":")
                fmt_tag_index[f] = keys.index(tag) if tag in keys else -1
            if snp.any():
                tag_pos = fmt[snp].map(fmt_tag_index).to_numpy()
                counts["tagged"] += int((tag_pos >= 0).sum())
            sel = snp.copy()
            if keep_row_mask is not None:
                line = offset + np.arange(n)
                inside = line < len(keep_row_mask)
                inside[inside] = keep_row_mask[line[inside]]
                sel &= inside
            ids = None
            if sel.any():
                sub = chunk[sel]
                ids = (sub[0] + "_" + sub[1] + "_" + sub[3] + "_" + sub[4])
                if keep is not None:
                    inside = keep.get_indexer(ids.to_numpy(dtype=object)) >= 0
                    idx = np.flatnonzero(sel)
                    sel[idx[~inside]] = False
                    ids = ids[inside]
            if sel.any():
                sub = chunk[sel]
                counts["kept"] += int(sel.sum())
                k_of_row = fmt[sel].map(fmt_tag_index).to_numpy()
                codes = np.empty((len(sub), len(samples)), dtype=object)
                for j, col in enumerate(sample_cols):
                    values = sub[col]
                    out = np.full(len(sub), ".", dtype=object)
                    for k in np.unique(k_of_row):
                        if k < 0:
                            continue
                        rows = np.flatnonzero(k_of_row == k)
                        part = _subfield(values.iloc[rows], int(k))
                        absent = part.isna().to_numpy()
                        if absent.any():
                            r = rows[np.flatnonzero(absent)[0]]
                            raise ValueError("%s: sample %s has no %s value in %r (FORMAT %s); vireo stops here"
                                             % (path, samples[j], tag, values.iloc[r], fmt[sel].iloc[r]))
                        out[rows] = part.to_numpy(dtype=object)
                    codes[:, j] = out
                parts["chrom"].append(sub[0].to_numpy(dtype=object))
                parts["pos"].append(sub[1].to_numpy(dtype=object))
                parts["ref"].append(sub[3].to_numpy(dtype=object))
                parts["alt"].append(sub[4].to_numpy(dtype=object))
                parts["info"].append(sub[7])
                parts["codes"].append(codes)
                parts["rows"].append(offset + np.flatnonzero(sel))
            offset += n
    finally:
        reader.close()

    if parts["chrom"]:
        chrom = np.concatenate(parts["chrom"])
        pos = np.concatenate(parts["pos"])
        ref = np.concatenate(parts["ref"])
        alt = np.concatenate(parts["alt"])
        codes = np.concatenate(parts["codes"], axis=0)
        rows = np.concatenate(parts["rows"])
        info_text = pd.concat(parts["info"], ignore_index=True)
    else:
        chrom = pos = ref = alt = np.zeros(0, dtype=object)
        codes = np.zeros((0, len(samples)), dtype=object)
        rows = np.zeros(0, dtype=np.int64)
        info_text = pd.Series([], dtype=object)
    gp, missing = genotype_probabilities(codes, tag)
    gp = gp.reshape(len(chrom), len(samples), 3)
    missing = missing.reshape(len(chrom), len(samples))
    info_values = {field: _info_values(info_text, field) for field in info}
    return Genotypes(variant_id(chrom, pos, ref, alt), chrom, pos, ref, alt, samples, gp, missing, info_values,
                     rows, header_samples, counts, tag)


def check_tag(genotypes, path=None):
    """Raise ValueError when the tag is in fewer than 10 per cent of the SNP records (vireo's check)."""
    c = genotypes.counts
    if c["tagged"] < 0.1 * c["snps"]:
        raise ValueError("%s: FORMAT %s is present in %d of %d SNP records; vireo stops when it is in fewer "
                         "than 10 per cent (try another --geno-tag)"
                         % (path or "donor VCF", genotypes.tag, c["tagged"], c["snps"]))


def chr_variants(ids):
    """The ids, plus each with "chr" added and with a leading "chr" removed (for vireo's fallbacks)."""
    out = set(ids)
    out.update("chr" + x for x in ids)
    out.update(x[3:] for x in ids if x.startswith("chr"))
    return out


def match_sites(cell_ids, donor_ids):
    """vireo's matching of pileup rows to donor records (io_utils.match_donor_VCF).

    Returns (idx_cell, idx_donor, how): see vireo_port.match_variants. When
    the pileup repeats an id, one copy is matched: the one vireo matches,
    which numpy's sort of the whole list decides. The same holds for a donor
    id that repeats, so `donor_ids` must be every SNP record vireo would read
    when any id among them repeats.
    """
    from demux_tools.vireo_port import match_variants

    return match_variants(cell_ids, donor_ids)


def read_matched(path, cell_ids, samples=None, tag="GT", info=(), match="auto", check=True):
    """Donor genotypes matched to pileup ids as vireo matches them; returns a dict.

    Reads the records whose id is a pileup id (or, when none matches, differs
    from one by a "chr" prefix) and matches them with `match_sites`. When an
    id among those records repeats, which copy vireo pairs depends on numpy's
    sort of every SNP id it read, so the match is then made against every SNP
    record of the file (`match="all"` does this always; without repeats both
    give the same rows). check: apply vireo's 10 per cent tag rule.
    Keys: genotypes (Genotypes of the records read), idx_cell (matched pileup
    rows, increasing), idx_donor (the record of each, an index into
    genotypes), how (see vireo_port.match_variants), repeated (repeated ids
    among the records at piled sites), matched_on ("records at piled sites"
    or "every SNP record").
    """
    if match not in ("auto", "all"):
        raise ValueError("match must be auto or all, not %r" % (match,))
    geno = read_genotypes(path, samples=samples, tag=tag, keep_ids=cell_ids, info=info)
    if check:
        check_tag(geno, path)
    idx_cell, idx_donor, how = match_sites(cell_ids, geno.variants)
    if len(idx_cell) == 0:
        # vireo's "chr" retries: read the records whose id differs by a "chr" prefix too
        geno = read_genotypes(path, samples=samples, tag=tag, keep_ids=chr_variants(cell_ids), info=info)
        idx_cell, idx_donor, how = match_sites(cell_ids, geno.variants)
    repeated = int(len(geno.variants) - len(set(geno.variants)))
    matched_on = "records at piled sites"
    if repeated or match == "all":
        all_ids, all_rows = snp_record_ids(path)
        idx_cell, idx_all, how = match_sites(cell_ids, all_ids)
        wanted = all_rows[idx_all]
        geno = read_genotypes(path, samples=samples, tag=tag, keep_rows=wanted, info=info)
        position = {int(r): k for k, r in enumerate(geno.rows)}
        idx_donor = np.array([position[int(r)] for r in wanted], dtype=np.int64)
        matched_on = "every SNP record"
    return {"genotypes": geno, "idx_cell": idx_cell, "idx_donor": np.asarray(idx_donor, dtype=np.int64),
            "how": how, "repeated": repeated, "matched_on": matched_on}
