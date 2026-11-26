"""Site lists, masks, thinning, annotation and site groups for the Scorer; the `sites` command.

A site is a pileup row, named by vireo's id "CHROM_POS_REF_ALT". The Scorer
scores a subset of sites by keeping those rows of the pileup and the genotype
array, so every function here returns a boolean mask over a capture's sites,
or a label per site.

Site lists are tables with a header: a `variant` column of ids, or `chrom`
and `pos` columns (with `ref` and `alt` when the alleles should match too). A
VCF works as a site list too: its records' CHROM, POS, REF and ALT.

This module needs numpy and pandas only.
"""
import bisect
import os
import sys

import numpy as np

from demux_tools._io import UsageError, is_gzip, open_text, read_table, start_clock, write_provenance, write_table
from demux_tools.cellsnp import read_base, variant_id

REGIONS = ("exonic", "intronic", "intergenic")


# ---------------------------------------------------------------------------
# Site lists and masks


def _is_vcf(path):
    with open_text(path) as f:
        first = f.readline()
    return first.startswith("##fileformat=VCF") or first.startswith("#CHROM")


def read_site_list(path):
    """Keys of a site list: ("variant", set of ids) or ("position", set of "CHROM:POS").

    A table with a `variant` column gives ids. A table with `chrom`, `pos`,
    `ref` and `alt` gives ids built from them; with `chrom` and `pos` only it
    gives positions. A VCF gives the ids of its records.
    """
    path = os.fspath(path)
    if _is_vcf(path):
        chrom, pos, ref, alt = read_base(path)
        return "variant", set(variant_id(chrom, pos, ref, alt))
    table = read_table(path)
    cols = set(table.columns)
    if "variant" in cols:
        return "variant", set(table["variant"])
    if {"chrom", "pos", "ref", "alt"} <= cols:
        return "variant", set(variant_id(table["chrom"], table["pos"], table["ref"], table["alt"]))
    if {"chrom", "pos"} <= cols:
        return "position", set(c + ":" + p for c, p in zip(table["chrom"], table["pos"]))
    raise UsageError("%s: a site list needs a variant column, or chrom and pos" % path)


def _positions(variants):
    """"CHROM:POS" of each vireo id (CHROM may contain "_"; the last three fields are POS, REF, ALT)."""
    out = []
    for v in variants:
        parts = v.rsplit("_", 3)
        out.append(parts[0] + ":" + parts[1])
    return out


def split_ids(variants):
    """(chrom, pos int64, ref, alt) arrays from vireo ids."""
    parts = [v.rsplit("_", 3) for v in variants]
    chrom = np.array([p[0] for p in parts], dtype=object)
    pos = np.array([int(p[1]) for p in parts], dtype=np.int64)
    ref = np.array([p[2] for p in parts], dtype=object)
    alt = np.array([p[3] for p in parts], dtype=object)
    return chrom, pos, ref, alt


def in_list(variants, site_list):
    """Boolean mask: which of `variants` are in a site list (path or the tuple from read_site_list)."""
    kind, keys = read_site_list(site_list) if not isinstance(site_list, tuple) else site_list
    names = list(variants) if kind == "variant" else _positions(variants)
    return np.array([n in keys for n in names], dtype=bool)


def mask(variants, keep=None, exclude=None):
    """Sites to score: in `keep` (all when None) and not in `exclude` (site lists or masks)."""
    n = len(variants)
    out = np.ones(n, dtype=bool)
    if keep is not None:
        out &= keep if isinstance(keep, np.ndarray) else in_list(variants, keep)
    if exclude is not None:
        out &= ~(exclude if isinstance(exclude, np.ndarray) else in_list(variants, exclude))
    return out


def thin(variants, window_bp, keep="first", depth=None):
    """Boolean mask that keeps no two sites within `window_bp` of each other on a chromosome.

    keep="first": walk each chromosome left to right and keep a site when it
    lies more than window_bp after the last kept one. keep="deepest": take
    sites by decreasing `depth` (UMIs over cells; ties by position) and keep a
    site when no kept site lies within window_bp. Sites closer than a read
    length often count the same molecules, so a cell's evidence there is
    counted more than once.
    """
    window_bp = int(window_bp)
    chrom, pos, _, _ = split_ids(variants)
    n = len(chrom)
    out = np.zeros(n, dtype=bool)
    if keep not in ("first", "deepest"):
        raise UsageError("thin keep must be first or deepest, not %r" % keep)
    if keep == "deepest" and depth is None:
        raise UsageError("thin keep=deepest needs a depth per site")
    for c in dict.fromkeys(chrom):
        idx = np.flatnonzero(chrom == c)
        if keep == "first":
            order = idx[np.argsort(pos[idx], kind="stable")]
            last = None
            for i in order:
                if last is None or pos[i] - last > window_bp:
                    out[i] = True
                    last = pos[i]
        else:
            d = np.asarray(depth)[idx]
            order = idx[np.lexsort((pos[idx], -d))]
            kept = []
            for i in order:
                p = pos[i]
                k = bisect.bisect_left(kept, p)
                if (k < len(kept) and kept[k] - p <= window_bp) or (k > 0 and p - kept[k - 1] <= window_bp):
                    continue
                kept.insert(k, p)
                out[i] = True
    return out


# ---------------------------------------------------------------------------
# GTF annotation


def read_gtf(path, features=("gene", "exon")):
    """Gene and exon records of a GTF as a DataFrame: feature, chrom, start, end, gene_id, gene_name."""
    import pandas as pd

    path = os.fspath(path)
    df = pd.read_csv(path, sep="\t", header=None, comment="#", usecols=[0, 2, 3, 4, 8],
                     names=["chrom", "feature", "start", "end", "attributes"],
                     dtype={"chrom": str, "feature": str, "start": np.int64, "end": np.int64, "attributes": str},
                     compression="gzip" if is_gzip(path) else None, engine="c")
    df = df[df["feature"].isin(features)].copy()
    df["gene_id"] = df["attributes"].str.extract(r'gene_id "([^"]*)"')[0].fillna("")
    df["gene_name"] = df["attributes"].str.extract(r'gene_name "([^"]*)"')[0]
    df["gene_name"] = df["gene_name"].fillna(df["gene_id"])
    return df.drop(columns=["attributes"]).reset_index(drop=True)


def annotate_regions(variants, gtf):
    """Region and genes of each site: DataFrame variant, region, gene_id, gene_name, n_genes.

    region: exonic (inside an exon of any gene), intronic (inside a gene but
    no exon), intergenic. gene_id, gene_name: of the gene whose exon holds the
    site, else of the overlapping gene; when several qualify, the one that
    starts first (then the smaller gene_id). n_genes: genes overlapping the
    site. GTF coordinates are 1-based and inclusive. Chromosome names must
    match the sites' exactly; when none does, a warning is printed.
    `gtf` is a path or the DataFrame of `read_gtf`.
    """
    import pandas as pd

    g = read_gtf(gtf) if not hasattr(gtf, "columns") else gtf
    genes = g[g["feature"] == "gene"]
    gene_start = dict(zip(genes["gene_id"], genes["start"].astype(np.int64)))
    name_of = dict(zip(g["gene_id"], g["gene_name"]))
    chrom, pos, _, _ = split_ids(variants)
    n = len(chrom)
    region = np.array(["intergenic"] * n, dtype=object)
    gene_id = np.array([""] * n, dtype=object)
    gene_name = np.array([""] * n, dtype=object)
    n_genes = np.zeros(n, dtype=np.int64)
    by_chrom = {c: sub for c, sub in g.groupby("chrom", sort=False)}
    if n and not set(chrom) & set(by_chrom):
        sys.stderr.write("demux_tools: warning: no site's chromosome (e.g. %s) is named in the GTF (e.g. %s); every "
                         "site is labelled intergenic\n" % (chrom[0], next(iter(by_chrom), "none")))
    for c in dict.fromkeys(chrom):
        sub = by_chrom.get(c)
        if sub is None:
            continue
        # Sweep along the chromosome. At one position, interval starts come first,
        # then sites, then interval ends, so both ends are inclusive.
        events = []
        for feature, start, end, gid in zip(sub["feature"], sub["start"], sub["end"], sub["gene_id"]):
            kind = 0 if feature == "gene" else 1
            events.append((int(start), 0, kind, gid))
            events.append((int(end), 2, kind, gid))
        for i in np.flatnonzero(chrom == c):
            events.append((int(pos[i]), 1, -1, int(i)))
        events.sort(key=lambda e: (e[0], e[1]))
        open_sets = ({}, {})
        for _, order, kind, key in events:
            if order == 0:
                open_sets[kind][key] = open_sets[kind].get(key, 0) + 1
            elif order == 2:
                left = open_sets[kind][key] - 1
                if left:
                    open_sets[kind][key] = left
                else:
                    del open_sets[kind][key]
            else:
                in_genes, in_exons = open_sets
                n_genes[key] = len(in_genes)
                pool = in_exons or in_genes
                if pool:
                    region[key] = "exonic" if in_exons else "intronic"
                    gid = min(pool, key=lambda k: (gene_start.get(k, 0), k))
                    gene_id[key], gene_name[key] = gid, name_of.get(gid, gid)
    return pd.DataFrame({"variant": list(variants), "region": region, "gene_id": gene_id, "gene_name": gene_name,
                         "n_genes": n_genes})


# ---------------------------------------------------------------------------
# Site groups


def info_bins(values, edges):
    """Labels "[a,b)" for each value by the bin edges (the last bin closed), "NA" for NaN or outside."""
    edges = [float(e) for e in edges]
    labels = ["[%g,%g)" % (a, b) for a, b in zip(edges[:-1], edges[1:])]
    labels[-1] = labels[-1][:-1] + "]"
    out = np.array(["NA"] * len(values), dtype=object)
    v = np.asarray(values, dtype=np.float64)
    for k, (a, b) in enumerate(zip(edges[:-1], edges[1:])):
        inside = (v >= a) & ((v < b) if k < len(labels) - 1 else (v <= b))
        out[inside] = labels[k]
    return out


def missing_pattern(missing, donors, top=10):
    """Per site, which donors are blank: "none", or donor names joined by "+".

    Patterns beyond the `top` most common are labelled "other".
    """
    import pandas as pd

    names = np.array(donors, dtype=object)
    labels = np.array(["+".join(names[row]) if row.any() else "none" for row in missing], dtype=object)
    if top is not None and len(labels):
        common = set(pd.Series(labels).value_counts().index[:top])
        common.add("none")
        labels = np.array([x if x in common else "other" for x in labels], dtype=object)
    return labels


def groups_from(spec, variants, missing=None, donors=None, info=None, site_table=None, top=10):
    """(name, labels) for one grouping of the sites.

    spec: "@missing" (which donors are blank), "@chrom", "@info:FIELD:e0,e1,...,ek"
    (bins of an INFO field), or a column of `site_table` (a DataFrame with a
    `variant` column; sites not in it get "NA").
    """
    if spec == "@missing":
        if missing is None:
            raise UsageError("@missing needs the genotype blanks")
        return "missing", missing_pattern(missing, donors, top=top)
    if spec == "@chrom":
        chrom, _, _, _ = split_ids(variants)
        return "chrom", chrom
    if spec.startswith("@info:"):
        parts = spec.split(":")
        if len(parts) != 3:
            raise UsageError("@info grouping is @info:FIELD:e0,e1,...; not %r" % spec)
        field, edges = parts[1], parts[2].split(",")
        if info is None or field not in info:
            raise UsageError("INFO field %s was not loaded" % field)
        return field, info_bins(info[field], edges)
    if site_table is None:
        raise UsageError("grouping by column %r needs --site-groups" % spec)
    if spec not in site_table.columns:
        raise UsageError("the site table has no column %r" % spec)
    lookup = dict(zip(site_table["variant"], site_table[spec]))
    return spec, np.array([lookup.get(v, "NA") for v in variants], dtype=object)


# ---------------------------------------------------------------------------
# Command: sites


def configure_sites(parser):
    parser.description = (
        "Make and describe site lists for the Scorer. annotate: one row per site with its region "
        "(exonic, intronic, intergenic) and gene from a GTF, and INFO fields from a panel. thin: keep no two "
        "sites within a window. from-panel: the sites of a VCF, optionally filtered on INFO fields. Each "
        "writes a table whose first column is the vireo id `variant`, usable as --sites or --exclude-sites.")
    sub = parser.add_subparsers(dest="action")
    sub.required = True
    a = sub.add_parser("annotate", help="region, gene and INFO fields per site")
    _site_source(a)
    a.add_argument("--gtf", help="GTF with gene and exon records (gzip or plain)")
    a.add_argument("--panel", help="VCF to take INFO fields from, matched by variant id")
    a.add_argument("--info", default="", help="comma-separated INFO fields to add (from --panel)")
    a.add_argument("--out", required=True)
    t = sub.add_parser("thin", help="keep no two sites within --window-bp")
    _site_source(t)
    t.add_argument("--window-bp", type=int, required=True, help="sites closer than this are thinned")
    t.add_argument("--keep", default="first", choices=("first", "deepest"),
                   help="first: left to right; deepest: by UMIs over cells (needs the DP matrix)")
    t.add_argument("--dp", help="cellSNP.tag.DP.mtx, for --keep deepest (default: the folder's)")
    t.add_argument("--out", required=True)
    f = sub.add_parser("from-panel", help="the sites of a VCF, filtered on INFO fields")
    f.add_argument("--vcf", required=True)
    f.add_argument("--min", action="append", default=[], metavar="FIELD=X", help="keep sites with FIELD >= X")
    f.add_argument("--max", action="append", default=[], metavar="FIELD=X", help="keep sites with FIELD <= X")
    f.add_argument("--snps-only", action="store_true", help="drop records vireo drops (REF or ALT longer than 1)")
    f.add_argument("--out", required=True)


def _site_source(p):
    p.add_argument("--cellsnp", help="cellsnp-lite folder: its base VCF's sites")
    p.add_argument("--base", help="cellSNP.base.vcf.gz (or any VCF) to take sites from")


def _source_sites(args):
    from demux_tools.cellsnp import resolve_cellsnp

    if args.base:
        path = args.base
    elif args.cellsnp:
        path = resolve_cellsnp(cellsnp=args.cellsnp)["base"]
    else:
        raise UsageError("give --cellsnp or --base")
    chrom, pos, ref, alt = read_base(path)
    return path, variant_id(chrom, pos, ref, alt), (chrom, pos, ref, alt)


def _site_frame(variants, parts):
    import pandas as pd

    chrom, pos, ref, alt = parts
    return pd.DataFrame({"variant": list(variants), "chrom": chrom, "pos": pos, "ref": ref, "alt": alt})


def _bounds(items):
    out = []
    for item in items:
        if "=" not in item:
            raise UsageError("expected FIELD=X, not %r" % item)
        field, value = item.split("=", 1)
        out.append((field, float(value)))
    return out


def run_sites(args):
    from demux_tools import vcf as vcf_mod

    clock = start_clock()
    inputs = []
    if args.action == "annotate":
        path, variants, parts = _source_sites(args)
        inputs.append(path)
        table = _site_frame(variants, parts)
        if args.gtf:
            inputs.append(args.gtf)
            regions = annotate_regions(variants, args.gtf)
            for col in ("region", "gene_id", "gene_name", "n_genes"):
                table[col] = regions[col].to_numpy()
        fields = [f for f in args.info.split(",") if f]
        if fields:
            if not args.panel:
                raise UsageError("--info needs --panel")
            inputs.append(args.panel)
            g = vcf_mod.read_genotypes(args.panel, samples=[], keep_ids=variants, info=fields)
            index = {v: k for k, v in enumerate(g.variants)}
            rows = np.array([index.get(v, -1) for v in variants])
            for f in fields:
                values = np.where(rows >= 0, g.info[f][np.maximum(rows, 0)], np.nan) if len(rows) else np.zeros(0)
                table[f] = values
            table["in_panel"] = np.where(rows >= 0, "true", "false")
        write_table(table, args.out)
        summary = {"sites": len(table)}
        if "region" in table.columns:
            summary["regions"] = {r: int((table["region"] == r).sum()) for r in REGIONS}
    elif args.action == "thin":
        path, variants, parts = _source_sites(args)
        inputs.append(path)
        depth = None
        if args.keep == "deepest":
            from demux_tools.cellsnp import read_mtx, resolve_cellsnp

            dp = args.dp or (resolve_cellsnp(cellsnp=args.cellsnp)["DP"] if args.cellsnp else None)
            if not dp:
                raise UsageError("--keep deepest needs --dp or --cellsnp")
            inputs.append(dp)
            shape, i, _, v = read_mtx(dp)
            depth = np.bincount(i, weights=v, minlength=shape[0])
        keep = thin(variants, args.window_bp, keep=args.keep, depth=depth)
        table = _site_frame(variants, parts)[keep]
        write_table(table, args.out)
        summary = {"sites": int(len(variants)), "kept": int(keep.sum()), "window_bp": args.window_bp}
    else:
        inputs.append(args.vcf)
        mins, maxs = _bounds(args.min), _bounds(args.max)
        fields = sorted(set(f for f, _ in mins + maxs))
        if args.snps_only:
            g = vcf_mod.read_genotypes(args.vcf, samples=[], info=fields)
            variants, info, parts = g.variants, g.info, (g.chrom, g.pos, g.ref, g.alt)
        else:
            chrom, pos, ref, alt = read_base(args.vcf)
            variants, parts = variant_id(chrom, pos, ref, alt), (chrom, pos, ref, alt)
            info = _all_records_info(args.vcf, fields) if fields else {}
        keep = np.ones(len(variants), dtype=bool)
        for f, x in mins:
            keep &= info[f] >= x
        for f, x in maxs:
            keep &= info[f] <= x
        table = _site_frame(variants, parts)
        for f in fields:
            table[f] = info[f]
        table = table[keep]
        write_table(table, args.out)
        summary = {"records": int(len(variants)), "kept": int(keep.sum())}
    write_provenance(args.out, args.argv, inputs, clock=clock, extra=summary)
    print("sites %s: %s" % (args.action, ", ".join("%s %s" % kv for kv in sorted(summary.items()))))
    return 0


def _all_records_info(path, fields):
    """INFO fields of every record of a VCF, in file order (floats, NaN where absent)."""
    import csv

    import pandas as pd

    from demux_tools.vcf import _info_values, read_header

    n_header, _ = read_header(path)
    df = pd.read_csv(path, sep="\t", header=None, skiprows=n_header, usecols=[7], dtype=str, keep_default_na=False,
                     na_filter=False, quoting=csv.QUOTE_NONE, compression="gzip" if is_gzip(path) else None)
    return {f: _info_values(df[7], f) for f in fields}
