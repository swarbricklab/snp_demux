"""The Scorer: vireo's fixed-mode calls as arithmetic, with each cell's evidence kept.

In fixed mode (as many donors as genotyped, genotypes not learnt) vireo
scores a cell against a donor by adding, over the sites where the cell has
reads, the expected log-likelihood of its ALT and REF UMIs under the donor's
genotype there. Three allele ratios (the ALT share expected at a 0/0, 0/1 and
1/1 site) are fitted per capture first. A doublet is scored the same way with
the two donors' genotypes mixed 50:50. The Scorer:

1. loads a capture: the pileup rows that match a SNP record of the donor
   VCF, with the donors' genotype probabilities, exactly as vireo builds
   them (`load_capture`; a numpy reader and an optional cache by default,
   vireoSNP's own readers with loader="vireo");
2. fits the allele ratios with vireoSNP's own code, called as vireo calls it
   in fixed mode (`fit`, backend "vireo");
3. recomputes the singlet and doublet log-likelihood matrices with the same
   arithmetic as vireo (`loglik_singlet`, `loglik_doublet`, `posteriors`),
   checks them against vireo's probabilities, and writes vireo's eight
   donor_ids.tsv columns plus nine columns of evidence per cell (`score`).

Under backend "vireo" the eight columns are vireo's own numbers, and the
internal check (`maxdiff_*` in fit.json) shows the evidence columns come
from the same arithmetic. scipy and vireoSNP are imported inside the
functions that need them, so the commands' --help needs neither.

A variant (`apply_variant`) scores a capture again with changed sites
(`subset_sites`), changed genotypes (`fill_missing`, `harden`, `add_donors`)
or held allele ratios, and adds the columns of `diagnostics`. backend
"native" fits without vireoSNP. Without these options the output is the
same as before they existed.
"""
import contextlib
import io
import os
import sys
import time

import numpy as np

from demux_tools import cache as cache_mod
from demux_tools._io import (UsageError, check_output_path, peak_rss_mb, provenance, read_lines, read_table,
                             start_clock, write_json, write_table)
from demux_tools.cellsnp import read_base, read_mtx, read_samples, resolve_cellsnp, variant_id
from demux_tools.vireo_files import DONOR_IDS_COLUMNS, pair_key, read_donor_ids, read_prob_header

DIAGNOSTIC_COLUMNS = ["best_pair", "second_donor", "umis", "ll_best", "ll_second", "margin", "ll_best_pair",
                      "doublet_margin", "fit_per_umi"]
CELL_COLUMNS = DONOR_IDS_COLUMNS + DIAGNOSTIC_COLUMNS
FLOAT_FORMATS = {"prob_max": "%.2e", "prob_doublet": "%.2e", "doublet_logLikRatio": "%.3f", "ll_best": "%.4f",
                 "ll_second": "%.4f", "margin": "%.4f", "ll_best_pair": "%.4f", "doublet_margin": "%.4f",
                 "fit_per_umi": "%.5f"}
VIREO_FORMATS = ("0.5.7", "0.5.9", "0.5.10")
GT_TEXT = ("0/0", "0/1", "1/1")
_UNIFORM_CODE = 3


def _pair_naming(vireo_format):
    return "permutations" if vireo_format == "0.5.7" else "combinations"


def _prob_doublet_column(vireo_format):
    return "prob_max_doublet" if vireo_format == "0.5.10" else "prob_doublet"


# ---------------------------------------------------------------------------
# Capture


class Capture(object):
    """One capture as vireo sees it in fixed mode.

    Attributes:
        name: the capture's name (labels only).
        barcodes: list of str, one per column.
        variants: list of vireo ids "CHROM_POS_REF_ALT", one per site used.
        AD, DP: scipy.sparse.csc_matrix int64 (sites, cells), built as vireo
            builds them (rows matched to the donor VCF, in pileup order).
        donors: donor names in the order used (vireo's order).
        GPb: float64 (sites, donors, 3), genotype probabilities as vireo
            parses them.
        missing: bool (sites, donors), True where the genotype was ".",
            "./." or ".|." (a third each).
        pileup_rows: int64, each site's 0-based row in the pileup files.
        vcf_rows: int64, each site's 0-based record in the donor VCF.
        info: {INFO field: float64 array}, one value per site.
        counts: dict of site, record and cell counts from loading.
        inputs: input file paths, for provenance.
        load_info: {"loader", "cache" ("hit", "miss" or "off"), "seconds"}.
    """

    def __init__(self, name, barcodes, variants, AD, DP, donors, GPb, missing, pileup_rows, vcf_rows, info,
                 counts, inputs, load_info):
        self.name = name
        self.barcodes = list(barcodes)
        self.variants = list(variants)
        self.AD = AD
        self.DP = DP
        self.donors = list(donors)
        self.GPb = GPb
        self.missing = missing
        self.pileup_rows = np.asarray(pileup_rows, dtype=np.int64)
        self.vcf_rows = np.asarray(vcf_rows, dtype=np.int64)
        self.info = info
        self.counts = counts
        self.inputs = list(inputs)
        self.load_info = load_info

    @property
    def n_cells(self):
        return self.AD.shape[1]

    @property
    def n_sites(self):
        return self.AD.shape[0]

    @property
    def n_donors(self):
        return len(self.donors)


def _donor_list(header_samples, donors, donor_order):
    """The donors to use, in the order to use them.

    donors: None (every sample of the VCF) or a list of names.
    donor_order: "vcf" (the VCF header's order, as vireo uses it), "list"
    (the order of `donors`), "sorted" (string order) or "prob-header:PATH"
    (the order of vireo's prob_singlet.tsv.gz header; it also sets the donors
    when `donors` is None).
    """
    order = donor_order or "vcf"
    from_header = None
    if order.startswith("prob-header:"):
        from_header = read_prob_header(order[len("prob-header:"):])
        if donors is None:
            donors = from_header
    chosen = list(header_samples) if donors is None else list(donors)
    if len(set(chosen)) != len(chosen):
        raise UsageError("a donor is named twice: %s" % ",".join(chosen))
    absent = [d for d in chosen if d not in header_samples]
    if absent:
        raise UsageError("donor(s) not in the VCF: %s" % ",".join(absent))
    if order == "vcf":
        return [s for s in header_samples if s in set(chosen)]
    if order == "list":
        return chosen
    if order == "sorted":
        return sorted(chosen)
    if from_header is not None:
        if sorted(from_header) != sorted(chosen):
            raise UsageError("the prob_singlet header lists %s, the donors chosen are %s"
                             % (",".join(from_header), ",".join(chosen)))
        return list(from_header)
    raise UsageError("donor order must be vcf, list, sorted or prob-header:PATH, not %r" % donor_order)


def _barcode_list(barcodes):
    if barcodes is None:
        return None
    if isinstance(barcodes, str):
        return [b.strip() for b in read_lines(barcodes) if b.strip()]
    return list(barcodes)


def _column_selection(all_barcodes, wanted):
    """Column indices of `wanted` barcodes, in pileup order (None for all)."""
    if wanted is None:
        return None
    position = {b: i for i, b in enumerate(all_barcodes)}
    absent = [b for b in wanted if b not in position]
    if absent:
        raise UsageError("%d barcode(s) not in the pileup, e.g. %s" % (len(absent), absent[0]))
    return np.array(sorted(set(position[b] for b in wanted)), dtype=np.int64)


def _csc(rows, cols, vals, shape):
    from scipy import sparse

    return sparse.csc_matrix((np.asarray(vals, dtype=np.int64), (rows, cols)), shape=shape)


def _encode_gp(gp, missing, tag):
    """Compact form of GPb for the cache: int8 codes for GT (exact), else the floats."""
    if tag != "GT":
        return {"gp": gp}
    codes = np.where(missing, _UNIFORM_CODE, gp.argmax(axis=2)).astype(np.int8)
    if not np.array_equal(_decode_gp({"gt_codes": codes}, "GT"), gp):
        return {"gp": gp}
    return {"gt_codes": codes}


def _decode_gp(arrays, tag):
    from demux_tools.vireo_port import parse_gt_code

    if "gp" in arrays:
        return np.asarray(arrays["gp"], dtype=np.float64)
    codes = arrays["gt_codes"]
    table = np.zeros((4, 3), dtype=np.float64)
    for g in range(3):
        table[g, g] = 1.0
    table[_UNIFORM_CODE] = parse_gt_code(".", "GT")
    return table[codes.astype(np.int64)]


def _store(capture, tag):
    arrays = {
        "barcodes": np.array(capture.barcodes, dtype=str),
        "variants": np.array(capture.variants, dtype=str),
        "donors": np.array(capture.donors, dtype=str),
        "missing": capture.missing,
        "pileup_rows": capture.pileup_rows,
        "vcf_rows": capture.vcf_rows,
    }
    for name in ("AD", "DP"):
        m = getattr(capture, name)
        arrays[name + "_data"] = m.data
        arrays[name + "_indices"] = m.indices
        arrays[name + "_indptr"] = m.indptr
    arrays.update(_encode_gp(capture.GPb, capture.missing, tag))
    for field, values in capture.info.items():
        arrays["info__" + field] = values
    return arrays


def _restore(name, arrays, meta, inputs, tag):
    from scipy import sparse

    shape = tuple(meta["shape"])
    mats = {}
    for m in ("AD", "DP"):
        mats[m] = sparse.csc_matrix((arrays[m + "_data"].astype(np.int64), arrays[m + "_indices"],
                                     arrays[m + "_indptr"]), shape=shape)
    info = {k[len("info__"):]: v for k, v in arrays.items() if k.startswith("info__")}
    return Capture(name, arrays["barcodes"].tolist(), arrays["variants"].tolist(), mats["AD"], mats["DP"],
                   arrays["donors"].tolist(), _decode_gp(arrays, tag), arrays["missing"].astype(bool),
                   arrays["pileup_rows"], arrays["vcf_rows"], info, meta["counts"], inputs, {})


def load_capture(name, cellsnp=None, base=None, samples=None, ad=None, dp=None, vcf=None, donors=None,
                 donor_order="vcf", geno_tag="GT", barcodes=None, info=(), cache_dir=None, loader="native",
                 match="auto"):
    """Load one capture for scoring; returns a Capture.

    cellsnp: a cellsnp-lite folder; base, samples, ad, dp: its files, which
        override the folder's (so DVC objects can be read in place).
    vcf: the donor VCF (any size; only records at piled sites are parsed).
    donors, donor_order: see `_donor_list`; vireo uses the VCF's order.
    geno_tag: GT, GP or PL. barcodes: a file of barcodes, or a list, to keep
        (pileup order is kept). info: INFO fields to keep per site.
    cache_dir: keep the loaded arrays there; a later load with the same
        inputs (path, size, modification time) and options reads them back.
        A cache that cannot be written gives a warning, not an error.
    loader: "native" (pandas readers) or "vireo" (vireoSNP's read_cellSNP,
        load_VCF, match_donor_VCF and parse_donor_GPb; no cache).
    match: "auto" matches pileup rows against the VCF records at piled sites,
        and against every SNP record of the VCF when an id among the former
        repeats (vireo pairs repeated ids by numpy's sort of the whole list);
        "all" always matches against every SNP record. Both give vireo's rows.
    Raises UsageError or ValueError when no site matches or an input is bad.
    """
    t0 = time.time()
    if vcf is None:
        raise UsageError("a donor VCF is needed")
    paths = resolve_cellsnp(cellsnp=cellsnp, base=base, samples=samples, ad=ad, dp=dp)
    for layer in ("AD", "DP"):
        if paths[layer] is None:
            raise UsageError("cellsnp %s matrix not given and not found" % layer)
    order_file = donor_order[len("prob-header:"):] if (donor_order or "").startswith("prob-header:") else None
    barcode_file = barcodes if isinstance(barcodes, str) else None
    inputs = [paths["base"], paths["samples"], paths["AD"], paths["DP"], os.fspath(vcf)]
    inputs += [p for p in (barcode_file, order_file) if p]
    info = tuple(info or ())
    if loader == "vireo":
        capture = _load_with_vireo(name, paths, vcf, donors, donor_order, geno_tag, barcodes, info, inputs)
        capture.load_info = {"loader": "vireo", "cache": "off", "seconds": round(time.time() - t0, 3)}
        return capture
    if loader != "native":
        raise UsageError("loader must be native or vireo, not %r" % loader)
    if match not in ("auto", "all"):
        raise UsageError("match must be auto or all, not %r" % match)

    options = {"tag": geno_tag, "donors": donors, "donor_order": "prob-header" if order_file else donor_order,
               "barcodes": None if barcodes is None or barcode_file else list(barcodes), "info": sorted(info),
               "match": match}
    key = None
    if cache_dir:
        key = cache_mod.key("capture", inputs, options)
        hit = cache_mod.load(cache_dir, key)
        if hit is not None:
            capture = _restore(name, hit[0], hit[1], inputs, geno_tag)
            capture.load_info = {"loader": "native", "cache": "hit", "cache_key": key,
                                 "seconds": round(time.time() - t0, 3)}
            return capture

    capture = _load_native(name, paths, vcf, donors, donor_order, geno_tag, barcodes, info, inputs, match)
    state = "miss" if cache_dir else "off"
    if cache_dir:
        meta = {"shape": list(capture.AD.shape), "counts": capture.counts, "name": name}
        if cache_mod.save(cache_dir, key, _store(capture, geno_tag), meta) is None:
            state = "miss, not written"
    capture.load_info = {"loader": "native", "cache": state, "cache_key": key,
                         "seconds": round(time.time() - t0, 3)}
    return capture


def _load_native(name, paths, vcf_path, donors, donor_order, tag, barcodes, info, inputs, match="auto"):
    from demux_tools import vcf as vcf_mod

    chrom, pos, ref, alt = read_base(paths["base"])
    cell_ids = variant_id(chrom, pos, ref, alt)
    all_barcodes = read_samples(paths["samples"])
    _, header_samples = vcf_mod.read_header(vcf_path)
    donor_list = _donor_list(header_samples, donors, donor_order)
    m = vcf_mod.read_matched(vcf_path, cell_ids, samples=donor_list, tag=tag, info=info, match=match)
    geno, idx_cell, idx_donor, how = m["genotypes"], m["idx_cell"], m["idx_donor"], m["how"]
    repeated, matched_on = m["repeated"], m["matched_on"]
    if len(idx_cell) == 0:
        raise ValueError("no pileup site matches a SNP record of %s" % vcf_path)
    n_rows = len(cell_ids)
    keep = np.zeros(n_rows, dtype=bool)
    keep[idx_cell] = True
    new_row = np.full(n_rows, -1, dtype=np.int64)
    new_row[idx_cell] = np.arange(len(idx_cell), dtype=np.int64)
    cols = _column_selection(all_barcodes, _barcode_list(barcodes))
    n_cols = len(all_barcodes)
    if cols is not None:
        new_col = np.full(n_cols, -1, dtype=np.int64)
        new_col[cols] = np.arange(len(cols), dtype=np.int64)
    mats = {}
    for layer in ("AD", "DP"):
        shape, i, j, v = read_mtx(paths[layer], keep_rows=keep)
        if shape != (n_rows, n_cols):
            raise ValueError("%s is %d x %d; the base VCF has %d records and the barcode list %d"
                             % (paths[layer], shape[0], shape[1], n_rows, n_cols))
        r = new_row[i]
        c = j.astype(np.int64)
        if cols is not None:
            c = new_col[c]
            inside = c >= 0
            r, c, v = r[inside], c[inside], v[inside]
        mats[layer] = _csc(r, c, v, (len(idx_cell), n_cols if cols is None else len(cols)))
    kept_barcodes = all_barcodes if cols is None else [all_barcodes[k] for k in cols]
    import pandas as pd

    ids = np.asarray(cell_ids, dtype=object)
    unmatched = np.flatnonzero(~keep)
    duplicate = int(pd.Series(ids[unmatched]).isin(pd.unique(ids[idx_cell])).sum())
    counts = {
        "n_sites_pileup": int(n_rows),
        "n_sites_matched": int(len(idx_cell)),
        "sites_dropped": {"not_in_panel": int(len(unmatched) - duplicate), "duplicate_row": int(duplicate)},
        "n_cells_pileup": int(n_cols),
        "n_cells": int(len(kept_barcodes)),
        "vcf_records": int(geno.counts["records"]),
        "vcf_snp_records": int(geno.counts["snps"]),
        "vcf_records_with_tag": int(geno.counts["tagged"]),
        "vcf_records_at_piled_sites": int(geno.counts["kept"]),
        "vcf_ids_repeated_at_piled_sites": repeated,
        "matched_on": matched_on,
        "match": how,
    }
    return Capture(name, kept_barcodes, [cell_ids[k] for k in idx_cell], mats["AD"], mats["DP"], donor_list,
                   geno.gp[idx_donor], geno.missing[idx_donor], idx_cell, geno.rows[idx_donor],
                   {f: v[idx_donor] for f, v in geno.info.items()}, counts, inputs, {})


def _load_with_vireo(name, paths, vcf_path, donors, donor_order, tag, barcodes, info, inputs):
    """The capture as vireoSNP's own readers build it (reference loader, slow on large panels)."""
    import shutil
    import tempfile

    from vireoSNP.utils.io_utils import match_donor_VCF, read_cellSNP
    from vireoSNP.utils.vcf_utils import load_VCF, parse_donor_GPb

    from demux_tools import vcf as vcf_mod
    from demux_tools._io import is_gzip
    from demux_tools.cellsnp import FILES
    from demux_tools.vireo_port import UNIFORM_CODES

    if info:
        raise UsageError("--info needs the native loader")
    tmp = tempfile.mkdtemp(prefix="demux_tools_vireo_")
    log = io.StringIO()
    try:
        for key in ("base", "samples", "AD", "DP"):
            os.symlink(os.path.realpath(paths[key]), os.path.join(tmp, FILES[key]))
        link = os.path.join(tmp, "donors.vcf.gz" if is_gzip(vcf_path) else "donors.vcf")
        os.symlink(os.path.realpath(vcf_path), link)
        with contextlib.redirect_stdout(log):
            cell_dat = read_cellSNP(tmp)
            all_ids = list(cell_dat["variants"])
            donor_vcf = load_VCF(link, biallelic_only=True, sparse=False, format_list=[tag])
            if donor_vcf["n_SNP_tagged"][0] < (0.1 * len(donor_vcf["GenoINFO"][tag])):
                raise ValueError("%s: too few records with FORMAT %s (vireo stops)" % (vcf_path, tag))
            n_records = len(donor_vcf["variants"])
            donor_ids_all = list(donor_vcf["variants"])
            cell_dat, donor_vcf = match_donor_VCF(cell_dat, donor_vcf)
            gp_all = parse_donor_GPb(donor_vcf["GenoINFO"][tag], tag) if len(donor_vcf["variants"]) else None
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    if gp_all is None:
        raise ValueError("no pileup site matches a SNP record of %s" % vcf_path)
    donor_list = _donor_list(donor_vcf["samples"], donors, donor_order)
    order = [donor_vcf["samples"].index(d) for d in donor_list]
    codes = np.array(donor_vcf["GenoINFO"][tag], dtype=object)[:, order]
    missing = np.isin(codes, list(UNIFORM_CODES))
    all_barcodes = [str(b) for b in np.atleast_1d(cell_dat["samples"])]
    AD, DP = cell_dat["AD"], cell_dat["DP"]
    cols = _column_selection(all_barcodes, _barcode_list(barcodes))
    if cols is not None:
        AD, DP = AD[:, cols], DP[:, cols]
        all_barcodes = [all_barcodes[k] for k in cols]
    idx_cell, _, how = vcf_mod.match_sites(all_ids, donor_ids_all)
    counts = {"n_sites_pileup": len(all_ids), "n_sites_matched": int(AD.shape[0]),
              "n_cells": len(all_barcodes), "vcf_snp_records": n_records, "match": how}
    return Capture(name, all_barcodes, cell_dat["variants"], AD.astype(np.int64), DP.astype(np.int64), donor_list,
                   gp_all[:, order, :], missing, idx_cell, np.full(AD.shape[0], -1), {}, counts, inputs,
                   {"vireo_log": log.getvalue().splitlines()})


# ---------------------------------------------------------------------------
# Fit


class Fit(object):
    """Fitted allele ratios and how they were obtained.

    theta_mean, theta_sum: float64 (1, 3), the Beta mean and concentration of
        the ALT share for genotypes 0/0, 0/1, 1/1 (vireo's beta_mu, beta_sum).
    backend: "vireo", "native", "held" (ratios held fixed by fit()), or
        "given" for ratios set by the caller.
    n_iter: iterations of vireo's two fitting passes, or None.
    elbo: vireo's final lower bound, or None.
    doublet_prior: the prior probability of a doublet, spread over the pairs.
    rv: vireo_wrap's result dict (backend vireo), else None.
    log: lines vireo printed.
    prior_is_vireos: False when doublet_prior differs from vireo's own, so
        vireo's probabilities in `rv` do not apply.
    """

    def __init__(self, theta_mean, theta_sum, backend="given", n_iter=None, elbo=None, doublet_prior=None,
                 rv=None, log=None, prior_is_vireos=True):
        self.theta_mean = np.asarray(theta_mean, dtype=np.float64).reshape(1, -1)
        self.theta_sum = np.asarray(theta_sum, dtype=np.float64).reshape(1, -1)
        self.backend = backend
        self.n_iter = n_iter
        self.elbo = elbo
        self.doublet_prior = doublet_prior
        self.rv = rv
        self.log = list(log or [])
        self.prior_is_vireos = prior_is_vireos

    def to_dict(self):
        return {"backend": self.backend, "theta_mean": self.theta_mean.ravel().tolist(),
                "theta_sum": self.theta_sum.ravel().tolist(), "n_iter": self.n_iter,
                "elbo": None if self.elbo is None else float(self.elbo), "doublet_prior": self.doublet_prior}


def default_doublet_prior(n_cells):
    """vireo's doublet prior: cells / 100,000, at most 0.5 (predict_doublet)."""
    return min(0.5, n_cells / 100000)


def _resolve_prior(capture, doublet_prior):
    if doublet_prior in (None, "auto"):
        return default_doublet_prior(capture.n_cells), True
    value = float(doublet_prior)
    if not 0.0 <= value < 1.0:
        raise UsageError("--doublet-prior must be in [0, 1), not %r" % doublet_prior)
    return value, False


def fit(capture, backend="vireo", theta=None, doublet_prior="auto"):
    """Fit the allele ratios of `capture`; returns a Fit.

    theta: None to fit, or (theta_mean, theta_sum) to hold the ratios fixed.
    backend "vireo": vireoSNP's vireo_wrap with the arguments vireo.py passes
    in fixed mode (learn_GT False, one start, no extra donors, doublets
    checked, no ambient step, one process). Its printed lines go to Fit.log.
    With one donor there are no pairs, so doublets are not checked.
    backend "native": vireo_port.fit_fixed_native, the same steps on
    (cells x donors) matrices, without vireoSNP; the probabilities then come
    from the Scorer's own matrices.
    """
    prior, auto = _resolve_prior(capture, doublet_prior)
    if theta is not None:
        return Fit(theta[0], theta[1], backend="held", doublet_prior=prior)
    if backend == "native":
        from demux_tools.vireo_port import fit_fixed_native

        res = fit_fixed_native(capture.AD, capture.DP, capture.GPb)
        return Fit(res["theta_mean"], res["theta_sum"], backend="native", n_iter=res["n_iter"], elbo=res["elbo"],
                   doublet_prior=prior)
    if backend != "vireo":
        raise UsageError("backend must be vireo or native, not %r" % backend)
    from vireoSNP.utils import vireo_wrap as vw

    lengths = []
    original = vw.Vireo

    class _Recording(original):
        # Counts the iterations of each fitting pass; the arithmetic is vireo's.
        def fit(self, *args, **kwargs):
            n0 = len(self.ELBO_)
            out = original.fit(self, *args, **kwargs)
            lengths.append(int(len(self.ELBO_) - n0 + 1))
            return out

    log = io.StringIO()
    vw.Vireo = _Recording
    try:
        with contextlib.redirect_stdout(log):
            RV = vw.vireo_wrap(capture.AD, capture.DP, n_donor=capture.n_donors, GT_prior=capture.GPb,
                               learn_GT=False, n_init=1, n_extra_donor=0, extra_donor_mode="distance",
                               check_doublet=capture.n_donors > 1, random_seed=None, ASE_mode=False,
                               check_ambient=False, nproc=1)
    finally:
        vw.Vireo = original
    return Fit(RV["theta_mean"], RV["theta_sum"], backend="vireo", n_iter=lengths, elbo=float(RV["LB_doublet"]),
               doublet_prior=prior, rv=RV, log=log.getvalue().splitlines(), prior_is_vireos=auto)


# ---------------------------------------------------------------------------
# Log-likelihoods and posteriors


def _terms(theta_mean, theta_sum):
    from demux_tools.vireo_port import digamma_terms

    return digamma_terms(np.asarray(theta_mean, dtype=np.float64), np.asarray(theta_sum, dtype=np.float64))


def loglik_singlet(AD, DP, GPb, theta_mean, theta_sum):
    """Each cell's log-likelihood under each donor (cells, donors), as vireo's predict_doublet sums it.

    GPb is normalised first, as vireo's model does with the genotypes it is
    given. A site subset is a row subset of AD, DP and GPb.
    """
    from demux_tools.vireo_port import loglik_columns, normalize

    BD = DP - AD
    return loglik_columns(AD, BD, DP, normalize(GPb), _terms(theta_mean, theta_sum))


def loglik_doublet(AD, DP, GPb, theta_mean, theta_sum, pair_chunk=16):
    """Each cell's log-likelihood under each donor pair as a 50:50 mix (cells, pairs).

    Pairs are in itertools.combinations order of the donors, as in vireo.
    The pair genotypes are built `pair_chunk` pairs at a time; the result is
    the same for any chunk size.
    """
    from demux_tools.vireo_port import add_doublet_gt, add_doublet_theta, donor_pairs, loglik_columns, normalize

    G = normalize(GPb)
    n_pairs = len(donor_pairs(G.shape[1]))
    out = np.zeros((AD.shape[1], n_pairs))
    if n_pairs == 0:
        return out
    mu_both, sum_both = add_doublet_theta(np.asarray(theta_mean, dtype=np.float64).reshape(1, -1),
                                          np.asarray(theta_sum, dtype=np.float64).reshape(1, -1))
    terms = _terms(mu_both, sum_both)
    BD = DP - AD
    step = max(1, int(pair_chunk))
    for start in range(0, n_pairs, step):
        pairs = np.arange(start, min(n_pairs, start + step))
        out[:, pairs] = loglik_columns(AD, BD, DP, add_doublet_gt(G, pairs), terms)
    return out


def posteriors(ll_s, ll_d, doublet_prior):
    """(prob_singlet, prob_doublet, doublet_LLR) from the two log-likelihood matrices (vireo's predict_doublet)."""
    from demux_tools.vireo_port import posteriors as _posteriors

    return _posteriors(ll_s, ll_d, doublet_prior)


def cell_site_terms(capture, fit_result, barcode):
    """One cell's log-likelihood under every donor, site by site (a DataFrame).

    One row per site where the cell has reads: variant, alt, ref, then for
    each donor `gt__<donor>` (its genotype: 0/0, 0/1, 1/1, ./. when blank, or
    the three probabilities when soft) and `term__<donor>`, the expected
    log-probability of the cell's UMIs there under that donor. Each term
    column sums to the cell's `loglik_singlet` value for that donor.
    """
    import pandas as pd

    from demux_tools.vireo_port import normalize

    if barcode not in capture.barcodes:
        raise UsageError("barcode %s is not in the capture" % barcode)
    j = capture.barcodes.index(barcode)
    a = np.asarray(capture.AD[:, j].todense()).ravel()
    d = np.asarray(capture.DP[:, j].todense()).ravel()
    idx = np.flatnonzero(d)
    dg1, dg2, dgs = _terms(fit_result.theta_mean, fit_result.theta_sum)
    e_alt = (dg1 - dgs).ravel()
    e_ref = (dg2 - dgs).ravel()
    per_gt = a[idx, None] * e_alt[None, :] + (d[idx] - a[idx])[:, None] * e_ref[None, :]
    G = normalize(capture.GPb[idx])
    terms = np.einsum("vdg,vg->vd", G, per_gt)
    table = {"cell": [barcode] * len(idx), "variant": [capture.variants[k] for k in idx],
             "alt": a[idx], "ref": d[idx] - a[idx]}
    for k, donor in enumerate(capture.donors):
        table["gt__" + donor] = genotype_text(capture.GPb[idx, k], capture.missing[idx, k])
    for k, donor in enumerate(capture.donors):
        table["term__" + donor] = terms[:, k]
    return pd.DataFrame(table)


def genotype_text(gp, missing):
    """Genotype labels: 0/0, 0/1, 1/1 for a hard genotype, ./. when blank, else the three probabilities."""
    out = []
    for p, m in zip(gp, missing):
        if m:
            out.append("./.")
        elif p.max() == 1.0:
            out.append(GT_TEXT[int(p.argmax())])
        else:
            out.append(",".join("%.3g" % x for x in p))
    return out


# ---------------------------------------------------------------------------
# Scores


class Scores(object):
    """The result of `score`.

    table: pandas DataFrame with CELL_COLUMNS (values; `write_cells` formats them).
    prob_singlet, prob_doublet, doublet_llr: the posteriors behind the calls.
    ll_singlet, ll_doublet: the log-likelihood matrices.
    n_vars, umis: sites with reads and UMIs per cell.
    source: "vireo" when the eight vireo columns are vireo's own numbers,
        "scorer" when they come from the recomputed matrices.
    check: largest differences between vireo's probabilities and the
        recomputed ones (backend vireo only).
    calls: counts of donor_id values.
    vireo_format: the vireo version whose naming the table follows.
    """

    def __init__(self, **kw):
        self.__dict__.update(kw)


def score(capture, fit_result, prob_max=0.9, prob_doublet=0.9, min_vars=10, vireo_format="0.5.7", pair_chunk=16):
    """Score every cell of `capture` with the allele ratios of `fit_result`; returns Scores.

    The calls follow vireo's rules in vireo's order: unassigned when prob_max
    is below `prob_max`, doublet when the best pair's probability is at least
    `prob_doublet`, unassigned when fewer than `min_vars` sites have reads.
    vireo_format sets the pair naming of best_doublet: "0.5.7" writes
    vireo 0.5.7's permutation name at the pair's index (wrong for most
    cells, as vireo writes it); "0.5.9" and "0.5.10" write the pair.
    """
    import pandas as pd

    from demux_tools.vireo_port import donor_calls, pair_names

    if vireo_format not in VIREO_FORMATS:
        raise UsageError("vireo format must be one of %s" % ", ".join(VIREO_FORMATS))
    AD, DP = capture.AD, capture.DP
    donors = capture.donors
    n_cell, n = capture.n_cells, capture.n_donors
    prior = fit_result.doublet_prior
    if prior is None:
        prior = default_doublet_prior(n_cell)
    ll_s = loglik_singlet(AD, DP, capture.GPb, fit_result.theta_mean, fit_result.theta_sum)
    ll_d = loglik_doublet(AD, DP, capture.GPb, fit_result.theta_mean, fit_result.theta_sum, pair_chunk=pair_chunk)
    p_s, p_d, llr = posteriors(ll_s, ll_d, prior)
    source, check = "scorer", {}
    rv = fit_result.rv
    if rv is not None and n > 1:
        check = {
            "maxdiff_prob_singlet": float(np.max(np.abs(p_s - rv["ID_prob"]))) if n_cell else 0.0,
            "maxdiff_prob_doublet": float(np.max(np.abs(p_d - rv["doublet_prob"]))) if n_cell else 0.0,
            "maxdiff_doublet_llr": float(np.max(np.abs(llr - rv["doublet_LLR"]))) if n_cell else 0.0,
        }
        if fit_result.prior_is_vireos:
            p_s, p_d, llr = rv["ID_prob"], rv["doublet_prob"], rv["doublet_LLR"]
            source = "vireo"
    n_vars = np.array(np.sum(DP > 0, axis=0)).reshape(-1)
    umis = np.array(DP.sum(axis=0)).reshape(-1).astype(np.int64)
    donor_id, pmax, pdbl, best_singlet, k = donor_calls(donors, p_s, p_d, n_vars, prob_max, prob_doublet, min_vars)
    names = np.array(pair_names(donors, _pair_naming(vireo_format)), "U100")    # vireo cuts names at 100 characters
    best_doublet = [str(names[x]) if x >= 0 else "" for x in k]

    rows = np.arange(n_cell)
    order = np.argsort(-ll_s, axis=1)
    ll_best = ll_s[rows, order[:, 0]] if n_cell else np.zeros(0)
    if n > 1:
        ll_second = ll_s[rows, order[:, 1]]
        second = np.array(donors, dtype=object)[order[:, 1]]
    else:
        ll_second = np.full(n_cell, np.nan)
        second = np.array([""] * n_cell, dtype=object)
    if ll_d.shape[1]:
        combos = [pair_key(donors[i], donors[j]) for i, j in _combinations(n)]
        best_pair = np.array(combos, dtype=object)[ll_d.argmax(1)]
        ll_best_pair = ll_d.max(1)
    else:
        best_pair = np.array([""] * n_cell, dtype=object)
        ll_best_pair = np.full(n_cell, np.nan)
    with np.errstate(divide="ignore", invalid="ignore"):
        fit_per_umi = np.where(umis > 0, ll_best / np.maximum(umis, 1), np.nan)
    table = pd.DataFrame({
        "cell": capture.barcodes,
        "donor_id": donor_id.astype(object),
        "prob_max": pmax,
        "prob_doublet": pdbl,
        "n_vars": n_vars.astype(np.int64),
        "best_singlet": best_singlet.astype(object),
        "best_doublet": best_doublet,
        "doublet_logLikRatio": llr,
        "best_pair": best_pair,
        "second_donor": second,
        "umis": umis,
        "ll_best": ll_best,
        "ll_second": ll_second,
        "margin": ll_best - ll_second,
        "ll_best_pair": ll_best_pair,
        "doublet_margin": ll_best - ll_best_pair,
        "fit_per_umi": fit_per_umi,
    }, columns=CELL_COLUMNS)
    values, counts = np.unique(donor_id, return_counts=True)
    return Scores(table=table, prob_singlet=p_s, prob_doublet=p_d, doublet_llr=llr, ll_singlet=ll_s,
                  ll_doublet=ll_d, n_vars=n_vars, umis=umis, source=source, check=check,
                  calls={str(v): int(c) for v, c in zip(values, counts)}, vireo_format=vireo_format,
                  doublet_prior=prior)


def _combinations(n):
    from demux_tools.vireo_port import donor_pairs

    return [tuple(p) for p in donor_pairs(n)]


def write_cells(scores, path, columns=None, extra=None):
    """Write the per-cell table (tab-separated, gzip for ".gz", deterministic); returns the path.

    Numbers use FLOAT_FORMATS; a missing value is written empty. With
    vireo_format 0.5.10 the prob_doublet column is named prob_max_doublet.
    extra: a DataFrame of further columns (one row per cell), appended after
    the base columns.
    """
    table = scores.table if columns is None else scores.table.loc[:, list(columns)]
    formats = {c: f for c, f in FLOAT_FORMATS.items() if c in table.columns}
    if extra is not None:
        import pandas as pd

        clash = [c for c in extra.columns if c in table.columns]
        if clash:
            raise ValueError("extra columns repeat %s" % ", ".join(clash))
        table = pd.concat([table.reset_index(drop=True), extra.reset_index(drop=True)], axis=1)
        formats.update(extra_formats(extra))
    renamed = _prob_doublet_column(scores.vireo_format)
    if renamed != "prob_doublet":
        table = table.rename(columns={"prob_doublet": renamed})
        if "prob_doublet" in formats:
            formats[renamed] = formats.pop("prob_doublet")
    return write_table(table, path, float_formats=formats)


def write_vireo_dir(capture, scores, out_dir):
    """vireo's donor_ids.tsv, summary.tsv and probability tables for these scores."""
    from demux_tools.vireo_port import write_vireo_files

    return write_vireo_files(out_dir, capture.donors, capture.barcodes, scores.n_vars, scores.prob_singlet,
                             scores.prob_doublet, scores.doublet_llr, naming=_pair_naming(scores.vireo_format),
                             prob_doublet_column=_prob_doublet_column(scores.vireo_format))


# ---------------------------------------------------------------------------
# Validation against vireo's own output


TEXT_COLUMNS = ["donor_id", "n_vars", "best_singlet", "best_doublet"]
NUMERIC_COLUMNS = ["prob_max", "prob_doublet", "doublet_logLikRatio"]


def validate_against_vireo(cells, vireo_dir, tolerance=0.0):
    """Compare a cells table with vireo's donor_ids.tsv, cell by cell.

    cells: a DataFrame of text (as read back from a cells table) or a path.
    Returns a dict: cells, same_barcodes, cells_identical (all eight columns
    identical as text), differ_<column> per column, maxabs_<column> for the
    three numeric columns, naming (of vireo's prob_doublet header, when
    present) and ok: identical, or, with a tolerance, the text columns
    identical and every numeric difference at most `tolerance`.
    """
    import pandas as pd

    from demux_tools.vireo_files import pair_naming

    if not isinstance(cells, pd.DataFrame):
        cells = read_donor_ids(cells)
    ref = read_donor_ids(os.path.join(vireo_dir, "donor_ids.tsv"))
    out = {"cells": int(len(ref)), "cells_scored": int(len(cells))}
    same = len(cells) == len(ref) and bool((cells["cell"].to_numpy() == ref["cell"].to_numpy()).all())
    out["same_barcodes"] = same
    prob_doublet_path = os.path.join(vireo_dir, "prob_doublet.tsv.gz")
    if os.path.exists(prob_doublet_path):
        try:
            out["naming"] = pair_naming(prob_doublet_path)["kind"]
        except ValueError:
            out["naming"] = "unknown"
    if not same:
        out["cells_identical"] = 0
        out["ok"] = False
        return out
    identical = np.ones(len(ref), dtype=bool)
    for c in DONOR_IDS_COLUMNS[1:]:
        differ = cells[c].to_numpy() != ref[c].to_numpy()
        out["differ_" + c] = int(differ.sum())
        identical &= ~differ
    out["cells_identical"] = int(identical.sum())
    for c in NUMERIC_COLUMNS:
        a = pd.to_numeric(cells[c], errors="coerce").to_numpy(dtype=np.float64)
        b = pd.to_numeric(ref[c], errors="coerce").to_numpy(dtype=np.float64)
        diff = np.abs(a - b)
        both_nan = np.isnan(a) & np.isnan(b)
        diff[both_nan] = 0.0
        out["maxabs_" + c] = float(np.nanmax(diff)) if len(diff) else 0.0
    if out["cells_identical"] == len(ref):
        out["ok"] = True
    elif tolerance > 0:
        text_ok = all(out["differ_" + c] == 0 for c in TEXT_COLUMNS)
        out["ok"] = bool(text_ok and all(out["maxabs_" + c] <= tolerance for c in NUMERIC_COLUMNS))
    else:
        out["ok"] = False
    out["tolerance"] = tolerance
    return out


def compare_vireo_files(ours, theirs):
    """Text comparison of donor_ids.tsv, summary.tsv and both probability tables (decompressed).

    Returns {file: "identical" | "differ" | "missing"}.
    """
    from demux_tools._io import open_text

    out = {}
    for name in ("donor_ids.tsv", "summary.tsv", "prob_singlet.tsv.gz", "prob_doublet.tsv.gz"):
        a, b = os.path.join(ours, name), os.path.join(theirs, name)
        if not (os.path.exists(a) and os.path.exists(b)):
            out[name] = "missing"
            continue
        with open_text(a) as fa, open_text(b) as fb:
            out[name] = "identical" if fa.read() == fb.read() else "differ"
    return out



# ---------------------------------------------------------------------------
# Changing the sites and the genotypes (each returns a new Capture)


def subset_sites(capture, keep):
    """A Capture with only the sites where `keep` (a boolean mask over sites) is True."""
    keep = np.asarray(keep, dtype=bool)
    if len(keep) != capture.n_sites:
        raise ValueError("site mask has %d entries for %d sites" % (len(keep), capture.n_sites))
    idx = np.flatnonzero(keep)
    counts = dict(capture.counts)
    counts["n_sites_masked"] = int(capture.n_sites - len(idx)) + int(capture.counts.get("n_sites_masked", 0))
    return Capture(capture.name, capture.barcodes, [capture.variants[k] for k in idx], capture.AD[idx],
                   capture.DP[idx], capture.donors, capture.GPb[idx], capture.missing[idx], capture.pileup_rows[idx],
                   capture.vcf_rows[idx], {f: v[idx] for f, v in capture.info.items()}, counts, capture.inputs,
                   capture.load_info)


def _with_genotypes(capture, GPb, missing, donors=None, counts=None):
    return Capture(capture.name, capture.barcodes, capture.variants, capture.AD, capture.DP,
                   capture.donors if donors is None else donors, GPb, missing, capture.pileup_rows, capture.vcf_rows,
                   capture.info, counts or capture.counts, capture.inputs, capture.load_info)


def hwe_probabilities(q):
    """Genotype probabilities ((1-q)^2, 2q(1-q), q^2) for ALT allele frequencies q (NaN rows stay NaN)."""
    q = np.asarray(q, dtype=np.float64)
    return np.stack([(1 - q) ** 2, 2 * q * (1 - q), q ** 2], axis=-1)


def fill_missing(capture, how="uniform", af=None, donors=None):
    """Replace blank genotypes; returns (Capture, counts).

    how: "uniform" (vireo's a third each: no change), "hwe" (Hardy-Weinberg
    proportions of `af`, one ALT allele frequency per site; a site with no
    frequency keeps a third each) or "drop-site" (drop every site where a
    donor is blank). donors: names of the donors to fill (default all).
    counts: blanks filled, blanks left (no frequency), sites dropped.
    """
    cols = list(range(capture.n_donors)) if donors is None else [capture.donors.index(d) for d in donors]
    target = np.zeros_like(capture.missing)
    target[:, cols] = capture.missing[:, cols]
    if how == "uniform":
        return capture, {"blanks": int(target.sum()), "filled": 0}
    if how == "drop-site":
        keep = ~target.any(axis=1)
        return subset_sites(capture, keep), {"blanks": int(target.sum()), "sites_dropped": int((~keep).sum())}
    if how != "hwe":
        raise UsageError("blank fill must be uniform, hwe or drop-site, not %r" % how)
    if af is None:
        raise UsageError("hwe needs an allele frequency per site")
    q = np.asarray(af, dtype=np.float64)
    have = ~np.isnan(q)
    fill = target & have[:, None]
    GPb = capture.GPb.copy()
    probs = hwe_probabilities(np.where(have, q, 0.0))
    v, d = np.nonzero(fill)
    GPb[v, d, :] = probs[v]
    return (_with_genotypes(capture, GPb, capture.missing),
            {"blanks": int(target.sum()), "filled": int(fill.sum()), "left_uniform": int((target & ~fill).sum())})


def harden(capture):
    """Soft genotypes made hard (one-hot at their most likely genotype); blanks stay blank."""
    G = capture.GPb
    hard = np.zeros_like(G)
    np.put_along_axis(hard, G.argmax(axis=2)[:, :, None], 1.0, axis=2)
    GPb = np.where(capture.missing[:, :, None], G, hard)
    changed = int(((GPb != G).any(axis=2)).sum())
    return _with_genotypes(capture, GPb, capture.missing), {"hardened": changed}


def add_donors(capture, vcf, donors=None, tag="GT", missing="uniform", af=None):
    """Append donors from another VCF (decoys), matched to the capture's sites by id; returns (Capture, counts).

    A site the other VCF lacks is blank for the new donors, then `missing`
    applies to those blanks: "uniform", "hwe" (with `af` per site) or
    "drop-site". Names already in the capture raise UsageError.
    """
    from demux_tools import vcf as vcf_mod

    _, header = vcf_mod.read_header(vcf)
    names = list(header) if donors is None else list(donors)
    clash = [d for d in names if d in capture.donors]
    if clash:
        raise UsageError("extra donor(s) already in the capture: %s" % ",".join(clash))
    m = vcf_mod.read_matched(vcf, capture.variants, samples=names, tag=tag, check=False)
    geno, idx_cell, idx_donor, how = m["genotypes"], m["idx_cell"], m["idx_donor"], m["how"]
    n = capture.n_sites
    gp = np.tile(np.array([1 / 3, 1 / 3, 1 / 3]), (n, len(names), 1))
    miss = np.ones((n, len(names)), dtype=bool)
    gp[idx_cell] = geno.gp[idx_donor]
    miss[idx_cell] = geno.missing[idx_donor]
    counts = dict(capture.counts)
    counts["extra_donors"] = names
    out = _with_genotypes(capture, np.concatenate([capture.GPb, gp], axis=1),
                          np.concatenate([capture.missing, miss], axis=1), donors=capture.donors + names,
                          counts=counts)
    out, filled = fill_missing(out, how=missing, af=af, donors=names)
    return out, {"extra_donors": names, "sites_in_extra_vcf": int(len(idx_cell)), "match": how, "fill": filled}


# ---------------------------------------------------------------------------
# Allele ratios held from elsewhere


def parse_theta(spec):
    """("refit", None), ("base", None), ("hold", path) or ("fixed", (theta_mean, theta_sum)) from a --theta value."""
    spec = (spec or "refit").strip()
    if spec in ("refit", "base"):
        return spec, None
    if spec.startswith("hold:"):
        return "hold", spec[len("hold:"):]
    if spec.startswith("fixed:"):
        parts = spec[len("fixed:"):].split(":")
        try:
            mean = [float(x) for x in parts[0].split(",")]
            total = [float(x) for x in parts[1].split(",")] if len(parts) > 1 else [1e6] * 3
        except (ValueError, IndexError):
            mean = total = []
        if len(mean) != 3 or len(total) != 3:
            raise UsageError("--theta fixed:m0,m1,m2[:s0,s1,s2] needs three means (and three concentrations)")
        if not all(0.0 < m < 1.0 for m in mean) or not all(0.0 < t < float("inf") for t in total):
            raise UsageError("--theta fixed: means must lie strictly between 0 and 1, and concentrations be "
                             "positive and finite (got %s)" % spec)
        return "fixed", (mean, total)
    raise UsageError("--theta must be refit, base, hold:FIT.json or fixed:m0,m1,m2:s0,s1,s2, not %r" % spec)


# ---------------------------------------------------------------------------
# Variants: one capture scored with changed sites, genotypes or allele ratios

VARIANT_COLUMNS = ("variant", "sites", "exclude_sites", "thin_bp", "thin_keep", "group_by", "fill_blanks", "harden",
                   "extra_vcf", "extra_donors", "extra_missing", "theta", "backend", "marginal", "pairwise_shared",
                   "top_sites", "extra_columns")
_FLAGS = ("harden", "marginal", "pairwise_shared", "extra_columns")


def default_variant(name="base"):
    """The variant that changes nothing: D2's scoring."""
    return {"variant": name, "sites": "", "exclude_sites": "", "thin_bp": "", "thin_keep": "first", "group_by": [],
            "fill_blanks": "keep", "harden": False, "extra_vcf": "", "extra_donors": "", "extra_missing": "uniform",
            "theta": "refit", "backend": "", "marginal": False, "pairwise_shared": False, "top_sites": "",
            "extra_columns": False}


def _truthy(value):
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in ("", "0", "false", "no", "off"):
        return False
    if text in ("1", "true", "yes", "on"):
        return True
    raise UsageError("expected true or false, not %r" % value)


def normalise_variant(row):
    """A variant dict with every key, text values cleaned and flags as bool; group_by as a list."""
    v = default_variant(str(row.get("variant") or "base"))
    for k in VARIANT_COLUMNS[1:]:
        if k not in row or row[k] is None:
            continue
        value = row[k]
        if k in _FLAGS:
            v[k] = _truthy(value)
        elif k == "group_by":
            v[k] = [g for g in (value if isinstance(value, list) else str(value).split(";")) if g]
        else:
            v[k] = str(value).strip()
    if v["fill_blanks"] in ("", "uniform"):
        v["fill_blanks"] = "keep"
    if not v["theta"]:
        v["theta"] = "refit"
    if not v["extra_missing"]:
        v["extra_missing"] = "uniform"
    if not v["thin_keep"]:
        v["thin_keep"] = "first"
    return v


def is_default_variant(v):
    d = default_variant(v["variant"])
    return all(v[k] == d[k] for k in VARIANT_COLUMNS[1:])


def read_variants(path):
    """The variants table: one row per variant, columns from VARIANT_COLUMNS (`variant` required).

    group_by holds groupings separated by ";". Flags take true or false.
    """
    table = read_table(path)
    if "variant" not in table.columns:
        raise UsageError("%s needs a variant column" % path)
    unknown = [c for c in table.columns if c not in VARIANT_COLUMNS]
    if unknown:
        raise UsageError("%s: unknown column(s) %s (known: %s)"
                         % (path, ", ".join(unknown), ", ".join(VARIANT_COLUMNS)))
    rows = [normalise_variant(r) for r in table.to_dict(orient="records")]
    names = [check_name(r["variant"], "variant name") for r in rows]
    if len(set(names)) != len(names):
        raise UsageError("%s names a variant twice" % path)
    for k, r in enumerate(rows):
        if r["theta"] == "base" and (k == 0 or parse_theta(rows[0]["theta"])[0] != "refit"):
            raise UsageError("theta=base holds the first variant's fit; the first variant must refit")
    return rows


def _af_spec(spec):
    """For "hwe:X": ("file", X) when X is an existing file, else ("info", X)."""
    if not spec.startswith("hwe:"):
        return None
    source = spec[len("hwe:"):]
    return ("file", source) if os.path.isfile(source) else ("info", source)


def info_fields_needed(variants):
    """INFO fields the variants read: @info groupings and hwe:FIELD fills."""
    fields = set()
    for v in variants:
        for g in v["group_by"]:
            if g.startswith("@info:"):
                fields.add(g.split(":")[1])
        for spec in (v["fill_blanks"], v["extra_missing"]):
            af = _af_spec(spec)
            if af and af[0] == "info":
                fields.add(af[1])
    return sorted(fields)


def _allele_frequencies(capture, spec):
    """ALT allele frequency per site from "hwe:FIELD" (INFO) or "hwe:FILE" (a table with variant and af)."""
    kind, source = _af_spec(spec)
    if kind == "info":
        if source not in capture.info:
            raise UsageError("INFO field %s was not loaded" % source)
        return capture.info[source]
    table = read_table(source)
    if "variant" not in table.columns or "af" not in table.columns:
        raise UsageError("%s needs columns variant and af" % source)
    lookup = dict(zip(table["variant"], table["af"]))
    return np.array([float(lookup[v]) if lookup.get(v, "") not in ("", "NA", ".") else np.nan
                     for v in capture.variants], dtype=np.float64)


def _fill(spec):
    if spec in ("keep", "uniform"):
        return "uniform"
    if spec == "drop-site":
        return "drop-site"
    if spec.startswith("hwe:"):
        return "hwe"
    raise UsageError("blank fill must be keep, hwe:FIELD, hwe:FILE or drop-site, not %r" % spec)


def held_theta(path):
    """(theta_mean, theta_sum) stored in a fit.json."""
    import json

    with open(path) as f:
        rec = json.load(f)
    return rec["theta_mean"], rec["theta_sum"]


class VariantResult(object):
    """What `apply_variant` returns: capture, fit, scores, extra (DataFrame of extra columns or None),
    sites (long top-sites table or None), record (dict for fit.json) and seconds (by step)."""

    def __init__(self, **kw):
        self.__dict__.update(kw)


def apply_variant(capture, variant, backend="vireo", doublet_prior="auto", prob_max=0.9, prob_doublet=0.9,
                  min_vars=10, vireo_format="0.5.7", base_fit=None, site_table=None, top_groups=10,
                  gene_column="gene_name", geno_tag="GT"):
    """Score `capture` under one variant (a dict from `normalise_variant`); returns a VariantResult.

    The steps run in this order: keep or drop sites (sites, exclude_sites,
    thin_bp), add donors (extra_vcf), fill blanks, harden, then fit or hold
    the allele ratios (theta: refit, base = `base_fit`, hold:FIT.json or
    fixed:...), score, and the diagnostics asked for. The variant's own
    backend, when set, overrides `backend`. geno_tag is the FORMAT field read
    from extra_vcf.
    """
    from demux_tools import diagnostics as dg
    from demux_tools import sites as sites_mod

    v = normalise_variant(variant)
    seconds = {}
    record = {}
    t = time.time()
    keep = None
    if v["sites"] or v["exclude_sites"]:
        keep = sites_mod.mask(capture.variants, keep=v["sites"] or None, exclude=v["exclude_sites"] or None)
    if v["thin_bp"]:
        base_keep = np.ones(capture.n_sites, dtype=bool) if keep is None else keep
        idx = np.flatnonzero(base_keep)
        depth = np.asarray(capture.DP.sum(axis=1)).ravel()[idx]
        thinned = sites_mod.thin([capture.variants[k] for k in idx], int(v["thin_bp"]), keep=v["thin_keep"],
                                 depth=depth)
        keep = np.zeros(capture.n_sites, dtype=bool)
        keep[idx[thinned]] = True
    work = capture
    if keep is not None:
        record["sites_kept"] = int(keep.sum())
        work = subset_sites(capture, keep)
    if v["extra_vcf"]:
        missing = _fill(v["extra_missing"])
        af = _allele_frequencies(work, v["extra_missing"]) if missing == "hwe" else None
        donors = [d for d in v["extra_donors"].split(",") if d] or None
        work, record["extra_donors"] = add_donors(work, v["extra_vcf"], donors=donors,
                                                  tag=geno_tag, missing=missing, af=af)
    if v["fill_blanks"] != "keep":
        how = _fill(v["fill_blanks"])
        af = _allele_frequencies(work, v["fill_blanks"]) if how == "hwe" else None
        work, record["fill_blanks"] = fill_missing(work, how=how, af=af)
    if v["harden"]:
        work, record["harden"] = harden(work)
    if work.n_sites == 0:
        raise ValueError("variant %s leaves no site" % v["variant"])
    seconds["overrides"] = round(time.time() - t, 3)

    t = time.time()
    mode, value = parse_theta(v["theta"])
    if mode == "refit":
        fit_result = fit(work, backend=v["backend"] or backend, doublet_prior=doublet_prior)
    elif mode == "base":
        if base_fit is None:
            raise UsageError("theta=base needs the first variant's fit")
        fit_result = fit(work, theta=(base_fit.theta_mean, base_fit.theta_sum), doublet_prior=doublet_prior)
    elif mode == "hold":
        fit_result = fit(work, theta=held_theta(value), doublet_prior=doublet_prior)
    else:
        fit_result = fit(work, theta=value, doublet_prior=doublet_prior)
    seconds["fit"] = round(time.time() - t, 3)
    t = time.time()
    scores = score(work, fit_result, prob_max=prob_max, prob_doublet=prob_doublet, min_vars=min_vars,
                   vireo_format=vireo_format)
    seconds["score"] = round(time.time() - t, 3)

    t = time.time()
    extra, long_sites = [], None
    wants = v["extra_columns"] or v["group_by"] or v["marginal"] or v["top_sites"]
    if wants or v["pairwise_shared"]:
        best, second = dg.best_and_second(scores)
        entries = dg.Entries(work) if wants else None
        labels = {}
        if v["extra_columns"]:
            extra.append(dg.cheap_columns(work, fit_result, best, second, entries=entries))
        for spec in v["group_by"]:
            name, lab = sites_mod.groups_from(spec, work.variants, missing=work.missing, donors=work.donors,
                                              info=work.info, site_table=site_table, top=top_groups)
            labels[name] = lab
            extra.append(dg.group_margins(work, fit_result, best, second, lab, name, entries=entries))
        if v["marginal"]:
            extra.append(dg.marginal_singlet(work, fit_result, scores.ll_singlet, best, entries=entries)[0])
        if v["pairwise_shared"]:
            extra.append(dg.pairwise_shared(work, fit_result))
        if v["top_sites"]:
            gene = None
            if site_table is not None and gene_column in site_table.columns:
                lookup = dict(zip(site_table["variant"], site_table[gene_column]))
                gene = np.array([lookup.get(x, "") for x in work.variants], dtype=object)
            long_sites, per_cell = dg.top_sites(work, fit_result, best, second, n=int(v["top_sites"]),
                                                labels=labels, gene=gene, entries=entries)
            extra.append(per_cell)
    seconds["diagnostics"] = round(time.time() - t, 3)
    if extra:
        import pandas as pd

        extra = pd.concat([e.reset_index(drop=True) for e in extra], axis=1)
    else:
        extra = None
    if not is_default_variant(v):
        record["options"] = {k: v[k] for k in VARIANT_COLUMNS}
    return VariantResult(capture=work, fit=fit_result, scores=scores, extra=extra, sites=long_sites, record=record,
                         seconds=seconds)


def extra_formats(extra):
    """%-formats for the extra columns: integers as they are, fit_expected_per_umi %.5f, other numbers %.4f."""
    formats = {}
    if extra is None:
        return formats
    for c in extra.columns:
        if extra[c].dtype.kind == "f":
            formats[c] = "%.5f" if c == "fit_expected_per_umi" else "%.4f"
    return formats


# ---------------------------------------------------------------------------
# Running one capture (shared by score and score-many)


def _fit_record(capture, fit_result, scores, timings, vireo_format):
    c = capture.counts
    rec = {
        "capture": capture.name,
        "backend": fit_result.backend,
        "vireo_format": vireo_format,
        "values_from": scores.source,
        "donors": capture.donors,
        "n_cells": capture.n_cells,
        "n_sites_pileup": c.get("n_sites_pileup"),
        "n_sites_matched": c.get("n_sites_matched"),
        "n_sites_used": capture.n_sites,
        "sites_dropped": c.get("sites_dropped"),
        "load_counts": c,
        "theta_mean": fit_result.theta_mean.ravel().tolist(),
        "theta_sum": fit_result.theta_sum.ravel().tolist(),
        "doublet_prior": scores.doublet_prior,
        "n_iter": fit_result.n_iter,
        "elbo": fit_result.elbo,
        "check": scores.check,
        "calls": scores.calls,
        "load": capture.load_info,
        "seconds": timings,
        "peak_rss_mb": peak_rss_mb(),
    }
    if fit_result.log:
        rec["vireo_log"] = fit_result.log
    return rec


def _capture_kwargs(args_or_row):
    get = args_or_row.get if isinstance(args_or_row, dict) else (lambda k, d=None: getattr(args_or_row, k, d))
    donors = get("donors")
    if isinstance(donors, str):
        donors = [d for d in donors.split(",") if d] or None
    donors_file = get("donors_file")
    if donors is None and donors_file:
        donors = [line.strip() for line in read_lines(donors_file) if line.strip()]
    return dict(cellsnp=get("cellsnp") or None, base=get("base") or None, samples=get("samples") or None,
                ad=get("ad") or None, dp=get("dp") or None, vcf=get("vcf") or None, donors=donors,
                donor_order=get("donor_order") or "vcf", geno_tag=get("geno_tag") or "GT",
                barcodes=get("barcodes") or None)


def _load_for(name, capture_kwargs, loader, cache_dir, n_donor, info):
    capture = load_capture(name, loader=loader, cache_dir=cache_dir, info=info, **capture_kwargs)
    if n_donor is not None and n_donor != capture.n_donors:
        raise UsageError("--n-donor %d differs from the %d donors used; vireo would not run in fixed mode "
                         "(learning mode, or its donor-subset path), which the Scorer does not reproduce"
                         % (n_donor, capture.n_donors))
    return capture


def _write_result(result, out_prefix, vireo_format, timings, vireo_dir=None, argv=None, clock=None):
    """Write PREFIX.cells.tsv.gz (and PREFIX.sites.tsv.gz, vireo's files) and PREFIX.fit.json; returns the record."""
    t = time.time()
    cells_path = out_prefix + ".cells.tsv.gz"
    write_cells(result.scores, cells_path, extra=result.extra)
    outputs = {"cells": cells_path, "vireo_dir": vireo_dir}
    if result.sites is not None:
        outputs["sites"] = out_prefix + ".sites.tsv.gz"
        formats = {c: "%.4f" for c in ("term_best", "term_second", "delta", "cum_share")}
        write_table(result.sites, outputs["sites"], float_formats=formats)
    if vireo_dir:
        write_vireo_dir(result.capture, result.scores, vireo_dir)
    timings = dict(timings)
    timings["write"] = round(time.time() - t, 3)
    record = _fit_record(result.capture, result.fit, result.scores, timings, vireo_format)
    if result.record:
        record["variant"] = result.record
    record["outputs"] = outputs
    record["provenance"] = provenance(argv or [], result.capture.inputs, clock=clock)
    write_json(record, out_prefix + ".fit.json")
    return record


def run_one(name, capture_kwargs, out_prefix, backend="vireo", loader="native", cache_dir=None,
            doublet_prior="auto", prob_max=0.9, prob_doublet=0.9, min_vars=10, vireo_format="0.5.7",
            vireo_dir=None, argv=None, clock=None, n_donor=None, variant=None, site_table=None, info=(),
            top_groups=10, gene_column="gene_name"):
    """Load, fit and score one capture; write PREFIX.cells.tsv.gz and PREFIX.fit.json.

    variant: a dict of variant options (default: none, D2's scoring).
    vireo_dir: also write vireo's files there. n_donor: refuse unless the
    capture has this many donors. Returns the fit.json record.
    """
    clock = clock or start_clock()
    variant = normalise_variant(variant or {})
    fields = sorted(set(info or ()) | set(info_fields_needed([variant])))
    t = time.time()
    capture = _load_for(name, capture_kwargs, loader, cache_dir, n_donor, fields)
    timings = {"load": round(time.time() - t, 3)}
    result = apply_variant(capture, variant, backend=backend, doublet_prior=doublet_prior, prob_max=prob_max,
                           prob_doublet=prob_doublet, min_vars=min_vars, vireo_format=vireo_format,
                           site_table=site_table, top_groups=top_groups, gene_column=gene_column,
                           geno_tag=capture_kwargs.get("geno_tag") or "GT")
    timings.update(_variant_seconds(result))
    return _write_result(result, out_prefix, vireo_format, timings, vireo_dir=vireo_dir, argv=argv, clock=clock)


def _variant_seconds(result):
    """fit and score as D2 reports them; the variant steps only when they ran."""
    out = {"fit": result.seconds["fit"], "score": result.seconds["score"]}
    if result.record or result.extra is not None:
        out["overrides"] = result.seconds["overrides"]
        out["diagnostics"] = result.seconds["diagnostics"]
    return out


def run_variants(name, capture_kwargs, out_dir, variants, backend="vireo", loader="native", cache_dir=None,
                 doublet_prior="auto", prob_max=0.9, prob_doublet=0.9, min_vars=10, vireo_format="0.5.7",
                 argv=None, site_table=None, info=(), top_groups=10, gene_column="gene_name"):
    """Load one capture once and score every variant; write OUT_DIR/<variant>.cells.tsv.gz and .fit.json each.

    Returns the list of fit.json records, in the order of `variants`.
    """
    variants = [normalise_variant(v) for v in variants]
    fields = sorted(set(info or ()) | set(info_fields_needed(variants)))
    t = time.time()
    capture = _load_for(name, capture_kwargs, loader, cache_dir, None, fields)
    load_seconds = round(time.time() - t, 3)
    records, base_fit = [], None
    for k, v in enumerate(variants):
        clock = start_clock()
        result = apply_variant(capture, v, backend=backend, doublet_prior=doublet_prior, prob_max=prob_max,
                               prob_doublet=prob_doublet, min_vars=min_vars, vireo_format=vireo_format,
                               base_fit=base_fit, site_table=site_table, top_groups=top_groups,
                               gene_column=gene_column, geno_tag=capture_kwargs.get("geno_tag") or "GT")
        if k == 0:
            base_fit = result.fit
        timings = {"load": load_seconds if k == 0 else 0.0}
        timings.update(_variant_seconds(result))
        timings["overrides"] = result.seconds["overrides"]
        timings["diagnostics"] = result.seconds["diagnostics"]
        records.append(_write_result(result, os.path.join(out_dir, v["variant"]), vireo_format, timings,
                                     argv=argv, clock=clock))
    return records


# ---------------------------------------------------------------------------
# Command-line options shared by the commands


def _add_capture_args(parser, manifest=False):
    g = parser.add_argument_group("capture (a cellsnp-lite folder, or its files)")
    if not manifest:
        g.add_argument("--cellsnp", help="cellsnp-lite folder (cellSNP.base.vcf.gz, cellSNP.samples.tsv, "
                                         "cellSNP.tag.AD.mtx, cellSNP.tag.DP.mtx)")
        g.add_argument("--base", help="cellSNP.base.vcf.gz (overrides the folder's)")
        g.add_argument("--samples", help="cellSNP.samples.tsv (overrides the folder's)")
        g.add_argument("--ad", help="cellSNP.tag.AD.mtx (overrides the folder's)")
        g.add_argument("--dp", help="cellSNP.tag.DP.mtx (overrides the folder's)")
        g.add_argument("--vcf", required=True,
                       help="donor genotype VCF (plain or gzip; only records at piled sites are parsed, so a whole "
                            "panel is fine)")
        g.add_argument("--donors", help="comma-separated donors to use (default: every sample of the VCF)")
        g.add_argument("--donors-file", help="file of donors, one per line")
        g.add_argument("--donor-order", default="vcf",
                       help="vcf (the VCF header's order, as vireo uses it; default), list (as given), sorted, "
                            "or prob-header:PATH (the order of a vireo prob_singlet.tsv.gz header)")
        g.add_argument("--barcodes", help="file of barcodes to keep, one per line (default: all)")
    g.add_argument("--geno-tag", default="GT", choices=("GT", "GP", "PL"), help="genotype FORMAT field (default GT)")


def _add_model_args(parser):
    g = parser.add_argument_group("model and calls")
    g.add_argument("--n-donor", type=int, default=None,
                   help="refuse to run unless this equals the number of donors used (vireo runs another "
                        "mode when --nDonor differs from the donors in the file)")
    g.add_argument("--backend", default="vireo", choices=("vireo", "native"),
                   help="how the allele ratios are fitted: vireo (vireoSNP's own code; default) or native (the "
                        "same steps without vireoSNP; calls then come from the Scorer's matrices)")
    g.add_argument("--loader", default="native", choices=("native", "vireo"),
                   help="native: pandas readers and the cache (default); vireo: vireoSNP's readers (slow, "
                        "for checking)")
    g.add_argument("--cache-dir", help="keep loaded inputs here and reuse them on the next run")
    g.add_argument("--doublet-prior", default="auto",
                   help="prior probability of a doublet: auto (vireo's: cells / 100,000, at most 0.5) or a number")
    g.add_argument("--prob-max", type=float, default=0.9, help="unassigned below this prob_max (vireo: 0.9)")
    g.add_argument("--prob-doublet", type=float, default=0.9,
                   help="doublet at or above this best-pair probability (vireo: 0.9)")
    g.add_argument("--min-vars", type=int, default=10, help="unassigned with fewer sites with reads (vireo: 10)")
    g.add_argument("--vireo-format", default="0.5.7", choices=VIREO_FORMATS,
                   help="whose output to mimic: 0.5.7 (permutation pair names, as the vireo image; default), "
                        "0.5.9 (pair names fixed) or 0.5.10 (prob_doublet renamed prob_max_doublet)")


def _model_kwargs(args):
    return dict(backend=args.backend, loader=args.loader, cache_dir=args.cache_dir,
                doublet_prior=args.doublet_prior, prob_max=args.prob_max, prob_doublet=args.prob_doublet,
                min_vars=args.min_vars, vireo_format=args.vireo_format)


def _add_shared_diagnostic_args(parser):
    g = parser.add_argument_group("site tables and INFO fields")
    g.add_argument("--site-groups", help="table of sites with a variant column and labels (e.g. from sites "
                                         "annotate), for --group-by COLUMN and --gene-column")
    g.add_argument("--info", default="", help="comma-separated INFO fields of the donor VCF to load per site")
    g.add_argument("--top-groups", type=int, default=10,
                   help="@missing keeps this many blank patterns and calls the rest 'other' (default 10)")
    g.add_argument("--gene-column", default="gene_name",
                   help="column of --site-groups that names the gene, for the top-sites gene shares")


def _add_variant_args(parser):
    g = parser.add_argument_group("a changed site list, genotypes or allele ratios (default: none)")
    g.add_argument("--sites", help="score only these sites (a site list or VCF)")
    g.add_argument("--exclude-sites", help="leave these sites out")
    g.add_argument("--thin-bp", type=int, help="keep no two sites within this many bp")
    g.add_argument("--thin-keep", default="first", choices=("first", "deepest"),
                   help="which site of a close group to keep (default first)")
    g.add_argument("--fill-blanks", default="keep",
                   help="blank genotypes: keep (a third each, as vireo), hwe:FIELD (Hardy-Weinberg from an INFO "
                        "allele frequency), hwe:FILE (a table with variant and af) or drop-site")
    g.add_argument("--harden", action="store_true", help="make soft genotypes hard (blanks stay blank)")
    g.add_argument("--extra-vcf", help="add donors from this VCF (decoys), matched by site")
    g.add_argument("--extra-donors", default="", help="comma-separated donors of --extra-vcf (default all)")
    g.add_argument("--extra-missing", default="uniform",
                   help="the extra donors' blanks and absent sites: uniform, hwe:FIELD, hwe:FILE or drop-site")
    g.add_argument("--theta", default="refit",
                   help="allele ratios: refit (default), hold:FIT.json (from an earlier score), or "
                        "fixed:m0,m1,m2:s0,s1,s2 (means and concentrations)")
    d = parser.add_argument_group("diagnostics (extra columns; default: none)")
    d.add_argument("--extra-columns", action="store_true",
                   help="alt_umis, n_sites_diff, umis_diff, fit_expected_per_umi")
    d.add_argument("--group-by", action="append", default=[],
                   help="margins by site group: a --site-groups column, @missing, @chrom or "
                        "@info:FIELD:e0,e1,...; repeatable")
    d.add_argument("--marginal", action="store_true",
                   help="exact marginal at soft and blank genotypes beside vireo's mean-field value")
    d.add_argument("--pairwise-shared", action="store_true", help="each donor pair compared on shared sites")
    d.add_argument("--top-sites", type=int, help="write PREFIX.sites.tsv.gz with each cell's N top sites")


def variant_from_args(args):
    """The variant a score command line asks for."""
    return normalise_variant({
        "variant": "base", "sites": args.sites, "exclude_sites": args.exclude_sites,
        "thin_bp": args.thin_bp or "", "thin_keep": args.thin_keep, "group_by": list(args.group_by),
        "fill_blanks": args.fill_blanks, "harden": args.harden, "extra_vcf": args.extra_vcf,
        "extra_donors": args.extra_donors, "extra_missing": args.extra_missing, "theta": args.theta,
        "backend": "", "marginal": args.marginal, "pairwise_shared": args.pairwise_shared,
        "top_sites": args.top_sites or "", "extra_columns": args.extra_columns})


def _diagnostic_kwargs(args):
    site_table = read_table(args.site_groups) if args.site_groups else None
    if site_table is not None and "variant" not in site_table.columns:
        raise UsageError("%s needs a variant column" % args.site_groups)
    return dict(site_table=site_table, info=[f for f in args.info.split(",") if f], top_groups=args.top_groups,
                gene_column=args.gene_column)


def _wrap_errors(fn, *a, **kw):
    try:
        return fn(*a, **kw)
    except (FileNotFoundError, ValueError) as e:
        raise UsageError(str(e))


def require_libraries(backends=("vireo",), loader="native"):
    """Raise UsageError, naming what to do, when scipy or (for vireo) vireoSNP cannot be imported."""
    import importlib

    needs = ["scipy"]
    if "vireo" in backends or loader == "vireo":
        needs.append("vireoSNP")
    for name in needs:
        try:
            importlib.import_module(name)
        except ImportError:
            if name == "scipy":
                raise UsageError("the Scorer needs scipy; run it in the vireo image (snp_demux-vireo:0.5.7) or "
                                 "an environment with numpy, pandas and scipy")
            raise UsageError("the vireo backend and loader need vireoSNP; run in the vireo image "
                             "(snp_demux-vireo:0.5.7), install vireoSNP==0.5.9, or pass --backend native")


def check_name(name, what):
    """Raise UsageError unless `name` can be a file or folder name (no "/", not empty, "." or "..")."""
    if not name or name in (".", "..") or "/" in name or "\x00" in name:
        raise UsageError("%s %r cannot be used as a file name" % (what, name))
    return name


# ---------------------------------------------------------------------------
# Command: score


def configure_score(parser):
    parser.description = (
        "Reproduce vireo's fixed-mode calls for one capture and write each cell's evidence: "
        "PREFIX.cells.tsv.gz (vireo's eight donor_ids.tsv columns, then best_pair, second_donor, umis, "
        "ll_best, ll_second, margin, ll_best_pair, doublet_margin, fit_per_umi) and PREFIX.fit.json "
        "(allele ratios, counts, the check against vireo's probabilities, timings, provenance).")
    _add_capture_args(parser)
    _add_model_args(parser)
    _add_shared_diagnostic_args(parser)
    _add_variant_args(parser)
    parser.add_argument("--vireo-out", help="also write donor_ids.tsv, summary.tsv, prob_singlet.tsv.gz and "
                                            "prob_doublet.tsv.gz to this folder, as vireo writes them (an output: "
                                            "a folder that already holds donor_ids.tsv is refused without --force)")
    parser.add_argument("--force", action="store_true", help="overwrite vireo's files in --vireo-out")
    parser.add_argument("--name", help="capture name for the records (default: the prefix's file name)")
    parser.add_argument("--out", required=True, help="output prefix")


def run_score(args):
    clock = start_clock()
    kwargs = _capture_kwargs(args)
    name = args.name or os.path.basename(args.out)
    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    for p in [args.vcf, args.cellsnp, args.base, args.samples, args.ad, args.dp]:
        if p:
            check_output_path(args.out + ".cells.tsv.gz", [p])
    if args.vireo_out:
        if os.path.exists(os.path.join(args.vireo_out, "donor_ids.tsv")) and not args.force:
            raise UsageError("%s already holds a donor_ids.tsv; --vireo-out is an output folder (pass --force to "
                             "overwrite it)" % args.vireo_out)
        if args.cellsnp and os.path.realpath(args.vireo_out) == os.path.realpath(args.cellsnp):
            raise UsageError("--vireo-out is the cellsnp input folder")
    if args.n_donor is not None:
        from demux_tools.vcf import read_header

        donor_list = _donor_list(read_header(args.vcf)[1], kwargs["donors"], kwargs["donor_order"])
        if args.n_donor != len(donor_list):
            raise UsageError("--n-donor %d differs from the %d donors used; vireo would not run in fixed mode "
                             "(learning mode, or its donor-subset path), which the Scorer does not reproduce"
                             % (args.n_donor, len(donor_list)))
    variant = variant_from_args(args)
    require_libraries((args.backend,) if parse_theta(variant["theta"])[0] == "refit" else (), args.loader)
    options = _model_kwargs(args)
    options.update(_diagnostic_kwargs(args))
    record = _wrap_errors(run_one, name, kwargs, args.out, vireo_dir=args.vireo_out, argv=args.argv, clock=clock,
                          n_donor=args.n_donor, variant=variant, **options)
    calls = record["calls"]
    print("score %s: %d cells, %d donors, %d sites used of %s piled; theta %s; %d doublet, %d unassigned; "
          "values from %s; cache %s"
          % (name, record["n_cells"], len(record["donors"]), record["n_sites_used"], record["n_sites_pileup"],
             ",".join("%.4g" % x for x in record["theta_mean"]), calls.get("doublet", 0),
             calls.get("unassigned", 0), record["values_from"], record["load"].get("cache")))
    return 0


# ---------------------------------------------------------------------------
# Command: score-many

MANIFEST_COLUMNS = ("capture", "cellsnp", "base", "samples", "ad", "dp", "vcf", "donors", "donor_order",
                    "barcodes", "vireo_dir", "geno_tag")


def read_manifest(path):
    """The manifest as a list of dicts (text). Requires capture and vcf, and cellsnp or base+samples+ad+dp."""
    table = read_table(path, comment=None)
    if "capture" not in table.columns or "vcf" not in table.columns:
        raise UsageError("%s needs columns capture and vcf" % path)
    unknown = [c for c in table.columns if c not in MANIFEST_COLUMNS]
    if unknown:
        raise UsageError("%s: unknown column(s) %s (known: %s)"
                         % (path, ", ".join(unknown), ", ".join(MANIFEST_COLUMNS)))
    rows = table.to_dict(orient="records")
    seen = set()
    for r in rows:
        check_name(r["capture"], "capture name")
        if r["capture"] in seen:
            raise UsageError("%s lists capture %s twice" % (path, r["capture"]))
        seen.add(r["capture"])
        if not r.get("cellsnp") and not all(r.get(k) for k in ("base", "samples", "ad", "dp")):
            raise UsageError("%s: capture %s needs cellsnp, or base, samples, ad and dp" % (path, r["capture"]))
    return rows


def _dp_size(row):
    try:
        paths = resolve_cellsnp(cellsnp=row.get("cellsnp") or None, base=row.get("base") or None,
                                samples=row.get("samples") or None, ad=row.get("ad") or None,
                                dp=row.get("dp") or None)
        return os.path.getsize(paths["DP"])
    except (OSError, TypeError):
        return 0


def _limit_threads():
    try:
        from threadpoolctl import threadpool_limits
        return threadpool_limits(limits=1)
    except Exception:
        return None


def _score_many_worker(job):
    """Score one manifest row in a worker process; returns (summary rows, validation row or None).

    Without variants there is one summary row (variant "base"); with
    variants, one per variant. A capture is validated against its vireo_dir
    on the variant named "base" when that variant changes nothing.
    """
    row, options, out_dir, argv, variants = job
    limiter = _limit_threads()
    cap = row["capture"]
    validation = None
    t0 = time.time()
    names = [v["variant"] for v in variants] if variants else ["base"]
    summaries = []
    try:
        cap_dir = os.path.join(out_dir, cap)
        os.makedirs(cap_dir, exist_ok=True)
        kwargs = _capture_kwargs(row)
        if variants:
            records = run_variants(cap, kwargs, cap_dir, variants, argv=argv, **options)
        else:
            records = [run_one(cap, kwargs, os.path.join(cap_dir, "base"), argv=argv, **options)]
        for name, record in zip(names, records):
            summary = {"capture": cap, "variant": name, "status": "ok"}
            summary.update(_summary_fields(record, with_variant_seconds=bool(variants)))
            summaries.append(summary)
        check = [v for v in (variants or [default_variant()]) if v["variant"] == "base" and is_default_variant(v)]
        if row.get("vireo_dir") and check:
            validation = {"capture": cap}
            validation.update(validate_against_vireo(os.path.join(cap_dir, "base.cells.tsv.gz"), row["vireo_dir"]))
    except Exception as e:  # recorded in summary.tsv; the run carries on with the other captures
        error = "error: %s: %s" % (type(e).__name__, str(e).replace("\t", " ").replace("\n", " "))
        done = set(s["variant"] for s in summaries)
        summaries += [{"capture": cap, "variant": n, "status": error} for n in names if n not in done]
    for summary in summaries:
        summary["seconds_total"] = round(time.time() - t0, 3)
        summary["peak_rss_mb"] = peak_rss_mb()
    if limiter is not None:
        limiter.restore_original_limits()
    return summaries, validation


def _summary_fields(record, with_variant_seconds=False):
    calls = record["calls"]
    out = {
        "cells": record["n_cells"],
        "donors": len(record["donors"]),
        "sites_pileup": record["n_sites_pileup"],
        "sites_used": record["n_sites_used"],
        "theta_mean": ",".join(repr(x) for x in record["theta_mean"]),
        "theta_sum": ",".join(repr(x) for x in record["theta_sum"]),
        "n_iter": ",".join(str(x) for x in (record["n_iter"] or [])),
        "singlets": sum(v for k, v in calls.items() if k not in ("doublet", "unassigned")),
        "doublets": calls.get("doublet", 0),
        "unassigned": calls.get("unassigned", 0),
        "values_from": record["values_from"],
        "cache": record["load"].get("cache"),
        "seconds_load": record["seconds"]["load"],
        "seconds_fit": record["seconds"]["fit"],
        "seconds_score": record["seconds"]["score"],
    }
    for k in ("maxdiff_prob_singlet", "maxdiff_prob_doublet", "maxdiff_doublet_llr"):
        out[k] = record["check"].get(k, "")
    if with_variant_seconds:
        out["backend"] = record["backend"]
        out["seconds_overrides"] = record["seconds"].get("overrides", "")
        out["seconds_diagnostics"] = record["seconds"].get("diagnostics", "")
    return out


def configure_score_many(parser):
    parser.description = (
        "Score many captures listed in a manifest, largest DP matrix first, in a pool of worker processes "
        "(one process per capture, BLAS threads 1). Writes OUT_DIR/<capture>/base.cells.tsv.gz and "
        "base.fit.json, OUT_DIR/summary.tsv (one row per capture), and OUT_DIR/validation.tsv for captures "
        "with a vireo_dir. Manifest columns (tab-separated, header): capture, vcf, and cellsnp or "
        "base+samples+ad+dp; optional donors (comma-separated), donor_order, barcodes, geno_tag, vireo_dir "
        "(vireo's output for that capture, to compare with). Exits 1 if a capture fails or differs from "
        "its vireo_dir. With --variants, each capture is loaded once and scored under every variant: "
        "OUT_DIR/<capture>/<variant>.cells.tsv.gz and .fit.json, one summary row per capture and variant.")
    parser.add_argument("--manifest", required=True, help="manifest table (see above)")
    parser.add_argument("--out-dir", required=True, help="output folder")
    parser.add_argument("--jobs", type=int, default=1, help="worker processes (default 1)")
    parser.add_argument("--variants", help="variants table: one row per variant; columns variant, sites, "
                                           "exclude_sites, thin_bp, thin_keep, group_by (';'-separated), "
                                           "fill_blanks, harden, extra_vcf, extra_donors, extra_missing, theta "
                                           "(refit, base = the first variant's fit, hold:FIT.json, fixed:...), "
                                           "backend, marginal, pairwise_shared, top_sites, extra_columns")
    _add_capture_args(parser, manifest=True)
    _add_model_args(parser)
    _add_shared_diagnostic_args(parser)


def run_score_many(args):
    import multiprocessing

    import pandas as pd

    clock = start_clock()
    rows = read_manifest(args.manifest)
    if args.n_donor is not None:
        raise UsageError("--n-donor applies to score only")
    variants = read_variants(args.variants) if args.variants else None
    backends = set(v["backend"] or args.backend for v in variants if parse_theta(v["theta"])[0] == "refit") \
        if variants else {args.backend}
    require_libraries(tuple(backends), args.loader)
    for r in rows:
        if not r.get("geno_tag"):
            r["geno_tag"] = args.geno_tag
    rows.sort(key=lambda r: -_dp_size(r))
    os.makedirs(args.out_dir, exist_ok=True)
    options = _model_kwargs(args)
    options.update(_diagnostic_kwargs(args))
    jobs = [(r, options, args.out_dir, args.argv, variants) for r in rows]
    summaries, validations = [], []
    ctx = multiprocessing.get_context("fork") if hasattr(os, "fork") else multiprocessing
    pool = ctx.Pool(processes=max(1, args.jobs), maxtasksperchild=1)
    try:
        for rows_done, validation in pool.imap_unordered(_score_many_worker, jobs):
            summaries.extend(rows_done)
            if validation is not None:
                validations.append(validation)
            for summary in rows_done:
                print("score-many: %s %s %s (%s s)" % (summary["capture"], summary["variant"], summary["status"],
                                                       summary["seconds_total"]))
            sys.stdout.flush()
    finally:
        pool.close()
        pool.join()
    order = {r["capture"]: i for i, r in enumerate(rows)}
    variant_order = {v["variant"]: i for i, v in enumerate(variants or [default_variant()])}
    summaries.sort(key=lambda s: (order[s["capture"]], variant_order.get(s["variant"], 0)))
    validations.sort(key=lambda s: order[s["capture"]])
    write_table(pd.DataFrame(summaries), os.path.join(args.out_dir, "summary.tsv"))
    failed = sorted(set(s["capture"] for s in summaries if s["status"] != "ok"), key=order.get)
    differ = []
    if validations:
        vt = pd.DataFrame(validations)
        write_table(vt, os.path.join(args.out_dir, "validation.tsv"))
        differ = [v["capture"] for v in validations if not v.get("ok")]
    extra = {"captures": len(rows), "failed": failed, "differ_from_vireo": differ,
             "validated": len(validations)}
    from demux_tools._io import write_provenance

    inputs = [args.manifest] + ([args.variants] if args.variants else [])
    write_provenance(os.path.join(args.out_dir, "summary.tsv"), args.argv, inputs, clock=clock, extra=extra)
    print("score-many: %d captures, %d failed, %d validated against vireo, %d differ"
          % (len(rows), len(failed), len(validations), len(differ)))
    return 1 if failed or differ else 0


# ---------------------------------------------------------------------------
# Command: validate-vireo


def configure_validate_vireo(parser):
    parser.description = (
        "Compare a Scorer cells table with vireo's donor_ids.tsv cell by cell: the barcodes in order and "
        "vireo's eight columns as text, plus the largest numeric difference in prob_max, prob_doublet and "
        "doublet_logLikRatio. Exits 0 when every cell is identical (or, with --tolerance, the text columns "
        "are and the numbers are within it), else 1.")
    parser.add_argument("--scored", required=True, help="Scorer output prefix, or a .cells.tsv.gz path")
    parser.add_argument("--vireo-dir", required=True, help="vireo's output folder (donor_ids.tsv)")
    parser.add_argument("--tolerance", type=float, default=0.0,
                        help="largest numeric difference accepted when text differs (default 0: text identical)")
    parser.add_argument("--files", help="also compare a folder written by score --vireo-out with vireo's "
                                        "donor_ids.tsv, summary.tsv and probability tables, as text")
    parser.add_argument("--out", help="write the result as JSON here")


def run_validate_vireo(args):
    clock = start_clock()
    cells = args.scored if args.scored.endswith(".tsv.gz") or args.scored.endswith(".tsv") \
        else args.scored + ".cells.tsv.gz"
    if not os.path.isfile(cells):
        raise UsageError("not found: %s" % cells)
    if not os.path.isfile(os.path.join(args.vireo_dir, "donor_ids.tsv")):
        raise UsageError("not found: %s" % os.path.join(args.vireo_dir, "donor_ids.tsv"))
    result = validate_against_vireo(cells, args.vireo_dir, tolerance=args.tolerance)
    if args.files:
        result["files"] = compare_vireo_files(args.files, args.vireo_dir)
        if any(v != "identical" for v in result["files"].values()):
            result["ok"] = False
    if args.out:
        rec = dict(result)
        rec["provenance"] = provenance(args.argv, [cells, os.path.join(args.vireo_dir, "donor_ids.tsv")], clock=clock)
        write_json(rec, args.out)
    differ = ", ".join("%s %d" % (c, result.get("differ_" + c, 0)) for c in DONOR_IDS_COLUMNS[1:]
                       if result.get("differ_" + c))
    print("validate-vireo: %s of %d cells identical in vireo's eight columns%s; largest differences: %s%s"
          % (result.get("cells_identical", 0), result["cells"],
             "" if result["same_barcodes"] else " (barcodes differ)",
             ", ".join("%s %.3g" % (c, result.get("maxabs_" + c, float("nan"))) for c in NUMERIC_COLUMNS),
             "; differing: " + differ if differ else ""))
    if args.files:
        print("validate-vireo: files %s" % ", ".join("%s %s" % kv for kv in sorted(result["files"].items())))
    return 0 if result["ok"] else 1


# ---------------------------------------------------------------------------
# Command: site-terms


def configure_site_terms(parser):
    parser.description = (
        "One cell's evidence, site by site: for each site where the cell has reads, its ALT and REF UMIs, "
        "each donor's genotype and the log-likelihood term under each donor. Each donor's terms sum to "
        "the cell's log-likelihood for that donor (ll_best for its best donor).")
    _add_capture_args(parser)
    parser.add_argument("--cell", action="append", required=True, help="barcode; repeatable")
    parser.add_argument("--fit", help="take the allele ratios from this fit.json (default: refit with vireo)")
    parser.add_argument("--cache-dir", help="as for score")
    parser.add_argument("--loader", default="native", choices=("native", "vireo"))
    parser.add_argument("--out", required=True, help="output table (.tsv or .tsv.gz)")


def fit_from_json(path, capture):
    """A Fit with the allele ratios stored in a fit.json; the doublet prior is vireo's for `capture`."""
    import json

    with open(path) as f:
        rec = json.load(f)
    out = Fit(rec["theta_mean"], rec["theta_sum"], backend="held:" + os.path.basename(path),
              doublet_prior=default_doublet_prior(capture.n_cells))
    return out


def run_site_terms(args):
    import pandas as pd

    clock = start_clock()
    name = os.path.basename(args.out)
    require_libraries(() if args.fit else ("vireo",), args.loader)
    capture = _wrap_errors(load_capture, name, loader=args.loader, cache_dir=args.cache_dir, **_capture_kwargs(args))
    fit_result = fit_from_json(args.fit, capture) if args.fit else fit(capture)
    tables = []
    totals = []
    for barcode in args.cell:
        t = cell_site_terms(capture, fit_result, barcode)
        tables.append(t)
        sums = {d: float(t["term__" + d].sum()) for d in capture.donors}
        best = max(sums, key=sums.get)
        totals.append({"cell": barcode, "sites": len(t), "umis": int(t["alt"].sum() + t["ref"].sum()),
                       "loglik": sums, "best": best})
        print("site-terms %s: %d sites, best donor %s %.4f" % (barcode, len(t), best, sums[best]))
    table = pd.concat(tables, ignore_index=True)
    formats = {c: "%.6f" for c in table.columns if c.startswith("term__")}
    write_table(table, args.out, float_formats=formats)
    from demux_tools._io import write_provenance

    write_provenance(args.out, args.argv, capture.inputs + ([args.fit] if args.fit else []), clock=clock,
                     extra={"theta_mean": fit_result.theta_mean.ravel().tolist(),
                            "theta_sum": fit_result.theta_sum.ravel().tolist(), "cells": totals})
    return 0
