# Portions of this file are adapted from vireoSNP 0.5.7 and 0.5.9
# (https://github.com/single-cell-genetics/vireo), by Yuanhua Huang,
# licensed under the Apache License, Version 2.0. The licence text is in
# LICENSE-vireoSNP.txt beside this file. You may not use this file except in
# compliance with that licence. Distributed on an "AS IS" BASIS, WITHOUT
# WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#
# Changes from the vireoSNP sources:
# - `normalize`, `loglik_amplify` and `add_doublet_theta` are copied unchanged
#   (vireo_base.py, vireo_doublet.py).
# - `match` gives the same pairing as vireo_base.match, computed with numpy
#   (argsort and searchsorted) instead of a Python loop, and returns -1 instead
#   of None for an id without a partner.
# - `match_variants` is the matching step of io_utils.match_donor_VCF with its
#   two "chr" fallbacks as vcf_utils.match_SNPs (0.5.9) writes them, returning
#   row indices instead of subsetting the inputs. 0.5.7 has the same first
#   pass, but its fallbacks stop with a TypeError unless every pileup id
#   matches once "chr" is added.
# - `parse_gt_code` is parse_GT_code of vcf_utils.parse_donor_GPb (0.5.7) with
#   min_prob = 0. `genotype_probabilities` applies it once per distinct code
#   (GT, PL) or vectorised over the distinct codes (GP), and raises ValueError
#   where vireoSNP would stop with an exception.
# - `add_doublet_gt` is vireo_doublet.add_doublet_GT restricted to a chosen set
#   of donor pairs, returning the pair columns only, so the doublet matrix can
#   be built a few pairs at a time.
# - `loglik_columns` and `posteriors` are the arithmetic of
#   vireo_doublet.predict_doublet, split so that the singlet and doublet
#   columns can be computed separately.
# - `write_vireo_files` writes the text of io_utils.write_donor_id with the pair
#   naming of 0.5.7 (permutations) or 0.5.9 (combinations) chosen by the
#   caller; gzip output has mtime 0, files are written atomically, and no
#   _log.txt is written.
"""Code ported from vireoSNP, kept in this one module (see the notice above).

Everything else in demux_tools is original. Where a function here returns
what a vireoSNP function returns, it gives the same numbers bit for bit on the
same numpy and scipy, because it performs the same floating-point operations
in the same order. The tests compare them with vireoSNP where it is installed.
scipy is imported inside the functions that need it.
"""
import itertools
import os

import numpy as np

UNIFORM_CODES = (".", "./.", ".|.")


# ---------------------------------------------------------------------------
# vireo_base.py


def normalize(X, axis=-1):
    """X divided by its sum along `axis` (vireo_base.normalize)."""
    X_sum = np.sum(X, axis=axis, keepdims=True)
    return X / X_sum


def loglik_amplify(X, axis=-1):
    """X minus its maximum along `axis` (vireo_base.loglik_amplify)."""
    X_max = np.max(X, axis=axis, keepdims=True)
    return X - X_max


def match(ref_ids, new_ids, uniq_ref_only=True):
    """Index into `new_ids` of each entry of `ref_ids`, or -1 (vireo_base.match).

    vireo sorts both lists with numpy's default argsort and walks them
    together. Among equal ids, the k-th one of `ref_ids` in that sorted order
    is paired with the k-th equal one of `new_ids` (with uniq_ref_only), so a
    repeated id in `ref_ids` beyond the copies in `new_ids` gets -1. Without
    uniq_ref_only every copy is paired with the first. Returns int64.
    """
    ref = np.asarray(ref_ids)
    new = np.asarray(new_ids)
    out = np.full(len(ref), -1, dtype=np.int64)
    if len(ref) == 0 or len(new) == 0:
        return out
    idx1 = np.argsort(ref)
    idx2 = np.argsort(new)
    sorted_ref = ref[idx1]
    sorted_new = new[idx2]
    left = np.searchsorted(sorted_new, sorted_ref, side="left")
    right = np.searchsorted(sorted_new, sorted_ref, side="right")
    if uniq_ref_only:
        start = np.ones(len(sorted_ref), dtype=bool)
        start[1:] = sorted_ref[1:] != sorted_ref[:-1]
        first = np.flatnonzero(start)
        rank = np.arange(len(sorted_ref)) - first[np.cumsum(start) - 1]
    else:
        rank = 0
    pos = left + rank
    ok = pos < right
    out[idx1[ok]] = idx2[pos[ok]]
    return out


# ---------------------------------------------------------------------------
# io_utils.match_donor_VCF (0.5.7)


def match_variants(cell_ids, donor_ids):
    """Pair pileup rows with donor VCF records as vireo does (match_donor_VCF, match_SNPs).

    Ids are "CHROM_POS_REF_ALT". When no id matches, vireo 0.5.9 retries with
    "chr" added to the pileup ids, then with "chr" added to the donor ids
    (0.5.7 stops with an error in these retries unless every pileup id
    matches).
    Returns (idx_cell, idx_donor, how): the matched pileup rows in increasing
    order, the donor record of each, and "as is", "chr added to pileup",
    "chr added to donors" or "none".
    """
    cell_ids = list(cell_ids)
    donor_ids = list(donor_ids)
    how = "as is"
    mm = match(cell_ids, donor_ids)
    if not (mm >= 0).any():
        how = "chr added to pileup"
        mm = match(["chr" + x for x in cell_ids], donor_ids)
    if not (mm >= 0).any():
        how = "chr added to donors"
        mm = match(cell_ids, ["chr" + x for x in donor_ids])
    if not (mm >= 0).any():
        how = "none"
    idx_cell = np.flatnonzero(mm >= 0).astype(np.int64)
    return idx_cell, mm[idx_cell], how


# ---------------------------------------------------------------------------
# vcf_utils.parse_donor_GPb (0.5.7)


def parse_gt_code(code, tag):
    """One genotype field as vireo reads it: probabilities of 0, 1 and 2 ALT alleles.

    ".", "./." and ".|." are a third each. GT adds the first and last
    characters (so "0|1" is 1, haploid "1" is 2); GP is three numbers; PL is
    vireo's 10^(-0.1 (PL - min PL) - 0.025). The result is divided by its sum.
    """
    if code in UNIFORM_CODES:
        return np.array([1 / 3, 1 / 3, 1 / 3])
    if tag == "GT":
        _prob = np.array([0, 0, 0], float)
        _prob[int(float(code[0]) + float(code[-1]))] = 1
    elif tag == "GP":
        _prob = np.array(code.split(","), float)
    elif tag == "PL":
        _Phred = np.array(code.split(","), float)
        _prob = 10 ** (-0.1 * (_Phred - min(_Phred)) - 0.025)
    else:
        raise ValueError("unsupported genotype tag %r (GT, GP or PL)" % (tag,))
    _prob += 0
    return _prob / np.sum(_prob)


def genotype_probabilities(codes, tag):
    """(probabilities, missing) for an array of genotype field strings.

    codes: array of str of any shape S. Returns float64 of shape S + (3,) with
    vireo's values for every code, and a bool array of shape S that is True for
    ".", "./." and ".|.". Raises ValueError for a code vireo cannot read (it
    would stop with an exception), naming the code.
    """
    import pandas as pd

    codes = np.asarray(codes, dtype=object)
    shape = codes.shape
    flat = codes.ravel()
    labels, uniques = pd.factorize(flat, sort=False)
    uniques = np.asarray(uniques, dtype=object)
    table = np.zeros((len(uniques), 3), dtype=np.float64)
    uniform = np.array([u in UNIFORM_CODES for u in uniques], dtype=bool)
    if uniform.any():
        table[uniform] = parse_gt_code(".", tag)
    rest = np.flatnonzero(~uniform)
    if tag == "GP" and len(rest):
        parts = [u.split(",") for u in uniques[rest]]
        bad = [i for i, p in enumerate(parts) if len(p) != 3]
        if bad:
            raise ValueError("GP value %r does not have three numbers" % (uniques[rest[bad[0]]],))
        try:
            values = np.array(parts, dtype=float)
        except ValueError:
            raise ValueError("a GP value is not three numbers, e.g. one of %r" % (list(uniques[rest[:5]]),))
        values += 0
        table[rest] = values / np.sum(values, axis=-1, keepdims=True)
    else:
        for i in rest:
            try:
                p = parse_gt_code(uniques[i], tag)
                if p.shape != (3,):
                    raise ValueError
            except (ValueError, IndexError):
                raise ValueError("vireo cannot read %s value %r" % (tag, uniques[i]))
            table[i] = p
    probs = table[labels].reshape(shape + (3,))
    missing = uniform[labels].reshape(shape)
    return probs, missing


# ---------------------------------------------------------------------------
# vireo_doublet.py


def add_doublet_theta(beta_mu, beta_sum):
    """Allele ratios of the doublet genotypes 0&1, 0&2 and 1&2 appended to the singlet ones."""
    combn_iter = itertools.combinations(range(beta_mu.shape[1]), 2)
    db_idx = np.array([x for x in combn_iter])

    beta_mu_db = (beta_mu[:, db_idx[:, 0]] + beta_mu[:, db_idx[:, 1]]) / 2.0
    beta_sum_db = np.sqrt(beta_sum[:, db_idx[:, 0]] * beta_sum[:, db_idx[:, 1]])

    return (np.append(beta_mu, beta_mu_db, axis=-1),
            np.append(beta_sum, beta_sum_db, axis=-1))


def donor_pairs(n_donor):
    """The donor pairs (i, j) in vireo's order, itertools.combinations; int array (n_pairs, 2)."""
    pairs = np.array(list(itertools.combinations(range(n_donor), 2)), dtype=np.int64)
    return pairs.reshape(-1, 2)


def add_doublet_gt(GT_prob, pairs=None):
    """Genotype probabilities of donor pairs, as add_doublet_GT builds them.

    GT_prob: (n_var, n_donor, 3). pairs: indices into `donor_pairs(n_donor)`
    (default all). Returns (n_var, len(pairs), 6): genotypes 0, 1, 2 of both
    donors, then the mixed genotypes 0&1, 0&2, 1&2, divided by their sum.
    Equal, value for value, to the pair columns of add_doublet_GT(GT_prob).
    """
    gt_idx = np.array(list(itertools.combinations(range(GT_prob.shape[2]), 2)))
    g_idx1 = gt_idx[:, 0]
    g_idx2 = gt_idx[:, 1]

    sp_idx = donor_pairs(GT_prob.shape[1])
    if pairs is not None:
        sp_idx = sp_idx[np.asarray(pairs, dtype=np.int64)]
    s_idx1 = sp_idx[:, 0]
    s_idx2 = sp_idx[:, 1]

    n_gt = GT_prob.shape[2]
    GT_prob2 = np.zeros((GT_prob.shape[0], sp_idx.shape[0], n_gt + gt_idx.shape[0]))

    GT_prob2[:, :, :n_gt] = (GT_prob[:, s_idx1, :] *
                             GT_prob[:, s_idx2, :])
    GT_prob2[:, :, n_gt:] = (GT_prob[:, s_idx1, :][:, :, g_idx1] *
                             GT_prob[:, s_idx2, :][:, :, g_idx2] +
                             GT_prob[:, s_idx1, :][:, :, g_idx2] *
                             GT_prob[:, s_idx2, :][:, :, g_idx1])
    return normalize(GT_prob2, axis=2)


def digamma_terms(theta_mean, theta_sum):
    """(dg1, dg2, dgs), each (1, 1, n_gt): digamma of the Beta shapes, as predict_doublet forms them."""
    from scipy.special import digamma

    _digamma1 = np.expand_dims(digamma(theta_sum * theta_mean), 1)
    _digamma2 = np.expand_dims(digamma(theta_sum * (1 - theta_mean)), 1)
    _digammas = np.expand_dims(digamma(theta_sum), 1)
    return _digamma1, _digamma2, _digammas


def loglik_columns(AD, BD, DP, GT, terms):
    """Sum over genotypes ig of AD.T @ (GT_ig * dg1_ig) + BD.T @ (GT_ig * dg2_ig) - DP.T @ (GT_ig * dgs_ig).

    AD, BD, DP: sparse (n_var, n_cell); GT: (n_var, n_col, n_gt); terms from
    `digamma_terms` for the same n_gt. Returns (n_cell, n_col). Each output
    column depends only on its own GT column, so computing a subset of
    columns gives the same values as computing all of them.
    """
    _digamma1, _digamma2, _digammas = terms
    logLik_ID = np.zeros((AD.shape[1], GT.shape[1]))
    for ig in range(GT.shape[2]):
        S1 = AD.T @ (GT[:, :, ig] * _digamma1[:, :, ig])
        S2 = BD.T @ (GT[:, :, ig] * _digamma2[:, :, ig])
        SS = DP.T @ (GT[:, :, ig] * _digammas[:, :, ig])
        logLik_ID += (S1 + S2 - SS)
    return logLik_ID


def posteriors(ll_singlet, ll_doublet, doublet_rate_prior):
    """(prob_singlet, prob_doublet, doublet_LLR) as predict_doublet computes them.

    The prior gives each donor (1 - rate) / n_donor and each pair
    rate / n_pairs. doublet_LLR is the best pair's log-likelihood minus the
    best donor's (NaN when there is no pair).
    """
    n_cell, n_donor = ll_singlet.shape
    n_doublet_pair = ll_doublet.shape[1]
    ID_prior = normalize(np.ones((n_cell, n_donor)))
    if n_doublet_pair:
        ID_prior_both = np.append(
            ID_prior * (1 - doublet_rate_prior),
            np.ones((n_cell, n_doublet_pair)) / n_doublet_pair * doublet_rate_prior, axis=1)
        logLik_ID = np.append(ll_singlet, ll_doublet, axis=1)
        logLik_ratio = logLik_ID[:, n_donor:].max(1) - logLik_ID[:, :n_donor].max(1)
    else:
        ID_prior_both = ID_prior
        logLik_ID = ll_singlet
        logLik_ratio = np.full(n_cell, np.nan)
    ID_prob_both = normalize(np.exp(loglik_amplify(logLik_ID + np.log(ID_prior_both))))
    return ID_prob_both[:, :n_donor], ID_prob_both[:, n_donor:], logLik_ratio


# ---------------------------------------------------------------------------
# io_utils.write_donor_id


def pair_names(donors, naming):
    """vireo's names for the pair columns: "permutations" (0.5.7) or "combinations" (0.5.9)."""
    gen = itertools.permutations if naming == "permutations" else itertools.combinations
    return [",".join(x) for x in gen(donors, 2)]


def donor_calls(donors, prob_singlet, prob_doublet, n_vars, prob_max_cut=0.9, prob_doublet_cut=0.9, min_vars=10):
    """write_donor_id's calls: (donor_id, prob_max, prob_doublet_max, best_singlet, best_pair_index).

    unassigned when prob_max < 0.9, doublet when the best pair's probability
    is at least 0.9, unassigned when fewer than 10 sites have reads; the later
    rule wins. best_pair_index is the argmax over pairs (combinations order),
    -1 when there are no pairs.
    """
    prob_max = np.max(prob_singlet, axis=1)
    if prob_doublet.shape[1]:
        prob_doublet_out = np.max(prob_doublet, axis=1)
        k = np.argmax(prob_doublet, axis=1)
    else:
        prob_doublet_out = np.zeros(prob_singlet.shape[0])
        k = np.full(prob_singlet.shape[0], -1, dtype=np.int64)
    donor_singlet = np.array(donors, "U100")[np.argmax(prob_singlet, axis=1)]
    donor_ids = donor_singlet.copy()
    donor_ids[prob_max < prob_max_cut] = "unassigned"
    donor_ids[prob_doublet_out >= prob_doublet_cut] = "doublet"
    donor_ids[np.asarray(n_vars) < min_vars] = "unassigned"
    return donor_ids, prob_max, prob_doublet_out, donor_singlet, k


def write_vireo_files(out_dir, donors, cells, n_vars, prob_singlet, prob_doublet, doublet_llr,
                      naming="permutations", prob_doublet_column="prob_doublet"):
    """Write donor_ids.tsv, summary.tsv, prob_singlet.tsv.gz and prob_doublet.tsv.gz as vireo does.

    naming: "permutations" writes 0.5.7's pair names (the k-th permutation for
    the k-th pair, in best_doublet and the prob_doublet header);
    "combinations" writes 0.5.9's. prob_doublet_column: "prob_max_doublet"
    for 0.5.10's column name. The text is vireo's, decompressed; the gzip
    bytes differ from vireo's (header and compression level). Names in
    best_singlet and best_doublet are cut at 100 characters, as vireo's
    "U100" arrays cut them. Returns the paths.
    """
    from demux_tools._io import open_text

    os.makedirs(out_dir, exist_ok=True)
    donor_ids, prob_max, prob_doublet_out, donor_singlet, k = donor_calls(donors, prob_singlet, prob_doublet, n_vars)
    doublet_names = pair_names(donors, naming)
    cut = np.array(doublet_names, "U100")
    donor_doublet = [str(cut[x]) if x >= 0 else "" for x in k]
    paths = {}

    paths["summary"] = os.path.join(out_dir, "summary.tsv")
    donor_ids_uniq, donor_ids_count = np.unique(donor_ids, return_counts=True)
    with open_text(paths["summary"], "wt") as fid:
        fid.write("Var1\tFreq\n")
        for i in range(len(donor_ids_uniq)):
            fid.write("%s\t%d\n" % (donor_ids_uniq[i], donor_ids_count[i]))

    paths["donor_ids"] = os.path.join(out_dir, "donor_ids.tsv")
    header = ["cell", "donor_id", "prob_max", prob_doublet_column, "n_vars",
              "best_singlet", "best_doublet", "doublet_logLikRatio"]
    with open_text(paths["donor_ids"], "wt") as fid:
        fid.write("\t".join(header) + "\n")
        for i in range(len(cells)):
            line = [cells[i], donor_ids[i], "%.2e" % prob_max[i],
                    "%.2e" % prob_doublet_out[i], "%d" % n_vars[i],
                    donor_singlet[i], donor_doublet[i],
                    "%.3f" % doublet_llr[i]]
            fid.write("\t".join(line) + "\n")

    paths["prob_singlet"] = os.path.join(out_dir, "prob_singlet.tsv.gz")
    with open_text(paths["prob_singlet"], "wt") as fid:
        fid.write("\t".join(["cell"] + list(donors)) + "\n")
        for i in range(len(cells)):
            line = ["%.2e" % x for x in prob_singlet[i, :]]
            fid.write("\t".join([cells[i]] + line) + "\n")

    paths["prob_doublet"] = os.path.join(out_dir, "prob_doublet.tsv.gz")
    with open_text(paths["prob_doublet"], "wt") as fid:
        fid.write("\t".join(["cell"] + doublet_names) + "\n")
        for i in range(len(cells)):
            line = ["%.2e" % x for x in prob_doublet[i, :]]
            fid.write("\t".join([cells[i]] + line) + "\n")
    return paths
