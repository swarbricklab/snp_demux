"""Per-cell diagnostics on top of the Scorer: where a cell's margin comes from, and what vireo's
mean-field average hides.

Every function takes a loaded capture (scorer.Capture), the allele ratios
(scorer.Fit) and, where it applies, each cell's best and runner-up donor as
indices into capture.donors. A cell's log-likelihood under a donor is a sum
of one term per site where the cell has reads (`entry_terms`), so a margin
splits exactly into the sites' contributions.

- `group_margins`: the margin between the best and runner-up donor, summed
  within each group of sites (region, INFO bin, which donors are blank...).
  The groups' margins add up to the margin.
- `pairwise_shared`: each pair of donors compared only on the sites where
  both are genotyped, so a blank genotype cannot decide between them.
- `marginal_singlet`: the exact marginal log-likelihood log sum_g GP_g
  exp(t_g) at sites with a soft or blank genotype, beside vireo's mean-field
  sum_g GP_g t_g. The two are equal at hard genotypes, and the exact value is
  never below the mean-field one.
- `top_sites`: the sites that carry most of each cell's margin.
- `expected_fit`: the score per UMI a true singlet of the best donor would
  reach on the same sites, to set beside `fit_per_umi`.
- `cheap_columns`: ALT UMIs, and the sites and UMIs where the best and
  runner-up donors' genotypes differ.

All use vireo's arithmetic (the genotypes normalised as vireo's model holds
them, the digamma expectations of the doublet step).
"""
import numpy as np

from demux_tools.vireo_port import normalize


class Entries(object):
    """The (site, cell) entries with reads: site, cell, alt and ref UMIs, in cell-major order."""

    def __init__(self, capture):
        dp = capture.DP.tocoo()
        order = np.lexsort((dp.row, dp.col))
        self.site = dp.row[order].astype(np.int64)
        self.cell = dp.col[order].astype(np.int64)
        depth = dp.data[order].astype(np.int64)
        ad = capture.AD.tocoo()
        n_sites = capture.n_sites
        key_d = self.cell * n_sites + self.site
        key_a = ad.col.astype(np.int64) * n_sites + ad.row.astype(np.int64)
        alt = np.zeros(len(key_d), dtype=np.int64)
        if len(key_a) and len(key_d):
            # AD and DP are in canonical form (duplicates summed), so each AD entry has one DP entry
            pos = np.searchsorted(key_d, key_a)
            ok = (pos < len(key_d)) & (key_d[np.minimum(pos, len(key_d) - 1)] == key_a)
            alt[pos[ok]] = ad.data[ok].astype(np.int64)
        keep = depth > 0
        self.site, self.cell, self.alt = self.site[keep], self.cell[keep], alt[keep]
        self.ref = depth[keep] - self.alt
        self.n_cells = capture.n_cells

    def __len__(self):
        return len(self.site)


def _expectations(fit):
    from demux_tools.vireo_port import digamma_terms

    dg1, dg2, dgs = digamma_terms(fit.theta_mean, fit.theta_sum)
    return (dg1 - dgs).ravel(), (dg2 - dgs).ravel()


def site_donor_terms(capture, fit):
    """(alt_term, ref_term), each (sites, donors): the term of one ALT and of one REF UMI at a site under a donor.

    A cell's term at a site under donor d is alt * alt_term + ref * ref_term,
    which equals sum_g GP_g (alt E[log theta_g] + ref E[log(1 - theta_g)]).
    """
    e_alt, e_ref = _expectations(fit)
    G = normalize(capture.GPb)
    return G @ e_alt, G @ e_ref


def entry_terms(capture, fit, entries, donor_of_entry, per_site=None):
    """Each entry's log-likelihood term under one donor per entry (an index array, one per entry)."""
    alt_term, ref_term = per_site or site_donor_terms(capture, fit)
    return entries.alt * alt_term[entries.site, donor_of_entry] + entries.ref * ref_term[entries.site, donor_of_entry]


def best_and_second(scores):
    """Index of each cell's best and runner-up donor, by log-likelihood, as score() orders them."""
    order = np.argsort(-scores.ll_singlet, axis=1)
    second = order[:, 1] if order.shape[1] > 1 else order[:, 0]
    return order[:, 0], second


def deltas(capture, fit, entries, best, second):
    """(term_best, term_second) per entry; their difference summed over a cell's entries is its margin."""
    per_site = site_donor_terms(capture, fit)
    tb = entry_terms(capture, fit, entries, best[entries.cell], per_site)
    ts = entry_terms(capture, fit, entries, second[entries.cell], per_site)
    return tb, ts


def group_margins(capture, fit, best, second, labels, name, entries=None):
    """Margin, UMIs and sites with reads per group of sites, per cell (a DataFrame).

    labels: one label per site. Columns margin__<name>__<label>,
    umis__<name>__<label>, nsites__<name>__<label>, labels in sorted order.
    """
    import pandas as pd

    entries = entries or Entries(capture)
    tb, ts = deltas(capture, fit, entries, best, second)
    values, codes = np.unique(np.asarray(labels, dtype=str), return_inverse=True)
    n_cells, k = capture.n_cells, len(values)
    flat = entries.cell * k + codes[entries.site]
    margin = np.bincount(flat, weights=tb - ts, minlength=n_cells * k).reshape(n_cells, k)
    umis = np.bincount(flat, weights=entries.alt + entries.ref, minlength=n_cells * k).reshape(n_cells, k)
    nsites = np.bincount(flat, minlength=n_cells * k).reshape(n_cells, k)
    cols = {}
    for j, v in enumerate(values):
        cols["margin__%s__%s" % (name, v)] = margin[:, j]
    for j, v in enumerate(values):
        cols["umis__%s__%s" % (name, v)] = umis[:, j].astype(np.int64)
    for j, v in enumerate(values):
        cols["nsites__%s__%s" % (name, v)] = nsites[:, j].astype(np.int64)
    return pd.DataFrame(cols)


def pairwise_shared(capture, fit):
    """Each pair of donors compared on the sites where both are genotyped (a DataFrame).

    For donors d and e, M[c, d, e] sums the cell's terms under d minus under e
    over the sites where neither is blank. A donor beats another when M > 0.
    Columns: pairwise_best (most wins; ties go to the higher total of M),
    pairwise_wins, pairwise_margin (the best donor's smallest M over the
    others) and pairwise_shared_umis (the cell's UMIs at the sites the best
    donor shares with the donor giving that smallest M). With no blank
    genotype, pairwise_best is the best singlet and pairwise_margin its margin.
    Sites are grouped by their pattern of blank donors, one log-likelihood
    matrix per pattern. M and the shared UMIs take 16 x cells x donors^2
    bytes (about 2.4 GB for 30,000 cells and 70 donors).
    """
    import pandas as pd

    from demux_tools.scorer import loglik_singlet

    missing = np.asarray(capture.missing, dtype=bool)
    n_cells, n = capture.n_cells, capture.n_donors
    M = np.zeros((n_cells, n, n))
    U = np.zeros((n_cells, n, n))
    if capture.n_sites:
        patterns, which = np.unique(missing, axis=0, return_inverse=True)
        which = np.asarray(which).reshape(-1)
    else:
        patterns, which = np.zeros((0, n), dtype=bool), np.zeros(0, dtype=np.int64)
    for p, pattern in enumerate(patterns):
        rows = np.flatnonzero(which == p)
        ll = loglik_singlet(capture.AD[rows], capture.DP[rows], capture.GPb[rows], fit.theta_mean, fit.theta_sum)
        umis = np.asarray(capture.DP[rows].sum(axis=0)).ravel()
        ok = ~pattern
        both = np.outer(ok, ok)
        M += (ll[:, :, None] - ll[:, None, :]) * both[None, :, :]
        U += umis[:, None, None] * both[None, :, :]
    eye = np.eye(n, dtype=bool)
    wins = ((M > 0) & ~eye[None, :, :]).sum(axis=2)
    total = np.where(eye[None, :, :], 0.0, M).sum(axis=2)
    best = np.array([max(range(n), key=lambda d: (wins[c, d], total[c, d], -d)) for c in range(n_cells)],
                    dtype=np.int64) if n_cells else np.zeros(0, dtype=np.int64)
    rows = np.arange(n_cells)
    if n > 1:
        others = np.where(eye[best], np.inf, M[rows, best, :])
        rival = np.argmin(others, axis=1)
        margin = others[rows, rival]
        shared = U[rows, best, rival].astype(np.int64)
    else:
        margin = np.full(n_cells, np.nan)
        shared = np.zeros(n_cells, dtype=np.int64)
    return pd.DataFrame({"pairwise_best": np.array(capture.donors, dtype=object)[best] if n_cells else [],
                         "pairwise_wins": wins[rows, best] if n_cells else [],
                         "pairwise_margin": margin, "pairwise_shared_umis": shared})


def marginal_singlet(capture, fit, ll_singlet, best, entries=None):
    """vireo's mean-field log-likelihoods with the exact marginal at soft and blank genotypes (a DataFrame).

    At a site where donor d's genotype probabilities GP are not one-hot, vireo
    adds sum_g GP_g t_g, where t_g is the term of the cell's UMIs under
    genotype g. The exact marginal is log sum_g GP_g exp(t_g). Columns:
    best_marginal (donor with the highest exact value), ll_best_marginal,
    ll_second_marginal, margin_marginal, and blank_umis (the cell's UMIs at
    sites where its best donor, by vireo's score, is blank).
    """
    import pandas as pd

    entries = entries or Entries(capture)
    e_alt, e_ref = _expectations(fit)
    G = normalize(capture.GPb)
    soft = G.max(axis=2) < 1.0
    ll = np.array(ll_singlet, dtype=np.float64, copy=True)
    for d in range(capture.n_donors):
        rows = soft[:, d]
        if not rows.any():
            continue
        sel = np.flatnonzero(rows[entries.site])
        if not len(sel):
            continue
        t = entries.alt[sel, None] * e_alt[None, :] + entries.ref[sel, None] * e_ref[None, :]
        g = G[entries.site[sel], d, :]
        mean_field = np.sum(g * t, axis=1)
        with np.errstate(divide="ignore"):
            logs = np.log(g) + t
        top = np.max(logs, axis=1)
        exact = top + np.log(np.sum(np.exp(logs - top[:, None]), axis=1))
        ll[:, d] += np.bincount(entries.cell[sel], weights=exact - mean_field, minlength=capture.n_cells)
    order = np.argsort(-ll, axis=1)
    rows = np.arange(capture.n_cells)
    best_m = order[:, 0]
    llb = ll[rows, best_m]
    lls = ll[rows, order[:, 1]] if capture.n_donors > 1 else np.full(capture.n_cells, np.nan)
    blank_entry = capture.missing[entries.site, best[entries.cell]]
    blank_umis = np.bincount(entries.cell[blank_entry], weights=(entries.alt + entries.ref)[blank_entry],
                             minlength=capture.n_cells).astype(np.int64)
    return pd.DataFrame({"best_marginal": np.array(capture.donors, dtype=object)[best_m],
                         "ll_best_marginal": llb, "ll_second_marginal": lls, "margin_marginal": llb - lls,
                         "blank_umis": blank_umis}), ll


def top_sites(capture, fit, best, second, n=10, labels=None, gene=None, entries=None):
    """The sites that carry each cell's margin.

    Returns (long, per_cell). long: one row per cell and rank (the n largest
    contributions, best donor's term minus runner-up's): cell, rank, variant,
    group labels, alt, ref, gt_best, gt_second, term_best, term_second, delta,
    cum_share (cumulative delta over the cell's margin). per_cell:
    top1_site_share and top10_sites_share (the largest and the ten largest
    deltas over the margin), and with `gene` labels top_gene, top_gene_share
    and top10_genes_share (deltas summed per gene first; sites whose gene is
    "" are left out, and a cell without such sites gets top_gene "" and shares
    of 0). Shares can exceed 1 when other sites pull the other way; they and
    cum_share are NaN when the margin is 0.
    """
    import pandas as pd

    from demux_tools.scorer import genotype_text

    entries = entries or Entries(capture)
    tb, ts = deltas(capture, fit, entries, best, second)
    delta = tb - ts
    n_cells = capture.n_cells
    margin = np.bincount(entries.cell, weights=delta, minlength=n_cells)
    order = np.lexsort((-delta, entries.cell))
    cell_sorted = entries.cell[order]
    start = np.searchsorted(cell_sorted, np.arange(n_cells))
    rank = np.arange(len(order)) - start[cell_sorted]
    d_sorted = delta[order]
    with np.errstate(divide="ignore", invalid="ignore"):
        top1 = np.where(margin != 0, np.bincount(cell_sorted[rank < 1], weights=d_sorted[rank < 1],
                                                 minlength=n_cells) / margin, np.nan)
        top10 = np.where(margin != 0, np.bincount(cell_sorted[rank < 10], weights=d_sorted[rank < 10],
                                                  minlength=n_cells) / margin, np.nan)
    per_cell = {"top1_site_share": top1, "top10_sites_share": top10}
    if gene is not None:
        gene = np.asarray(gene, dtype=str)
        named = gene[entries.site] != ""
        names, gcode = np.unique(gene[entries.site][named], return_inverse=True)
        key = entries.cell[named] * max(len(names), 1) + np.asarray(gcode).reshape(-1)
        uk, inv = np.unique(key, return_inverse=True)
        sums = np.bincount(np.asarray(inv).reshape(-1), weights=delta[named], minlength=len(uk))
        gcell = uk // max(len(names), 1)
        gname = names[uk % max(len(names), 1)] if len(names) else np.zeros(0, dtype=str)
        go = np.lexsort((-sums, gcell))
        gcell_s = gcell[go]
        gstart = np.searchsorted(gcell_s, np.arange(n_cells))
        grank = np.arange(len(go)) - gstart[gcell_s]
        top_gene = np.array([""] * n_cells, dtype=object)
        first = go[grank == 0]
        top_gene[gcell[first]] = gname[first]
        with np.errstate(divide="ignore", invalid="ignore"):
            per_cell["top_gene"] = top_gene
            per_cell["top_gene_share"] = np.where(
                margin != 0, np.bincount(gcell_s[grank < 1], weights=sums[go][grank < 1], minlength=n_cells) / margin,
                np.nan)
            per_cell["top10_genes_share"] = np.where(
                margin != 0, np.bincount(gcell_s[grank < 10], weights=sums[go][grank < 10], minlength=n_cells) / margin,
                np.nan)
    keep = order[rank < n]
    kd = delta[keep]
    kc = entries.cell[keep]
    with np.errstate(divide="ignore", invalid="ignore"):
        # cumulative within each cell (rows are grouped by cell, in rank order)
        csum = np.cumsum(kd)
        cstart = np.searchsorted(kc, kc, side="left")
        base = np.concatenate([[0.0], csum])[cstart]
        cum = np.where(margin[kc] != 0, (csum - base) / margin[kc], np.nan)
    sites = entries.site[keep]
    long = {"cell": np.array(capture.barcodes, dtype=object)[kc], "rank": rank[rank < n] + 1,
            "variant": np.array(capture.variants, dtype=object)[sites]}
    for lname, lab in (labels or {}).items():
        long[lname] = np.asarray(lab, dtype=object)[sites]
    long.update({"alt": entries.alt[keep], "ref": entries.ref[keep],
                 "gt_best": genotype_text(capture.GPb[sites, best[kc]], capture.missing[sites, best[kc]]),
                 "gt_second": genotype_text(capture.GPb[sites, second[kc]], capture.missing[sites, second[kc]]),
                 "term_best": tb[keep], "term_second": ts[keep], "delta": kd, "cum_share": cum})
    return pd.DataFrame(long), pd.DataFrame(per_cell)


def expected_fit(capture, fit, best, entries=None):
    """Score per UMI expected of a true singlet of the best donor on the cell's own sites and depths.

    At a site with n UMIs, a singlet of donor b has an expected ALT share
    theta_bar = sum_g GP_g theta_g (theta_mean of the fit), so its expected
    term is n * sum_g GP_g (theta_bar e_alt_g + (1 - theta_bar) e_ref_g).
    Summed over the cell's sites and divided by its UMIs; NaN without UMIs.
    """
    entries = entries or Entries(capture)
    alt_term, ref_term = site_donor_terms(capture, fit)
    theta_bar = normalize(capture.GPb) @ np.asarray(fit.theta_mean, dtype=np.float64).ravel()
    per_umi = theta_bar * alt_term + (1 - theta_bar) * ref_term          # (sites, donors)
    depth = entries.alt + entries.ref
    term = depth * per_umi[entries.site, best[entries.cell]]
    total = np.bincount(entries.cell, weights=term, minlength=capture.n_cells)
    umis = np.bincount(entries.cell, weights=depth, minlength=capture.n_cells)
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(umis > 0, total / umis, np.nan)


def cheap_columns(capture, fit, best, second, entries=None):
    """alt_umis, n_sites_diff, umis_diff and fit_expected_per_umi (a DataFrame).

    n_sites_diff and umis_diff count the cell's sites (and UMIs there) where
    the best and runner-up donors' genotype probabilities differ.
    """
    import pandas as pd

    entries = entries or Entries(capture)
    G = normalize(capture.GPb)
    _, code = np.unique(G.reshape(-1, G.shape[2]), axis=0, return_inverse=True)
    code = np.asarray(code).reshape(G.shape[:2])                       # one code per distinct genotype vector
    differ = code[entries.site, best[entries.cell]] != code[entries.site, second[entries.cell]]
    n_cells = capture.n_cells
    depth = entries.alt + entries.ref
    return pd.DataFrame({
        "alt_umis": np.bincount(entries.cell, weights=entries.alt, minlength=n_cells).astype(np.int64),
        "n_sites_diff": np.bincount(entries.cell[differ], minlength=n_cells).astype(np.int64),
        "umis_diff": np.bincount(entries.cell[differ], weights=depth[differ], minlength=n_cells).astype(np.int64),
        "fit_expected_per_umi": expected_fit(capture, fit, best, entries=entries),
    })
