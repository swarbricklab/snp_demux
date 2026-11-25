"""Doublet arithmetic: expected doublets, involvement of the largest donor, UMI ratio, calls with no male cell.

Definitions (worked numbers in docs/experiment_tools/scorecard.md):

- Expected doublets: a capture of n cells holds rate * n / 1000 doublets as a
  share of its cells (rate 0.008, the 0.8 per cent per 1,000 cells rule of
  thumb). Of those, a share 1 - sum(p_i^2) pairs two different donors, where
  p_i are the donors' shares of the capture's singlets.
- Involvement of the largest donor: over doublet calls with a named pair, how
  often the capture's largest donor (share p) is in the pair, over the rate
  expected if doublets paired donors in proportion to their shares,
  2 p (1 - p) / (1 - sum(p_i^2)). Needs 20 named pairs; NaN when one donor
  holds every singlet (no expectation).
- UMI ratio: a cell's UMIs over the median UMIs of the same run's singlets of
  its cell type in its capture. Two cells together read near 2.
- Calls with no male cell: a doublet call that names one male and one female
  donor and holds a male cell lacks a chrY UMI with the chance p that a male
  singlet at half its UMIs lacks one. A call with no male cell lacks one almost
  always. With n calls, O of them without a chrY UMI and E = sum of p over the
  calls, the share of calls with no male cell is (O - E) / (n - E): each call
  with no male cell adds (1 - p) to O - E, not 1. p is read per bin of 250
  UMIs from the run's male singlets; a call whose bin holds no male singlet
  takes the nearest bin's p (the lower one on a tie), and such calls are
  counted.
"""
import numpy as np
import pandas as pd

RATE_PER_1000 = 0.008
BIN_WIDTH = 250
MAX_BIN = 80


def expected_total(n_cells, rate_per_1000=RATE_PER_1000):
    """Expected doublet share of a capture of n_cells cells."""
    return rate_per_1000 * n_cells / 1000.0


def expected_between(n_cells, shares, rate_per_1000=RATE_PER_1000):
    """Expected share of cells that are doublets of two different donors; shares: donor shares of singlets."""
    shares = np.asarray(shares, dtype=float)
    return expected_total(n_cells, rate_per_1000) * (1.0 - (shares ** 2).sum())


def involvement(pairs, shares, min_pairs=20):
    """(observed over expected rate of the largest donor in the pairs, named pairs).

    pairs: iterable of "a,b" texts (others ignored). shares: pandas Series of
    donor shares, largest first (value_counts(normalize=True)). NaN when there
    are fewer than min_pairs named pairs or no shares.
    """
    pairs = [p.split(",") for p in pairs if p.count(",") == 1]
    if len(pairs) < min_pairs or len(shares) == 0:
        return float("nan"), len(pairs)
    top = shares.index[0]
    observed = np.mean([top in p for p in pairs])
    p = float(shares.iloc[0])
    spread = 1 - float((shares ** 2).sum())
    if spread <= 0 or p >= 1:
        return float("nan"), len(pairs)
    expected = 2 * p * (1 - p) / spread
    return observed / expected, len(pairs)


def umi_ratio(umis, singlet, capture, cell_type):
    """Each cell's UMIs over the median UMIs of the singlets of its capture and cell type.

    umis: float Series (NaN when unknown); singlet: bool Series (the run's
    singlets); capture, cell_type: text Series ("" is a cell type of its own).
    Returns a float Series aligned with umis.
    """
    frame = pd.DataFrame({"capture": capture, "cell_type": cell_type, "umis": umis})
    ref = frame[np.asarray(singlet, dtype=bool)].groupby(["capture", "cell_type"]).umis.median().rename("ref")
    return frame.umis / frame.join(ref, on=["capture", "cell_type"]).ref


def p_no_chry(umis, chry, keep, bin_width=BIN_WIDTH, max_bin=MAX_BIN):
    """P(no chrY UMI) by UMI bin (bins of bin_width, capped at max_bin), over the cells in keep.

    keep: bool Series selecting male singlets with a sex call and known UMIs.
    Returns a Series indexed by bin (float), sorted.
    """
    m = np.asarray(keep, dtype=bool)
    return (chry[m] == 0).groupby((umis[m] // bin_width).clip(upper=max_bin)).mean()


def mixed_pairs(pairs, males, all_males):
    """True for a pair "a,b" with exactly one donor in `males` and not both in `all_males`."""
    males, all_males = set(males), set(all_males)
    out = []
    for p in pairs:
        s = p.split(",")
        out.append(len(s) == 2 and ((s[0] in males) != (s[1] in males))
                   and not (s[0] in all_males and s[1] in all_males))
    return np.array(out, dtype=bool)


FILLS = ("nearest", "highest")


def nearest_bin_values(bins, p0, fill="nearest"):
    """(values, filled): p0 at each bin, a bin p0 lacks taking the nearest one's value (lower on a tie).

    bins: float array of bin numbers; p0: Series indexed by bin, sorted.
    fill: "nearest", or "highest" to give every missing bin the highest bin's
    value (the rule of the first scorecard, which undercounts E for calls
    below the lowest observed bin). filled: how many bins were not in p0.
    NaN values when p0 is empty.
    """
    if fill not in FILLS:
        raise ValueError("fill must be one of %s, not %r" % (", ".join(FILLS), fill))
    b = np.asarray(bins, dtype=float)
    if len(p0) == 0:
        return np.full(len(b), np.nan), int(len(b))
    idx = np.asarray(p0.index, dtype=float)
    vals = p0.to_numpy(dtype=float)
    pos = np.searchsorted(idx, b)
    lo = np.clip(pos - 1, 0, len(idx) - 1)
    hi = np.clip(pos, 0, len(idx) - 1)
    take_hi = np.abs(idx[hi] - b) < np.abs(b - idx[lo])
    near = np.where(take_hi, hi, lo)
    missing = idx[near] != b
    if fill == "highest":
        near = np.where(missing, len(idx) - 1, near)
    return vals[near], int(missing.sum())


def no_male_share(pairs, sel, has_sex, umis, chry, p0, males, all_males=None,
                  bin_width=BIN_WIDTH, max_bin=MAX_BIN, return_filled=False, fill="nearest"):
    """(calls, observed, expected, share) for doublet calls in sel naming one of `males` and a female donor.

    pairs: text Series of donor pairs; sel: bool Series of the doublet calls
    considered; has_sex: bool Series of cells with a sex call; umis, chry:
    float Series; p0: from p_no_chry. all_males: every male donor (a pair of
    two males is never counted), default `males`. A call whose bin p0 lacks
    takes the nearest bin's value (nearest_bin_values, `fill`); with return_filled a
    fifth value counts those calls. share is NaN when there are no more calls
    than expected.
    """
    all_males = males if all_males is None else all_males
    mixed = mixed_pairs(pairs, males, all_males)
    s = np.asarray(sel, dtype=bool) & mixed & np.asarray(has_sex, dtype=bool) & umis.notna().to_numpy()
    bins = ((umis[s] / 2) // bin_width).clip(upper=max_bin)
    values, filled = nearest_bin_values(bins.to_numpy(), p0, fill)
    expected = float(values.sum()) if len(values) else 0.0
    observed = int((chry[s] == 0).sum())
    n = int(s.sum())
    share = (observed - expected) / (n - expected) if n > expected else float("nan")
    if return_filled:
        return n, observed, expected, share, filled
    return n, observed, expected, share


def tool_groups(called, labels):
    """Doublet tool group per cell from per-tool doublet calls.

    called: list of bool arrays, one per tool; labels: their names. One tool
    gives "<a>" or "not <a>"; two give "both tools", "<a> only", "<b> only" or
    "neither <a> nor <b>".
    """
    if len(called) == 1:
        return np.where(called[0], labels[0], "not " + labels[0]).astype(object)
    if len(called) != 2:
        raise ValueError("doublet tool groups need one or two tools, not %d" % len(called))
    a, b = (np.asarray(x, dtype=bool) for x in called)
    return np.select([a & b, a, b], ["both tools", labels[0] + " only", labels[1] + " only"],
                     "neither %s nor %s" % (labels[0], labels[1])).astype(object)


def group_names(labels):
    """The doublet groups `tool_groups` writes, in the order reports list them (without "neither")."""
    if len(labels) == 1:
        return [labels[0]]
    return ["both tools", labels[0] + " only", labels[1] + " only"]


def pick_pair(called, pairs):
    """Per cell: the pair of the first tool that called a doublet; the first tool's pair when none did.

    called: list of bool arrays; pairs: list of text arrays, one per tool.
    """
    out = np.asarray(pairs[0], dtype=object).copy()
    done = np.asarray(called[0], dtype=bool).copy()
    for c, p in zip(called[1:], pairs[1:]):
        take = np.asarray(c, dtype=bool) & ~done
        out[take] = np.asarray(p, dtype=object)[take]
        done |= take
    return out
