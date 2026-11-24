# Data dictionary

Column definitions for snp_demux outputs that other pipelines consume.

General conventions:
- Missing values are written as `NA`.
- Derived floats are rounded to 6 decimal places.
- Sample ids are often bare digits, so read every id column as a string. A reader that infers types will read `007` as `7`.

## `doublet_summary.csv`

**Paths:**
- `assignment/<capture_id>/doublet_summary.csv`: one row for one capture (rule `doublet_summary`).
- `assignment/report/doublet_summary.csv`: the per-capture rows of every capture, concatenated (rule `combine_doublet_summary`).

**Purpose:** the inputs needed to turn the doublets that genotype demultiplexing *can* see into an expected *total* doublet rate for the capture. Downstream expression-based doublet detection (e.g. scDblFinder in mkobj) uses that rate.

**Inputs, per capture:**
- the consensus call per barcode, from `cell_assignment.tsv` (`status`, `assignment`);
- vireo's per-barcode doublet posterior, from `sample_barcode_assignments.tsv` (`Vireo_prob_doublet`);
- the samples declared for the capture, from the samples sheet (`config["deps"]["samples"]`).

### Background

Genotype demultiplexing only sees **cross-donor** doublets, meaning two cells from different samples. If droplets pair up cells at random, a doublet is cross-donor with probability

    1 − Σ p_i²

where `p_i` is sample *i*'s share of the capture's cells. Σ p_i² is the chance that two randomly paired cells come from the same donor, so the cross-donor share is its complement. The observed cross-donor doublet rate (written *f_x* in the formulas below) and the sample proportions `p_i` together determine the expected total doublet rate. This file records both.

**Symbols used in the formulas:** *f_x* = `doublet_rate_cross_donor`; Σ p_i² = `same_donor_pair_prob`.

### Columns

| column | type | definition |
|---|---|---|
| `capture_id` | string | The capture this row describes. |
| `n_barcodes` | integer | Number of barcodes in the capture's `cell_assignment.tsv`, i.e. every barcode the demultiplexers were run on. This is the denominator of both `doublet_rate_cross_donor*` columns. |
| `n_singlet` | integer | Barcodes whose consensus `status` is `singlet`. |
| `n_doublet_hard` | integer | Barcodes whose consensus `status` is `doublet`: the hard (thresholded) doublet calls. |
| `n_unassigned` | integer | Barcodes whose consensus `status` is `unassigned`: the consensus couldn't call them. This includes a singlet call that couldn't be resolved to a sample of this capture. |
| `sum_prob_doublet` | float | Sum of `Vireo_prob_doublet` over the capture's barcodes, skipping barcodes vireo didn't score (`NA`). |
| `doublet_rate_cross_donor` | float in [0, 1] | `n_doublet_hard / n_barcodes`: the observed cross-donor doublet rate from the hard calls (*f_x*). |
| `doublet_rate_cross_donor_posterior` | float in [0, 1] | `sum_prob_doublet / n_barcodes`: the expected cross-donor doublet rate from vireo's posteriors. It counts the uncertain barcodes that the hard calls leave out. |
| `n_samples` | integer | Number of `sample_id`s declared for the capture in the samples sheet, **including passengers** (declared samples that contributed no cells). |
| `same_donor_pair_prob` | float in (0, 1] | `Σ p_i²` over all declared samples: the probability that two randomly paired cells come from the same sample. `p_i` = consensus singlets assigned to sample *i* ÷ consensus singlets assigned to any declared sample. `1` for a single-sample capture. `NA` if the capture has no singlets assigned to a declared sample. |
| `sample_props` | string | Every declared sample's `p_i`, as `sample:p;sample:p;…`, in samples-sheet order. Passengers appear as `sample:0.0`. `NA` when `same_donor_pair_prob` is `NA`. |

### Rules behind the columns

- **The consensus is whichever one is in force.** `status` and `assignment` come from `cell_assignment.tsv`, so the counts follow the consensus that produced it.
- **The status counts don't have to add up to `n_barcodes`.** `n_singlet + n_doublet_hard + n_unassigned` equals `n_barcodes` only when every barcode has one of those three statuses. A missing `status` is counted in `n_barcodes` but in none of the three.
- **`p_i` uses singlets only.** Doublets and unassigned barcodes have no single sample. Singlets assigned to a sample *not* declared for the capture are ignored, so `p_i` always sums to 1 over the declared samples.
- **Passengers stay in `n_samples` and `same_donor_pair_prob`.** A declared sample with no cells contributes `p_i = 0`. Keeping it in `sample_props` makes the dropout visible.
- **`doublet_rate_cross_donor_posterior` is not a bound in either direction.** vireo's `prob_doublet` is the posterior of the single most likely donor **pair**, not the sum over all pairs. Where several pairs are plausible, a barcode's total doublet probability is higher than its `prob_doublet`. On the other hand, the posterior counts partial doublet probability on many barcodes that no hard call includes. On the OneK1K test fixture, `doublet_rate_cross_donor_posterior` (0.071) is above `doublet_rate_cross_donor` (0.047).
- **Within-donor doublets are invisible to both `doublet_rate_cross_donor*` columns.** That includes two cells from one sample, and same-sample aggregates from incomplete dissociation.

### Worked example

Take a capture with 4 declared samples, A–D, and 10,000 barcodes:
- 600 consensus doublets;
- 2,350 singlets from each sample.

Then:
- `n_barcodes = 10000`, `n_singlet = 9400`, `n_doublet_hard = 600`, `n_unassigned = 0`.
- `doublet_rate_cross_donor = 600 / 10000 = 0.06`.
- `p_A = p_B = p_C = p_D = 2350 / 9400 = 0.25`, so `same_donor_pair_prob = 4 × 0.25² = 0.25`.
- `sample_props = A:0.25;B:0.25;C:0.25;D:0.25`.
- Under random pairing, `1 − 0.25 = 0.75` of doublets are cross-donor. The 6% observed are therefore three quarters of an expected 8% total.
