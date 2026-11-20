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

**Purpose:** the doublets that genotype demultiplexing *can* see, the inputs needed to turn them into an expected *total* doublet rate for the capture, and that rate. Downstream expression-based doublet detection (e.g. scDblFinder in mkobj) uses that rate.

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
| `n_recovered` | integer | Barcodes whose consensus `status` is `recovered`: exactly one method called them, as a singlet of a named donor, and that method is in `params.consensus.recover_from`. They are not counted in `n_singlet` and do not contribute to `p_i`. |
| `n_unassigned` | integer | Barcodes whose consensus `status` is `unassigned`: the consensus couldn't call them. This includes a singlet call that couldn't be resolved to a sample of this capture. |
| `sum_prob_doublet` | float | Sum of `Vireo_prob_doublet` over the capture's barcodes, skipping barcodes vireo didn't score (`NA`). |
| `doublet_rate_cross_donor` | float in [0, 1] | `n_doublet_hard / n_barcodes`: the observed cross-donor doublet rate from the hard calls (*f_x*). |
| `doublet_rate_cross_donor_posterior` | float in [0, 1] | `sum_prob_doublet / n_barcodes`: the expected cross-donor doublet rate from vireo's posteriors. It counts the uncertain barcodes that the hard calls leave out. |
| `n_samples` | integer | Number of `sample_id`s declared for the capture in the samples sheet, **including passengers** (declared samples that contributed no cells). |
| `same_donor_pair_prob` | float in (0, 1] | `Σ p_i²` over all declared samples: the probability that two randomly paired cells come from the same sample. `p_i` = consensus singlets assigned to sample *i* ÷ consensus singlets assigned to any declared sample. `1` for a single-sample capture. `NA` if the capture has no singlets assigned to a declared sample. |
| `sample_props` | string | Every declared sample's `p_i`, as `sample:p;sample:p;…`, in samples-sheet order. Passengers appear as `sample:0.0`. `NA` when `same_donor_pair_prob` is `NA`. |
| `loading_rate_per_1k_cells` | float | Doublet rate per 1,000 cells **recovered** (here, barcodes demultiplexed), from 10x's multiplet-rate table for the capture's chip: e.g. `0.008` means 0.8% doublets per 1,000 cells. The chip comes from an optional `chip` column in `captures.csv` (blank → `params.doublet.default_chip`), and its rate from `params.doublet.rate_per_1k`. Defaults: `standard: 0.008`, `HT: 0.004`. Chip names match case-insensitively; a chip with no rate is an error. |
| `doublet_rate_total` | float | Expected **total** doublet rate: *f_x* / (1 − Σ p_i²), i.e. `doublet_rate_cross_donor / (1 − same_donor_pair_prob)`. This is the hard cross-donor rate scaled up by the share of doublets that are cross-donor under random pairing. `NA` when `doublet_rate_fallback` is set or `n_barcodes` is 0. |
| `doublet_rate_total_sd` | float | Standard error of `doublet_rate_total`: the binomial SE of *f_x*, √(*f_x* (1 − *f_x*) / `n_barcodes`), divided by (1 − Σ p_i²). Uncertainty in `same_donor_pair_prob` is ignored. `NA` whenever `doublet_rate_total` is. |
| `doublet_rate_within_donor` | float | Expected **within-donor** doublet rate, invisible to genotype methods: `doublet_rate_total − doublet_rate_cross_donor`. `NA` whenever `doublet_rate_total` is. |
| `doublet_rate_loading` | float | What the chip's loading curve predicts: `n_barcodes / 1000 × loading_rate_per_1k_cells`. Always computed. `n_barcodes` stands in for cells recovered, so it understates them when the barcode list is a subset of the capture or was filtered before demultiplexing, and `doublet_rate_vs_loading` is then too high. The OneK1K test fixture is such a subset: 1,400 of a capture's 20,982 barcodes, so its loading rate is 0.011 rather than 0.168, and the ratio reads 12.2 rather than 0.81. |
| `doublet_rate_vs_loading` | float | `doublet_rate_total / doublet_rate_loading`, as a QC: about 1 means the genotype-based estimate agrees with the loading curve; well above 1 suggests overloading, cell clumping or a demultiplexing problem. `NA` when `doublet_rate_total` is `NA` or `doublet_rate_loading` is 0. |
| `doublet_rate_fallback` | string | Empty when the rates were computed. Otherwise why `doublet_rate_total`, `doublet_rate_total_sd` and `doublet_rate_within_donor` are `NA`: `single_sample` (`n_samples < 2`, so no doublet can be cross-donor); `no_singlets` (`same_donor_pair_prob` is `NA`); `one_sample_only` (every singlet came from one declared sample, so `1 − same_donor_pair_prob = 0`). Written as an empty field, which most CSV readers return as missing. |

### Rules behind the columns

- **The consensus is whichever one is in force.** `status` and `assignment` come from `cell_assignment.tsv`, so the counts follow the consensus that produced it.
- **The status counts don't have to add up to `n_barcodes`.** `n_singlet + n_doublet_hard + n_recovered + n_unassigned` equals `n_barcodes` only when every barcode has one of those four statuses. A missing `status` is counted in `n_barcodes` but in none of the four.
- **`p_i` uses singlets only.** Doublets and unassigned barcodes have no single sample, and recovered barcodes rest on one method's call. Singlets assigned to a sample *not* declared for the capture are ignored, so `p_i` always sums to 1 over the declared samples. A sample whose cells are more often left unassigned or called doublets is under-counted, which raises `same_donor_pair_prob`. This can be the case for a sample without a genotype. For example, two samples of 500 cells each give 0.5; if half of one sample's cells are left unassigned, the singlets are 500 and 250, `same_donor_pair_prob` is 0.556, and the factor 1 / (1 − Σ p_i²) rises from 2.0 to 2.25.
- **Passengers stay in `n_samples` and `same_donor_pair_prob`.** A declared sample with no cells contributes `p_i = 0`. Keeping it in `sample_props` makes the dropout visible.
- **`doublet_rate_cross_donor_posterior` is not a bound in either direction.** vireo's `prob_doublet` is the posterior of the single most likely donor **pair**, not the sum over all pairs. Where several pairs are plausible, a barcode's total doublet probability is higher than its `prob_doublet`. On the other hand, the posterior counts partial doublet probability on many barcodes that no hard call includes. On the OneK1K test fixture, `doublet_rate_cross_donor_posterior` is 0.071 and vireo's own hard doublets are 0.065 of barcodes. `doublet_rate_cross_donor` is not directly comparable: it counts the consensus of four methods, while the posterior column is vireo alone. It depends on the consensus in force: 0.047 on the fixture under the informative-majority consensus (the default), 0.126 under `at_least_half_singlet`. The posterior also inherits vireo's doublet prior, which is the number of barcodes ÷ 100,000 (capped at 0.5), so it moves with the size of the barcode list.
- **A false cross-donor doublet is scaled up like a true one.** `n_doublet_hard` includes any singlet wrongly called a doublet, for example one carrying ambient RNA from another donor, which from allele counts alone can look like an unequal doublet. Dividing by `1 − same_donor_pair_prob` to reach a total rate multiplies those as well: by 1.14 for eight equal samples (two points of false doublets add 2.3 points) and by 5.56 for two samples at 90 and 10 per cent (11.1 points).
- **The expected rates are biased low, with one exception.** `doublet_rate_total` and `doublet_rate_within_donor` scale up `doublet_rate_cross_donor`, which undercounts: genotype methods miss some cross-donor doublets, and same-sample aggregates from incomplete dissociation are always within-donor. False cross-donor doublets (above) push the other way, so the rates are not strict lower bounds. mkobj passes `doublet_rate_total` / `doublet_rate_total_sd` to scDblFinder as `dbr` / `dbr.sd`.
- **Within-donor doublets are invisible to both `doublet_rate_cross_donor*` columns.** That includes two cells from one sample, and same-sample aggregates from incomplete dissociation.

### Worked example

Take a capture with 4 declared samples, A–D, and 10,000 barcodes:
- 600 consensus doublets;
- 2,350 singlets from each sample.

Then:
- `n_barcodes = 10000`, `n_singlet = 9400`, `n_doublet_hard = 600`, `n_recovered = 0`, `n_unassigned = 0`.
- `doublet_rate_cross_donor = 600 / 10000 = 0.06`.
- `p_A = p_B = p_C = p_D = 2350 / 9400 = 0.25`, so `same_donor_pair_prob = 4 × 0.25² = 0.25`.
- `sample_props = A:0.25;B:0.25;C:0.25;D:0.25`.
- Under random pairing, `1 − 0.25 = 0.75` of doublets are cross-donor. The 6% observed are therefore three quarters of an expected 8% total.
- `doublet_rate_total = 0.06 / 0.75 = 0.08`, and `doublet_rate_within_donor = 0.08 − 0.06 = 0.02`.
- `doublet_rate_total_sd = sqrt(0.06 × 0.94 / 10000) / 0.75 ≈ 0.003166`.
- On a standard chip (`loading_rate_per_1k_cells = 0.008`): `doublet_rate_loading = 10000 / 1000 × 0.008 = 0.08`, so `doublet_rate_vs_loading = 1.0` and `doublet_rate_fallback` is empty.
- On an HT chip (`0.004`), `doublet_rate_loading = 0.04` and `doublet_rate_vs_loading = 2.0`.
