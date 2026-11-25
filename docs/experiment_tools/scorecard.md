# Scorecard: `cell-table`, `scorecard`, `calibrate` and `doublet-excess`

These tools score sets of calls against named truth: per-cell sex calls, cell
types, donor attributes and the loading of each capture. They need numpy and
pandas only (`truth.py`, `doublets.py`, `scorecard.py` and
`data/cell_classes.tsv`).

Every number they write names its truth: a count of cells, sex, cell type,
decoy donors or UMIs. Study-specific inputs (paths, which donors are male,
which captures form a pool) come in as tables. The repository holds none.

## Inputs

| Table | Columns | Used by |
|---|---|---|
| Runs manifest | `capture`, `run`, `path`; optionally `kind`, `prob_doublet`, `label_map`, `label_pattern`, `cluster_map`, `sba` (as in [Compare runs](compare-runs.md)) | `cell-table` |
| Truth manifest | `capture`, `path`; optionally `columns` (comma-separated; default every column but `barcode`) | `cell-table`, `calibrate` |
| Captures | `capture`; optionally `pool`, `capture_set`, `donors`, `sister_key`, `genotyped_<run>` | all |
| Donors | `donor`; optionally `sex`, `immune_only`, `decoy` | `scorecard`, `calibrate` |
| Cell classes | `label`, `class` (`immune` or `non-immune`) | `scorecard`, `doublet-excess` |
| Ambient | `capture`, `vireo_ambient` (vireo's `ambient_summary.csv`), `souporcell_ambient` (`ambient_rna.txt`) | `doublet-excess` |

- Each truth file is a per-cell table with a `barcode` column, such as a sex
  call table (`n_umi`, `chrY_counts`, `xist_counts`, `sex_call`) or a cell
  type table (`coarse_cell_type`). Files of one capture are joined by
  barcode.
- A capture with no `pool` gets `donors:` plus its sorted donors, so captures
  without a pool that share their donors share a pool. `capture_set` defaults
  to `main`. `sister_key` defaults to the pool plus the sorted donors.
- `immune_only` marks a donor whose cells are all immune, such as a blood
  control. `decoy` marks a donor added to the genotypes that is not in the
  capture, such as a decoy scored with the Scorer.
- A donor not in the donor table has the sex `--default-sex` (default
  `unknown`, which is not scored).
- `data/cell_classes.tsv` holds 33 generic coarse cell-type labels: 12
  immune and 21 non-immune. Any other label is `unclear` and an empty label
  `no annotation`; neither enters the cell-type truth. Labels such as `stem
  cell` and `abnormal cell` can hold tumour cells and immune cells alike, so
  they are left out.

## `cell-table`

```bash
python -I tools/demux_tools cell-table --runs runs.tsv --truth truth.tsv --captures captures.tsv \
    --out cells.tsv.gz
```

One row per barcode of each capture: `capture`, `pool`, `capture_set`,
`barcode`, the truth columns, then `<run>_status`, `<run>_donor` and
`<run>_pair` for each run, in the manifest's order. A run without a cell has
empty fields there.

- The rows are the barcodes of every run, or of the runs named by
  `--universe`; the other runs are joined to them.
- A cluster map `majority:RUN` names a souporcell run's clusters by the donor
  most of their cells carry in RUN, counting cells both call singlet.

## `scorecard`

```bash
python -I tools/demux_tools scorecard --cells cells.tsv.gz --run cons --run legacy \
    --tools cons=vireo:cons_vireo,souporcell:cons_souporcell \
    --donors donors.tsv --captures captures.tsv --out scorecard/
```

`--tools RUN=LABEL:STATUS_RUN[/PAIR_RUN],...` names one or two tool runs of
the cell table for a scored run. Their doublet calls split its doublet tier
into groups, and they give the pair of each doublet call:

| vireo calls a doublet | souporcell calls a doublet | Group | Pair taken from |
|---|---|---|---|
| yes | yes | `both tools` | vireo |
| yes | no | `vireo only` | vireo |
| no | yes | `souporcell only` | souporcell |
| no | no | `neither vireo nor souporcell` | vireo |

`/PAIR_RUN` takes a tool's pair from another run than its calls, for example
Demuxafy's vireo call with the pair from vireo's own file.

### Measures

Each row of `scorecard.tsv` is one share per capture set and run:
`capture_set`, `measure`, `run`, `truth`, `value`, `low`, `high`,
`numerator`, `denominator`, `observed_no_chrY`, `expected_no_chrY`,
`calls_in_bins_without_data`.

| Measure | Truth | Numerator over denominator |
|---|---|---|
| Singlets kept, recovered tier, doublet tier, unassigned | count | That tier over all cells of the table |
| Wrong donor by sex: female share of singlets labelled a male donor (also per male donor, and for recovered cells) | sex | Female calls over cells with a male or female call |
| Wrong donor by sex: male share of singlets labelled a female donor, captures with a male donor | sex | Male calls over the same |
| Floor: male share of singlets, captures with no male donor | sex | The sex call's own false rate |
| Wrong donor by cell type: non-immune share of singlets (and recovered cells) labelled an immune-only donor | cell type | Non-immune over immune and non-immune cells |
| Wrong donor by decoy: share of singlets labelled a decoy donor | decoy | Written when the donor table has a decoy |
| Ambiguous-sex cells called doublet, or singlet, captures with a male donor | sex | Cells whose sex call is `ambiguous` (both chrY and XIST signal) |
| Ambiguous-sex share of cells with a sex call value, with and without a male donor | sex | Cells whose sex call is `ambiguous` over cells with any sex call value; run `any` |
| UMI ratio of doublet calls (and of singlets): share at 1.5 or more | UMIs | Cells with a UMI ratio |
| Male-female doublet calls: estimated share with no male cell | sex and UMIs | See below |
| UMI ratio of doublet calls: median | UMIs | Over all doublet calls |

Rows with no denominator are left out, and so are recovered-tier rows of a
run that has no recovered cell. "Captures with a male donor" are those whose
`donors` in the capture table include a male donor. Without a capture table,
or when it gives no donors for a capture, it is unknown which captures hold a
male donor: the six measures that need it (the floors, the male share of
singlets labelled a female donor and the ambiguous-sex rows) are not written,
and the command says so on standard error and in the provenance file.
Sex calls are compared in lower case.

The **UMI ratio** of a cell is its UMIs over the median UMIs of the same
run's singlets of its coarse cell type in its capture. Two cells together
read near 2.

The **share of male-female doublet calls with no male cell** rests on chrY.
A real male-female doublet gets about half its UMIs from the male cell, so it
lacks a chrY UMI with the chance p that a male singlet at half its UMIs lacks
one. A call with no male cell lacks one almost always. A worked example: 100
calls, each at 4,000 UMIs, where a male singlet at 2,000 UMIs lacks chrY with
p = 0.25.

| | Calls | Without a chrY UMI |
|---|---|---|
| Real doublets | 40 | 10 (40 × 0.25) |
| No male cell | 60 | 60 |
| Observed, O | 100 | 70 |
| Expected if all were real, E | | 25 (100 × 0.25) |

The share with no male cell is (O - E) / (n - E) = 45 / 75 = 0.6. Each call
with no male cell adds (1 - p) to O - E, not 1, so dividing by n would give
0.45. p is read per bin of 250 UMIs (capped at bin 80) from the run's male
singlets in the capture set. A call whose bin holds no male singlet takes the
nearest bin's p (the lower one on a tie), and `calls_in_bins_without_data`
counts those calls. `--no-chry-fill highest` gives every such call the
highest bin's p instead, as the first scorecard did; for calls below the
lowest observed bin that p is far too small, so E is too small and the share
too large. A pair of two male donors is never counted. The share is
empty when there are no more calls than expected, for example when every
male singlet lacks a chrY UMI.

### Intervals

An interval is the 2.5th and 97.5th percentile of the share over 1,000
resamples of pools (`--n-boot`). A resample draws, with replacement, as many
pools as the capture set has, from numpy's `default_rng(--seed)` started
afresh for each capture set. Captures of one pool share donors and loading,
so they are not independent; the pool is the unit. A resample whose pools
hold none of the denominator's cells is left out. There is no interval when
the denominator's cells sit in fewer than five pools (`--min-pools`) or fewer
than 95 per cent of the resamples are left. `--interval-sets` names the sets
that get intervals.

### Other outputs

- `lost_cells.tsv`: per capture set and run, the composition of each tier
  (`singlet`, `recovered`, `doublet`, `unassigned`) by cell class and by UMI
  bin (`--umi-bins`, default under 1,000, 1,000 to 2,999, 3,000 to 9,999,
  10,000 or more). Each tier's column sums to 1.
- `donor_recall.tsv`: per capture set, capture and sample-sheet donor, each
  run's singlets, recovered cells and share of the capture's singlets, and
  the median share of the same donor in sister captures (same `sister_key`).
  A run is flagged when its share is under a fifth of its sisters' median or
  of another run's share, where that reference is 0.02 or more. A donor that a
  run had no genotype for (`genotyped_<run>`) gets "no genotype in the
  <run> run" in place of that run's flags.

## `calibrate`

The wrong-donor rate by sex in bins of a per-cell score, such as the Scorer's
margin.

```bash
python -I tools/demux_tools calibrate --cells scorer_cells.tsv --truth truth.tsv --donors donors.tsv \
    --score-column margin --bins 0,1,2.2,5,10,20,50,100,inf --out calibration.tsv
```

- `--cells` lists per-capture tables (`capture`, `path`) with a barcode
  (`cell`), the score, the donor each cell is scored as (`best_singlet`) and
  the run's call (`donor_id`). Calls in `--exclude` (default `doublet`) are
  left out; `--unassigned` names the call of a cell the run did not assign.
- Captures without a truth file, and captures where no cell is labelled a
  male donor, are left out. A truth file without the `--sex-column` column
  stops the command.
- Per bin: cells labelled a male donor and those with a sex call, the female
  calls among them and their share with a Wilson 95 per cent interval; the
  same for cells the run assigned; the same for cells labelled a female
  donor (male calls); the unassigned share; the median UMIs.
- A margin of 2.2 is about where a 0.9 probability cut falls when two donors
  compete, since log(0.9 / 0.1) = 2.2.

## `doublet-excess`

Per capture, the doublet calls of each run split by the tool that made them,
against the loading expectation.

```bash
python -I tools/demux_tools doublet-excess --cells cells.tsv.gz --run cons \
    --tools cons=vireo:cons_vireo,souporcell:cons_souporcell --ambient ambient.tsv --out excess.tsv
```

- Expected doublets: 0.8 per cent per 1,000 cells, times 1 minus the sum of
  the squared donor shares of the capture's singlets for doublets of two
  donors. A capture of 10,000 cells from two donors in equal shares expects
  8 per cent doublets, 4 per cent between donors.
- Excess: the reference run's doublet share minus the expectation between
  donors.
- Largest-donor involvement: how often the capture's largest donor (share p)
  is in a named doublet pair, over 2p(1 - p) / (1 - sum of squared shares).
  It needs 20 named pairs, and is empty when one donor holds every singlet.
- `--reference-run` must be one of the `--run` values. A missing ambient file
  stops the command.
- Per run R: `R_doublets`, `R_doublet`, per tool group `R_doublet_<group>`
  and `R_umi_ratio_<group>` (median), `R_largest_donor_in_pairs`,
  `R_named_pairs`. For the reference run also `expected_total`,
  `expected_between_donors`, `expected_between_donors_cells`, `excess`,
  `R_largest_donor_in_pairs_<group>`, and with `--ambient` the median vireo
  `ambient_other` of its singlets and souporcell's ambient estimate.
