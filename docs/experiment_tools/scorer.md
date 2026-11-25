# Scorer: vireo's fixed-mode calls as arithmetic

The Scorer reproduces vireo's per-cell calls from a cellsnp-lite pileup and a
donor genotype VCF, without a vireo run, and writes the evidence behind each
call. It covers vireo's fixed mode only: as many donors as the VCF genotypes,
genotypes taken as given. It needs numpy, pandas and scipy, and vireoSNP
0.5.7 or 0.5.9 for its fit.

## What vireo computes in fixed mode

vireo first fits three allele ratios per capture: the expected share of ALT
UMIs at a 0/0, a 0/1 and a 1/1 site. It then scores each cell against each
donor by adding, over the sites where the cell has reads, the expected
log-probability of the cell's ALT and REF UMIs under the donor's genotype
there. A doublet is scored the same way, with the two donors' genotypes mixed
50:50. The scores, a prior that puts cells / 100,000 (at most 0.5) on
doublets, and three thresholds give the call:

1. `unassigned` when the best donor's probability (`prob_max`) is below 0.9;
2. `doublet` when the best pair's probability is at least 0.9;
3. `unassigned` when fewer than 10 sites have reads (`n_vars`).

A later rule overrides an earlier one. A blank genotype (`./.`) counts as a
third each for 0/0, 0/1 and 1/1.

The Scorer runs vireoSNP's own fitting code, called with the arguments vireo
uses in fixed mode, then recomputes the singlet and doublet score matrices
with the same floating-point operations as vireo. So it gives vireo's calls
and, for each cell, how far apart the candidates are.

## Commands

### `score`

```bash
python -I tools/demux_tools score --cellsnp PILEUP_DIR --vcf DONORS.vcf.gz --out PREFIX
```

| Option | Meaning |
|---|---|
| `--cellsnp DIR` | A cellsnp-lite folder |
| `--base`, `--samples`, `--ad`, `--dp` | Its files one by one; they override the folder's, so files named by their md5 can be read in place |
| `--vcf F` | The donor VCF, plain or gzip. Only records at piled sites are parsed, so a whole panel can be given |
| `--donors a,b` or `--donors-file F` | Donors to use (default: every sample of the VCF). The result is vireo's on a VCF that holds only these donors, not vireo's when given the full VCF and a smaller `--nDonor` (another mode) |
| `--donor-order` | `vcf` (the VCF header's order, which vireo uses; default), `list`, `sorted`, or `prob-header:PATH` to take the order from a vireo `prob_singlet.tsv.gz`. The workflow's `subset_vcf` rule puts the donors in the sample sheet's order, so reproducing a workflow run needs `list` (with the sheet's order) or `prob-header:` |
| `--geno-tag` | `GT` (default), `GP` or `PL` |
| `--barcodes F` | Keep only these barcodes (pileup order is kept). The allele ratios and the doublet prior are then fitted on these cells, so the result is vireo's on a pileup that holds only them, not vireo's calls for the full capture |
| `--n-donor N` | Stop unless N equals the donors used. vireo runs another mode when `--nDonor` differs from the donors in the file |
| `--loader` | `native` (default) or `vireo` (vireoSNP's readers; slow, for checking) |
| `--cache-dir D` | Keep the loaded arrays and reuse them when the inputs and options are unchanged (see below) |
| `--doublet-prior` | `auto` (vireo's) or a number. With a number, the calls come from the recomputed matrices |
| `--prob-max`, `--prob-doublet`, `--min-vars` | The three thresholds (vireo: 0.9, 0.9, 10) |
| `--vireo-format` | `0.5.7` (default), `0.5.9` or `0.5.10`: which vireo's output to mimic (see below) |
| `--vireo-out D` | Also write `donor_ids.tsv`, `summary.tsv`, `prob_singlet.tsv.gz` and `prob_doublet.tsv.gz` to this folder, as vireo writes them. A folder that already holds a `donor_ids.tsv` is refused unless `--force` is given |

It writes two files:

- `PREFIX.cells.tsv.gz`, one row per cell (columns below);
- `PREFIX.fit.json`, with the donors, the sites piled, matched and used, the
  sites dropped by reason, the allele ratios (`theta_mean`, `theta_sum`) at
  full precision, the doublet prior, vireo's iterations and lower bound, the
  internal check, the call counts, timings, whether the cache was hit, peak
  memory and the provenance record.

### `score-many`

```bash
python -I tools/demux_tools score-many --manifest captures.tsv --out-dir OUT --jobs 4 --cache-dir CACHE
```

The manifest is tab-separated with a header:

| Column | Holds |
|---|---|
| `capture` | A name, unique in the manifest; it becomes a folder name, so it holds no "/" and is not "." or ".." |
| `cellsnp`, or `base`, `samples`, `ad` and `dp` | The pileup |
| `vcf` | The donor VCF |
| `donors`, `donor_order`, `barcodes`, `geno_tag` | Optional, as for `score` |
| `vireo_dir` | Optional: vireo's output for this capture, to compare with |

Captures run largest DP matrix first, one worker process each, with BLAS
limited to one thread. Each writes `OUT/<capture>/base.cells.tsv.gz` and
`base.fit.json`. `OUT/summary.tsv` has one row per capture (status, cells,
donors, sites, allele ratios, call counts, load, fit and score seconds, cache
hit, peak memory, internal check). `OUT/validation.tsv` has one row per
capture with a `vireo_dir` (the columns of `validate-vireo`). The command
exits 1 when a capture fails or differs from its `vireo_dir`; the tables are
written first.

### `validate-vireo`

```bash
python -I tools/demux_tools validate-vireo --scored PREFIX --vireo-dir VIREO_OUTPUT [--files SCORER_VIREO_OUT]
```

Compares the barcodes in order and vireo's eight columns as text, and prints
the largest numeric difference in `prob_max`, `prob_doublet` and
`doublet_logLikRatio`. `--files` also compares a folder written by
`score --vireo-out` with vireo's four files, as decompressed text. It exits 0
when every cell is identical. With `--tolerance X` it also exits 0 when the
text columns are identical and every numeric difference is at most X.

### `site-terms`

```bash
python -I tools/demux_tools site-terms --cellsnp DIR --vcf F --cell BARCODE [--fit PREFIX.fit.json] --out terms.tsv
```

One row per site where the cell has reads: the site, its ALT and REF UMIs,
each donor's genotype (`gt__<donor>`) and the cell's log-likelihood term at
that site under each donor (`term__<donor>`). Each donor's terms add up to the
cell's score for that donor. `--fit` takes the allele ratios from an earlier
`score` run instead of refitting.

## The cells table

The first eight columns are vireo's `donor_ids.tsv`, with vireo's text.

| Column | Format | Meaning |
|---|---|---|
| `cell` | | Barcode |
| `donor_id` | | The call, by vireo's three rules |
| `prob_max` | `%.2e` | Probability of the best donor |
| `prob_doublet` | `%.2e` | Probability of the best pair (`prob_max_doublet` with `--vireo-format 0.5.10`) |
| `n_vars` | `%d` | Sites with reads among the sites used |
| `best_singlet` | | Donor with the highest probability |
| `best_doublet` | | The best pair as the chosen vireo version writes it (see below) |
| `doublet_logLikRatio` | `%.3f` | Best pair's score minus best donor's |
| `best_pair` | | The best pair, its two donors in string order |
| `second_donor` | | The runner-up donor |
| `umis` | `%d` | UMIs at the sites used |
| `ll_best`, `ll_second` | `%.4f` | Score (log-likelihood) under the best and runner-up donor |
| `margin` | `%.4f` | `ll_best` minus `ll_second` |
| `ll_best_pair` | `%.4f` | Score under the best pair |
| `doublet_margin` | `%.4f` | `ll_best` minus `ll_best_pair`: negative when a doublet fits better |
| `fit_per_umi` | `%.5f` | `ll_best` divided by `umis`: how well even the best donor fits |

A margin of 2.2 is about where vireo's 0.9 threshold falls when two donors
compete, since log(0.9 / 0.1) = 2.2.

vireoSNP 0.5.7 names the pair of a doublet wrongly. It scores pairs in
`itertools.combinations` order and names them in `itertools.permutations`
order, so `best_doublet` names the right pair only when it includes the first
donor (see `vireo_files.py`). `--vireo-format 0.5.7` writes these names, so
the table matches the vireo image's output. `0.5.9` and `0.5.10` write the
pair itself. `best_pair` is always the pair that was scored.

## One cell, worked through

These numbers come from a synthetic capture, made and scored with:

```bash
python -I tools/demux_tools synth --out cap --donors 3 --sites 80 --cells 30 --depth 1:2 \
    --coverage 0.15 --seed 5 --theta 0.02,0.5,0.98
python -I tools/demux_tools score --cellsnp cap/cellsnp --vcf cap/donors.vcf.gz --out s
python -I tools/demux_tools site-terms --cellsnp cap/cellsnp --vcf cap/donors.vcf.gz \
    --cell GCAAGAGCGGCGTGGT-1 --fit s.fit.json --out terms.tsv
```

The fitted allele ratios are 0.027 at 0/0 sites, 0.49 at 0/1 and 0.97 at 1/1.
One UMI then costs:

| UMI | 0/0 site | 0/1 site | 1/1 site |
|---|---|---|---|
| ALT | 3.71 | 0.71 | 0.03 |
| REF | 0.03 | 0.68 | 3.75 |

Cell `GCAAGAGCGGCGTGGT-1` has 14 UMIs at 9 sites. Its best donor, d2, scores
-3.09 and the runner-up, d1, scores -5.72, so the margin is 2.63. The two
donors share a genotype at one site, which adds the same to both. The other
eight sites carry the margin:

| Site | ALT, REF UMIs | d2 | d1 | d2 gains |
|---|---|---|---|---|
| chr1_17073_A_C | 2, 0 | 1/1 | 0/1 | 1.36 |
| chr1_48041_T_C | 0, 2 | 0/0 | 0/1 | 1.32 |
| chr2_35245_A_G | 0, 2 | 0/0 | 0/1 | 1.32 |
| chr2_31467_T_A | 1, 0 | 1/1 | 0/1 | 0.68 |
| chr2_49091_G_C | 1, 0 | 1/1 | 0/1 | 0.68 |
| chr1_14486_C_A | 1, 0 | 0/1 | 1/1 | -0.68 |
| chr2_41190_C_G | 1, 0 | 0/1 | 1/1 | -0.68 |
| chr1_32080_G_A | 2, 0 | 0/1 | 1/1 | -1.36 |

vireo gives d2 a probability of 0.932, above its 0.9 threshold, but calls the
cell `unassigned` because only 9 sites have reads.

## Identity with vireo

- **Inputs.** Pileup rows are matched to the VCF's SNP records by
  `CHROM_POS_REF_ALT`. When no id matches, the Scorer retries with "chr"
  added to either side, as vireoSNP 0.5.9 does; vireoSNP 0.5.7 stops with an
  error there instead. When the pileup or the VCF repeats an id, vireo uses
  the copy that numpy's sort of the whole list puts first, and so does the
  Scorer: it then matches against every SNP record of the VCF. Genotypes are
  parsed with vireo's rules: `vcf.py` lists them. The matrices are built with
  the same structure as vireo's, so the sparse sums run in the same order.
- **Names.** vireo cuts donor and pair names to 100 characters in
  `best_singlet` and `best_doublet`; the Scorer cuts them the same way there,
  and keeps them whole in `best_pair`, `second_donor` and the file headers.
- **Fit.** `--backend vireo` calls vireoSNP's `vireo_wrap` as vireo's command
  line does in fixed mode: one start, no extra donors, doublets checked, the
  ambient step off, one process. Fixed mode needs no seed: the random start is
  replaced at the first iteration.
- **Columns.** Under backend vireo the eight vireo columns are vireo's own
  numbers. `fit.json` records the largest difference between them and the
  Scorer's recomputed probabilities (`check`). In the vireo image this
  difference is 0.
- **Environments.** Exact reproduction is checked in the vireo image
  (`snp_demux-vireo:0.5.7`). Elsewhere, a different BLAS or numpy can move a
  probability in its last bits, which can flip the last digit of the text in a
  rare cell.

The tests run vireoSNP's own `vireo` command in fixed mode on synthetic
captures and require `donor_ids.tsv`, `summary.tsv` and both probability
tables to be identical as text.

## Speed and the cache

The native loader reads the matrices with pandas' C reader and keeps only the
rows that match the VCF while it reads. It parses only the VCF records at
piled sites. With `--cache-dir`, the matched arrays are saved as one `.npz`
file. The key covers each input's real path, size and modification time and
the options that shape the arrays, so a replaced or touched input is read
again. File contents are not read for the key: a file rewritten with the
same size within the same second (Lustre keeps modification times to the
second) would be served from the cache. An entry that cannot be read counts
as absent, and a cache folder that cannot be written gives a warning while
the run carries on. vireo's fit then takes most of a warm run.

## Library use

```python
import sys
sys.dont_write_bytecode = True
sys.path.insert(0, "<checkout>/tools")
from demux_tools import scorer as S

capture = S.load_capture("pool", cellsnp="pileup/", vcf="donors.vcf.gz", cache_dir="cache/")
fit = S.fit(capture)                                  # vireo's fit
scores = S.score(capture, fit)                        # scores.table, scores.ll_singlet, ...
held = S.Fit(fit.theta_mean, fit.theta_sum, doublet_prior=fit.doublet_prior)
```

`loglik_singlet(AD, DP, GPb, theta_mean, theta_sum)` and `loglik_doublet(...)`
take any row subset of the sites (a subset of the rows of `AD`, `DP` and
`GPb`) and any genotype array. `Fit(theta_mean, theta_sum)` holds allele
ratios fixed, for example those of an earlier fit.
