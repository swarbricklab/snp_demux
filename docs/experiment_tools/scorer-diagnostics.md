# Scorer diagnostics: changed sites, genotypes and allele ratios, and where a margin comes from

This layer extends the [Scorer](scorer.md). It scores a capture again with a
changed site list, changed genotypes or held allele ratios, without reloading
it, and adds columns that show which sites carry each cell's margin. It also
adds a fit that needs no vireoSNP. Every addition is off by default: without
the new options, `score` and `score-many` write what the Scorer writes, byte
for byte.

## The native fit

`--backend native` fits the allele ratios with `vireo_port.fit_fixed_native`
instead of vireoSNP. With the genotypes fixed, each sparse product vireo forms
in an iteration is a fixed (cells x donors) matrix times a number, so the six
products are formed once and each iteration works on (cells x donors)
matrices. The loop, its two passes, the delay of the ratio update and the
stopping rule are vireo's, and the lower bound includes vireo's constant
terms. It therefore stops at the iteration where vireo stops, unless a
floating-point difference moves a step of the lower bound across vireo's
0.01 threshold.

The ratios then differ from vireo's only by the order of floating-point sums.
The calls and probabilities come from the Scorer's own matrices (`values_from`
is `scorer` in `fit.json`). A probability can differ from vireo's in its last
bits, so the last printed digit of `prob_max` or `prob_doublet` can differ for
a rare cell.

## Changing the sites

| Option | Effect |
|---|---|
| `--sites F` | Score only the sites in this list |
| `--exclude-sites F` | Leave these sites out |
| `--thin-bp N` | Keep no two sites within N bp: `--thin-keep first` walks each chromosome left to right; `deepest` keeps the sites with most UMIs over cells first |

A site list is a table with a `variant` column (vireo's `CHROM_POS_REF_ALT`),
or `chrom` and `pos` (and `ref`, `alt`), or a VCF. The `sites` command makes
them:

```bash
python -I tools/demux_tools sites annotate --cellsnp PILEUP --gtf genes.gtf.gz --panel panel.vcf.gz --info R2,AF --out sites.tsv.gz
python -I tools/demux_tools sites thin --cellsnp PILEUP --window-bp 100 --keep deepest --out thinned.tsv
python -I tools/demux_tools sites from-panel --vcf panel.vcf.gz --min R2=0.8 --snps-only --out r2.tsv
```

`annotate` gives each site its region (`exonic` inside an exon of any gene,
`intronic` inside a gene but no exon, `intergenic`), its gene (the gene whose
exon holds it, else the overlapping gene that starts first), the number of
genes overlapping it, and INFO fields from a panel (a VCF with or without
sample columns). Chromosome names must match the sites' exactly ("chr1" is
not "1"); when no site's chromosome is in the GTF, it prints a warning. Its
table also serves as
`--site-groups`.

## Changing the genotypes

| Option | Effect |
|---|---|
| `--fill-blanks keep` | vireo's a third each for a blank genotype (default) |
| `--fill-blanks hwe:FIELD` | Hardy-Weinberg proportions ((1-q)^2, 2q(1-q), q^2) of the ALT frequency q in an INFO field of the donor VCF (loaded with `--info FIELD`) |
| `--fill-blanks hwe:FILE` | The same, with q from a table with columns `variant` and `af` |
| `--fill-blanks drop-site` | Drop every site where a donor is blank |
| `--harden` | Make soft genotypes hard at their most likely genotype; blanks stay blank |
| `--extra-vcf F --extra-donors a,b` | Add donors from another VCF, matched by site (decoys). `--extra-missing` sets what their blanks and absent sites become: `uniform`, `hwe:FIELD`, `hwe:FILE` or `drop-site` |

A site keeps its record as blank after a fill, so `@missing`, `blank_umis` and
the pairwise columns still see it as blank.

## Holding the allele ratios

`--theta refit` fits the ratios (default). `--theta hold:PREFIX.fit.json`
takes them from an earlier run. `--theta fixed:m0,m1,m2:s0,s1,s2` sets the
means (strictly between 0 and 1) and concentrations (positive) for 0/0,
0/1 and 1/1. A changed site list or
genotype array moves the ratios when they are refitted; holding them shows
the change in the evidence alone.

## Diagnostic columns

Each option appends its columns after the Scorer's seventeen.

| Option | Columns |
|---|---|
| `--extra-columns` | `alt_umis`; `n_sites_diff` and `umis_diff`, the sites with reads (and UMIs there) where the best and runner-up donors' genotypes differ; `fit_expected_per_umi`, the score per UMI a true singlet of the best donor would reach at the same sites and depths, to set beside `fit_per_umi` |
| `--group-by SPEC` | `margin__<group>__<label>`, `umis__<group>__<label>`, `nsites__<group>__<label>`: the margin between the best and runner-up donor summed over each group's sites. The groups' margins add up to `margin`. SPEC is a column of `--site-groups`, `@missing` (which donors are blank at the site; patterns beyond `--top-groups` are `other`), `@chrom`, or `@info:FIELD:e0,e1,...` (bins of an INFO field; the last bin is closed) |
| `--marginal` | `best_marginal`, `ll_best_marginal`, `ll_second_marginal`, `margin_marginal`: scores with the exact marginal at soft and blank genotypes (below); `blank_umis`, the cell's UMIs at sites where its best donor is blank |
| `--pairwise-shared` | `pairwise_best`, `pairwise_wins`, `pairwise_margin`, `pairwise_shared_umis`: each pair of donors compared only on the sites where both are genotyped (below) |
| `--top-sites N` | `top1_site_share`, `top10_sites_share` (the largest and the ten largest site contributions over the margin) and, with a gene column in `--site-groups`, `top_gene`, `top_gene_share`, `top10_genes_share` (sites with no gene are left out of these three). Writes `PREFIX.sites.tsv.gz` with each cell's N largest contributions: site, group labels, UMIs, both donors' genotypes and terms, their difference and the cumulative share |

A share is a contribution divided by the margin. It can exceed 1 when other
sites pull the other way, and it is empty when the margin is 0 (so is
`cum_share`).

`--pairwise-shared` keeps two (cells x donors x donors) arrays of 8-byte
numbers, so its memory grows with the square of the donors: about 2.4 GB for
30,000 cells and 70 donors.

### Mean-field against exact at a blank genotype

vireo scores a soft or blank genotype by averaging the log-likelihoods over
its genotype probabilities. That is a mean-field value. The exact marginal
averages the likelihoods instead. Both use vireo's expectations: an ALT UMI
under genotype g counts E[log theta_g] (a digamma difference of the fitted
Beta), so the exact marginal averages exp(E[log theta_g]), which is close to
theta_g when the fitted concentrations are large, as they are for a capture.
Take one ALT UMI at a site where the donor is blank (a third each), with
allele ratios 0.01, 0.5 and 0.99 and large concentrations:

| | Value | Cost of the UMI |
|---|---|---|
| Mean-field (vireo) | (log 0.01 + log 0.5 + log 0.99) / 3 | 1.77 |
| Exact | log((0.01 + 0.5 + 0.99) / 3) | 0.69 |

The exact value is never below the mean-field one, and the two are equal at
a hard genotype. With `--fill-blanks hwe:AF` at an ALT frequency of 0.05 the
genotype is (0.9025, 0.095, 0.0025), and vireo's mean-field cost is 0.09 for
a REF UMI and 4.22 for an ALT UMI.

### Pairwise scores on shared sites

For donors d and e, M(d, e) sums the cell's terms under d minus under e over
the sites where neither is blank. A blank genotype then cannot decide between
the two donors it affects. `pairwise_best` is the donor that beats the most
others (ties go to the higher total of M), `pairwise_margin` its smallest M,
and `pairwise_shared_umis` the cell's UMIs at the sites it shares with that
closest donor. Without blanks, `pairwise_best` is the best singlet and
`pairwise_margin` its margin.

## Many variants per load

```bash
python -I tools/demux_tools score-many --manifest captures.tsv --variants variants.tsv --out-dir OUT --cache-dir CACHE
```

`variants.tsv` has one row per variant. Each capture is loaded once and scored
under every variant, in the order given, writing
`OUT/<capture>/<variant>.cells.tsv.gz`, `.fit.json` and, with `top_sites`,
`.sites.tsv.gz`. `summary.tsv` has one row per capture and variant, with the
backend and the seconds spent on overrides, fit, score and diagnostics.

| Column | Holds |
|---|---|
| `variant` | A name, unique in the table; it becomes a file name, so it holds no "/" and is not "." or ".." |
| `sites`, `exclude_sites`, `thin_bp`, `thin_keep` | As the `score` options |
| `group_by` | Groupings separated by `;` |
| `fill_blanks`, `harden`, `extra_vcf`, `extra_donors`, `extra_missing` | As the `score` options |
| `theta` | `refit`, `base` (hold the first variant's ratios; the first variant must refit), `hold:FIT.json` or `fixed:...` |
| `backend` | `vireo` or `native` for this variant (default: the command's `--backend`) |
| `marginal`, `pairwise_shared`, `extra_columns` | `true` or `false` |
| `top_sites` | N, or empty |

A capture with a `vireo_dir` is checked against it on the variant named
`base`, when that variant changes nothing.

## Library use

```python
from demux_tools import scorer as S, diagnostics as D, sites

capture = S.load_capture("pool", cellsnp="pileup/", vcf="donors.vcf.gz", info=("AF",))
base = S.fit(capture)                                          # or backend="native"
filled, counts = S.fill_missing(capture, how="hwe", af=capture.info["AF"])
held = S.fit(filled, theta=(base.theta_mean, base.theta_sum))
scores = S.score(filled, held)
best, second = D.best_and_second(scores)
margins = D.group_margins(filled, held, best, second, sites.missing_pattern(filled.missing, filled.donors), "missing")
```

`S.apply_variant(capture, {...})` runs one variant as `score` does, and
`S.run_variants(...)` many on one load.
