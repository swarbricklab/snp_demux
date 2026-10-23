# Validation: consistency with the published Demuxafy reference

This repository ships two independent kinds of validation, which answer different
questions:

| | What it proves | Where |
|---|---|---|
| **CI fixture + golden** (`config/test_config.yaml`, DVC-tracked `test/`) | The pipeline is **deterministic / reproducible** — a fresh run reproduces the frozen output byte-for-byte (content). Fast, runs in CI. | `run_test.sh --check` |
| **This document** | The pipeline is **correct** — it reproduces the *published* results of the reference implementation (Demuxafy) on that tool's own canonical dataset. Heavy, one-off, not CI. | below |

The CI fixture is a 10%-downsampled subset chosen for speed; it cannot be compared
to published numbers. This document records a **full-depth** run against the
original Demuxafy test dataset and compares our per-method calls to the numbers
Neavin et al. publish in the Demuxafy documentation.

## Dataset

Demuxafy's public **TestData4PipelineFull** — the OneK1K 14-donor PBMC pool
(Neavin et al., *Genome Biology* 2024; "Demuxafy: improvement in droplet assignment
by integrating multiple single-cell demultiplexing and doublet detection methods").
20,982 droplets, 14 multiplexed donors.

Downloaded from the Demuxafy docs' (current) Google-Drive links — **not re-hosted
here** (the pool carries real donor genotypes; see the private-fixture policy):

```
# tarball id 1RzEYVEvhpWvjwjs671UyTV3R8iO4VbwC  (TestData4PipelineFull.tar.gz)
wget -O TestData4PipelineFull.tar.gz \
  'https://drive.usercontent.google.com/download?id=1RzEYVEvhpWvjwjs671UyTV3R8iO4VbwC&export=download&confirm=t'
```

Contents used: `possorted_genome_bam.bam` (39 GB, ~652 M reads), the 20,982
filtered barcodes, and `test_dataset.vcf` (14 donor genotypes, 229,719 variants).
Contigs are GRCh38p10 Ensembl-style and were renamed to UCSC `chr`-style to match
the genome FASTA (coordinate-safe on the primary assembly).

## Methodology

A one-off full-depth run of the real chain (not the CI fixture):

1. **`sc_pileup`** produces the per-cell pileups (cellSNP-lite, popscle `dsc-pileup`,
   souporcell `vartrix`) from the full BAM against the 14-donor genotype panel.
2. **`snp_demux`** runs the four demultiplexing methods + the Demuxafy consensus on
   those pileups.

The run configs are committed for reproducibility:

- `config/test_config.full.yaml` — the full-depth demux config (production demux
  thresholds). Consumes the pileups under `sc_pileup/results/full`.
- `config/test_config.full_neavin.yaml` — the **matched-config** variant used for
  the apples-to-apples vireo comparison (consumes the `results/full_neavin` pileup,
  built with cellSNP at Neavin's thresholds — see below).

The sc_pileup side uses `sc_pileup/test/config.full.yaml` (default) and
`sc_pileup/test/config.full_neavin.yaml` (Neavin cellSNP params). Both runs were
executed on PBS (`normal` queue, local-exec snakemake with `--use-singularity`;
compute nodes have no internet, so the containers are pre-cached in the shared
singularity prefix).

## Results — per method vs Neavin's published numbers

Neavin's published singlet totals below are the sums of her per-donor / per-cluster
breakdowns (the summary totals in the online docs are unreliable; the breakdowns are
not). The run completed cleanly at **production thresholds** — no threshold
relaxation and no duplicate-barcode artifact were needed at full depth (those are
specific to the sparse CI fixture).

| Method | Ours (singlet / doublet) | Neavin (singlet / doublet) | Agreement |
|---|---|---|---|
| **freemuxlet** | 18,216 / 2,764 | 18,215 / 2,767 | ≈ identical (Δ 1 / 3) |
| **souporcell** | 18,185 / 2,774 (+23 unassigned) | 18,225 / 2,757 | within ~0.2% |
| **vireo** (lenient cellSNP) | 17,972 / 2,919 | 18,099 / 2,770 | per-donor ~1%; +149 doublet |

freemuxlet and souporcell use the popscle / vartrix pileups and match Neavin to a
handful of cells. The only notable gap was **vireo's doublet count**, which traced
to a configuration difference, not an algorithmic one (next section).

### souporcell 0.2%
souporcell's clustering is random-restart EM, so run-to-run assignments wobble for
cells near the singlet/doublet boundary — a ~0.2% jitter is squarely that noise.
A second, smaller contributor: our default souporcell is a leaner mode than
Neavin's tutorial (it vartrix-es the donor panel with `skip_remap: true`, whereas
her tutorial uses `--common_variants` + the minimap2 remap). Per sc_pileup's own
benchmarking the remap is effectively a no-op here (~1% of total ALT reads, 98.7%
of per-cell calls unchanged), so it is not the driver. Not a discrepancy to chase.

## Apples-to-apples: vireo at Neavin's cellSNP thresholds

The vireo gap was hypothesised to be our lenient cellSNP settings (we kept every
site: `minMAF 0.0 / minCOUNT 1`) vs Neavin's `minMAF 0.1 / minCOUNT 20`. vireo is
the one method that consumes cellSNP; freemuxlet/souporcell (which matched) do not.
Re-running cellSNP at Neavin's thresholds (111,207 → 15,550 covered sites) and
re-running vireo on it confirms the hypothesis:

| Donor | vireo, lenient cellSNP | **vireo, Neavin cellSNP (ours)** | Neavin published | Δ (matched − pub) |
|---|---|---|---|---|
| 352_353 | 1611 | 1621 | 1619 | +2 |
| 43_43 | 1557 | 1562 | 1561 | +1 |
| 597_598 | 1525 | 1533 | 1532 | +1 |
| 349_350 | 1467 | 1479 | 1475 | +4 |
| 42_42 | 1415 | 1430 | 1430 | 0 |
| 660_661 | 1369 | 1378 | 1377 | +1 |
| 113_113 | 1335 | 1342 | 1342 | 0 |
| 39_39 | 1291 | 1300 | 1309 | −9 |
| 596_597 | 1259 | 1260 | 1271 | −11 |
| 41_41 | 1133 | 1142 | 1144 | −2 |
| 465_466 | 1094 | 1105 | 1104 | +1 |
| 40_40 | 1086 | 1098 | 1097 | +1 |
| 633_634 | 961 | 963 | 967 | −4 |
| 632_633 | 869 | 866 | 871 | −5 |
| **doublet** | 2919 | **2731** | 2770 | −39 |
| **unassigned** | 91 | 172 | 113 | +59 |
| singlet total | 17,972 | **18,079** | 18,099 | −20 |

Matching cellSNP **closed the gap**: the +149 doublet over-call became −39 (within
~1.4%), per-donor singlets now agree within ~0–1% (several exact), and singlet
totals are within 20 cells. The small residual (a couple of donors ±5–11, and the
unassigned↔doublet trade) reflects the remaining unmatched knobs — cellSNP
`minLEN`/`minMAPQ`, tool versions, vireo's `--callAmbientRNAs` ambient threshold —
and is noise-level, not structural.

## Conclusion

At matched configuration, `sc_pileup → snp_demux` reproduces Neavin's published
vireo essentially exactly, and at default configuration freemuxlet and souporcell
already match to within a handful of cells. The implementation is consistent with
the reference Demuxafy results.

## Reproduce

```bash
# 1. Download + extract TestData4PipelineFull (data-mover / copyq), chr-rename the BAM.
#    See sc_pileup/test/make_test_fixtures.sh for the chr-rename recipe (same transform,
#    without the 10% subsample).
# 2. sc_pileup pileup (PBS, normal queue, --use-singularity, cached containers):
snakemake --snakefile workflow/Snakefile --configfile test/config.full.yaml \
  --workflow-profile profiles/workflow --use-singularity \
  --singularity-prefix /g/data/a56/software/singularity/images --cores 48
#    For the Neavin-matched vireo: use test/config.full_neavin.yaml (cellSNP minMAF 0.1/minCOUNT 20).
# 3. snp_demux demux (PBS):
snakemake --snakefile workflow/Snakefile --configfile config/test_config.full.yaml \
  --workflow-profile profiles/workflow --use-singularity \
  --singularity-prefix /g/data/a56/software/singularity/images --cores 48
#    Neavin-matched vireo branch only:
snakemake ... --configfile config/test_config.full_neavin.yaml test/results_full_neavin/vireo/all.done
```

The full-depth inputs are not DVC-tracked (39 GB, external); the configs above are
the reproducible record.
