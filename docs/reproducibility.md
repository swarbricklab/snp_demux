# Reproducibility

snp_demux is built to be reproducible: pinned containers, a frozen upstream
pileup, and demultiplexing that is deterministic to the bit for every **call**.
This note records the container pinning, what is and isn't bit-reproducible, and
how the container split was verified.

## Containers

Every rule runs in a pinned image (config `containers:`). The images that wrap
the Demuxafy toolchain are pinned to the **exact versions the retired
`Demuxafy.sif` shipped**, so results are unchanged from the original monolith:

| Image | Tool(s) | Pin / provenance |
|-------|---------|------------------|
| `quay.io/biocontainers/bcftools:1.13--h3a49de5_0` | bcftools (subset/reheader/query) | public biocontainer |
| `ghcr.io/swarbricklab/sc_pileup-popscle:da70fc78` | popscle (demuxlet, freemuxlet) | built from `statgen/popscle@da70fc78`, the commit `Demuxafy.sif` shipped; shared with sc_pileup |
| `ghcr.io/swarbricklab/snp_demux-vireo:0.5.7` | vireo (demux CLI) | conda env pinned to `Demuxafy.sif`'s exact vireo stack: vireoSNP 0.5.7, numpy 1.21.5, scipy 1.7.3, MKL 2021.4, threadpoolctl 3.1.0 |
| `swarbricklab/vireo_snp:0.5.6` | vireoSNP genotyping (`annotate.py`) + cell counts | the pre-split image for those rules (vireoSNP 0.5.6, pandas 1.3.0) |
| `cumulusprod/souporcell:2.5` | souporcell | published image |
| `ghcr.io/swarbricklab/snp_demux-r:r4.2` | all R/Rmd rules | conda env = union of the per-rule R deps |
| `ghcr.io/swarbricklab/snp_demux-demuxafy:63fb228` | Demuxafy consensus/summary/assign scripts | R + the Demuxafy scripts pinned to repo commit `63fb228` (the sif's commit) |

The image recipes live under [`containers/`](../containers). The per-rule conda
envs under [`workflow/envs/`](../workflow/envs) are retained and referenced as
commented `conda:` directives next to each `container:`, so the R/Python rules can
run via `--use-conda` on a platform without singularity.

## What is bit-reproducible

All of the **demultiplexing results** are deterministic and reproduce byte-for-byte
given the same frozen pileups + panel:

- per-engine outputs — `demuxlet.best`, `freemuxlet.clust1.samples.gz`,
  `ref_clust_pearson_correlations.tsv`, vireo `donor_ids.tsv`,
  souporcell `clusters.tsv`;
- the Demuxafy consensus (`combined_results.tsv`);
- all downstream assignment/stats tables (`sample_barcode_assignments.tsv`,
  `cell_assignment.tsv`, `all_cells_stats.tsv`, `all_singlets_stats.tsv`,
  `all_*_stats_full.tsv`).

souporcell seeds its RNG from a fixed default (`--seed` defaults to 4) and derives
each restart's seed from it, so its **clustering** is deterministic — `clusters.tsv`
reproduces exactly given the same inputs and thread count. Its per-cluster
**genotype inference** is not bit-stable, though (see below); that affects only
`cluster_genotypes.vcf`, never the cluster assignments or anything downstream.

## What is **not** bit-reproducible

vireo's **ambient-RNA proportion** estimate (`--callAmbientRNAs`) is **not**
bit-reproducible. It is computed with multiple threads (`--nproc`), and the
parallel floating-point reductions are not associative, so the per-barcode ambient
fractions vary in the ~3rd–4th decimal place **run to run with the same image**.
This affects only:

- `vireo/<capture>/prop_ambient.tsv`
- `vireo/<capture>/ambient_summary.csv`
- `vireo/<capture>/ambient_stats.csv`

It does **not** affect any donor assignment — the hard calls are argmax-stable, so
`donor_ids.tsv` and everything downstream stay byte-identical. Treat the ambient
figures as approximate QC, not as reproducible quantities. (This is a property of
vireo, not of the container split: it was equally non-reproducible in the original
`Demuxafy.sif` workflow.)

souporcell is **genotype-free** and its per-cluster inference is likewise not
bit-reproducible (same non-associative parallel float reductions). This affects:

- `souporcell/<capture>/ambient_rna.txt` — an ambient estimate, not a call;
- `souporcell/<capture>/cluster_genotypes.vcf` — the inferred genotypes of the
  clusters; the per-cluster genotype *values* wobble run to run.

It does **not** affect the cluster *assignments* (`clusters.tsv`), their named
annotation (`clusters_annotated.tsv`) or any downstream donor assignment: the
annotation step re-matches clusters to the **known** donor genotypes, which is
robust to this wobble. As with vireo, this is a property of souporcell, not of the
container split.

## How the container split was verified

The focused containers replaced the monolithic `Demuxafy.sif` with no change to
results. This was checked end-to-end on internal pooled captures by running the
workflow with the new images and diffing against the `Demuxafy.sif` baseline:
every per-engine output, the consensus, and every R-processed assignment/stats
table were **byte-identical**; the only difference was the ambient QC above, which
was confirmed to vary run-to-run with a fixed image (i.e. inherent, not a
regression).

## Reproducing a run

1. Freeze the upstream [sc_pileup](https://github.com/swarbricklab/sc_pileup)
   outputs and the per-pool `informative.vcf.gz` panel.
2. `snakemake check_inputs --workflow-profile profiles/workflow` — preflight the
   panels (see [`data_prep_requirements.md`](data_prep_requirements.md)).
3. Run with `--use-singularity` and the pinned `containers:` map. Because the
   images are tag-pinned (and the GHCR ones are content-addressed by digest), the
   same pileups + panel give the same calls.
