# SNP Demux

This workflow demultiplexes "pooled captures".
A "pooled capture" is when multiple tissues samples are dissociated into single cells together in the same test tube, and then captured together using the 10X Chromium platform to generate scRNA-seq libraries.

TODO: insert figure

When these libraries are sequenced, it is not obvious which cell barcode came from which sample.
But by utilizing SNP profiles for each sample in the pool, we can attempt to demultiplex the capture by assigning each barcode to an individual samples.
We can also use SNP information to detect "heterotypic doublets", where a barcode represents a doublet of two cells, each from a different sample.

## Contributions

This workflow is a simplification and refactoring of the [demuxafy](https://github.com/swarbricklab/demuxafy) workflow developed by [Daniel Roden](https://github.com/dlroden) and [Beata Kiedik](https://github.com/BeataKiedik). 
The demuxafy workflow is, in turn, a Snakemake wrapper around the [Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/index.html) package and tools developed by [Drew Neavin](https://github.com/drneavin) at the Powell Lab.

The main difference between this workflow and the [demuxafy](https://github.com/swarbricklab/demuxafy) workflow is that the [genotyping](https://github.com/swarbricklab/genotyping) and [imputation](https://github.com/swarbricklab/imputation) phases have been factored out into separate workflows, which are prerequisites for this workflow. 

## Method

Briefly, this workflow runs four well-established SNP demux methods in parralel and then uses [Demuxafy](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/index.html) to call doublets and find a consensus sample-barcode assignment.
The four methods are
1. vireo [docs](https://vireosnp.readthedocs.io/en/latest/)| [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1865-2) | [GitHub](https://github.com/single-cell-genetics/vireo)
2. souporcell [paper](https://www.nature.com/articles/s41592-020-0820-1) | [GitHub](https://github.com/wheaton5/souporcell)
3. demuxlet [docs](https://github.com/statgen/demuxlet/blob/master/tutorial/README.MD) | [paper](https://www.nature.com/articles/nbt.4042) | [GitHub](https://github.com/statgen/demuxlet)
4. freemuxlet [GitHub](https://github.com/statgen/popscle)

See the links above for details on the algorithms used by each of these methods.

The results of each method are combined using Demuxafy's [AtLeastHalfSinglet](https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/CombineResults.html) option, but there is scope to extend this workflow to include other options for combining results.

At a high level of abstraction, the overall processing flow is as follows:

TODO: insert figure

## Requirements

This workflow requires
1. `cellranger` "outs" for each capture (barcode and bam files)
2. SNP profiles for each donor (*) as a VCF file

### `cellranger` outs

This workflow requires the following `cellranger` "outs":
1. `barcodes.tsv.gz`
2. `possorted_genome_bam.bam`

In particular, the `barcodes.tsv.gz` file under `flitered_bc_matrix` should be preferred over the `raw_bc_matrix`, because `cellranger` filters out empty droplets and other low-quality barcodes.
Excluding these barcodes should improve the performance and accuracy of this workflow.

> Theoretically, it may be possible to further improve performance and accuracy by further filtering the barcode list, based on criteria such as gene counts, mitochondrial content or other doublet detection methods.
> However this theoretical benefit has not been benchmarked.

### VCF file

This workflow has been tested with SNP profiles based on microarray data.
Raw microarray data (`.CEL` files) for each donor sample can be aggregated into an appropriate VCF file via the [genotyping](https://github.com/swarbricklab/genotyping) workflow. 
The resulting VCF can be used directly as an input to this workflow.
Alternatively, additional SNPs can be imputed based on the 1000 Genomes project via the [imputation](https://github.com/swarbricklab/imputation) workflow.
Imputation provides more SNPs and therefore a better chance of finding matching reads in sparse scRNA-seq data, which [apparently](https://demultiplexing-doublet-detecting-docs.readthedocs.io/_/downloads/en/v0.0.3/pdf/) yields better results but we have not benchmarked this.


> (*) Note that SNP profiles do not have to be generated from the same tissue as the sample being demultiplexed.
> We get good results with other tissues from the same donor, including blood samples.

Our experience so far is with SNP profiles based on SNP microarrays, but theoretically VCF files derived from whole genome sequencing should also work, if appropriately filtered to exclude non-SNP variants.

> Note: make sure that the chromosome notation (`chr1`, `chr2`, .. versus `1`, `,2`, ...) in the VCF file matches that used by the reference transcriptome specified with `cellranger`.

## Processing

This module is defined as [Snakemake](https://snakemake.readthedocs.io/en/v7.32.3/) workflow.
The rule graph is as follows:

TODO: insert rule graph
