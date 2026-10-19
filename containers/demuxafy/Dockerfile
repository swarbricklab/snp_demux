# Slimmed replacement for the monolithic Demuxafy.sif (~6.9 GB), carrying ONLY the
# Demuxafy consensus/summary/assign scripts this workflow invokes plus their R
# dependencies:
#   * Combine_Results.R         (the demuxafy rule)         -> consensus
#   * Assign_Indiv_by_Geno.R    (the freemuxlet_assign rule) -> name freemuxlet clusters
#   * Demuxlet_summary.sh / Freemuxlet_summary.sh (summaries)
# The demux ENGINES live in their own focused images (popscle for demuxlet/
# freemuxlet, vireo, souporcell, bcftools) — see config `containers:`. This image
# has NO popscle: freemuxlet clustering runs in the popscle image, and only the
# genotype-naming step (Assign_Indiv_by_Geno.R) runs here.
# Scripts are pinned to the exact commit the retired sif shipped.
# Demuxafy is MIT-licensed (Drew Neavin, Powell Lab):
#   https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs
# Built via absconda --dockerfile -> GHCR.
FROM mambaorg/micromamba:1.5.8
ARG DDD_COMMIT=63fb228ab9db437386e70bdb2798edcef309445a
USER root

# R + the union of Combine_Results.R and Assign_Indiv_by_Geno.R dependencies.
RUN micromamba install -y -n base -c conda-forge -c bioconda \
        r-base=4.2 r-data.table r-r.utils r-argparse r-tidyverse r-future.apply \
        r-complexupset r-rcolorbrewer r-vcfr r-lsa bioconductor-complexheatmap \
        git \
 && micromamba clean -afy
ENV PATH=/opt/conda/bin:$PATH

# Demuxafy scripts pinned to the exact commit the retired sif shipped.
RUN git clone https://github.com/drneavin/Demultiplexing_Doublet_Detecting_Docs.git /opt/ddd \
 && git -C /opt/ddd checkout "${DDD_COMMIT}" \
 && chmod +x /opt/ddd/scripts/*.R /opt/ddd/scripts/*.sh
ENV PATH=/opt/ddd/scripts:$PATH
