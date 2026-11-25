#!/bin/bash
# Touch empty placeholder inputs so `snakemake -n` can build the full DAG in CI
# without any real data. Paths mirror the sc_pileup output layout + the panel.
set -eu
base=.tests/ci
for cap in Pool_A; do
  for plat in short; do
    mkdir -p $base/pileup/cellsnp/$plat/$cap $base/pileup/popscle/$plat/$cap $base/pileup/souporcell/$plat/$cap $base/panel/$cap
    touch $base/pileup/cellsnp/$plat/$cap/cellSNP.base.vcf.gz \
          $base/pileup/cellsnp/$plat/$cap/cellSNP.samples.tsv \
          $base/pileup/cellsnp/$plat/$cap/cellSNP.tag.AD.mtx \
          $base/pileup/cellsnp/$plat/$cap/cellSNP.tag.DP.mtx \
          $base/pileup/popscle/$plat/$cap/pileup.var.gz \
          $base/pileup/popscle/$plat/$cap/pileup.cel.gz \
          $base/pileup/popscle/$plat/$cap/pileup.plp.gz \
          $base/pileup/popscle/$plat/$cap/pileup.umi.gz \
          $base/pileup/souporcell/$plat/$cap/alt.mtx \
          $base/pileup/souporcell/$plat/$cap/ref.mtx \
          $base/pileup/souporcell/$plat/$cap/barcodes.tsv \
          $base/pileup/souporcell/$plat/$cap/variants.vcf \
          $base/panel/$cap/informative.vcf.gz \
          $base/panel/$cap/informative.vcf.gz.csi
  done
done
touch $base/barcodes.tsv $base/genome.fa
echo "placeholders touched"
