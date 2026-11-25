# Input data-prep requirements

snp_demux consumes per-cell pileups (from [`sc_pileup`](https://github.com/swarbricklab/sc_pileup))
and a per-pool donor-genotype panel (`<panel_root>/<capture>/informative.vcf.gz`).
The demultiplexing tools (vireo / demuxlet / freemuxlet / souporcell) impose a few
requirements on the **BAM** and the **panel/sites VCF**. cellSNP-lite and souporcell
tolerate all of these; **popscle (dsc-pileup / demuxlet) and the BAM↔VCF contract do
not** — getting them wrong produces hard FATAL errors mid-run. They are properties
of the *inputs*, not workflow bugs, so the panel-building stage should guarantee
them.

## 1. The sites/panel VCF must span every autosome (chr1–chr22)

`popscle dsc-pileup` runs once per autosome and **fatals on any chromosome whose
(split) VCF has zero variants**:

```
[E:...cmdCramDscPileup] Cannot read any single variant from .../chrN.vcf.gz
```

A genome-wide panel never triggers this, but a region-subset one (e.g. a
downsampled test fixture) will. If you subset the panel by region, keep at least
one site on **every** autosome (reads may be absent on a chromosome — dsc-pileup
tolerates sites with zero reads, just not zero sites).

## 2. The VCF `##contig` order must match the BAM `@SQ` order

`popscle dsc-pileup` requires the VCF's contig ordering to match the BAM's:

```
[E:...cmdCramDscPileup] Your VCF/BCF files and SAM/BAM/CRAM files have different
ordering of chromosomes. SAM/BAM/CRAM file has chr19 before chr2, but VCF/BCF
file has chr19 after chr2
```

Cell Ranger's GRCh38 reference is **lexicographically** ordered
(`chr1, chr10, chr11, …, chr19, chr2, …`), whereas many panels are **karyotypic**
(`chr1, chr2, chr3, …`). Note `bcftools sort` re-orders to the header's `##contig`
order, so it can silently *introduce* a mismatch. Reheader the panel's `##contig`
block to the BAM's `@SQ` order, then sort to match:

```bash
# ##contig block in BAM @SQ order
samtools view -H "$BAM" | awk -F'\t' '/^@SQ/{n="";l="";
  for(i=1;i<=NF;i++){if($i~/^SN:/)n=substr($i,4); if($i~/^LN:/)l=substr($i,4)}
  print "##contig=<ID="n",length="l">"}' > contigs.txt
{ zcat panel.vcf.gz | grep '^##' | grep -v '^##contig'; cat contigs.txt;
  zcat panel.vcf.gz | grep -m1 '^#CHROM'; } > new_header.txt
bcftools reheader -h new_header.txt panel.vcf.gz | bcftools sort -Oz -o panel.fixed.vcf.gz
bcftools index panel.fixed.vcf.gz
```

## 3. demuxlet needs a populated `R2` INFO field

The demuxlet rule uses `--geno-error-coeff`/`--geno-error-offset`, which model the
per-site genotype error from the imputation **`R2`** (imputation accuracy) INFO
field. demuxlet **fatals on the first site that lacks an `R2` value**:

```
Cannot extract R2 (1 float value) from INFO field at chrN:POS. Cannot use --geno-error-coeff
```

Imputed sites carry `R2`, but **genotyped/`TYPED` sites carry `ER2` (empirical R2),
not `R2`** — so a panel mixing typed + imputed variants may have `R2` on only a
fraction of sites. Ensure every panel site has an `R2` value. Genotyped sites are
high-confidence, so `R2=1` is a reasonable fill:

```bash
zcat panel.vcf.gz \
  | awk -F'\t' 'BEGIN{OFS="\t"} /^#/{print;next}
                { if($8 !~ /(^|;)R2=/) $8=$8";R2=1"; print }' \
  | bgzip > panel.r2.vcf.gz
bcftools index panel.r2.vcf.gz
```

(Alternatively, drop `--geno-error-coeff`/`--geno-error-offset` to use demuxlet's
fixed `--geno-error` instead — but that changes the error model.)

---

All three were observed validating the sc_pileup→snp_demux chain on real pooled
captures; cellSNP-lite (vireo) and souporcell (vartrix) were unaffected in every
case. The upstream panel-building stage (e.g. `sc_pileup`'s informative-panel
builder, or `bcarc_chromium build_informative_panels`) is the right place to
enforce them.
