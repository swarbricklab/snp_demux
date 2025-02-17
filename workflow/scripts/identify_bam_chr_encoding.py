# ===================
# title: "annotation"
# description: Idnetifies the Chr encoding used in bam files
# project: "workflow-demuxafy"
# author: "Daniel Roden"
# email: "d.roden@garvan.org.au"
# company: "Garvan Institute"
# created: 30/05/2023
# -------------------
# NOTE: this is a snakmake script that uses the "snakemake" object internals
# ===================

import logging as log
import pysam

log.basicConfig(filename = snakemake.log[0], 
                    level = log.INFO, 
                    format = '%(asctime)s:%(levelname)s:%(message)s')                    
log.info("LOG FILE: " + snakemake.log[0])

log.info("----------------------------------") 
log.info("CONFIG")

log.info("INPUT")

# ----------
# DEBUG
#file_bam = "resources/testing/input/int_chr_encoding/cellranger/Atlas_Pool_12a/sample_alignments.bam"
#file_chr_encoding = "temp/Atlas_Pool_12a.chr_encoding.txt"
#file_bam = "resources/testing/input/souporcell/BrCa_Pilot2/input/cellranger/Mixed_1/sample_alignments.bam"
#file_chr_encoding = "temp/Mixed_1.chr_encoding.txt"
# ----------

file_bam = snakemake.input["bam"]
print("DEBUG: file_bam: ")
print(file_bam)
log.info("file_bam: ".join(file_bam))
log.info("")

log.info("OUTPUT")
file_chr_encoding = snakemake.output["chr_encoding"]
log.info("file_chr_encoding: " + file_chr_encoding)
log.info("")

log.info("WILDCARDS")
capture_id = snakemake.wildcards["capture_id"]
log.info("capture_id: " + capture_id)
log.info("")

log.info("PARAMS")
log.info("")

log.info("----------------------------------") 

log.info("get the bam header ")
bam_header = pysam.head("-h", "2", file_bam)

log.info("is the encoding chr or (numeric/non-chr encoded prefix)?")
# TODO: function for parsing bam header
bam_header = bam_header.splitlines()[1]
chr = bam_header.split(sep="\t")[1]
chr = chr.split(":")[1]
chr_encoding = "int_encoded"
if(chr.startswith("chr")):
    chr_encoding = "chr_encoded"

log.info("write to file")
file_chr_encoding = open(file_chr_encoding, "wt")
n = file_chr_encoding.write(chr_encoding)
file_chr_encoding.close()

log.info("")
log.info("PROCESSING COMPLETE")
