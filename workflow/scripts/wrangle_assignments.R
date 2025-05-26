library(tibble)
library(dplyr)
library(stringr)
library(readr)
library(VariantAnnotation)

# Get input, output, log and wildcard from snakemake
file_sample_map                <- snakemake@input[["sample_map"]]
file_vcf                       <- snakemake@input[["vcf"]]
file_demuxafy_assignments      <- snakemake@input[["demuxafy_assignments"]]
combined_edited                <- snakemake@output[["sample_barcode_assignments"]]
summary_assignments            <- snakemake@output[["cell_assignment"]]
file_log                       <- snakemake@log[[1]]
capture                        <- snakemake@wildcards[["capture_id"]]
assignment_column              <- snakemake@params[["assignment_column"]]
droplettype_column             <- snakemake@params[["droplettype_column"]]

# Set up logging
log   <- file(snakemake@log[[1]], open="wt")
sink(log); sink(log, type="message")

# Load inputs
df_map     <- read_csv(file_sample_map)
df_assign  <- read_tsv(file_demuxafy_assignments)
vcf_header <- scanVcfHeader(file_vcf)

# Identify samples lacking genotype by comparing sample sheet and VCF header
donors_in_capture <- df_map %>%
  dplyr::filter(capture_id == capture) %>%
  dplyr::pull(donor_id)
print("Donors in capture: ")
print(paste(donors_in_capture, sep = " "))

known_genotype_donors <- samples(vcf_header)
print("Donors with known genotypes: ")
print(paste(known_genotype_donors, sep = " "))

donors_without_SNP <- setdiff(donors_in_capture, known_genotype_donors)
print("Donors without known genotypes")
print(paste(donors_without_SNP, sep = " "))

# Defensive check: only one sample should lack genotypes per capture
if(length(donors_without_SNP) > 1) {
  stop(
    paste("Error: multiple samples without genotype for", capture,":", donors_without_SNP)
  )
}

# If there is a donor without a known genotype, overwrite the assignment with the souporcell+vireosnp call
if(length(donors_without_SNP) > 0) {
  message("Found donor/sample without known genotype.")
  message("Overriding assignments with souporcell+vireosnp calls")
  glimpse(df_assign)
  df_assign <- df_assign %>%
    # Preserve original columns via rename()
    dplyr::mutate(
      original_assignment   = !!assignment_column,
      original_droplettype  = !!droplettype_column
    ) 
    # Override for samples without SNP genotypes
    df_assign <- df_assign %>% dplyr::mutate(
      !!assignment_column := if_else(
        (Souporcell_vireosnp_Individual_Assignment %in% donors_without_SNP) & (!!droplettype_column == "singlet"),
        Souporcell_vireosnp_Individual_Assignment,
        !!sym(assignment_column)
      ),
      !!droplettype_column  := if_else(
        Souporcell_vireosnp_DropletType %in% donors_without_SNP,
        Souporcell_vireosnp_DropletType,
        !!sym(droplettype_column)
      )
    )
}

# Build donor→sample lookup
map_vec <- df_map %>%
  dplyr::filter(capture_id == capture) %>%
  dplyr::select(donor_id, sample_id) %>%
  distinct() %>%
  deframe()
print(map_vec)

# Recode *_Individual_Assignment columns from donors to samples
# Prefix barcode with capture id
df_assign <- df_assign %>%
  dplyr::mutate(
    across(
      ends_with("_Individual_Assignment"),
      ~ recode(.x, !!!map_vec, .default = .x)
    )
  ) %>%
  dplyr::mutate(
    barcode = paste0(capture,"_",Barcode)
  ) %>%
  dplyr::select(-Barcode) %>%
  dplyr::select(
    barcode, everything()
  )

message("Creating assignment summary")
df_summary <- df_assign %>%
  dplyr::select(
    barcode,
    assignment = !!sym(assignment_column),
    status     = !!sym(droplettype_column)
  ) %>%
  dplyr::mutate(
    status = if_else(
        str_detect(assignment, "singlet"), 
        "unassigned", 
        status
    )
  )

message("Writing sample barcode assignments")
write_tsv(df_assign, combined_edited)
message("Writing summary assignments")
write_tsv(df_summary, summary_assignments)
