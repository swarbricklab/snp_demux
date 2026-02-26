# relabel_vireo.R
# This script relabels Vireo output files using donor information from the samples CSV and writes relabelled outputs to the desired locations.
# Inputs and outputs are provided by Snakemake via snakemake@input and snakemake@output.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
  library(pdftools)
})  

# Snakemake input/output
samples_csv <- snakemake@input[["samples_csv"]]
vireo_tsv <- snakemake@input[["donor_ids"]]
prob_doub_tsv <- snakemake@input[["prob_doublet"]]
prob_sing_tsv <- snakemake@input[["prob_singlet"]]
ambient_tsv <- snakemake@input[["prop_ambient"]]
est_dist_pdf <- snakemake@input[["fig_est"]]
input_dist_pdf <- snakemake@input[["fig_input"]]
summary_file <- snakemake@input[["summary"]]

capture <- snakemake@wildcards[["capture_id"]]

# Outputs
vireo_out <- snakemake@output[["donor_ids"]]
sample_map <- snakemake@output[["sample_map"]]
prob_doublet_out <- snakemake@output[["prob_doublet"]]
prob_singlet_out <- snakemake@output[["prob_singlet"]]
ambient_out <- snakemake@output[["prop_ambient"]]
input_dist_png <- snakemake@output[["input_dist_png"]]
est_dist_png <- snakemake@output[["fig_est_png"]]
summary_out <- snakemake@output[["summary"]]

# Redirect messages and warnings to the Snakemake log file
log_file <- snakemake@log[[1]]
if (!is.null(log_file) && nzchar(log_file)) {
    log_con <- file(log_file, open = "wt")
    sink(log_con, type = "message")
    sink(log_con, type = "output")
    on.exit({
        sink(type = "message")
        sink(type = "output")
        close(log_con)
    }, add = TRUE)
}

# -------------------------------------------------

# Create a map between the donor IDs in the Vireo results and the actual donor IDs from the samples CSV file
# Use this map to relabel the Vireo results and save both the map and the relabelled outputs

# Read the samples CSV file
samples <- read_csv(samples_csv, col_types = cols(.default = "c"))
# Read the Vireo TSV file
vireo <- read_tsv(vireo_tsv, col_types = cols(.default = "c"))

# Construct correct donor-sample mapping 
all_donors <- samples %>%
    filter(capture_id == capture) %>%
    pull(donor_id) %>%
    sort()
message("All donors in the samples CSV: ", paste(all_donors, collapse = ", "))

# Construct list of unique donor IDs from the Vireo results
# This will include the placeholder 'donorX' for any missing donors
# as well as "doublet" and "unassigned"
vireo_donors <- vireo %>%
    select(donor_id) %>%
    distinct() %>%
    pull(donor_id) %>%
    sort()
message("Donors in Vireo results: ", paste(vireo_donors, collapse = ", "))

# Check for the presence of 'donorX' placeholders in the Vireo results
# If not present, it means all donors are accounted for and
# we can skip the relabeling step
if (!any(str_starts(vireo_donors, "donor"))) {
    message("No missing donors in Vireo results. Skipping relabeling.")
    write_tsv(vireo, vireo_out)
    # Create an empty sample map since no relabeling is needed
    file.create(sample_map)
    # Copy the other outputs without relabeling
    file.copy(prob_doub_tsv, prob_doublet_out)
    file.copy(prob_sing_tsv, prob_singlet_out)
    file.copy(ambient_tsv, ambient_out)
    # Convert PDFs to PNGs (handle empty placeholders)
    for (pair in list(
        list(src = est_dist_pdf, dst = est_dist_png),
        list(src = input_dist_pdf, dst = input_dist_png)
    )) {
        if (file.size(pair$src) == 0) {
            file.create(pair$dst)
        } else {
            pdf_convert(pair$src, format = "png", pages = 1, dpi = 300, filenames = pair$dst)
        }
    }
    file.copy(summary_file, summary_out)
    # Exit the script early
    quit("no", status = 0)
}

# Find which (if any) donor is not in the vireo results
missing_donors <- setdiff(all_donors, vireo_donors)
message("Missing donors in Vireo results: ", paste(missing_donors, collapse = ", "))

# There should be exactly one missing donor
if (length(missing_donors) > 1) {
  stop("More than one missing donor in Vireo results. Please check the input files.")
} else if (length(missing_donors) == 0) {
  stop("donorX found but no matches found. Please check the input files.")
}

# Construct map for relabeling vireo results
# Map the 'donorX' placeholder to the actual missing donor ID
message("Creating relabel map for Vireo results.")
relabel_map <- tibble(
        old = vireo_donors
    ) %>%
    mutate(
        new = if_else(
            str_starts(old, "donor"),
            missing_donors[[1]],
            old)
    ) %>% deframe()

# Save the relabel map to a file for use with bcftools
message("Saving relabel map to ", sample_map)
write_tsv(
    tibble(donor_id = names(relabel_map), sample_id = relabel_map),
    sample_map,
    col_names = FALSE
)

# First split the comma-separated 'best_doublet' column into two columns
# Then relabel the donor_id, best_singlet, and doublet columns using the relabel map
# before recombining 'best_doublet' back into a single column
message("Relabelling Vireo results.")
vireo_relabelled <- vireo %>%
    separate(
        best_doublet, 
        into = c("doublet_1", "doublet_2"), 
        sep = ",", 
        remove = FALSE
    ) %>% 
    mutate(
        across(
            c(donor_id, best_singlet, doublet_1, doublet_2),
            ~ relabel_map[.x]
        )
    ) %>%
    mutate(
        best_doublet = paste(doublet_1, doublet_2, sep = ",")
    ) %>%
    select(-doublet_1, -doublet_2)

# Save the relabelled Vireo results to a new TSV file
message("Saving relabelled Vireo results to ", vireo_out)
write_tsv(vireo_relabelled, vireo_out)  

# -------------------------------------------------

# Relabel the column names in the doublet probabilities TSV file, and save to a new file

# Read and relabel only the header of the doublet probabilities TSV file
message("Reading header from ", prob_doub_tsv)
header_line <- readLines(prob_doub_tsv, n = 1)
header_fields <- strsplit(header_line, "\t")[[1]]

# Relabel doublet columns (skip "cell")
columns <- sapply(header_fields[-1], function(x) {
    donors <- strsplit(x, ",")[[1]]
    paste(relabel_map[donors[1]], relabel_map[donors[2]], sep = ",")
})

new_header <-  paste(c("cell", columns), collapse = "\t")

# Write the relabelled header and the rest of the file to the output (gzipped)
message("Saving relabelled doublet probabilities to ", prob_doublet_out)
input_con <- file(prob_doub_tsv, open = "rt")
output_con <- gzfile(prob_doublet_out, open = "wt")
writeLines(new_header, output_con)
while (length(line <- readLines(input_con, n = 1)) > 0) {
    if (startsWith(line, header_line)) next  # skip original header
    writeLines(line, output_con)
}
close(input_con)
close(output_con)

# -------------------------------------------------

# Relabel donor column names in the singlet probabilities TSV file, and save to a new file

# Relabel donor column names in the singlet probabilities TSV file, and save to a new file

message("Relabelling singlet probabilities column names.")
header_line <- readLines(prob_sing_tsv, n = 1)
header_fields <- strsplit(header_line, "\t")[[1]]

# Relabel donor columns (skip "cell")
columns <- sapply(header_fields[-1], function(x) relabel_map[x])

new_header <- paste(c("cell", columns), collapse = "\t")

input_con <- file(prob_sing_tsv, open = "rt")
output_con <- gzfile(prob_singlet_out, open = "wt")
writeLines(new_header, output_con)
while (length(line <- readLines(input_con, n = 1)) > 0) {
    if (startsWith(line, header_line)) next  # skip original header
    writeLines(line, output_con)
}
close(input_con)
close(output_con)

# -------------------------------------------------

# Relabel column names in the ambient proportions TSV file, and save to a new file

message("Relabelling ambient proportions column names.")
header_line <- readLines(ambient_tsv, n = 1)
header_fields <- strsplit(header_line, "\t")[[1]]

# Relabel donor columns (skip "cell" and "logLik_ratio")
columns <- c(
    header_fields[1],  # "cell"
    sapply(header_fields[-c(1, length(header_fields))], function(x) relabel_map[x]),
    header_fields[length(header_fields)]  # "logLik_ratio"
)

new_header <- paste(columns, collapse = "\t")

input_con <- file(ambient_tsv, open = "rt")
output_con <- file(ambient_out, open = "wt")
writeLines(new_header, output_con)
while (length(line <- readLines(input_con, n = 1)) > 0) {
    if (startsWith(line, header_line)) next  # skip original header
    writeLines(line, output_con)
}
close(input_con)
close(output_con)

# -------------------------------------------

# Convert PDF files to PNG
message("Converting PDF figures to PNG format.")
images <- setNames(
    c(input_dist_png, est_dist_png),
    c(input_dist_pdf, est_dist_pdf)
)
for (source_image in names(images)) {
    target_image <- images[[source_image]]
    # Skip empty placeholder PDFs (vireo touches these when not produced)
    if (file.size(source_image) == 0) {
        message("Skipping empty PDF: ", source_image)
        file.create(target_image)
        next
    }
    message("Converting ", source_image, " to ", target_image)
    pdf_convert(
        pdf      = source_image,
        format   = "png",
        pages    = 1,
        dpi      = 300,
        filenames = target_image
    )
}


# -------------------------------------------

# Relabel the summary file and save to a new file
message("Relabelling summary file and saving to ", summary_out)
read_tsv(
    summary_file,
    col_types = cols(.default = "c")
) %>%
mutate(
    across(
        c(Var1), ~ relabel_map[.x]
    )
) %>%
write_tsv(summary_out)