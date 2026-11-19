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
# Helpers

# Look up donor IDs in relabel_map (defined below), dropping the names that
# indexing a named vector attaches.
relabel_donors <- function(x) unname(relabel_map[x])

# Copy a vireo TSV (plain or gzipped) to out_path, rewriting only its header.
# prob_doublet, prob_singlet and prop_ambient carry donor IDs in their column
# names but not in their bodies, so the body is streamed through unchanged.
# relabel_fields maps the split header fields to the new ones.
relabel_header <- function(in_path, out_path, relabel_fields, gzip = FALSE) {
    header_line <- readLines(in_path, n = 1)
    header_fields <- strsplit(header_line, "\t")[[1]]
    new_header <- paste(relabel_fields(header_fields), collapse = "\t")

    input_con <- file(in_path, open = "rt")
    output_con <- if (gzip) gzfile(out_path, open = "wt") else file(out_path, open = "wt")
    on.exit({
        close(input_con)
        close(output_con)
    }, add = TRUE)
    writeLines(new_header, output_con)
    # Stream in chunks rather than line by line: these matrices have one row
    # per barcode, and per-line readLines() is slow on large captures.
    while (length(chunk <- readLines(input_con, n = 10000)) > 0) {
        writeLines(chunk[!startsWith(chunk, header_line)], output_con)  # drop original header
    }
    invisible(out_path)
}

# Render page 1 of a vireo PDF figure as PNG. vireo only draws some figures
# (e.g. the input-genotype distance plot needs unknown donors), and the vireo
# rule touches empty placeholders for the rest; those become empty PNGs.
convert_pdf_to_png <- function(src, dst) {
    if (file.size(src) == 0) {
        message("Skipping empty PDF: ", src)
        file.create(dst)
        return(invisible(dst))
    }
    message("Converting ", src, " to ", dst)
    pdf_convert(pdf = src, format = "png", pages = 1, dpi = 300, filenames = dst)
    invisible(dst)
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
    convert_pdf_to_png(est_dist_pdf, est_dist_png)
    convert_pdf_to_png(input_dist_pdf, input_dist_png)
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

# Relabel the donor IDs in the column names of the per-cell probability
# matrices. The bodies are numeric and are copied unchanged.

# prob_doublet: "cell" + one column per donor pair, "<donorA>,<donorB>".
message("Saving relabelled doublet probabilities to ", prob_doublet_out)
relabel_header(prob_doub_tsv, prob_doublet_out, gzip = TRUE, function(fields) {
    pairs <- strsplit(fields[-1], ",")
    c("cell", vapply(pairs, function(d) paste(relabel_donors(d[1:2]), collapse = ","), character(1)))
})

# prob_singlet: "cell" + one column per donor.
message("Saving relabelled singlet probabilities to ", prob_singlet_out)
relabel_header(prob_sing_tsv, prob_singlet_out, gzip = TRUE, function(fields) {
    c("cell", relabel_donors(fields[-1]))
})

# prop_ambient: "cell", one column per donor, then "logLik_ratio".
message("Saving relabelled ambient proportions to ", ambient_out)
relabel_header(ambient_tsv, ambient_out, function(fields) {
    n <- length(fields)
    c(fields[1], relabel_donors(fields[-c(1, n)]), fields[n])
})

# -------------------------------------------

message("Converting PDF figures to PNG format.")
convert_pdf_to_png(input_dist_pdf, input_dist_png)
convert_pdf_to_png(est_dist_pdf, est_dist_png)


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