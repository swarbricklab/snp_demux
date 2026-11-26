# Per-capture roll-up consumed by both the capture page and the index.
#
# Split out of report_capture.Rmd because Snakemake's RMarkdown integration
# permits only a single output file, so the Rmd cannot emit this itself. Having
# it here is better anyway: the page and the index now read the same numbers
# rather than each deriving their own.

log_file <- file(snakemake@log[[1]], open = "wt")
sink(log_file); sink(log_file, type = "message")

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(readr)
})

capture  <- snakemake@wildcards[["capture_id"]]
pool     <- snakemake@params[["pool_id"]]
declared <- strsplit(snakemake@params[["declared_samples"]], ",")[[1]]
declared <- declared[nzchar(declared)]

assign_tbl <- read_tsv(snakemake@input[["assignments"]], show_col_types = FALSE)
cells      <- read_tsv(snakemake@input[["cell_assignment"]], show_col_types = FALSE)
amb_sum    <- read_csv(snakemake@input[["ambient_summary"]], show_col_types = FALSE)

# souporcell writes a bare sentence, e.g. "ambient RNA estimated as 2.011...%"
soup_amb <- suppressWarnings(as.numeric(sub(
  ".*?([0-9.]+)%.*", "\\1",
  paste(readLines(snakemake@input[["souporcell_ambient"]]), collapse = " ")
)))

methods   <- c("Demuxlet", "Freemuxlet", "Souporcell", "Vireo", "Souporcell_vireosnp")
consensus <- "Consensus"

# A call is the sample id for a singlet, or the droplet type otherwise, so that
# identity and doublet/unassigned are comparable on one axis.
as_call <- function(dtype, ind) {
  ifelse(dtype == "singlet" & !is.na(ind) & !ind %in% c("doublet", "unassigned"),
         ind, dtype)
}

long <- bind_rows(lapply(methods, function(m) tibble(
  Barcode = assign_tbl$Barcode,
  call    = as_call(assign_tbl[[paste0(m, "_DropletType")]],
                    assign_tbl[[paste0(m, "_Individual_Assignment")]])
)))
# a recovered barcode's call is its sample, like a singlet's
cons_dtype <- assign_tbl[[paste0(consensus, "_DropletType")]]
cons <- tibble(
  Barcode   = assign_tbl$Barcode,
  cons_call = as_call(ifelse(cons_dtype %in% "recovered", "singlet", cons_dtype),
                      assign_tbl[[paste0(consensus, "_Individual_Assignment")]])
)

mean_agree <- long %>%
  left_join(cons, by = "Barcode") %>%
  group_by(Barcode) %>%
  summarise(n_agree = sum(call == cons_call, na.rm = TRUE), .groups = "drop") %>%
  pull(n_agree) %>% mean(na.rm = TRUE)

n_cells   <- nrow(cells)
n_singlet <- sum(cells$status == "singlet", na.rm = TRUE)
n_doublet <- sum(cells$status == "doublet", na.rm = TRUE)
# "recovered": one method's lone singlet call (params.consensus.recover_from).
# Counted apart from singlets; n_recovered below is samples, not barcodes.
n_rescued <- sum(cells$status == "recovered", na.rm = TRUE)
recovered <- unique(cells$assignment[cells$status == "singlet"])
recovered <- recovered[!is.na(recovered)]

tibble(
  capture_id            = capture,
  pool_id               = pool,
  n_barcodes            = n_cells,
  n_singlet             = n_singlet,
  n_doublet             = n_doublet,
  n_recovered_barcodes  = n_rescued,
  n_unassigned          = n_cells - n_singlet - n_doublet - n_rescued,
  frac_singlet          = n_singlet / n_cells,
  n_declared            = length(declared),
  n_recovered           = length(intersect(declared, recovered)),
  missing_samples       = paste(setdiff(declared, recovered), collapse = ";"),
  ambient_souporcell    = soup_amb / 100,
  ambient_vireo_median  = median(amb_sum$ambient_other, na.rm = TRUE),
  mean_methods_agreeing = mean_agree
) %>% write_tsv(snakemake@output[["metrics"]])

message("wrote ", snakemake@output[["metrics"]])
