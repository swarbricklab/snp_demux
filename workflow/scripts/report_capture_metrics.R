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
# Keep only tools actually present in the combined table: with a configurable
# params.tools subset (e.g. [souporcell, vireo]), the disabled tools' *_DropletType
# columns are absent, and referencing them yields a length-0 `call` column that
# breaks the tibble() below.
methods   <- methods[paste0(methods, "_DropletType") %in% names(assign_tbl)]
consensus <- "AtLeastHalfSinglet"

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
cons <- tibble(
  Barcode   = assign_tbl$Barcode,
  cons_call = as_call(assign_tbl[[paste0(consensus, "_DropletType")]],
                      assign_tbl[[paste0(consensus, "_Individual_Assignment")]])
)

mean_agree <- long %>%
  left_join(cons, by = "Barcode") %>%
  group_by(Barcode) %>%
  summarise(n_agree = sum(call == cons_call, na.rm = TRUE), .groups = "drop") %>%
  pull(n_agree) %>% mean(na.rm = TRUE)

# Diagnostics roll-up (design doc §D gates). Robust to empty/absent files: an
# empty-but-valid TSV yields NA silhouette / 0 absorbed donors rather than error.
read_diag <- function(key) {
  path <- snakemake@input[[key]]
  if (is.null(path) || !file.exists(path) || file.info(path)$size == 0)
    return(NULL)
  out <- tryCatch(read_tsv(path, show_col_types = FALSE), error = function(e) NULL)
  if (is.null(out) || nrow(out) == 0) NULL else out
}

emb_metrics <- read_diag("embedding_metrics")
silhouette_overall <- if (!is.null(emb_metrics) && "cluster" %in% names(emb_metrics)) {
  ov <- emb_metrics$silhouette[emb_metrics$cluster == "overall"]
  if (length(ov)) ov[1] else NA_real_
} else NA_real_

loo <- read_diag("leave_one_out")
n_absorbed_donors <- if (!is.null(loo) && "verdict" %in% names(loo)) {
  sum(loo$verdict == "absorbed", na.rm = TRUE)
} else NA_integer_

# Pre-flight present/absent screen: roll the per-donor verdicts into counts so
# the landing page can flag captures with absent donors.
pa <- read_diag("present_absent")
if (!is.null(pa) && "verdict" %in% names(pa)) {
  n_screened_donors <- nrow(pa)
  n_present         <- sum(pa$verdict == "present", na.rm = TRUE)
  n_absent          <- sum(pa$verdict == "absent", na.rm = TRUE)
  n_untestable      <- sum(pa$verdict == "untestable", na.rm = TRUE)
} else {
  n_screened_donors <- NA_integer_
  n_present <- NA_integer_
  n_absent <- NA_integer_;  n_untestable <- NA_integer_
}

n_cells   <- nrow(cells)
n_singlet <- sum(cells$status == "singlet", na.rm = TRUE)
n_doublet <- sum(cells$status == "doublet", na.rm = TRUE)
recovered <- unique(cells$assignment[cells$status == "singlet"])
recovered <- recovered[!is.na(recovered)]

tibble(
  capture_id            = capture,
  pool_id               = pool,
  n_barcodes            = n_cells,
  n_singlet             = n_singlet,
  n_doublet             = n_doublet,
  n_unassigned          = n_cells - n_singlet - n_doublet,
  frac_singlet          = n_singlet / n_cells,
  n_declared            = length(declared),
  n_recovered           = length(intersect(declared, recovered)),
  missing_samples       = paste(setdiff(declared, recovered), collapse = ";"),
  ambient_souporcell    = soup_amb / 100,
  ambient_vireo_median  = median(amb_sum$ambient_other, na.rm = TRUE),
  mean_methods_agreeing = mean_agree,
  silhouette_overall    = silhouette_overall,
  n_absorbed_donors     = n_absorbed_donors,
  pa_n_donors           = n_screened_donors,
  pa_n_present          = n_present,
  pa_n_absent           = n_absent,
  pa_n_untestable       = n_untestable
) %>% write_tsv(snakemake@output[["metrics"]])

message("wrote ", snakemake@output[["metrics"]])
