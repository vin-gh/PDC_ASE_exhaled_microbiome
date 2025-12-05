# =============================================================================
# File:    WilcoxonGnusD.R
# Title:   Genus-level Wilcoxon test (WC vs PDC) with effect sizes
#
# Purpose
#   Compute two-sample Wilcoxon rank-sum tests between wet-wall cyclone (ASE)
#   and PDC-sampler groups at the Genus level, and export a table with:
#   Taxa, avg(ASE), sd(ASE), avg(PDC), sd(PDC), p.value, q.values,
#   interval lower, interval upper, log2FoldChange, Combine ES.
#   The result is intended for manual filtering/sorting in Excel, then
#   visualization by the companion script WilcoxonGnusF.R (see manuscript Fig. 7a).
#
# Input
#   ./data/featureTable.sample.g.relative.txt
#   (tab-delimited; columns: Taxonomy | <samples...> [| optional Tax_detail])
#
# How to run
#   R >= 4.0 with packages: tidyverse, readr, effsize
#   install.packages(c("tidyverse","readr","effsize"))
#   From repo root: Rscript scripts/WilcoxonGnusD.R
#
# Output
#   ./results/stats/Wilcoxon_Genus_ASE_vs_PDC_results.csv
#
# Notes
#   - log2FoldChange = log2( mean(ASE)+eps  /  mean(PDC)+eps )
#   - Combine ES = Cliff's delta (ASE vs PDC; >0 means ASE > PDC)
#   - interval lower/upper = 95% bootstrap CI for (PDC_mean - ASE_mean)
#
# Author : Gao Han
# Affil. : School of Biological Science and Medical Engineering, Southeast University
# Version: 1.0
# Date   : 2025-11-04
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(effsize)  # Cliff's delta
})

# ----------------------------
# Repo-relative paths
# ----------------------------
indir   <- file.path("data")
infile  <- file.path(indir, "featureTable.sample.g.relative.txt")

outdir  <- file.path("results", "stats")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
outfile <- file.path(outdir, "Wilcoxon_Genus_ASE_vs_PDC_results.csv")

# ----------------------------
# Parameters
# ----------------------------
set.seed(123)
B   <- 2000   # bootstrap iterations for CI of mean difference (PDC - ASE)
eps <- 1e-6   # small pseudocount for log2FC

# ----------------------------
# Load and reshape (Genus x Samples, relative abundance)
# ----------------------------
if (!file.exists(infile)) stop("Input not found: ", infile)

dat_raw <- read_tsv(infile, show_col_types = FALSE)

# Sample columns = everything except taxonomy columns (if present)
sample_cols <- setdiff(colnames(dat_raw), c("Taxonomy","Tax_detail"))

dat_long <- dat_raw %>%
  pivot_longer(all_of(sample_cols), names_to = "SampleID", values_to = "RelAbund") %>%
  mutate(
    Group    = ifelse(grepl("^ASE", SampleID), "ASE", "PDC"),
    RelAbund = as.numeric(RelAbund)
  )

# ----------------------------
# Bootstrap CI helper for mean difference (PDC - ASE)
# ----------------------------
boot_ci_diff <- function(x_ase, x_pdc, B = B) {
  x_ase <- as.numeric(x_ase); x_pdc <- as.numeric(x_pdc)
  x_ase <- x_ase[is.finite(x_ase)]; x_pdc <- x_pdc[is.finite(x_pdc)]
  if (length(x_ase) == 0 || length(x_pdc) == 0) return(c(NA_real_, NA_real_))
  diffs <- replicate(B, mean(sample(x_pdc, replace = TRUE)) - mean(sample(x_ase, replace = TRUE)))
  stats::quantile(diffs, probs = c(0.025, 0.975), na.rm = TRUE) %>% as.numeric()
}

# ----------------------------
# Per-genus statistics
# ----------------------------
res_tbl <- dat_long %>%
  group_by(Taxonomy) %>%
  summarise(
    ASE_vec = list(RelAbund[Group == "ASE"]),
    PDC_vec = list(RelAbund[Group == "PDC"]),
    avg_ASE = mean(RelAbund[Group == "ASE"], na.rm = TRUE),
    sd_ASE  = sd(RelAbund[Group == "ASE"], na.rm = TRUE),
    avg_PDC = mean(RelAbund[Group == "PDC"], na.rm = TRUE),
    sd_PDC  = sd(RelAbund[Group == "PDC"], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    p_value = tryCatch(
      wilcox.test(unlist(ASE_vec), unlist(PDC_vec),
                  alternative = "two.sided", exact = FALSE)$p.value,
      error = function(e) NA_real_
    ),
    combine_es = tryCatch(
      effsize::cliff.delta(unlist(ASE_vec), unlist(PDC_vec))$estimate,
      error = function(e) NA_real_
    ),
    diff_CI = list(boot_ci_diff(unlist(ASE_vec), unlist(PDC_vec), B = B)),
    log2FC  = log2((avg_ASE + eps) / (avg_PDC + eps))
  ) %>%
  ungroup() %>%
  mutate(
    q_values         = p.adjust(p_value, method = "BH"),
    `interval lower` = purrr::map_dbl(diff_CI, 1),
    `interval upper` = purrr::map_dbl(diff_CI, 2)
  ) %>%
  select(
    Taxa              = Taxonomy,
    `avg(ASE)`        = avg_ASE,
    `sd(ASE)`         = sd_ASE,
    `avg(PDC)`        = avg_PDC,
    `sd(PDC)`         = sd_PDC,
    `p.value`         = p_value,
    `q.values`        = q_values,
    `interval lower`,
    `interval upper`,
    `log2FoldChange`  = log2FC,
    `Combine ES`      = combine_es
  ) %>%
  arrange(`p.value`)

# ----------------------------
# Write out
# ----------------------------
write.csv(res_tbl, outfile, row.names = FALSE)
message("[WilcoxonGnusD] Saved -> ", outfile)
message("Note: Manually filter/sort by p.value and log2FoldChange in Excel, then visualize with WilcoxonGnusF.R.")
