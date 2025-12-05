# =====================================================================
# File:    ANCOMBC2.R
# Title:   ANCOM-BC2 differential abundance (PDC vs WC) + LFC heatmap
#
# Purpose
#   Perform ANCOM-BC2 at Genus level to compare PDC-sampler vs wet-wall
#   cyclone (ASE), using a paired random intercept for Volunteer. Falls
#   back to a GLMM (Poisson, offset by library size, random intercept
#   for Volunteer) if ANCOM-BC2 fails. Exports a CSV table and an LFC
#   heatmap matching Fig. 7b (and Fig. S3).
#
# How to use
#   1) Place the input table under ./data/
#        - featureTable.sample.total.absolute.txt
#        (tab-delimited; columns = ASV, <sample columns...>, Taxonomy)
#   2) In R (>= 4.0), install dependencies:
#        install.packages(c("tidyverse","svglite"))
#        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#        BiocManager::install(c("phyloseq","ANCOMBC"))
#   3) Run:
#        source("ANCOMBC2.R")
#
# Input assumptions
#   - Sample names encode group and subject as "ASE.V01.*" or "PDC.V01.*"
#     so that Platform = ASE/PDC and Volunteer = Vxx can be parsed.
#
# Output
#   - ./figures/differential/ANCOMBC2_PDC_vs_ASE_results_Genus.csv
#   - ./figures/differential/ANCOMBC2_LFC_heatmap_PDC_vs_ASE.(pdf|svg|tiff)
#     (Positive LFC means PDC > ASE)
#
# Author : Gao Han
# Affil. : School of Biological Science and Medical Engineering, Southeast University
# Version: 1.0
# Date   : 2025-11-04
# Notes  : Genus-level analysis; prevalence/depth filters mimic article settings.
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(phyloseq)
  library(ANCOMBC)
  library(svglite)
})

set.seed(1)

# ----------------------------
# Paths (repo-relative)
# ----------------------------
indir   <- file.path("data")
in_file <- file.path(indir, "featureTable.sample.total.absolute.txt")

outdir <- file.path("figures", "differential")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

f_pdf <- file.path(outdir, "ANCOMBC2_LFC_heatmap_PDC_vs_ASE.pdf")
f_svg <- file.path(outdir, "ANCOMBC2_LFC_heatmap_PDC_vs_ASE.svg")
f_tif <- file.path(outdir, "ANCOMBC2_LFC_heatmap_PDC_vs_ASE.tiff")
f_csv <- file.path(outdir, "ANCOMBC2_PDC_vs_ASE_results_Genus.csv")

# ----------------------------
# Read ASV table and parse taxonomy
#   Expect: first column = ASV, last column = Taxonomy,
#           sample columns in between.
# ----------------------------
if (!file.exists(in_file)) stop("Input file not found: ", in_file)

raw <- readr::read_tsv(in_file, show_col_types = FALSE, col_names = TRUE)
colnames(raw)[1] <- "ASV"
if (!"Taxonomy" %in% colnames(raw)) colnames(raw)[ncol(raw)] <- "Taxonomy"

sample_cols <- setdiff(colnames(raw), c("ASV","Taxonomy"))

otu_mat <- raw %>%
  select(all_of(c("ASV", sample_cols))) %>%
  column_to_rownames("ASV") %>%
  as.matrix()

extract_rank <- function(tx, prefix){
  x <- stringr::str_extract(tx, paste0(prefix, "[^;]*"))
  x <- gsub(paste0("^", prefix), "", x)
  ifelse(is.na(x) | x == "", NA, x)
}

tax_df <- tibble(
  ASV     = raw$ASV,
  Kingdom = extract_rank(raw$Taxonomy, "k__"),
  Phylum  = extract_rank(raw$Taxonomy, "p__"),
  Class   = extract_rank(raw$Taxonomy, "c__"),
  Order   = extract_rank(raw$Taxonomy, "o__"),
  Family  = extract_rank(raw$Taxonomy, "f__"),
  Genus   = extract_rank(raw$Taxonomy, "g__"),
  Species = extract_rank(raw$Taxonomy, "s__")
) %>%
  mutate(Genus = ifelse(is.na(Genus) | Genus == "", "Unclassified", Genus)) %>%
  column_to_rownames("ASV")

tax_mat <- as.matrix(tax_df)

# ----------------------------
# Sample metadata from sample names
#   Platform = ASE/PDC (reference = ASE)
#   Volunteer = Vxx (paired random effect)
#   Replicate = 3rd token if present
#   libsize = library size (for offsets / QC)
# ----------------------------
libsize <- colSums(otu_mat, na.rm = TRUE)

meta <- tibble(SampleID = sample_cols) %>%
  mutate(
    Platform  = ifelse(stringr::str_detect(SampleID, "^ASE"), "ASE", "PDC"),
    Volunteer = stringr::str_extract(SampleID, "V\\d{2}"),
    Replicate = stringr::str_split_fixed(SampleID, "\\.", 3)[,3],
    libsize   = as.numeric(libsize[SampleID])
  ) %>%
  mutate(
    Platform  = factor(Platform, levels = c("ASE","PDC")),
    Volunteer = factor(Volunteer),
    Replicate = factor(Replicate)
  ) %>%
  column_to_rownames("SampleID")

# ----------------------------
# Build phyloseq object
# ----------------------------
ps <- phyloseq(
  otu_table(otu_mat, taxa_are_rows = TRUE),
  tax_table(tax_mat),
  sample_data(meta)
)

# ----------------------------
# Heatmap helper (uses Genus + LFC)
#   - If no q<0.05, shows top-N by |LFC|
#   - Positive LFC => PDC > ASE
# ----------------------------
plot_and_save <- function(df_res, title_prefix){
  show_tbl <- df_res %>%
    filter(!is.na(lfc)) %>%
    { if (sum(.$q < 0.05, na.rm = TRUE) > 0) filter(., q < 0.05) else head(., 25) } %>%
    arrange(lfc)
  
  if (nrow(show_tbl) == 0) stop("No taxa to display (all LFC are NA).")
  
  plt_df <- show_tbl %>%
    mutate(Taxon = factor(Taxon, levels = Taxon),
           sig   = ifelse(q < 0.05, "sig","ns")) %>%
    transmute(Taxon, Contrast = "PDC_vs_ASE", LFC = lfc, sig)
  
  lim <- max(abs(plt_df$LFC), na.rm = TRUE); lim <- max(lim, 0.5)
  
  theme_heat <- ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid    = element_blank(),
      axis.title    = element_blank(),
      axis.text.x   = element_text(face = "bold"),
      axis.text.y   = element_text(size = 10),
      legend.position = "right",
      plot.title    = element_text(face = "bold", hjust = 0.5)
    )
  
  p <- ggplot(plt_df, aes(x = Contrast, y = Taxon, fill = LFC)) +
    geom_tile(color = "grey90", linewidth = 0.4) +
    geom_text(aes(label = sprintf("%.2f", LFC),
                  color = factor(sig, levels = c("sig","ns"))),
              size = 3.2) +
    scale_color_manual(values = c("sig" = "#1B9E77", "ns" = "grey20"), guide = "none") +
    scale_fill_gradient2(
      low = "#3B8EE3", mid = "white", high = "#D75445",
      midpoint = 0, limits = c(-lim, lim),
      name = "Log fold-change\n(PDC vs ASE)"
    ) +
    labs(title = paste0(title_prefix, " (Genus)\nPositive LFC = PDC higher than ASE")) +
    theme_heat
  
  W <- 6
  H <- max(4, 0.35 * nrow(plt_df) + 1.5)
  ggsave(f_pdf, p, width = W, height = H, dpi = 300)
  ggsave(f_svg, p, width = W, height = H, dpi = 300)
  ggsave(f_tif, p, width = W, height = H, dpi = 600, compression = "lzw", bg = "white")
  invisible(p)
}

# ----------------------------
# ANCOM-BC2 (primary)
#   - Genus level
#   - Fixed: Platform
#   - Random: ~ 1 | Volunteer
#   - Prevalence and depth cutoffs follow article-style heuristics
# ----------------------------
run_ancombc2 <- function(ps_obj){
  ancombc2(
    data         = ps_obj,
    tax_level    = "Genus",
    fix_formula  = "Platform",
    rand_formula = "~ 1 | Volunteer",
    group        = "Platform",
    struc_zero   = TRUE,
    neg_lb       = TRUE,
    prv_cut      = 0.10,
    lib_cut      = 1000,
    s0_perc      = 0.05,
    alpha        = 0.05,
    p_adj_method = "BY",
    global       = FALSE
  )
}

try_A <- try(run_ancombc2(ps), silent = TRUE)

if (!inherits(try_A, "try-error")) {
  lfc_tbl <- as.data.frame(try_A$res$lfc)
  q_tbl   <- as.data.frame(try_A$res$q_val)
  pick    <- which(grepl("Platform", colnames(lfc_tbl))); if (length(pick) == 0) pick <- 1
  colname <- colnames(lfc_tbl)[pick[1]]
  
  df_out <- tibble(
    Taxon = rownames(lfc_tbl),
    lfc   = lfc_tbl[[colname]],
    q     = q_tbl[[colname]]
  ) %>% arrange(desc(abs(lfc)))
  
  write.csv(df_out, f_csv, row.names = FALSE)
  plot_and_save(df_out, "ANCOM-BC2")
} else {
  # --------------------------
  # Retry without lmerTest attached (if any)
  # --------------------------
  if ("package:lmerTest" %in% search()) {
    detach("package:lmerTest", unload = TRUE, character.only = TRUE)
  }
  try_B <- try(run_ancombc2(ps), silent = TRUE)
  
  if (!inherits(try_B, "try-error")) {
    lfc_tbl <- as.data.frame(try_B$res$lfc)
    q_tbl   <- as.data.frame(try_B$res$q_val)
    pick    <- which(grepl("Platform", colnames(lfc_tbl))); if (length(pick) == 0) pick <- 1
    colname <- colnames(lfc_tbl)[pick[1]]
    
    df_out <- tibble(
      Taxon = rownames(lfc_tbl),
      lfc   = lfc_tbl[[colname]],
      q     = q_tbl[[colname]]
    ) %>% arrange(desc(abs(lfc)))
    
    write.csv(df_out, f_csv, row.names = FALSE)
    plot_and_save(df_out, "ANCOM-BC2 (no lmerTest)")
  } else {
    # --------------------------
    # Fallback: GLMM (Poisson + offset(log libsize) + (1|Volunteer))
    # --------------------------
    suppressPackageStartupMessages(library(lme4))
    
    genus_counts <- as.data.frame(otu_mat) %>%
      rownames_to_column("ASV") %>%
      left_join(
        tibble(ASV = rownames(tax_mat), Genus = as.character(tax_mat[,"Genus"])),
        by = "ASV"
      ) %>%
      select(-ASV) %>%
      relocate(Genus) %>%
      group_by(Genus) %>%
      summarise(across(everything(), ~ sum(.x, na.rm = TRUE)), .groups = "drop")
    
    fit_one_genus <- function(row_i){
      gname <- row_i$Genus
      long_df <- row_i %>%
        pivot_longer(-Genus, names_to = "SampleID", values_to = "count") %>%
        left_join(meta %>% rownames_to_column("SampleID"), by = "SampleID") %>%
        mutate(count = as.numeric(count))
      
      if (sum(long_df$count, na.rm = TRUE) == 0) return(tibble(Taxon = gname, lfc = NA, p = NA))
      if (length(unique(long_df$Platform[long_df$count > 0])) < 2) return(tibble(Taxon = gname, lfc = NA, p = NA))
      
      fm  <- count ~ Platform + offset(log(libsize)) + (1|Volunteer)
      fit <- try(
        glmer(fm, data = long_df, family = poisson(link = "log"),
              control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))),
        silent = TRUE
      )
      if (inherits(fit, "try-error")) return(tibble(Taxon = gname, lfc = NA, p = NA))
      
      cf <- summary(fit)$coef
      if (!"PlatformPDC" %in% rownames(cf)) return(tibble(Taxon = gname, lfc = NA, p = NA))
      
      tibble(Taxon = gname,
             lfc   = unname(cf["PlatformPDC","Estimate"]),
             p     = unname(cf["PlatformPDC","Pr(>|z|)"]))
    }
    
    res_glmm <- genus_counts %>%
      rowwise() %>% group_split() %>%
      purrr::map_dfr(~ fit_one_genus(.x)) %>%
      mutate(q = p.adjust(p, method = "BY")) %>%
      arrange(desc(abs(lfc)))
    
    write.csv(res_glmm, f_csv, row.names = FALSE)
    plot_and_save(res_glmm, "GLMM fallback (Poisson + offset)")
  }
}

message("Done. Results written to: ", f_csv)
message("Heatmap written to: ", f_pdf, " / ", f_svg, " / ", f_tif)
