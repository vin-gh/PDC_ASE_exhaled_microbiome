# =============================================================================
# File:    Shared_Taxson.R
# Title:   Pie chart of shared genera by combined abundance (WC + PDC)
#
# Purpose
#   Compute genus-level composition using the *shared* genera between
#   wet-wall cyclone (ASE) and PDC-sampler groups, define the sector size
#   as Combined = ASE + PDC, and export a publication-ready pie chart.
#   (Labels/indices are intentionally minimal; fine-tune in Illustrator.)
#
# How to use
#   1) Place inputs under ./data/
#        - featureTable.group.g.relative.txt
#          (tab-delimited; columns: Genus | WC | PDC [| optional Tax_detail])
#        - ASE-PDC_ASE_PDC.txt  (optional; Venn-derived shared ASVs for strict filtering)
#   2) In R (>= 4.0), install dependencies once:
#        install.packages(c("tidyverse","readr","ggpubr","ggsci","svglite"))
#   3) Run in the repo root:
#        source("scripts/Shared_Taxson.R")     # if this file is under ./scripts/
#      or from shell:
#        Rscript scripts/Shared_Taxson.R
#
# Outputs
#   - ./figures/composition/Pie_Shared_Genus_byCombinedAbundance.(pdf|svg|tiff)
#
# Notes
#   - “Shared genera” can be defined in two ways:
#       (a) Threshold-based (default): genus present in BOTH groups >= min_abund
#       (b) Strict (optional): restrict genera using Venn-derived shared ASVs
#   - After export, adjust label placement/numbering in a vector editor as needed.
#
# Author : Gao Han
# Affil. : School of Biological Science and Medical Engineering, Southeast University
# Version: 1.0
# Date   : 2025-11-04
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readr)
  library(ggpubr)
  library(ggsci)
  library(svglite)
})

# ----------------------------
# Repo-relative paths
# ----------------------------
indir  <- file.path("data")
in_rel <- file.path(indir, "featureTable.group.g.relative.txt")  # Genus x 2 groups (relative abundance)
in_shared_asv <- file.path(indir, "ASE-PDC_ASE_PDC.txt")         # Optional: Venn shared ASVs (for strict mode)

outdir <- file.path("figures", "composition")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile_pdf  <- file.path(outdir, "Pie_Shared_Genus_byCombinedAbundance.pdf")
outfile_svg  <- file.path(outdir, "Pie_Shared_Genus_byCombinedAbundance.svg")
outfile_tiff <- file.path(outdir, "Pie_Shared_Genus_byCombinedAbundance.tiff")

# ----------------------------
# Parameters
# ----------------------------
min_abund        <- 1e-4  # Presence threshold per group (relative abundance)
topN             <- 10    # Show top-N genera as sectors; others -> "Others"
use_shared_file  <- FALSE # If TRUE, restrict genera via Venn shared ASVs file

# ----------------------------
# Load main input (Genus-level relative abundance for ASE & PDC)
# Expected columns: Genus | ASE | PDC [| Tax_detail]
# ----------------------------
if (!file.exists(in_rel)) stop("Input file not found: ", in_rel)

rel <- read_tsv(in_rel, show_col_types = FALSE, col_names = TRUE) %>%
  rename(Genus = 1, ASE = 2, PDC = 3) %>%
  mutate(
    ASE = as.numeric(ASE),
    PDC = as.numeric(PDC)
  ) %>%
  filter(!is.na(Genus), (!is.na(ASE) | !is.na(PDC))) %>%
  group_by(Genus) %>%                                      # collapse duplicated genus rows if any
  summarise(
    ASE = sum(replace_na(ASE, 0)),
    PDC = sum(replace_na(PDC, 0)),
    .groups = "drop"
  )

# ----------------------------
# Define shared genera
#   A) strict (Venn-based) if use_shared_file = TRUE
#   B) threshold-based (default): ASE>=min_abund & PDC>=min_abund
# ----------------------------
if (use_shared_file) {
  if (!file.exists(in_shared_asv)) {
    stop("Shared-ASV file not found but use_shared_file=TRUE: ", in_shared_asv)
  }
  shared_raw <- read_tsv(in_shared_asv, show_col_types = FALSE, col_names = TRUE)
  
  extract_genus <- function(tx) {
    g <- stringr::str_extract(tx, "g__[^;]*")
    g <- gsub("^g__", "", g)
    ifelse(is.na(g) | g == "", "unclassified", g)
  }
  
  shared_genera <- shared_raw %>%
    mutate(Genus = extract_genus(.data[[names(shared_raw)[1]]])) %>% # first col contains taxonomy string
    distinct(Genus) %>%
    pull(Genus) %>%
    gsub("^g__", "", .)
  
  rel_shared <- rel %>% filter(Genus %in% shared_genera)
} else {
  rel_shared <- rel %>% filter(ASE >= min_abund & PDC >= min_abund)
}

if (nrow(rel_shared) == 0) {
  stop("Shared genus set is empty. Lower 'min_abund' or set use_shared_file=TRUE.")
}

# ----------------------------
# Combined abundance (ASE + PDC) and percent
# ----------------------------
pie_df <- rel_shared %>%
  transmute(
    Genus,
    Combined = ASE + PDC
  ) %>%
  arrange(desc(Combined)) %>%
  mutate(Percent = Combined / sum(Combined))

# ----------------------------
# Keep top-N genera, lump the rest to "Others"
# ----------------------------
pie_top <- pie_df %>%
  mutate(Genus_grp = forcats::fct_lump_n(Genus, n = topN, w = Combined, other_level = "Others")) %>%
  group_by(Genus_grp) %>%
  summarise(Combined = sum(Combined), .groups = "drop") %>%
  arrange(desc(Combined)) %>%
  mutate(
    Percent   = Combined / sum(Combined),
    Genus_grp = forcats::fct_reorder(Genus_grp, Combined, .desc = TRUE)
  )

# ----------------------------
# Palette (NPG). If k>10, extend by interpolation.
# ----------------------------
k <- nrow(pie_top)
pal_base <- ggsci::pal_npg("nrc")(10)
pal_vec  <- if (k <= 10) pal_base else grDevices::colorRampPalette(pal_base)(k)

# ----------------------------
# Labels (short; will be fine-tuned in Illustrator)
# ----------------------------
wrap_genus <- function(x, width = 16) stringr::str_wrap(x, width = width)
min_label_pct <- 0.00  # label all by default; raise (e.g., 0.02) to label only ≥2%

pie_top <- pie_top %>%
  mutate(
    Label = ifelse(Percent >= min_label_pct | Genus_grp == "Others",
                   paste0(wrap_genus(as.character(Genus_grp)), ": ",
                          scales::percent(Percent, accuracy = 0.01)),
                   "")
  )

# ----------------------------
# Plot
# ----------------------------
p_pie <- ggpubr::ggpie(
  data    = pie_top,
  x       = "Combined",
  label   = "Label",
  fill    = "Genus_grp",
  color   = "white",
  lab.pos = "out",
  repel   = TRUE,
  lab.font = c(4, "plain", 12),
  radius  = 0.82
) +
  scale_fill_manual(values = pal_vec) +
  ggpubr::theme_pubr(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
    plot.margin = margin(10, 80, 10, 10)  # extra right margin for outer labels
  ) +
  labs(title = "Shared genera between ASE and PDC (by combined abundance)")

# ----------------------------
# Export (large canvas; edit labels in AI if needed)
# ----------------------------
W <- 10
H <- 8
ggsave(outfile_pdf,  p_pie, width = W, height = H, dpi = 300)
ggsave(outfile_svg,  p_pie, width = W, height = H, dpi = 300)
ggsave(outfile_tiff, p_pie, width = W, height = H, dpi = 600, compression = "lzw", bg = "white")

# Console summary
cat("\n[Shared_Taxson] Done.\n")
cat(sprintf("  Shared genera (n) : %d\n", nrow(rel_shared)))
cat(sprintf("  Top-N displayed   : %d (Others lumped if applicable)\n", topN))
cat(sprintf("  Saved: %s\n         %s\n         %s\n", outfile_pdf, outfile_svg, outfile_tiff))
cat("  Note: adjust label positions/indices in Illustrator for final layout.\n")
