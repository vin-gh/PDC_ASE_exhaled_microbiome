# =====================================================================
# File:    AlphaDiversity.R
# Title:   Alpha-diversity boxplots (ASE vs PDC) for four indices
#
# Purpose: Reproduce Fig. 4a–d. Visualize alpha diversity metrics
#          (observed_features, chao1, Shannon, simpson) as 2×2
#          boxplots comparing ASE vs PDC with Wilcoxon significance.
#
# Usage:
#   - Place input Excel at ./data/Alpha.xlsx
#   - In R (R >= 4.0):
#       install.packages(c("tidyverse","readxl","ggpubr","svglite"))
#       source("AlphaDiversity.R")
#
# Input:
#   ./data/Alpha.xlsx
#   Required columns:
#     Sample_Name, observed_features, chao1, shannon (or Shannon), simpson
#   Sample_Name should encode platform as "ASE.*" or "PDC.*"
#
# Output:
#   ./figures/alpha/Alpha_Boxplots_4in1.(pdf|svg|tiff)
#
# Author : Gao Han
# Affil. : School of Biological Science and Medical Engineering, Southeast University
# Version: 1.0
# Date   : 2025-11-04
# Notes  : Significance via Wilcoxon rank-sum (unpaired). Methods/context
#          per main text; panel mapping matches Fig. 4a–d.
# =====================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(ggpubr)
  library(svglite)
})

# ----------------------------
# Paths (repo-relative)
# ----------------------------
infile <- file.path("data", "Alpha.xlsx")
outdir <- file.path("figures", "alpha")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

outfile_pdf  <- file.path(outdir, "Alpha_Boxplots_4in1.pdf")
outfile_svg  <- file.path(outdir, "Alpha_Boxplots_4in1.svg")
outfile_tiff <- file.path(outdir, "Alpha_Boxplots_4in1.tiff")

# ----------------------------
# Load & validate
# ----------------------------
alpha_raw <- read_excel(infile)

# Normalize metric column names (case-insensitive safety)
names(alpha_raw) <- tolower(trimws(names(alpha_raw)))

required_cols <- c("sample_name","observed_features","chao1","shannon","simpson")
missing_cols  <- setdiff(required_cols, names(alpha_raw))
if (length(missing_cols) > 0) {
  stop("Missing required columns in Alpha.xlsx: ", paste(missing_cols, collapse = ", "))
}

# ----------------------------
# Parse platform from Sample_Name
# ----------------------------
alpha_df <- alpha_raw %>%
  mutate(
    Platform = if_else(str_detect(sample_name, "^ASE\\b", ignore_case = TRUE), "ASE",
                       if_else(str_detect(sample_name, "^PDC\\b", ignore_case = TRUE), "PDC", NA_character_)),
    Volunteer = str_split_fixed(sample_name, "\\.", 3)[, 2],
    Platform  = factor(Platform, levels = c("ASE", "PDC"))
  ) %>%
  filter(!is.na(Platform))

# ----------------------------
# Long format for faceting (2×2)
# ----------------------------
alpha_long <- alpha_df %>%
  pivot_longer(
    cols = c(observed_features, chao1, shannon, simpson),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = c("observed_features","chao1","shannon","simpson"),
      labels = c("Observed features","Chao1","Shannon","Simpson")
    )
  )

# ----------------------------
# Aesthetics (edit to customize)
# ----------------------------
my_cols <- c("ASE" = "#E09848", "PDC" = "#48B880")

my_theme <- theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_line(linewidth = 0.3, linetype = "dashed", colour = "grey85"),
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "grey92", colour = NA),
    strip.text       = element_text(face = "bold"),
    axis.title.x     = element_blank(),
    axis.title.y     = element_text(size = 13),
    axis.text.x      = element_text(size = 12),
    legend.position  = "none",
    plot.margin      = margin(8, 12, 8, 12)
  )

# ----------------------------
# Statistics (Wilcoxon rank-sum)
# ----------------------------
compare_pairs   <- list(c("ASE","PDC"))
wilcox_args     <- list(exact = FALSE, correct = TRUE)  # robust to ties
sigmap <- list(cutpoints = c(0, .001, .01, .05, 1), symbols = c("***","**","*","ns"))

# ----------------------------
# Plot
# ----------------------------
p_alpha <- ggplot(alpha_long, aes(x = Platform, y = value, fill = Platform)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = 21,
    outlier.stroke = 0.3,
    outlier.alpha  = 0.7
  ) +
  geom_jitter(
    width = 0.12,
    alpha = 0.7,
    size  = 1.8,
    colour = "grey25"
  ) +
  scale_fill_manual(values = my_cols) +
  facet_wrap(~ metric, ncol = 2, nrow = 2, scales = "free_y") +
  ylab("Alpha diversity value") +
  my_theme +
  ggpubr::stat_compare_means(
    comparisons  = compare_pairs,
    method       = "wilcox.test",
    method.args  = wilcox_args,
    label        = "p.signif",          # use "p.format" to print numeric p
    symnum.args  = sigmap,
    label.y.npc  = "top",
    size         = 4
    # , hide.ns  = TRUE                 # uncomment to hide "ns"
    # , paired   = TRUE                 # uncomment if strictly paired by Volunteer
  )

# ----------------------------
# Save (vector + high-res raster)
# ----------------------------
ggsave(outfile_pdf,  plot = p_alpha, width = 6.5, height = 9.5, units = "in")          # PDF (vector)
ggsave(outfile_svg,  plot = p_alpha, width = 6.5, height = 9.5, units = "in", device = "svg")
ggsave(outfile_tiff, plot = p_alpha, width = 6.5, height = 9.5, units = "in", dpi = 600, compression = "lzw")

message("Alpha-diversity figure saved to: ", outdir)

# ----------------------------
# Minimal customization guide (edit-and-run)
# ----------------------------
# - Colors:     change hex in `my_cols`.
# - Facet order: edit `levels` in factor(metric) above.
# - Y scales:   set `scales = "fixed"` in facet_wrap for common y-axis.
# - Signif:     switch to numeric p by label = "p.format"; hide ns via `hide.ns = TRUE`.
# - Paired test: add `paired = TRUE` in stat_compare_means if samples are matched by Volunteer.
