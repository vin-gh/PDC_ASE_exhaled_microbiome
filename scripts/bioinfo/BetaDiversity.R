# =====================================================================
# File:    BetaDiversity.R
# Title:   Beta-diversity PCoA (Bray–Curtis, Weighted UniFrac, Unweighted UniFrac)
#
# Purpose:
#   Render three vertically stacked PCoA panels (Bray–Curtis → Weighted UniFrac
#   → Unweighted UniFrac) comparing ASE vs PDC, matching Fig. 4e–g.
#
# Usage:
#   - Place three Excel files under ./data/ :
#       bray_curtis_PCoA.xlsx
#       weighted_unifrac_PCoA.xlsx
#       unweighted_unifrac_PCoA.xlsx
#     Each must contain columns: sample, PC1, PC2
#     (Optional) If present, PC1_pct and PC2_pct (%) will be appended to axis labels.
#   - In R (>= 4.0):
#       install.packages(c("tidyverse","readxl","ggpubr","svglite"))
#       source("BetaDiversity.R")
#
# Output:
#   ./figures/beta/PCoA_3panels_vertical.(pdf|svg|tiff)
#
# Author : Gao Han
# Affil. : School of Biological Science and Medical Engineering, Southeast University
# Version: 1.0
# Date   : 2025-11-04
# Notes  : Shapes—ASE: square (22), PDC: circle (21). 95% normal ellipses per group.
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
indir  <- file.path("data")
outdir <- file.path("figures", "beta")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

file_bray <- file.path(indir, "bray_curtis_PCoA.xlsx")
file_wu   <- file.path(indir, "weighted_unifrac_PCoA.xlsx")
file_uwu  <- file.path(indir, "unweighted_unifrac_PCoA.xlsx")

outfile_pdf  <- file.path(outdir, "PCoA_3panels_vertical.pdf")
outfile_svg  <- file.path(outdir, "PCoA_3panels_vertical.svg")
outfile_tiff <- file.path(outdir, "PCoA_3panels_vertical.tiff")

# ----------------------------
# Read & parse helper
# ----------------------------
read_pcoa_tbl <- function(file_path){
  if (!file.exists(file_path)) stop("File not found: ", file_path)
  
  tb <- read_excel(file_path)
  names(tb) <- tolower(trimws(names(tb)))
  
  # allow either "sample" or first column as sample
  if (!("sample" %in% names(tb))) names(tb)[1] <- "sample"
  if (!all(c("pc1","pc2") %in% names(tb))) {
    # if second/third columns are PC1/PC2, rename defensively
    if (ncol(tb) >= 3) {
      names(tb)[2:3] <- c("pc1","pc2")
    } else {
      stop("Expect columns: sample, PC1, PC2 in ", basename(file_path))
    }
  }
  
  tb %>%
    transmute(
      sample   = as.character(sample),
      PC1      = as.numeric(pc1),
      PC2      = as.numeric(pc2),
      PC1_pct  = dplyr::coalesce(as.numeric(tb[["pc1_pct"]]), NA_real_),
      PC2_pct  = dplyr::coalesce(as.numeric(tb[["pc2_pct"]]), NA_real_),
      Platform = stringr::str_split_fixed(sample, "\\.", 3)[,1],
      Volunteer= stringr::str_split_fixed(sample, "\\.", 3)[,2]
    ) %>%
    mutate(
      Platform = dplyr::if_else(
        stringr::str_detect(Platform, "^ASE$", ignore_case = TRUE), "ASE",
        dplyr::if_else(stringr::str_detect(Platform, "^PDC$", ignore_case = TRUE), "PDC", NA_character_)
      ),
      Platform = factor(Platform, levels = c("ASE","PDC"))
    ) %>%
    filter(!is.na(Platform))
}

tb_bray <- read_pcoa_tbl(file_bray)
tb_wu   <- read_pcoa_tbl(file_wu)
tb_uwu  <- read_pcoa_tbl(file_uwu)

# ----------------------------
# Palette & theme
# ----------------------------
pal_grp <- c("ASE" = "#E09848", "PDC" = "#48B880")  # aligned with alpha plot colors

theme_pcoa <- theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_line(linewidth = 0.3, linetype = "dashed", colour = "grey85"),
    panel.grid.minor = element_blank(),
    panel.border     = element_rect(colour = "black", linewidth = 0.6),
    axis.title       = element_text(face = "bold", size = 15),
    axis.text        = element_text(size = 12),
    legend.position  = "bottom",
    legend.title     = element_blank(),
    legend.text      = element_text(size = 14),
    plot.title       = element_text(face = "bold", hjust = 0, size = 16),
    plot.margin      = margin(10, 14, 10, 14)
  )

# ----------------------------
# Axis label helper (optionally append % variance)
# ----------------------------
axis_lab <- function(x_default, pct){
  if (is.na(pct)) x_default else sprintf("%s (%.1f%%)", x_default, pct)
}

# ----------------------------
# Single-panel constructor
# ----------------------------
make_pcoa_plot <- function(df, title_lab, ellipse_level = 0.95) {
  xlab <- axis_lab("PCoA1", unique(df$PC1_pct)[!is.na(unique(df$PC1_pct))][1] %||% NA_real_)
  ylab <- axis_lab("PCoA2", unique(df$PC2_pct)[!is.na(unique(df$PC2_pct))][1] %||% NA_real_)
  
  ggplot(df, aes(x = PC1, y = PC2)) +
    geom_hline(yintercept = 0, linetype = "dotted", colour = "grey75", linewidth = 0.4) +
    geom_vline(xintercept = 0, linetype = "dotted", colour = "grey75", linewidth = 0.4) +
    stat_ellipse(
      aes(fill = Platform),
      type = "norm", level = ellipse_level,
      geom = "polygon", alpha = 0.18, show.legend = FALSE, segments = 400
    ) +
    stat_ellipse(
      aes(group = Platform),
      type = "norm", level = ellipse_level,
      color = "grey20", linetype = "dashed", linewidth = 0.9,
      show.legend = FALSE, segments = 400
    ) +
    geom_point(
      aes(fill = Platform, shape = Platform),
      size = 4.6, color = "grey20", stroke = 0.6
    ) +
    scale_shape_manual(values = c("ASE" = 22, "PDC" = 21)) +
    scale_fill_manual(values = pal_grp) +
    labs(title = title_lab, x = xlab, y = ylab) +
    scale_x_continuous(expand = expansion(mult = 0.06)) +
    scale_y_continuous(expand = expansion(mult = 0.06)) +
    theme_pcoa
}

# Panels (top→bottom)
p_bray <- make_pcoa_plot(tb_bray, "Bray–Curtis")
p_wu   <- make_pcoa_plot(tb_wu,   "Weighted UniFrac")
p_uwu  <- make_pcoa_plot(tb_uwu,  "Unweighted UniFrac")

# Shared legend (bottom)
legend_g <- ggpubr::get_legend(
  ggplot(tb_bray, aes(PC1, PC2, fill = Platform, shape = Platform)) +
    geom_point() +
    scale_fill_manual(values = pal_grp) +
    scale_shape_manual(values = c("ASE" = 22, "PDC" = 21)) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom", legend.title = element_blank())
)

# Arrange (one column, three rows + legend)
p_combined <- ggarrange(
  p_bray + theme(legend.position = "none"),
  p_wu   + theme(legend.position = "none"),
  p_uwu  + theme(legend.position = "none"),
  ncol = 1, nrow = 3, heights = c(1,1,1), align = "v"
)

final_plot <- ggarrange(p_combined, legend_g, ncol = 1, heights = c(1, 0.10))

# ----------------------------
# Export
# ----------------------------
W <- 6   # inches
H <- 15  # inches (tall page for three panels + legend)

ggsave(outfile_pdf,  final_plot, width = W, height = H, units = "in")               # vector
ggsave(outfile_svg,  final_plot, width = W, height = H, units = "in", device = "svg")
ggsave(outfile_tiff, final_plot, width = W, height = H, units = "in",
       dpi = 600, compression = "lzw", bg = "white")

message("PCoA figure saved to: ", outdir)
