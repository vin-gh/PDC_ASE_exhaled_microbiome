# =============================================================================
# File:    WilcoxonGnusF.R
# Title:   Visualization of curated Genus-level Wilcoxon results (WC vs PDC)
#
# Purpose
#   Plot five independent panels from the *curated* results produced by
#   WilcoxonGnusD.R after manual filtering/sorting in Excel:
#     (1) Group means ± SD (percent)
#     (2) Forest plot of mean difference (PDC − ASE) with 95% CI
#     (3) p-values with significance stars
#     (4) Combined effect size (Cliff’s delta) lollipop
#     (5) log2FoldChange (ASE vs PDC) lollipop
#   This script generates figure files for manuscript Fig. 7a.
#
# Input (place one of the following in ./results/stats/)
#   - Wilcoxon_Genus_ASE_vs_PDC_results.xlsx   (recommended for curated file), or
#   - Wilcoxon_Genus_ASE_vs_PDC_results.csv
#   Required columns (exact names, case-sensitive):
#     Taxa (or Taxonomy), avg(ASE), sd(ASE), avg(PDC), sd(PDC),
#     p.value, q.values, interval lower, interval upper,
#     log2FoldChange, Combine ES
#
# How to run
#   R >= 4.0, packages: tidyverse, readxl, readr, ggsci, svglite
#   install.packages(c("tidyverse","readxl","readr","ggsci","svglite"))
#   From repo root:
#     Rscript scripts/WilcoxonGnusF.R
#
# Output
#   ./results/figures/fig7a/Fig7a_Panel*.pdf / .svg / .tiff
#
# Notes
#   - This script does not perform statistical tests. It only visualizes the
#     curated table you prepared from WilcoxonGnusD.R output.
#   - Differences and CIs are interpreted as (PDC − ASE). log2FC is ASE vs PDC.
#
# Author : Gao Han
# Affil. : School of Biological Science and Medical Engineering, Southeast University
# Version: 1.0
# Date   : 2025-11-04
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(readr)
  library(ggsci)
  library(svglite)
})

# ----------------------------
# Repo-relative paths
# ----------------------------
in_dir  <- file.path("results", "stats")
out_dir <- file.path("results", "figures", "fig7a")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Prefer the curated Excel file; fall back to CSV if needed
in_file_xlsx <- file.path(in_dir, "Wilcoxon_Genus_ASE_vs_PDC_results.xlsx")
in_file_csv  <- file.path(in_dir, "Wilcoxon_Genus_ASE_vs_PDC_results.csv")

if (file.exists(in_file_xlsx)) {
  raw <- read_xlsx(in_file_xlsx)
} else if (file.exists(in_file_csv)) {
  raw <- read_csv(in_file_csv, show_col_types = FALSE)
} else {
  stop("Input not found. Provide either:\n  ",
       in_file_xlsx, "\n  or\n  ", in_file_csv)
}

# ----------------------------
# Standardize columns
# ----------------------------
tax_col <- dplyr::case_when(
  "Taxa" %in% names(raw)      ~ "Taxa",
  "Taxonomy" %in% names(raw)  ~ "Taxonomy",
  TRUE                        ~ ""
)
if (tax_col == "") stop("Missing required taxonomy column: 'Taxa' or 'Taxonomy'.")

raw <- raw %>%
  rename(Taxon = all_of(tax_col)) %>%
  mutate(across(everything(), ~ .x)) %>%  # noop to keep dplyr chain
  { names(.) <- trimws(names(.)); . } %>%
  mutate(Taxon = trimws(Taxon))

# Enforce numeric for all but Taxon
num_cols <- setdiff(names(raw), "Taxon")
df <- raw %>%
  mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.)))) %>%
  # Required renames
  rename(
    avg_ASE = `avg(ASE)`, sd_ASE = `sd(ASE)`,
    avg_PDC = `avg(PDC)`, sd_PDC = `sd(PDC)`,
    p_value = `p.value`,  q_values = `q.values`,
    int_low = `interval lower`, int_up = `interval upper`,
    log2FC  = `log2FoldChange`,
    combine_es = `Combine ES`
  )

need_cols <- c("avg_ASE","sd_ASE","avg_PDC","sd_PDC",
               "p_value","q_values","int_low","int_up","log2FC","combine_es")
if (!all(need_cols %in% names(df))) {
  missing <- setdiff(need_cols, names(df))
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}
if (any(!sapply(df[need_cols], is.numeric))) {
  bad <- names(df[need_cols])[!sapply(df[need_cols], is.numeric)]
  stop("Non-numeric values in numeric columns: ", paste(bad, collapse = ", "),
       ". Check your curated file (e.g., stray characters like '–' or blanks).")
}

# Preserve current order (top row shown on top)
df$Taxon <- factor(df$Taxon, levels = rev(unique(df$Taxon)))

# Derived display columns (percent units)
df <- df %>%
  mutate(
    avg_ASE_pct = avg_ASE * 100,
    avg_PDC_pct = avg_PDC * 100,
    sd_ASE_pct  = sd_ASE  * 100,
    sd_PDC_pct  = sd_PDC  * 100,
    diff        = avg_PDC - avg_ASE,
    diff_pct    = diff * 100,
    int_low_pct = int_low * 100,
    int_up_pct  = int_up  * 100
  )

# ----------------------------
# Plot helpers
# ----------------------------
save_fig <- function(plot_obj, filename_base,
                     width = 8,
                     height = max(6, 0.45 * nlevels(plot_obj$data$Taxon) + 2)) {
  ggsave(file.path(out_dir, paste0(filename_base, ".pdf")),  plot_obj, width = width, height = height, dpi = 300)
  ggsave(file.path(out_dir, paste0(filename_base, ".svg")),  plot_obj, width = width, height = height, dpi = 300)
  ggsave(file.path(out_dir, paste0(filename_base, ".tiff")), plot_obj, width = width, height = height, dpi = 600,
         compression = "lzw", bg = "white")
}

theme_base <- theme_bw(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    axis.title.y     = element_blank()
  )

p_stars <- function(p){
  case_when(
    is.na(p)   ~ "",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ ""
  )
}

pal_groups <- ggsci::pal_nejm()(2); names(pal_groups) <- c("ASE","PDC")
col_dir <- c("pos" = "#1B9E77", "neg" = "#D75445")

# ----------------------------
# Panel 1: Means ± SD (percent)
# ----------------------------
df_bar <- tibble(
  Taxon = rep(df$Taxon, each = 2),
  Group = rep(c("ASE","PDC"), times = nrow(df)),
  mean  = c(rbind(df$avg_ASE_pct, df$avg_PDC_pct)),
  sd    = c(rbind(df$sd_ASE_pct,  df$sd_PDC_pct))
)

p1 <- ggplot(df_bar, aes(x = Taxon, y = mean, fill = Group)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd),
                position = position_dodge(width = 0.7), width = 0.25, linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(values = pal_groups, name = NULL) +
  labs(y = "Proportions (%)", title = "Proportions (%)") +
  theme_base +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", hjust = 0.5))

save_fig(p1, "Fig7a_Panel1_Proportions")

# ----------------------------
# Panel 2: Difference (PDC − ASE) with 95% CI
# ----------------------------
df_diff <- df %>% mutate(dir = ifelse(diff_pct >= 0, "pos", "neg"))

p2 <- ggplot(df_diff, aes(x = Taxon, y = diff_pct, color = dir)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(ymin = int_low_pct, ymax = int_up_pct), size = 0.4) +
  coord_flip() +
  scale_color_manual(values = col_dir, guide = "none") +
  labs(y = "Difference between proportions (%)", title = "95% confidence intervals") +
  theme_base +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

save_fig(p2, "Fig7a_Panel2_DiffCI")

# ----------------------------
# Panel 3: p-values with significance stars
# ----------------------------
df_p <- df %>%
  mutate(
    stars = p_stars(p_value),
    p_lab = formatC(p_value, format = "e", digits = 2),
    label = paste0(stars, "  ", p_lab)
  )

p3 <- ggplot(df_p, aes(x = 1, y = Taxon, label = label)) +
  geom_text(color = "firebrick", size = 4.2, fontface = "bold") +
  coord_cartesian(xlim = c(0.5, 1.5)) +
  labs(title = "P value") +
  theme_base +
  theme(
    plot.title   = element_text(face = "bold", hjust = 0.5),
    panel.grid   = element_blank(),
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  )

save_fig(p3, "Fig7a_Panel3_Pvalue")

# ----------------------------
# Panel 4: Combined ES (Cliff’s delta) lollipop
# ----------------------------
p4 <- ggplot(df, aes(x = Taxon, y = combine_es,
                     color = ifelse(combine_es >= 0, "pos", "neg"))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_segment(aes(xend = Taxon, y = 0, yend = combine_es), linewidth = 0.8) +
  geom_point(size = 2.4) +
  coord_flip() +
  scale_color_manual(values = col_dir, guide = "none") +
  labs(title = "Combine ES (Cliff’s delta)", y = NULL) +
  theme_base +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

save_fig(p4, "Fig7a_Panel4_CombineES")

# ----------------------------
# Panel 5: log2FoldChange (ASE vs PDC) lollipop
# ----------------------------
p5 <- ggplot(df, aes(x = Taxon, y = log2FC,
                     color = ifelse(log2FC >= 0, "pos", "neg"))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_segment(aes(xend = Taxon, y = 0, yend = log2FC), linewidth = 0.8, color = "grey40") +
  geom_point(size = 2.4, color = "grey10") +
  coord_flip() +
  scale_color_manual(values = col_dir, guide = "none") +
  labs(title = "log2FC (ASE vs PDC)", y = NULL) +
  theme_base +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

save_fig(p5, "Fig7a_Panel5_Log2FC")

message("Done. Figures written to: ", out_dir)
