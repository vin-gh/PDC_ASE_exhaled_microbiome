# =============================================================================
# File:    WilcoxonPhylumF.R
# Title:   Visualization of curated Phylum-level Wilcoxon results (ASE vs PDC)
#
# Purpose
#   Visualize the *curated* Phylum-level table produced by the companion script
#   (WilcoxonGenusD.R adapted to phylum with input
#   featureTable.sample.p.relative.selected.txt). After you manually filter/sort
#   the results in Excel, this script renders publication-ready panels for
#   manuscript Fig. 7c.
#
# Input (place one of the following in ./results/stats/)
#   - Wilcoxon_Phylum_ASE_vs_PDC_results.xlsx   (preferred), or
#   - Wilcoxon_Phylum_ASE_vs_PDC_results.csv
#   Required columns (case-sensitive):
#     Taxa (or Taxonomy), avg(ASE), sd(ASE), avg(PDC), sd(PDC),
#     p.value, q.values, interval lower, interval upper,
#     log2FoldChange, Combine ES
#
# Output
#   ./results/figures/fig7c/Fig7c_Panel*.{pdf,svg,tiff}
#   Panels:
#     (1) Group means ± SD (percent)
#     (2) Mean difference (PDC − ASE) with 95% CI
#     (3) p-values with significance stars
#     (4) log2FoldChange (ASE vs PDC) lollipop
#
# Usage
#   R >= 4.0; packages: tidyverse, readxl, readr, ggsci, svglite
#   install.packages(c("tidyverse","readxl","readr","ggsci","svglite"))
#   From repo root: Rscript scripts/WilcoxonPhylumF.R
#
# Metadata
#   Version: 1.0
#   Author : Gao Han
#   Affil. : School of Biological Science and Medical Engineering, Southeast University
#   Date   : 2025-11-04
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
out_dir <- file.path("results", "figures", "fig7c")
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

in_xlsx <- file.path(in_dir, "Wilcoxon_Phylum_ASE_vs_PDC_results.xlsx")
in_csv  <- file.path(in_dir, "Wilcoxon_Phylum_ASE_vs_PDC_results.csv")

# Prefer curated Excel; fallback to CSV
if (file.exists(in_xlsx)) {
  raw <- read_xlsx(in_xlsx)
} else if (file.exists(in_csv)) {
  raw <- read_csv(in_csv, show_col_types = FALSE)
} else {
  stop("Input not found. Provide either:\n  ",
       in_xlsx, "\n  or\n  ", in_csv)
}

# ----------------------------
# Standardize columns and types
# ----------------------------
tax_col <- dplyr::case_when(
  "Taxa" %in% names(raw)     ~ "Taxa",
  "Taxonomy" %in% names(raw) ~ "Taxonomy",
  TRUE                       ~ ""
)
if (tax_col == "") stop("Missing taxonomy column: 'Taxa' or 'Taxonomy'.")

raw <- raw %>%
  rename(Taxon = all_of(tax_col)) %>%
  { names(.) <- trimws(names(.)); . } %>%
  mutate(Taxon = trimws(Taxon))

num_cols <- setdiff(names(raw), "Taxon")
df <- raw %>%
  mutate(across(all_of(num_cols), ~ suppressWarnings(as.numeric(.)))) %>%
  rename(
    avg_ASE = `avg(ASE)`, sd_ASE = `sd(ASE)`,
    avg_PDC = `avg(PDC)`, sd_PDC = `sd(PDC)`,
    p_value = `p.value`,  q_values = `q.values`,
    int_low = `interval lower`, int_up = `interval upper`,
    log2FC  = `log2FoldChange`,
    combine_es = `Combine ES`
  )

req <- c("avg_ASE","sd_ASE","avg_PDC","sd_PDC",
         "p_value","q_values","int_low","int_up","log2FC","combine_es")
if (!all(req %in% names(df))) {
  stop("Missing required columns: ", paste(setdiff(req, names(df)), collapse = ", "))
}
if (any(!sapply(df[req], is.numeric))) {
  bad <- names(df[req])[!sapply(df[req], is.numeric)]
  stop("Non-numeric values in columns: ", paste(bad, collapse = ", "),
       ". Clean the curated file (remove stray characters, dashes, blanks).")
}

# Keep the curated order (top row displayed on top)
df$Taxon <- factor(df$Taxon, levels = rev(unique(df$Taxon)))

# Display-friendly derivatives (percent units, differences)
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
save_fig <- function(p, base, width = 7.5,
                     height = max(6, 0.42 * nlevels(p$data$Taxon) + 2)) {
  ggsave(file.path(out_dir, paste0(base, ".pdf")),  p, width = width, height = height, dpi = 300)
  ggsave(file.path(out_dir, paste0(base, ".svg")),  p, width = width, height = height, dpi = 300)
  ggsave(file.path(out_dir, paste0(base, ".tiff")), p, width = width, height = height, dpi = 600,
         compression = "lzw", bg = "white")
}

theme_base <- theme_bw(base_size = 12) +
  theme(panel.grid.minor = element_blank(),
        axis.title.y = element_blank())

p_stars <- function(p) case_when(
  is.na(p)   ~ "",
  p < 0.001  ~ "***",
  p < 0.01   ~ "**",
  p < 0.05   ~ "*",
  TRUE       ~ ""
)

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
  geom_col(position = position_dodge(width = 0.7), width = 0.64) +
  geom_errorbar(aes(ymin = pmax(mean - sd, 0), ymax = mean + sd),
                position = position_dodge(width = 0.7), width = 0.25, linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(values = pal_groups, name = NULL) +
  labs(y = "Relative abundance (%)", title = "Proportions (%)") +
  theme_base +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", hjust = 0.5))

save_fig(p1, "Fig7c_Panel1_Proportions")

# ----------------------------
# Panel 2: Difference (PDC − ASE) with 95% CI
# ----------------------------
df_diff <- df %>% mutate(dir = ifelse(diff_pct >= 0, "pos", "neg"))

p2 <- ggplot(df_diff, aes(x = Taxon, y = diff_pct, color = dir)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_pointrange(aes(ymin = int_low_pct, ymax = int_up_pct), size = 0.4) +
  coord_flip() +
  scale_color_manual(values = col_dir, guide = "none") +
  labs(y = "Difference between proportions (%)", title = "95% confidence intervals") +
  theme_base +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

save_fig(p2, "Fig7c_Panel2_DiffCI")

# ----------------------------
# Panel 3: p-values with significance stars
# ----------------------------
df_p <- df %>% mutate(label = paste0(p_stars(p_value), "  ", formatC(p_value, format = "e", digits = 2)))

p3 <- ggplot(df_p, aes(x = 1, y = Taxon, label = label)) +
  geom_text(color = "firebrick", size = 4.0, fontface = "bold") +
  coord_cartesian(xlim = c(0.5, 1.5)) +
  labs(title = "P value") +
  theme_base +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

save_fig(p3, "Fig7c_Panel3_Pvalue")

# ----------------------------
# Panel 4: log2FoldChange (ASE vs PDC) lollipop
# ----------------------------
p4 <- ggplot(df, aes(x = Taxon, y = log2FC,
                     color = ifelse(log2FC >= 0, "pos", "neg"))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey55") +
  geom_segment(aes(xend = Taxon, y = 0, yend = log2FC), linewidth = 0.8, color = "grey40") +
  geom_point(size = 2.3, color = "grey10") +
  coord_flip() +
  scale_color_manual(values = col_dir, guide = "none") +
  labs(title = "log2FC (ASE vs PDC)", y = NULL) +
  theme_base +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

save_fig(p4, "Fig7c_Panel4_Log2FC")

message("Done. Figures written to: ", out_dir)
