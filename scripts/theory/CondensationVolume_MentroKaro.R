# =====================================================================
# Title:    Monte Carlo simulation of PDC-sampler condensate generation
# Purpose:  Reproduce Fig. 2a (0–30 min cumulative condensate + 95% CI)
# Method:   Simple humidity-mixing / saturation-approach model;
#           stochastic RH_exh and collection efficiency (eta).
# Theory:   Supplementary Text S1; parameters in Supplementary Table S2.
#
# Usage:
#   - Put this script in repo root (e.g., scripts/ or root).
#   - Outputs are written to relative folders: ./data/ and ./figures/
#   - Run in R: source("pdc_montecarlo_condensate.R")
#     (Optional) set.seed(20251104) for exact reproducibility.
#
# Dependencies: ggplot2, dplyr (install.packages(c("ggplot2","dplyr")))
#
# Outputs:
#   data/PDC_MonteCarlo_cumulative_time_series.csv
#   data/PDC_MonteCarlo_summary_by_time.csv
#   figures/PDC_MonteCarlo_CumVolume_vs_Time.{svg,pdf,tiff}
#
# Metadata:
#   Version: 1.0
#   Author:  Gao Han
#   Affil.:  School of Biological Science and Medical Engineering,
#            Southeast University
#   Date:    2025-11-04
#   Notes:   No personal data; simulation only.
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# --------------------------
# Fixed physical parameters
# --------------------------
p_kPa        <- 101.325   # total pressure (kPa)
R_kPaL_molK  <- 8.314     # gas constant (kPa·L/(mol·K))
rho_w_kg_L   <- 1.0       # water density (kg/L)
M_da_kg_mol  <- 0.02897   # dry air molar mass (kg/mol)

# --------------------------
# Operating conditions
# --------------------------
Vdot_exh_Lpm <- 10        # exhaled flow (L/min)
T_exh_C      <- 28        # exhaled gas temperature (°C)
Vdot_amb_Lpm <- 15        # make-up air flow (L/min)
T_amb_C      <- 20        # make-up air temperature (°C)
RH_amb       <- 0.50      # make-up air relative humidity (-)
Vdot_PDC_Lpm <- 25        # total draw by PDC (L/min) [Vdot_exh + Vdot_amb]
T_out_C      <- 4         # outlet gas temperature after condensation (°C)

# --------------------------
# Monte Carlo settings
# --------------------------
# set.seed(20251104)      # uncomment for reproducibility
n_runs       <- 15        # number of simulated trajectories
t_end_min    <- 30        # total simulated time (min)
dt_min       <- 1         # time step (min)
vary_per_min <- TRUE      # TRUE: RH_exh and eta vary every minute

# Stochastic ranges (uniform)
RH_exh_min <- 0.60        # exhaled RH lower bound (-)
RH_exh_max <- 1.00        # exhaled RH upper bound (-)
eta_min    <- 0.20        # collection efficiency lower bound (-)
eta_max    <- 0.90        # collection efficiency upper bound (-)

# --------------------------
# I/O (GitHub-friendly paths)
# --------------------------
dir_data <- file.path("data")
dir_fig  <- file.path("figures")
if (!dir.exists(dir_data)) dir.create(dir_data, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(dir_fig))  dir.create(dir_fig,  recursive = TRUE, showWarnings = FALSE)

# --------------------------
# Saturation vapor pressure (kPa), Magnus formula
# --------------------------
es_kPa <- function(T_C) {
  0.61094 * exp(17.625 * T_C / (T_C + 243.04))
}

# ------------------------------------------------------------
# One-step (1 min) condensate flow rate
# Returns:
#   Vdot_upper_Lpm : theoretical upper limit (no efficiency)
#   Vdot_eta_Lpm   : effective value with efficiency eta
# ------------------------------------------------------------
step_Vdot <- function(RH_exh, eta_effect) {
  # vapor partial pressures
  e_exh_kPa <- RH_exh * es_kPa(T_exh_C)
  e_amb_kPa <- RH_amb  * es_kPa(T_amb_C)
  
  # humidity ratios (kg water/kg dry air)
  w_exh <- 0.62198 * e_exh_kPa / (p_kPa - e_exh_kPa)
  w_amb <- 0.62198 * e_amb_kPa / (p_kPa - e_amb_kPa)
  
  # total molar flow rates (mol/min)
  n_tot_exh <- p_kPa * Vdot_exh_Lpm / (R_kPaL_molK * (T_exh_C + 273.15))
  n_tot_amb <- p_kPa * Vdot_amb_Lpm / (R_kPaL_molK * (T_amb_C + 273.15))
  
  # dry-air molar flows (mol/min)
  xw_exh <- e_exh_kPa / p_kPa
  xw_amb <- e_amb_kPa / p_kPa
  n_dry_exh <- n_tot_exh * (1 - xw_exh)
  n_dry_amb <- n_tot_amb * (1 - xw_amb)
  
  # dry-air mass flows (kg/min)
  mdot_dry_exh <- n_dry_exh * M_da_kg_mol
  mdot_dry_amb <- n_dry_amb * M_da_kg_mol
  mdot_dry_mix <- mdot_dry_exh + mdot_dry_amb
  
  # mixed humidity and saturated humidity at outlet temperature
  mdot_w_mix <- w_exh * mdot_dry_exh + w_amb * mdot_dry_amb
  w_mix      <- mdot_w_mix / mdot_dry_mix
  w_out_star <- 0.62198 * es_kPa(T_out_C) / (p_kPa - es_kPa(T_out_C))
  
  # condensate mass and volume rates
  delta_w <- max(0, w_mix - w_out_star)
  mdot_cond_upper <- delta_w * mdot_dry_mix            # kg/min
  Vdot_upper_Lpm  <- mdot_cond_upper / rho_w_kg_L      # L/min
  Vdot_eta_Lpm    <- eta_effect * Vdot_upper_Lpm       # L/min
  
  c(Vdot_upper_Lpm = Vdot_upper_Lpm, Vdot_eta_Lpm = Vdot_eta_Lpm)
}

# ------------------------------------------------------------
# Simulate one trajectory (returns cumulative volume by minute)
# ------------------------------------------------------------
simulate_one_run <- function(run_id) {
  times   <- seq(0, t_end_min, by = dt_min)
  n_steps <- length(times) - 1
  
  if (vary_per_min) {
    RH_series  <- runif(n_steps, RH_exh_min, RH_exh_max)
    eta_series <- runif(n_steps, eta_min,    eta_max)
  } else {
    RH_val  <- runif(1, RH_exh_min, RH_exh_max)
    eta_val <- runif(1, eta_min,    eta_max)
    RH_series  <- rep(RH_val,  n_steps)
    eta_series <- rep(eta_val, n_steps)
  }
  
  V_upper_cum <- numeric(n_steps + 1)  # mL
  V_eta_cum   <- numeric(n_steps + 1)
  
  for (k in 1:n_steps) {
    rates <- step_Vdot(RH_exh = RH_series[k], eta_effect = eta_series[k]) # L/min
    V_upper_cum[k + 1] <- V_upper_cum[k] + rates["Vdot_upper_Lpm"] * dt_min * 1000  # mL
    V_eta_cum[k + 1]   <- V_eta_cum[k]   + rates["Vdot_eta_Lpm"]   * dt_min * 1000  # mL
  }
  
  data.frame(
    run_id          = run_id,
    time_min        = times,
    V_cum_upper_mL  = V_upper_cum,
    V_cum_eta_mL    = V_eta_cum
  )
}

# --------------------------
# Monte Carlo execution
# --------------------------
mc_df <- dplyr::bind_rows(lapply(1:n_runs, simulate_one_run))

# Primary dataset for plotting (effective condensate only)
plot_df <- mc_df %>% select(run_id, time_min, V_cum_eta_mL)

# Time-wise summary statistics (for reference)
summary_df <- plot_df %>%
  group_by(time_min) %>%
  summarise(
    mean_mL = mean(V_cum_eta_mL),
    sd_mL   = sd(V_cum_eta_mL),
    q025    = quantile(V_cum_eta_mL, 0.025),
    q975    = quantile(V_cum_eta_mL, 0.975),
    .groups = "drop"
  )

# --------------------------
# Visualization
# --------------------------
p <- ggplot(plot_df, aes(x = time_min, y = V_cum_eta_mL)) +
  geom_point(alpha = 0.25, size = 1.0) +
  stat_smooth(method = "lm", se = TRUE, linewidth = 1.0) +
  labs(
    x = "Time (min)",
    y = "Cumulative condensate (mL)",
    title = "PDC sampler: Monte Carlo of condensate accumulation (0–30 min)",
    subtitle = sprintf("RH_exh ~ U[%.0f%%, %.0f%%]; eta ~ U[%.0f%%, %.0f%%]; n = %d; vary_per_min = %s",
                       RH_exh_min*100, RH_exh_max*100, eta_min*100, eta_max*100,
                       n_runs, ifelse(vary_per_min, "TRUE", "FALSE")),
    caption = "Fit: linear regression with 95% CI; points are per-minute cumulative values for each run."
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11)
  )

print(p)

# --------------------------
# Save data and figures
# --------------------------
file_data_detail  <- file.path(dir_data, "PDC_MonteCarlo_cumulative_time_series.csv")
file_data_summary <- file.path(dir_data, "PDC_MonteCarlo_summary_by_time.csv")
write.csv(mc_df,      file_data_detail,  row.names = FALSE)
write.csv(summary_df, file_data_summary, row.names = FALSE)

file_svg <- file.path(dir_fig, "PDC_MonteCarlo_CumVolume_vs_Time.svg")
file_pdf <- file.path(dir_fig, "PDC_MonteCarlo_CumVolume_vs_Time.pdf")   # bug fixed
file_tif <- file.path(dir_fig, "PDC_MonteCarlo_CumVolume_vs_Time.tiff")

ggsave(file_svg, p, width = 7, height = 5, units = "in")
ggsave(file_pdf, p, width = 7, height = 5, units = "in")
ggsave(file_tif, p, width = 7, height = 5, units = "in", dpi = 600, compression = "lzw")

# --------------------------
# Optional: add dashed mean of the theoretical upper bound
# --------------------------
# lines_df <- mc_df %>% group_by(time_min) %>%
#   summarise(mean_upper = mean(V_cum_upper_mL), .groups = "drop")
# p2 <- p + geom_line(data = lines_df, aes(x = time_min, y = mean_upper),
#                     linewidth = 0.9, linetype = "dashed")
# ggsave(file.path(dir_fig, "PDC_MonteCarlo_CumVolume_vs_Time_withUpper.pdf"),
#        p2, width = 7, height = 5, units = "in")
