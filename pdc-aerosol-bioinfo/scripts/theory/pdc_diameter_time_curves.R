# =====================================================================
# File:    pdc_diameter_time_curves.R
# Title:   Diameter–time growth trajectories at S = 1.60 in a PDC condenser
# Model:   Maxwell r^2 growth + Fuchs–Sutugin correction + κ-Köhler equilibrium
#
# Purpose: Reproduce Fig. 2e. Compute d(t) for multiple initial wet diameters
#          (at inlet RH) under the experimental PDC operating condition S = 1.60.
#          Save tidy results (Excel) and a publication-grade figure.
#
# Usage:   In the repository root:
#            install.packages(c("ggplot2","dplyr","tidyr","openxlsx"))
#            source("pdc_diameter_time_curves.R")
#
# Inputs:  None (all parameters defined below)
# Outputs: ./data/PDC_diameter_time_traj_S1p60.xlsx
#          ./figures/theory/Diameter_vs_Time_S1p60.(svg|pdf|tiff)
#
# Version: 1.0
# Author : Gao Han
# Affil. : School of Biological Science and Medical Engineering, Southeast University
# Date   : 2025-11-04
# Notes  : Theory summarized in Supplementary Text S1; parameter settings in
#          Supplementary Table S2. See LICENSE in the repository.
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
})

# ----------------------------
# Parameters (edit as needed)
# ----------------------------

## (A) Operating conditions → inlet state & residence time
p_Pa        <- 101325      # ambient pressure (Pa)
Vdot_exh_Lm <- 10          # exhaled flow (L/min)
T_exh_C     <- 37          # exhaled temperature (°C)
RH_exh      <- 1.00        # exhaled RH (-)

Vdot_amb_Lm <- 15          # make-up air flow (L/min)
T_amb_C     <- 20          # make-up air temperature (°C)
RH_amb      <- 0.50        # make-up air RH (-)

Vdot_total_Lm <- 25        # total PDC draw (L/min)
D_tube_m      <- 3.8e-3    # condenser tube inner diameter (m)
L_tube_m      <- 0.30      # condenser tube length (m)
T_out_C       <- 4         # outlet characteristic temperature (°C)

## (B) Growth/engineering parameters
S_eff     <- 1.60          # effective supersaturation (-) for this figure
kappa     <- 0.30          # κ-Köhler hygroscopicity (-)
sigma_Npm <- 0.072         # surface tension (N/m)
lambda_m  <- 65e-9         # mean free path (m)
alpha_m   <- 1.0           # mass accommodation (-)

## (C) Thermophysical properties (assumed constant at effective temperature)
rho_l_kgm3 <- 1000         # liquid water density (kg/m^3)
L_v_Jkg    <- 2.5e6        # latent heat (J/kg)
R_v_JkgK   <- 461.5        # water vapor gas constant (J/(kg·K))
k_air_WmK  <- 0.026        # air thermal conductivity (W/(m·K))
Dv_m2s     <- 2.3e-5       # vapor diffusivity in air (m^2/s)
Mw_kgmol   <- 0.018015     # water molar mass (kg/mol)

## (D) Time axis and initial wet diameters
t_end_ms <- 15             # max time for trajectories (ms)
n_t_pts  <- 30             # number of time points

# Initial wet diameters (m) at inlet RH; edit to any set you need
d0_list_m <- c(100e-9, 200e-9, 300e-9, 500e-9, 700e-9, 1e-6, 1.5e-6, 2e-6)

## (E) Repo-friendly output paths
dir_data <- file.path("data")
dir_fig  <- file.path("figures", "theory")
if (!dir.exists(dir_data)) dir.create(dir_data, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(dir_fig))  dir.create(dir_fig,  recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Psychrometrics & model pieces
# ----------------------------

es_Pa <- function(T_C) 0.61094 * exp(17.625 * T_C / (T_C + 243.04)) * 1000      # saturation vapor pressure (Pa)
w_from_e   <- function(e_Pa, p_Pa) 0.62198 * e_Pa / (p_Pa - e_Pa)               # humidity ratio (kg/kg_dry)
w_from_TRH <- function(T_C, RH, p_Pa) w_from_e(RH * es_Pa(T_C), p_Pa)
h_kJkgdry  <- function(T_C, w) 1.006 * T_C + w * (2501 + 1.86 * T_C)            # moist air enthalpy (kJ/kg_dry)

beta_FS <- function(dp_m, lambda_m, alpha_m = 1.0) {                            # Fuchs–Sutugin correction
  Kn <- 2 * lambda_m / dp_m
  (1 + Kn) / (1 + 1.71 * Kn + 1.33 * Kn^2)  # α = 1 form
}

S_eq_kohler <- function(r_m, kappa, rd_m, T_K, sigma_Npm, rho_l_kgm3, Mw_kgmol) {
  # κ-Köhler equilibrium supersaturation S_eq = a_w * exp(A/r)
  A  <- 2 * sigma_Npm * Mw_kgmol / (rho_l_kgm3 * 8.314462618 * T_K)  # Kelvin term
  aw <- 1 / (1 + kappa * (rd_m^3) / (r_m^3 - rd_m^3))                # activity term
  aw * exp(A / r_m)
}

solve_rd_from_r0 <- function(r0_m, S_in, kappa, T_K, sigma_Npm, rho_l_kgm3, Mw_kgmol) {
  # invert dry radius rd given initial wet radius r0 at inlet supersaturation S_in
  f <- function(rd) S_eq_kohler(r0_m, kappa, rd, T_K, sigma_Npm, rho_l_kgm3, Mw_kgmol) - S_in
  uniroot(f, lower = 1e-10, upper = r0_m * 0.999)$root
}

K_term <- function(L_v_Jkg, rho_l_kgm3, k_air_WmK, R_v_JkgK, T_K) {
  # thermal resistance term in Maxwell r^2 growth
  L_v_Jkg^2 * rho_l_kgm3 / (k_air_WmK * R_v_JkgK * T_K^2)
}
D_term <- function(rho_l_kgm3, R_v_JkgK, T_K, es_Pa_val, Dv_m2s) {
  # diffusive resistance term in Maxwell r^2 growth
  rho_l_kgm3 * R_v_JkgK * T_K / (es_Pa_val * Dv_m2s)
}

# ----------------------------
# Inlet mixed state & residence time
# ----------------------------

V_exh_m3s <- Vdot_exh_Lm   / 1000 / 60
V_amb_m3s <- Vdot_amb_Lm   / 1000 / 60
V_tot_m3s <- Vdot_total_Lm / 1000 / 60

n_tot_exh <- p_Pa * V_exh_m3s / (8.314462618 * (T_exh_C + 273.15))
n_tot_amb <- p_Pa * V_amb_m3s / (8.314462618 * (T_amb_C + 273.15))
xw_exh    <- (RH_exh * es_Pa(T_exh_C)) / p_Pa
xw_amb    <- (RH_amb * es_Pa(T_amb_C)) / p_Pa
n_dry_exh <- n_tot_exh * (1 - xw_exh)
n_dry_amb <- n_tot_amb * (1 - xw_amb)

M_da <- 0.02897  # kg/mol
mdot_dry_exh <- n_dry_exh * M_da
mdot_dry_amb <- n_dry_amb * M_da
mdot_dry_mix <- mdot_dry_exh + mdot_dry_amb

w_exh <- w_from_TRH(T_exh_C, RH_exh, p_Pa)
w_amb <- w_from_TRH(T_amb_C, RH_amb, p_Pa)
h_exh <- h_kJkgdry(T_exh_C, w_exh)
h_amb <- h_kJkgdry(T_amb_C, w_amb)

w_in <- (w_exh * mdot_dry_exh + w_amb * mdot_dry_amb) / mdot_dry_mix
h_in <- (h_exh * mdot_dry_exh + h_amb * mdot_dry_amb) / mdot_dry_mix
T_in_C <- uniroot(function(Tc) h_kJkgdry(Tc, w_in) - h_in, c(-20, 50))$root

e_in  <- w_in * p_Pa / (0.62198 + w_in)
RH_in <- e_in / es_Pa(T_in_C)

U_ms  <- V_tot_m3s / (pi * (D_tube_m/2)^2)
t_res <- L_tube_m / U_ms

T_eff_C <- (T_in_C + T_out_C) / 2
T_eff_K <- T_eff_C + 273.15

K_val <- K_term(L_v_Jkg, rho_l_kgm3, k_air_WmK, R_v_JkgK, T_eff_K)
D0    <- D_term(rho_l_kgm3, R_v_JkgK, T_eff_K, es_Pa(T_eff_C), Dv_m2s)

# ----------------------------
# Trajectory solver: d0 (wet, at RH_in) → rd → integrate to get d(t)
# ----------------------------

traj_one <- function(d0_m, S_eff, t_end_ms, n_t_pts) {
  r0 <- d0_m / 2
  rd <- solve_rd_from_r0(r0, RH_in, kappa, T_eff_K, sigma_Npm, rho_l_kgm3, Mw_kgmol)
  
  t_grid_s <- seq(0, t_end_ms/1000, length.out = n_t_pts)
  dt <- diff(t_grid_s); dt <- c(dt[1], dt)  # use first step for the initial interval
  
  r <- r0
  d_out <- numeric(length(t_grid_s))
  d_out[1] <- 2 * r
  
  for (i in 2:length(t_grid_s)) {
    dp   <- 2 * r
    bfs  <- beta_FS(dp, lambda_m, alpha_m)
    Def  <- D0 / bfs
    Seq  <- S_eq_kohler(r, kappa, rd, T_eff_K, sigma_Npm, rho_l_kgm3, Mw_kgmol)
    drive <- max(0, S_eff - Seq)                 # no growth if drive ≤ 0
    dr2dt <- 2 * drive / (K_val + Def)
    r2    <- r^2 + dr2dt * dt[i]
    r     <- sqrt(max(r2, rd^2 * 1.0000001))     # keep r > rd
    d_out[i] <- 2 * r
  }
  
  data.frame(
    time_ms = t_grid_s * 1000,
    d_um    = d_out * 1e6,
    d0_um   = d0_m * 1e6
  )
}

traj_list <- lapply(d0_list_m, function(d0) traj_one(d0, S_eff, t_end_ms, n_t_pts))
traj_df   <- bind_rows(traj_list) %>%
  mutate(d0_label = paste0(
    "d0 = ",
    ifelse(d0_um < 1, sprintf("%.0f nm", d0_um*1000), sprintf("%.2f µm", d0_um))
  ))

# ----------------------------
# Save data (Excel)
# ----------------------------
xlsx_path <- file.path(dir_data, "PDC_diameter_time_traj_S1p60.xlsx")
wb <- createWorkbook()
addWorksheet(wb, "d_vs_t")
writeData(wb, "d_vs_t", traj_df)
saveWorkbook(wb, xlsx_path, overwrite = TRUE)

# ----------------------------
# Plot (blue gradient; dashed residence-time line)
# ----------------------------

col_fun <- colorRampPalette(c("#c6dbef", "#6baed6", "#2171b5"))  # light→dark blue
cols <- col_fun(length(unique(traj_df$d0_label)))

p <- ggplot(traj_df, aes(x = time_ms, y = d_um, color = d0_label)) +
  geom_line(linewidth = 1.1) +
  geom_vline(xintercept = t_res*1000, linetype = "dashed") +
  annotate(
    "text", x = t_res*1000, y = max(traj_df$d_um)*0.97,
    label = sprintf("Residence time ≈ %.2f ms", t_res*1000),
    angle = 90, vjust = -0.5, size = 3.5
  ) +
  scale_color_manual(values = cols) +
  labs(
    title = "Hygroscopic growth trajectories in PDC (S = 1.60)",
    subtitle = sprintf("Inlet mix: T_in ≈ %.1f °C, RH_in ≈ %.0f%%; effective T ≈ %.1f °C",
                       T_in_C, RH_in*100, T_eff_C),
    x = "Time (ms)",
    y = "Diameter d (µm)",
    color = "Initial diameter"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title    = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    legend.position = "right",
    axis.title = element_text(face = "bold")
  )

fig_svg  <- file.path(dir_fig, "Diameter_vs_Time_S1p60.svg")
fig_pdf  <- file.path(dir_fig, "Diameter_vs_Time_S1p60.pdf")
fig_tif  <- file.path(dir_fig, "Diameter_vs_Time_S1p60.tiff")

ggsave(filename = fig_svg, plot = p, width = 7, height = 5, units = "in")
ggsave(filename = fig_pdf, plot = p, width = 7, height = 5, units = "in")
ggsave(filename = fig_tif, plot = p, width = 7, height = 5, units = "in", dpi = 600, compression = "lzw")

cat("\n== Done ==\n")
cat(sprintf("Data  : %s\n", xlsx_path))
cat(sprintf("Figure: %s\n       %s\n       %s\n", fig_svg, fig_pdf, fig_tif))
