# =====================================================================
# File:    pdc_gf_curve.R
# Title:   Growth-factor (GF) curves vs initial wet diameter (d0)
#          for particles traversing a PDC condenser under experimental
#          conditions (Maxwell r^2 + Fuchs–Sutugin + κ-Köhler).
#
# Purpose: Reproduce Fig. 2d in the paper. Computes GF(d0) for a set of
#          effective supersaturations S and saves:
#            - tidy table (Excel) to ./data/
#            - publication-ready figure (SVG/PDF/TIFF) to ./figures/theory/
#          Parameters can be edited in the block below to explore scenarios.
#
# Usage:   In the repo root:
#            install.packages(c("ggplot2","dplyr","tidyr","openxlsx","scales"))
#            source("pdc_gf_curve.R")
#
# Inputs:  No external input; all parameters defined below.
# Outputs: ./data/PDC_GF_vs_d0_Sweep.xlsx
#          ./figures/theory/GF_vs_d0_by_S.(svg|pdf|tiff)
#
# Version: 1.0
# Author : Gao Han
# Affil. : School of Biological Science and Medical Engineering,
#          Southeast University
# Date   : 2025-11-04
# Notes  : Theory summarized in Supplementary Text S1; parameter notes in
#          Supplementary Table S2. See LICENSE in the repository.
# =====================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(openxlsx)
  library(scales)
})

# ----------------------------
# Parameters (edit as needed)
# ----------------------------

## (A) Operating conditions (set the inlet state and residence time)
p_Pa        <- 101325       # ambient pressure (Pa)
Vdot_exh_Lm <- 10           # exhaled flow (L/min)
T_exh_C     <- 37           # exhaled temperature (°C)
RH_exh      <- 1.00         # exhaled RH (-)

Vdot_amb_Lm <- 15           # make-up air flow (L/min)
T_amb_C     <- 20           # make-up air temperature (°C)
RH_amb      <- 0.50         # make-up air RH (-)

Vdot_total_Lm <- 25         # total PDC draw (L/min)
D_tube_m      <- 3.8e-3     # condenser tube ID (m)
L_tube_m      <- 0.30       # condenser tube length (m)
T_out_C       <- 4          # outlet characteristic temperature (°C)

## (B) Model parameters (govern growth kinetics)
S_vec     <- c(1.03, 1.05, 1.10, 1.15, 1.30, 1.60, 1.80, 2.00)  # effective S
kappa     <- 0.30           # κ-Köhler hygroscopicity (-)
sigma_Npm <- 0.072          # surface tension (N/m)
lambda_m  <- 65e-9          # mean free path (m)
alpha_m   <- 1.0            # mass accommodation (-)

## (C) Thermophysical properties (assumed constant at an effective T)
rho_l_kgm3 <- 1000          # liquid water density (kg/m^3)
L_v_Jkg    <- 2.5e6         # latent heat (J/kg)
R_v_JkgK   <- 461.5         # water vapor gas constant (J/(kg·K))
k_air_WmK  <- 0.026         # air thermal conductivity (W/(m·K))
Dv_m2s     <- 2.3e-5        # vapor diffusivity in air (m^2/s)
Mw_kgmol   <- 0.018015      # water molar mass (kg/mol)

## (D) Diameter grid & time discretization
d0_min_um  <- 0.10          # x-axis start (μm)
d0_max_um  <- 5.00          # x-axis end (μm)
n_diam_pts <- 250           # number of x samples
n_steps    <- 400           # time steps (larger -> smoother, slower)

## (E) Repo-friendly output paths
dir_data <- file.path("data")
dir_fig  <- file.path("figures", "theory")
if (!dir.exists(dir_data)) dir.create(dir_data, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(dir_fig))  dir.create(dir_fig,  recursive = TRUE, showWarnings = FALSE)

# ----------------------------
# Psychrometrics & model pieces
# ----------------------------

es_Pa <- function(T_C) 0.61094 * exp(17.625 * T_C / (T_C + 243.04)) * 1000  # Pa

w_from_e   <- function(e_Pa, p_Pa) 0.62198 * e_Pa / (p_Pa - e_Pa)
w_from_TRH <- function(T_C, RH, p_Pa) w_from_e(RH * es_Pa(T_C), p_Pa)

h_kJkgdry <- function(T_C, w) 1.006 * T_C + w * (2501 + 1.86 * T_C)

beta_FS <- function(dp_m, lambda_m, alpha_m = 1.0) {
  Kn <- 2 * lambda_m / dp_m
  (1 + Kn) / (1 + 1.71 * Kn + 1.33 * Kn^2)     # α = 1 form
}

S_eq_kohler <- function(r_m, kappa, rd_m, T_K, sigma_Npm, rho_l_kgm3, Mw_kgmol) {
  A  <- 2 * sigma_Npm * Mw_kgmol / (rho_l_kgm3 * 8.314462618 * T_K)  # Kelvin term
  aw <- 1 / (1 + kappa * (rd_m^3) / (r_m^3 - rd_m^3))                # activity term
  aw * exp(A / r_m)
}

solve_rd_from_r0 <- function(r0_m, S_in, kappa, T_K, sigma_Npm, rho_l_kgm3, Mw_kgmol) {
  f <- function(rd) S_eq_kohler(r0_m, kappa, rd, T_K, sigma_Npm, rho_l_kgm3, Mw_kgmol) - S_in
  uniroot(f, lower = 1e-10, upper = r0_m * 0.999)$root
}

K_term <- function(L_v_Jkg, rho_l_kgm3, k_air_WmK, R_v_JkgK, T_K) {
  L_v_Jkg^2 * rho_l_kgm3 / (k_air_WmK * R_v_JkgK * T_K^2)
}
D_term <- function(rho_l_kgm3, R_v_JkgK, T_K, es_Pa_val, Dv_m2s) {
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
# Growth solver: d0 (wet, at RH_in) -> GF(d0; S)
# ----------------------------

grow_one <- function(d0_m, S_eff) {
  r0 <- d0_m / 2
  rd <- solve_rd_from_r0(r0, RH_in, kappa, T_eff_K, sigma_Npm, rho_l_kgm3, Mw_kgmol)
  
  dt <- t_res / n_steps
  r  <- r0
  for (i in 1:n_steps) {
    dp   <- 2 * r
    bfs  <- beta_FS(dp, lambda_m, alpha_m)
    Def  <- D0 / bfs
    Seq  <- S_eq_kohler(r, kappa, rd, T_eff_K, sigma_Npm, rho_l_kgm3, Mw_kgmol)
    drive <- max(0, S_eff - Seq)
    dr2dt <- 2 * drive / (K_val + Def)
    r2    <- r^2 + dr2dt * dt
    r     <- sqrt(max(r2, rd^2 * 1.0000001))
  }
  d_final <- 2 * r
  list(d0 = d0_m, d_final = d_final, rd = rd, S = S_eff)
}

d0_grid_m <- seq(d0_min_um, d0_max_um, length.out = n_diam_pts) * 1e-6

res_list <- vector("list", length(S_vec) * length(d0_grid_m))
idx <- 1
for (S_eff in S_vec) {
  for (d0 in d0_grid_m) {
    res_list[[idx]] <- grow_one(d0, S_eff)
    idx <- idx + 1
  }
}

res_df <- do.call(rbind, lapply(res_list, as.data.frame)) %>%
  mutate(
    d0_um      = d0 * 1e6,
    d_final_um = d_final * 1e6,
    GF         = d_final / d0,
    S_label    = paste0("S = ", sprintf("%.2f", S))
  ) %>%
  select(S_label, S, d0_um, d_final_um, GF)

# ----------------------------
# Save data (Excel)
# ----------------------------
xlsx_path <- file.path(dir_data, "PDC_GF_vs_d0_Sweep.xlsx")
wb <- createWorkbook()
addWorksheet(wb, "GF_vs_d0")
writeData(wb, "GF_vs_d0", res_df)
saveWorkbook(wb, xlsx_path, overwrite = TRUE)

# ----------------------------
# Plot (journal style; blue→green gradient)
# ----------------------------
col_fun <- colorRampPalette(c("#1f78b4", "#2ca25f"))
cols <- col_fun(length(unique(res_df$S_label)))

p <- ggplot(res_df, aes(x = d0_um, y = GF, color = S_label)) +
  geom_line(linewidth = 1.1) +
  labs(
    title = "Particle Growth Factor in PDC vs Initial Wet Diameter",
    subtitle = sprintf("Inlet mix: T_in ≈ %.1f °C, RH_in ≈ %.0f%%; residence time ≈ %.2f ms; κ = %.2f",
                       T_in_C, RH_in*100, t_res*1e3, kappa),
    x = "Initial wet diameter d0 (μm)",
    y = "Growth factor GF = d_final / d0",
    color = "Supersaturation"
  ) +
  scale_color_manual(values = cols) +
  coord_cartesian(xlim = c(d0_min_um, d0_max_um)) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(size = 11),
    legend.position = "right",
    axis.title = element_text(face = "bold")
  )

fig_svg  <- file.path(dir_fig, "GF_vs_d0_by_S.svg")
fig_pdf  <- file.path(dir_fig, "GF_vs_d0_by_S.pdf")
fig_tif  <- file.path(dir_fig, "GF_vs_d0_by_S.tiff")

ggsave(fig_svg, p, width = 7, height = 5, units = "in")
ggsave(fig_pdf, p, width = 7, height = 5, units = "in")
ggsave(fig_tif, p, width = 7, height = 5, units = "in", dpi = 600, compression = "lzw")

cat("\n== Done ==\n")
cat(sprintf("Data  : %s\n", xlsx_path))
cat(sprintf("Figure: %s\n       %s\n       %s\n", fig_svg, fig_pdf, fig_tif))
