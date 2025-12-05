# =====================================================================
# File:    pdc_particle_growth.R
# Title:   Particle-size growth inside a PDC sampler
# Purpose: Compute hygroscopic growth of aerosol particles traversing the
#          PDC condenser under experimental conditions using the main-text
#          theory (Maxwell r^2 growth with Fuchs–Sutugin correction and
#          κ-Köhler equilibrium). This script evaluates the growth for
#          user-defined operating/thermophysical parameters and prints
#          results. No plotting is produced.
#
# How to run:
#   - R (>= 4.0). Base R only (no extra packages required).
#   - Place this script in the repo root and run:
#       source("pdc_particle_growth.R")
#   - Edit the parameter block to test different scenarios.
#
# Optional output:
#   - Set SAVE_RESULTS = TRUE to write a CSV into ./data/
#
# Metadata:
#   Version : 1.0
#   Author  : Gao Han
#   Affil.  : School of Biological Science and Medical Engineering,
#             Southeast University
#   Date    : 2025-11-04
#   Links   : Theory in Supplementary Text S1; parameter notes in Suppl. Table S2
# =====================================================================

# ----------------------------
# Parameters (edit as needed)
# ----------------------------

## (A) Operating conditions
p_Pa        <- 101325      # ambient pressure (Pa)
Vdot_exh_Lm <- 10          # exhaled flow (L/min)
T_exh_C     <- 37          # exhaled-gas temperature (°C)
RH_exh      <- 1.00        # exhaled-gas relative humidity (-)

Vdot_amb_Lm <- 15          # make-up air flow (L/min)
T_amb_C     <- 20          # make-up air temperature (°C)
RH_amb      <- 0.50        # make-up air RH (-)

Vdot_total_Lm <- 25        # total PDC draw (L/min)
D_tube_m      <- 3.8e-3    # condenser tube inner diameter (m)
L_tube_m      <- 0.50      # condenser tube length (m)
T_out_C       <- 4         # outlet characteristic temperature (°C)

## (B) Growth-model parameters (govern the kinetics)
S_eff     <- 1.10          # effective supersaturation S (-), e.g., 1.05 = 5%
kappa     <- 0.30          # κ-Köhler hygroscopicity (-)
sigma_Npm <- 0.072         # surface tension (N/m)
lambda_m  <- 65e-9         # mean free path (m)
alpha_m   <- 1.0           # mass accommodation coefficient (-)

## (C) Thermophysical properties (treated as constants at an effective T)
T_eff_C <- ((T_exh_C*Vdot_exh_Lm + T_amb_C*Vdot_amb_Lm) /
              (Vdot_exh_Lm + Vdot_amb_Lm) + T_out_C)/2
T_eff_K <- T_eff_C + 273.15

rho_l_kgm3 <- 1000         # liquid-water density (kg/m^3)
L_v_Jkg    <- 2.5e6        # latent heat (J/kg)
R_v_JkgK   <- 461.5        # water vapor gas constant (J/(kg·K))
k_air_WmK  <- 0.026        # air thermal conductivity (W/(m·K))
Dv_m2s     <- 2.3e-5       # vapor diffusivity in air (m^2/s)
Mw_kgmol   <- 0.018015     # water molar mass (kg/mol)

## (D) Time discretization
n_steps <- 400             # number of steps; larger -> finer but slower

## (E) Initial wet diameters at the inlet RH (m)
d0_list_m <- c(100e-9, 500e-9, 1e-6)

## (F) Optional CSV output (GitHub-friendly path)
SAVE_RESULTS <- TRUE
out_dir      <- "data"
out_csv      <- file.path(out_dir, "PDC_ParticleGrowth_results.csv")

# ----------------------------
# Helper functions
# ----------------------------

# Saturation vapor pressure (Pa), Magnus formula
es_Pa <- function(T_C) {
  0.61094 * exp(17.625 * T_C / (T_C + 243.04)) * 1000
}

# Humidity ratio w (kg/kg dry air) from vapor pressure e (Pa) and total p (Pa)
w_from_e <- function(e_Pa, p_Pa) {
  0.62198 * e_Pa / (p_Pa - e_Pa)
}

# Humidity ratio from T (°C) and RH (-)
w_from_TRH <- function(T_C, RH, p_Pa) {
  e <- RH * es_Pa(T_C)
  w_from_e(e, p_Pa)
}

# Moist-air specific enthalpy (kJ/kg dry air): h = 1.006 T + w (2501 + 1.86 T)
h_kJkgdry <- function(T_C, w) {
  1.006 * T_C + w * (2501 + 1.86 * T_C)
}

# Fuchs–Sutugin correction factor (mass transfer correction)
beta_FS <- function(dp_m, lambda_m, alpha_m = 1.0) {
  Kn <- 2 * lambda_m / dp_m
  (1 + Kn) / (1 + 1.71 * Kn + 1.33 * Kn^2)  # α = 1 form
}

# κ-Köhler equilibrium saturation ratio S_eq(r)
S_eq_kohler <- function(r_m, kappa, rd_m, T_K, sigma_Npm, rho_l_kgm3, Mw_kgmol, p_Pa) {
  A  <- 2 * sigma_Npm * Mw_kgmol / (rho_l_kgm3 * 8.314462618 * T_K)  # Kelvin term
  aw <- 1 / (1 + kappa * (rd_m^3) / (r_m^3 - rd_m^3))
  aw * exp(A / r_m)
}

# Infer dry radius rd from inlet wet radius r0 and inlet saturation S_in
solve_rd_from_r0 <- function(r0_m, S_in, kappa, T_K, sigma_Npm, rho_l_kgm3, Mw_kgmol, p_Pa) {
  f <- function(rd) {
    S_eq_kohler(r0_m, kappa, rd, T_K, sigma_Npm, rho_l_kgm3, Mw_kgmol, p_Pa) - S_in
  }
  uniroot(f, lower = 1e-10, upper = r0_m * 0.999)$root
}

# Thermal and diffusive resistance terms (s/m)
K_term <- function(L_v_Jkg, rho_l_kgm3, k_air_WmK, R_v_JkgK, T_K) {
  L_v_Jkg^2 * rho_l_kgm3 / (k_air_WmK * R_v_JkgK * T_K^2)
}
D_term <- function(rho_l_kgm3, R_v_JkgK, T_K, es_Pa_val, Dv_m2s) {
  rho_l_kgm3 * R_v_JkgK * T_K / (es_Pa_val * Dv_m2s)
}

# ----------------------------
# Inlet mixed state and residence time
# ----------------------------

# Volumetric flows (m^3/s)
V_exh_m3s <- Vdot_exh_Lm   / 1000 / 60
V_amb_m3s <- Vdot_amb_Lm   / 1000 / 60
V_tot_m3s <- Vdot_total_Lm / 1000 / 60

# Total and dry-air molar flows (mol/s)
n_tot_exh <- p_Pa * V_exh_m3s / (8.314462618 * (T_exh_C + 273.15))
n_tot_amb <- p_Pa * V_amb_m3s / (8.314462618 * (T_amb_C + 273.15))
xw_exh    <- (RH_exh * es_Pa(T_exh_C)) / p_Pa
xw_amb    <- (RH_amb * es_Pa(T_amb_C)) / p_Pa
n_dry_exh <- n_tot_exh * (1 - xw_exh)
n_dry_amb <- n_tot_amb * (1 - xw_amb)

M_da_kgmol <- 0.02897
mdot_dry_exh <- n_dry_exh * M_da_kgmol
mdot_dry_amb <- n_dry_amb * M_da_kgmol
mdot_dry_mix <- mdot_dry_exh + mdot_dry_amb

# Inlet humidity ratio and enthalpy (dry-air-weighted mixing)
w_exh <- w_from_TRH(T_exh_C, RH_exh, p_Pa)
w_amb <- w_from_TRH(T_amb_C, RH_amb, p_Pa)
h_exh <- h_kJkgdry(T_exh_C, w_exh)
h_amb <- h_kJkgdry(T_amb_C, w_amb)

w_in <- (w_exh*mdot_dry_exh + w_amb*mdot_dry_amb) / mdot_dry_mix
h_in <- (h_exh*mdot_dry_exh + h_amb*mdot_dry_amb) / mdot_dry_mix

# Recover inlet temperature from h_in and w_in
T_in_C <- uniroot(function(Tc) h_kJkgdry(Tc, w_in) - h_in, c(-20, 50))$root

# Inlet RH
e_in  <- w_in * p_Pa / (0.62198 + w_in)
RH_in <- e_in / es_Pa(T_in_C)

# Residence time (s)
U_ms  <- V_tot_m3s / (pi * (D_tube_m/2)^2) # mean velocity
t_res <- L_tube_m / U_ms

# ----------------------------
# Assemble constant terms
# ----------------------------
K_val <- K_term(L_v_Jkg, rho_l_kgm3, k_air_WmK, R_v_JkgK, T_eff_K)
D0    <- D_term(rho_l_kgm3, R_v_JkgK, T_eff_K, es_Pa(T_eff_C), Dv_m2s)

# ----------------------------
# Time integration (Maxwell r^2 with FS + κ-Köhler)
# ----------------------------
grow_one <- function(d0_m) {
  r0 <- d0_m / 2
  rd <- solve_rd_from_r0(r0, RH_in, kappa, T_eff_K, sigma_Npm, rho_l_kgm3, Mw_kgmol, p_Pa)
  
  dt <- t_res / n_steps
  r  <- r0
  
  for (i in 1:n_steps) {
    dp  <- 2 * r
    bfs <- beta_FS(dp, lambda_m, alpha_m)     # 0 < bfs <= 1
    Def <- D0 / bfs                           # increased resistance when bfs < 1
    Seq <- S_eq_kohler(r, kappa, rd, T_eff_K, sigma_Npm, rho_l_kgm3, Mw_kgmol, p_Pa)
    drive <- max(0, S_eff - Seq)              # no shrinkage in this demonstration
    dr2dt <- 2 * drive / (K_val + Def)        # Maxwell r^2 form
    r2    <- r^2 + dr2dt * dt
    r     <- sqrt(max(r2, rd^2 * 1.0000001))  # enforce r > rd
  }
  
  d_final <- 2 * r
  list(rd = rd, d_final = d_final, d0 = d0_m)
}

# Batch evaluation
res_list <- lapply(d0_list_m, grow_one)
res_df <- data.frame(
  d0_nm      = sapply(res_list, function(x) x$d0)      * 1e9,
  rd_nm      = sapply(res_list, function(x) x$rd)      * 1e9,
  d_final_nm = sapply(res_list, function(x) x$d_final) * 1e9
)
res_df$growth_factor <- res_df$d_final_nm / res_df$d0_nm

# ----------------------------
# Print summary
# ----------------------------
cat("==== Inlet/Mixing Summary ====\n")
cat(sprintf("T_in (°C)     : %.2f\n", T_in_C))
cat(sprintf("RH_in (%%)     : %.1f\n", RH_in*100))
cat(sprintf("T_eff (°C)    : %.2f\n", T_eff_C))
cat(sprintf("t_res (ms)    : %.3f\n", t_res*1e3))
cat(sprintf("S_eff (-)     : %.2f  (~%d%%)\n", S_eff, round((S_eff-1)*100)))
cat(sprintf("kappa (-)     : %.2f\n", kappa))
cat(sprintf("lambda (nm)   : %.0f\n", lambda_m*1e9))
cat("==============================\n\n")

cat("---- Growth Results (wet diameters) ----\n")
print(res_df[, c("d0_nm","rd_nm","d_final_nm","growth_factor")], row.names = FALSE)

# ----------------------------
# Optional CSV output
# ----------------------------
if (isTRUE(SAVE_RESULTS)) {
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(res_df, out_csv, row.names = FALSE)
  cat(sprintf("\nResults written to: %s\n", out_csv))
}
