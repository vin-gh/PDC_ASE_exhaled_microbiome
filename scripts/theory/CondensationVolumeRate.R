# ===========================================
# PDC 冷凝水体积生成速率计算（R 版）
# 作者：<你的名字>（可替换）
# 说明：
# 1) 修改“参数区”的变量即可复现不同工况
# 2) 输出同时包含 不加η 的上限值 与 乘以η 的有效值
# 3) 采用常用 Magnus 公式计算饱和水汽压（kPa），准确、稳定
# ===========================================

# ------------------------
# 参数区（按需修改）
# ------------------------
p_kPa        <- 101.325        # 总压（kPa），默认标准大气压
R_kPaL_molK  <- 8.314          # 通用气体常数（kPa·L/(mol·K)）
rho_w_kg_L   <- 1.0            # 液态水密度（kg/L），常温近似 1
M_da_kg_mol  <- 0.02897        # 干空气摩尔质量（kg/mol）

# 三股空气：呼出气 + 环境新风 = PDC 总抽吸
Vdot_exh_Lpm <- 10             # 呼出气体积流量（L/min）
T_exh_C      <- 28             # 呼出气温度（°C）
RH_exh       <- 0.9           # 呼出气相对湿度（–，100% = 1.00）

Vdot_amb_Lpm <- 15             # 环境新风体积流量（L/min）
T_amb_C      <- 20             # 新风温度（°C）
RH_amb       <- 0.50           # 新风相对湿度（–，50% = 0.50）

Vdot_PDC_Lpm <- 25             # PDC 抽吸总流量（L/min）——建议与上面两股相加一致
auto_balance <- TRUE           # 若两股和 ≠ PDC，总流量是否自动用 PDC 值平衡新风？

T_out_C      <- 4              # 出口/冷凝后气体温度（°C），用气旋器外壁温度近似
eta_effect   <- 0.24           # η 效果系数（0~1），代表传热传质/收集等综合有效性

# ------------------------
# 工具函数：饱和水汽压（kPa），Magnus 公式（常用、可靠）
# 参考：e_s(T) = 0.61094 * exp(17.625*T / (T + 243.04))
# ------------------------
es_kPa <- function(T_C) {
  0.61094 * exp(17.625 * T_C / (T_C + 243.04))
}

# ------------------------
# 核心计算函数（返回列表 + 友好打印）
# ------------------------
compute_condensation <- function(
    p_kPa, R_kPaL_molK, rho_w_kg_L, M_da_kg_mol,
    Vdot_exh_Lpm, T_exh_C, RH_exh,
    Vdot_amb_Lpm, T_amb_C, RH_amb,
    Vdot_PDC_Lpm, auto_balance,
    T_out_C, eta_effect
) {
  # 1) 如需，自动平衡新风流量 = 总抽吸 - 呼出气；避免不一致导致质量不守恒
  if (auto_balance) {
    Vdot_amb_Lpm <- Vdot_PDC_Lpm - Vdot_exh_Lpm
    if (Vdot_amb_Lpm < 0) {
      stop("总抽吸 < 呼出气流量，无法平衡。请检查 Vdot_PDC_Lpm 与 Vdot_exh_Lpm。")
    }
  } else {
    # 如不自动平衡，给出提示
    if (abs((Vdot_exh_Lpm + Vdot_amb_Lpm) - Vdot_PDC_Lpm) > 1e-6) {
      warning("注意：呼出气 + 新风 ≠ PDC 总抽吸。建议开启 auto_balance 或手动调整。")
    }
  }
  
  # 2) 计算各股气流的水汽分压 e（kPa）
  e_exh_kPa <- RH_exh * es_kPa(T_exh_C)
  e_amb_kPa <- RH_amb * es_kPa(T_amb_C)
  
  # 3) 含湿量 w = 0.62198 * e / (p - e)（kg水/kg干空气）
  w_exh <- 0.62198 * e_exh_kPa / (p_kPa - e_exh_kPa)
  w_amb <- 0.62198 * e_amb_kPa / (p_kPa - e_amb_kPa)
  
  # 4) 各股总摩尔流率（mol/min）：n_tot = pV/RT（kPa·L / (kPa·L/mol/K * K)）
  T_exh_K <- T_exh_C + 273.15
  T_amb_K <- T_amb_C + 273.15
  n_tot_exh_mpm <- p_kPa * Vdot_exh_Lpm / (R_kPaL_molK * T_exh_K)
  n_tot_amb_mpm <- p_kPa * Vdot_amb_Lpm / (R_kPaL_molK * T_amb_K)
  
  # 5) 水汽摩尔分数 xw = e/p；干空气摩尔流率 n_dry = n_tot * (1 - xw)
  xw_exh <- e_exh_kPa / p_kPa
  xw_amb <- e_amb_kPa / p_kPa
  n_dry_exh_mpm <- n_tot_exh_mpm * (1 - xw_exh)
  n_dry_amb_mpm <- n_tot_amb_mpm * (1 - xw_amb)
  
  # 6) 干空气质量流率（kg/min）
  mdot_dry_exh_kgpm <- n_dry_exh_mpm * M_da_kg_mol
  mdot_dry_amb_kgpm <- n_dry_amb_mpm * M_da_kg_mol
  mdot_dry_mix_kgpm <- mdot_dry_exh_kgpm + mdot_dry_amb_kgpm
  
  # 7) 混合水汽质量流率（kg/min）与混合含湿量 w_mix
  mdot_w_exh_kgpm <- w_exh * mdot_dry_exh_kgpm
  mdot_w_amb_kgpm <- w_amb * mdot_dry_amb_kgpm
  mdot_w_mix_kgpm <- mdot_w_exh_kgpm + mdot_w_amb_kgpm
  w_mix <- mdot_w_mix_kgpm / mdot_dry_mix_kgpm
  
  # 8) 出口（T_out）饱和含湿量 w*（kg/kg）
  e_out_kPa <- es_kPa(T_out_C)
  w_out_star <- 0.62198 * e_out_kPa / (p_kPa - e_out_kPa)
  
  # 9) 冷凝质量流率（上限，kg/min）与体积流率（L/min）
  delta_w <- max(0, w_mix - w_out_star)  # 防止负值（不发生冷凝时为 0）
  mdot_cond_kgpm_upper <- delta_w * mdot_dry_mix_kgpm
  Vdot_cond_Lpm_upper  <- mdot_cond_kgpm_upper / rho_w_kg_L
  
  # 10) 乘以 η 的有效值（若 η=1 则等于上限）
  Vdot_cond_Lpm_eta <- eta_effect * Vdot_cond_Lpm_upper
  
  # 11) 便于阅读的列表与打印
  out <- list(
    inputs = list(
      p_kPa = p_kPa, Vdot_exh_Lpm = Vdot_exh_Lpm, T_exh_C = T_exh_C, RH_exh = RH_exh,
      Vdot_amb_Lpm = Vdot_amb_Lpm, T_amb_C = T_amb_C, RH_amb = RH_amb,
      Vdot_PDC_Lpm = Vdot_PDC_Lpm, auto_balance = auto_balance,
      T_out_C = T_out_C, eta_effect = eta_effect
    ),
    psychrometrics = list(
      es_exh_kPa = e_exh_kPa, es_amb_kPa = e_amb_kPa, es_out_kPa = e_out_kPa,
      w_exh = w_exh, w_amb = w_amb, w_mix = w_mix, w_out_star = w_out_star
    ),
    dry_air = list(
      n_dry_exh_mpm = n_dry_exh_mpm, n_dry_amb_mpm = n_dry_amb_mpm,
      mdot_dry_exh_kgpm = mdot_dry_exh_kgpm, mdot_dry_amb_kgpm = mdot_dry_amb_kgpm,
      mdot_dry_mix_kgpm = mdot_dry_mix_kgpm
    ),
    results = list(
      Vdot_cond_Lpm_upper = Vdot_cond_Lpm_upper,
      Vdot_cond_Lpm_eta   = Vdot_cond_Lpm_eta,
      mdot_cond_kgpm_upper = mdot_cond_kgpm_upper
    )
  )
  
  # 友好输出
  cat("====== PDC 冷凝水体积生成速率（理论）======\n")
  cat(sprintf("总压 p           : %.3f kPa\n", p_kPa))
  cat(sprintf("呼出气          : %5.2f L/min, %4.1f °C, RH=%3.0f%%\n",
              Vdot_exh_Lpm, T_exh_C, RH_exh*100))
  cat(sprintf("环境新风        : %5.2f L/min, %4.1f °C, RH=%3.0f%%\n",
              Vdot_amb_Lpm, T_amb_C, RH_amb*100))
  cat(sprintf("PDC 总抽吸       : %5.2f L/min (auto_balance=%s)\n",
              Vdot_PDC_Lpm, as.character(auto_balance)))
  cat(sprintf("出口/冷凝温度    : %4.1f °C\n", T_out_C))
  cat(sprintf("η（效果系数）    : %.2f\n", eta_effect))
  cat("-------------------------------------------------\n")
  cat(sprintf("混合含湿量 w_mix : %.5f kg水/kg干空气\n", w_mix))
  cat(sprintf("出口饱和 w*      : %.5f kg水/kg干空气 @ %.1f °C\n", w_out_star, T_out_C))
  cat(sprintf("干空气质量流率   : %.6f kg/min\n", mdot_dry_mix_kgpm))
  cat("-------------------------------------------------\n")
  cat(sprintf("冷凝 质量流率上限: %.6f kg/min\n", mdot_cond_kgpm_upper))
  cat(sprintf("冷凝 体积速率上限: %.3f mL/min (≈ %.2f mL/10min)\n",
              Vdot_cond_Lpm_upper*1000, Vdot_cond_Lpm_upper*10000))
  cat(sprintf("考虑 η 后体积速率: %.3f mL/min (≈ %.2f mL/10min)\n",
              Vdot_cond_Lpm_eta*1000,   Vdot_cond_Lpm_eta*10000))
  cat("===============================================\n")
  
  invisible(out)
}

# ------------------------
# 运行（用上面的默认参数）
# ------------------------
res <- compute_condensation(
  p_kPa, R_kPaL_molK, rho_w_kg_L, M_da_kg_mol,
  Vdot_exh_Lpm, T_exh_C, RH_exh,
  Vdot_amb_Lpm, T_amb_C, RH_amb,
  Vdot_PDC_Lpm, auto_balance,
  T_out_C, eta_effect
)

# 你可以通过 res$results 或 res$psychrometrics 查看细节，例如：
# str(res)
# res$results$Vdot_cond_Lpm_upper  # 不考虑 η 的上限 (L/min)
# res$results$Vdot_cond_Lpm_eta    # 考虑 η 的值 (L/min)
