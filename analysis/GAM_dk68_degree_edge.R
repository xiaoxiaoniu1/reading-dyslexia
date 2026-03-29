# ============================================================
# R Script: GAM Interaction Analysis (Diagnosis x Continuous Age)
# Purpose: Test if developmental trajectories differ between TD and DD
# Method:  Generalized Additive Models (mgcv) with Difference Smooths
# ============================================================

# --- 0. Load Packages ---
if(!require(pacman)) install.packages("pacman",repos = "https://cloud.r-project.org/")
pacman::p_load(readxl, dplyr, purrr, stringr, tools, mgcv, ggplot2)

# ============================================================
# 1) Paths & Settings
# ============================================================
# 请务必确认路径正确
demo_file    <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_any_2/all_data_cqt_any_2.xlsx"
mind_csv_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_any_2/MIND_DK68_combat/"

out_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_any_2/MIND_GAM_DK68"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_degree_csv <- file.path(out_dir, "GAM_DK68_degree_results.csv")
out_edge_csv   <- file.path(out_dir, "GAM_DK68_edge_results.csv")

# 专门存放显著结果图片的文件夹
plot_dir <- file.path(out_dir, "Significant_Plots") 
dir.create(plot_dir, showWarnings = FALSE)

# ============================================================
# 2) Data Loading & Preprocessing
# ============================================================
cat("\n[Step 1] Loading Demographics...\n")
df <- read_excel(demo_file, sheet = "Sheet1")

df <- df %>%
  mutate(
    subj_id = as.character(id),
    subj_id = trimws(subj_id),
    file_base = paste0(subj_id, "_combat"),
    mind_file = file.path(mind_csv_dir, paste0(file_base, ".csv")),
    
    # === 核心变量定义 ===
    Diagnosis_Raw = ifelse(group_d_or_c == 0, "TD", "DD"),
    Sex      = factor(ifelse(sex == 1, "Male", "Female")),
    Site     = as.factor(site),
    Age_Raw  = as.numeric(age_month), # 本身就是年龄
    
    has_file = file.exists(mind_file)
  ) %>%
  filter(!is.na(id), subj_id != "", !is.na(Age_Raw), has_file == TRUE)

# === GAM 关键设置 ===
# 1. 年龄中心化：让截距代表平均年龄时的组间差异，减少多重共线性
df$Age_Centered <- df$Age_Raw - mean(df$Age_Raw, na.rm = TRUE)

# 2. 因子设置：TD 为参考水平 (Reference Level)
df$Diagnosis_Fac <- factor(df$Diagnosis_Raw, levels = c("TD", "DD"))

# 3. 构造数值型哑变量 (Dummy Variable)
# 用于 by=DD，从而拟合差异曲线。TD=0, DD=1。
df$DD <- as.numeric(df$Diagnosis_Fac == "DD")

df2 <- df
cat(" -> Included Subjects:", nrow(df2), "\n")
cat(" -> Age Range:", round(min(df2$Age_Raw), 1), "-", round(max(df2$Age_Raw), 1), "years\n")

# ============================================================
# 3) Matrix Loading Functions
# ============================================================
read_mind_csv <- function(fp) {
  if (!file.exists(fp)) return(NULL)
  mat <- tryCatch(as.matrix(read.csv(fp, row.names = 1, check.names = FALSE)), error=function(e) NULL)
  if (!is.null(mat)) storage.mode(mat) <- "numeric"
  mat
}

calc_degree <- function(mat) {
  diag(mat) <- NA
  rowMeans(mat, na.rm = TRUE)
}

# --- Initialize Matrix ---
cat("\n[Step 2] Reading Matrices...\n")
tmp_mat <- read_mind_csv(df2$mind_file[1])
n_roi <- nrow(tmp_mat)
roi_names <- rownames(tmp_mat)
if (is.null(roi_names)) roi_names <- paste0("ROI_", 1:n_roi)

# 提取上三角索引
upper_idx <- which(upper.tri(tmp_mat), arr.ind = TRUE)
n_edge <- nrow(upper_idx)

degree_mat <- matrix(NA_real_, nrow = nrow(df2), ncol = n_roi, dimnames = list(df2$file_base, roi_names))
edge_mat   <- matrix(NA_real_, nrow = nrow(df2), ncol = n_edge, dimnames = list(df2$file_base, paste0("E_", 1:n_edge)))

pb <- txtProgressBar(min = 0, max = nrow(df2), style = 3)
for (i in 1:nrow(df2)) {
  mat <- read_mind_csv(df2$mind_file[i])
  if (!is.null(mat) && nrow(mat) == n_roi) {
    degree_mat[i, ] <- calc_degree(mat)
    edge_mat[i, ]   <- mat[upper.tri(mat)]
  }
  setTxtProgressBar(pb, i)
}
close(pb)
cat("\n")

# ============================================================
# 4) Visualization Function (Robust)
# ============================================================
plot_gam_trajectory <- function(gam_model, data, feature_name, out_path) {
  # 1. 创建预测网格
  age_seq <- seq(min(data$Age_Centered), max(data$Age_Centered), length.out = 100)
  age_raw_seq <- age_seq + mean(data$Age_Raw) # 还原真实年龄用于X轴
  
  # 2. 协变量控制：设为众数或第一行数据
  sex_ref  <- factor(levels(data$Sex)[1], levels = levels(data$Sex)) 
  site_ref <- factor(levels(data$Site)[1], levels = levels(data$Site))
  
  # 3. 构造预测数据框
  # Predict for TD (DD=0)
  new_td <- data.frame(
  Age_Centered = age_seq,
  Diagnosis_Fac = factor("TD", levels = levels(data$Diagnosis_Fac)),
  DD = 0, Sex = sex_ref, Site = site_ref
  )
  # Predict for DD (DD=1)
  new_dd <- data.frame(
  Age_Centered = age_seq,
  Diagnosis_Fac = factor("DD", levels = levels(data$Diagnosis_Fac)),
  DD = 1, Sex = sex_ref, Site = site_ref
  )
  
  # 4. 获取预测值和标准误
  pred_td <- predict(gam_model, newdata = new_td, se.fit = TRUE)
  pred_dd <- predict(gam_model, newdata = new_dd, se.fit = TRUE)
  
  plot_df <- rbind(
    data.frame(Age = age_raw_seq, Fit = pred_td$fit, SE = pred_td$se.fit, Group = "TD"),
    data.frame(Age = age_raw_seq, Fit = pred_dd$fit, SE = pred_dd$se.fit, Group = "DD")
  )
  
  # 5. 绘图 (ggplot2)
  p <- ggplot() +
    # 原始散点
    geom_point(data = data, aes(x = Age_Raw, y = y, color = Diagnosis_Fac), alpha = 0.2, size = 1) +
    # 置信区间带
    geom_ribbon(data = plot_df, aes(x = Age, ymin = Fit - 1.96*SE, ymax = Fit + 1.96*SE, fill = Group), alpha = 0.2, color = NA) +
    # 拟合线
    geom_line(data = plot_df, aes(x = Age, y = Fit, color = Group), size = 1.2) +
    
    scale_color_manual(values = c("TD" = "#d62728", "DD" = "#1f77b4")) +
    scale_fill_manual(values = c("TD" = "#d62728", "DD" = "#1f77b4")) +
    
    labs(title = paste0(feature_name), 
         subtitle = "GAM Trajectory (Shaded = 95% CI)",
         x = "Age (Years)", y = "MSN Degree") +
    theme_bw(base_size = 14) +
    theme(legend.position = "top", panel.grid.minor = element_blank())
  
  ggsave(out_path, p, width = 6, height = 4, dpi = 300)
}

plot_gam_difference <- function(gam_model, data, feature_name, out_path) {
  age_seq <- seq(min(data$Age_Centered), max(data$Age_Centered), length.out = 100)
  age_raw_seq <- age_seq + mean(data$Age_Raw)
  
  sex_ref  <- factor(levels(data$Sex)[1], levels = levels(data$Sex)) 
  site_ref <- factor(levels(data$Site)[1], levels = levels(data$Site))
  
  new_td <- data.frame(
  Age_Centered = age_seq,
  Diagnosis_Fac = factor("TD", levels = levels(data$Diagnosis_Fac)),
  DD = 0, Sex = sex_ref, Site = site_ref
  )
  new_dd <- data.frame(
  Age_Centered = age_seq,
  Diagnosis_Fac = factor("DD", levels = levels(data$Diagnosis_Fac)),
  DD = 1, Sex = sex_ref, Site = site_ref
  )
  
  pred_td <- predict(gam_model, newdata = new_td, se.fit = TRUE)
  pred_dd <- predict(gam_model, newdata = new_dd, se.fit = TRUE)
  
  diff_fit <- pred_dd$fit - pred_td$fit
  diff_se  <- sqrt(pred_dd$se.fit^2 + pred_td$se.fit^2)
  
  plot_df <- data.frame(Age = age_raw_seq, Diff = diff_fit, SE = diff_se)
  
  p <- ggplot(plot_df, aes(x = Age, y = Diff)) +
    geom_ribbon(aes(ymin = Diff - 1.96*SE, ymax = Diff + 1.96*SE), alpha = 0.2, fill = "#1f77b4", color = NA) +
    geom_line(size = 1.2, color = "#1f77b4") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    labs(title = paste0(feature_name, " (DD - TD)"),
         subtitle = "GAM Difference Curve (Shaded = 95% CI)",
         x = "Age (Years)", y = "DD - TD") +
    theme_bw(base_size = 14) +
    theme(panel.grid.minor = element_blank())
  
  ggsave(out_path, p, width = 6, height = 4, dpi = 300)
}

# ============================================================
# 5) GAM Core Engine
# ============================================================
run_gam_core <- function(y_vec, df_sub, feature_name, out_dir, save_plot = FALSE) {
  
  d_run <- df_sub
  d_run$y <- y_vec
  d_run <- d_run %>% filter(is.finite(y))
  
  # 样本量检查
  if (nrow(d_run) < 30 || var(d_run$y) < 1e-12) {
    return(list(stats = c(n=nrow(d_run), p_Diag=NA, p_Age=NA, p_Int=NA, R2=NA), plot=FALSE))
  }
  
  # === 模型定义 ===
  # y ~ 截距组别差异 + s(全体/TD趋势) + s(DD的偏差趋势) + 协变量
  fit <- tryCatch({
    gam(y ~ Diagnosis_Fac + s(Age_Centered, k=5) + s(Age_Centered, by=DD, k=5) + Sex,
        data = d_run, method = "REML")
  }, error = function(e) NULL)
  
  if (is.null(fit)) return(list(stats = c(n=nrow(d_run), p_Diag=NA, p_Age=NA, p_Int=NA, R2=NA), plot=FALSE))
  
  summ <- summary(fit)
  
  # === 提取 P 值 (更加稳健的 grep) ===
  # 1. 截距差异 (Group Difference at Mean Age)
  p_diag <- NA
  if ("Diagnosis_FacDD" %in% rownames(summ$p.table)) {
    p_diag <- summ$p.table["Diagnosis_FacDD", "Pr(>|t|)"]
  }
  
  # 2. 交互效应 (Difference Smooth): 查找含有 ":DD" 的平滑项
  # 这代表 DD 的发育轨迹形状是否显著偏离 TD
  idx_int <- grep(":DD", rownames(summ$s.table), fixed = TRUE)
  p_int   <- if(length(idx_int)>0) summ$s.table[idx_int[1], "p-value"] else NA
  
  # 3. 年龄主效应 (Common/Reference Smooth)
  # 查找完全匹配 "s(Age_Centered)" 的项
  idx_age <- grep("^s\\(Age_Centered\\)$", rownames(summ$s.table))
  p_age   <- if(length(idx_age)>0) summ$s.table[idx_age[1], "p-value"] else NA
  
  # === 自动绘图 (仅当交互显著 P < 0.05 时) ===
  plot_done <- FALSE
  if (save_plot && !is.na(p_int) && p_int < 0.05) {
    safe_fname <- gsub("[^a-zA-Z0-9_]", "_", feature_name) 
    out_path <- file.path(out_dir, paste0("GAM_SigInt_", safe_fname, ".png"))
    out_path_diff <- file.path(out_dir, paste0("GAM_SigInt_Diff_", safe_fname, ".png"))
    try({ plot_gam_trajectory(fit, d_run, feature_name, out_path); plot_done <- TRUE })
    try({ plot_gam_difference(fit, d_run, feature_name, out_path_diff); plot_done <- TRUE })
  }
  
  return(list(
    stats = c(n=nrow(d_run), p_Diag=p_diag, p_Age=p_age, p_Int=p_int, R2=summ$r.sq), 
    plot = plot_done
  ))
}

# ============================================================
# 6) Execution: ROI-wise
# ============================================================
cat("\n[Step 3] Running ROI-wise GAM...\n")
roi_res_list <- list()

for (k in 1:ncol(degree_mat)) {
  fname <- colnames(degree_mat)[k]
  # 对 Degree 我们总是尝试画图，因为特征数少
  res <- run_gam_core(degree_mat[,k], df2, fname, plot_dir, save_plot = TRUE)
  
  roi_res_list[[k]] <- data.frame(
    feature = fname,
    t(res$stats),
    Has_Plot = res$plot
  )
}

roi_out <- bind_rows(roi_res_list)
# FDR 校正
roi_out$p_Diag_FDR <- p.adjust(roi_out$p_Diag, "fdr")
roi_out$p_Int_FDR  <- p.adjust(roi_out$p_Int, "fdr") 
roi_out$p_Age_FDR  <- p.adjust(roi_out$p_Age, "fdr")

cat("Significant Interactions (ROI, FDR<0.05):", sum(roi_out$p_Int_FDR < 0.05, na.rm=TRUE), "\n")
write.csv(roi_out, out_degree_csv, row.names = FALSE)

# ============================================================
# 7) Execution: Edge-wise (Chunked)
# ============================================================
cat("\n[Step 4] Running Edge-wise GAM (Chunked)...\n")
# Edge 数量大，不自动画图，避免生成数万张图片，除非你特别想看
chunk_size <- 2000 
chunks <- split(seq_len(n_edge), ceiling(seq_len(n_edge) / chunk_size))
edge_out_list <- list()

pb <- txtProgressBar(min = 0, max = length(chunks), style = 3)

for (ci in seq_along(chunks)) {
  idxs <- chunks[[ci]]
  
  chunk_res <- lapply(idxs, function(k) {
    # 构造易读的 Edge 名称
    e_name <- paste0(roi_names[upper_idx[k,1]], "__", roi_names[upper_idx[k,2]])
    
    # 运行 GAM (save_plot = FALSE)
    res <- run_gam_core(edge_mat[, k], df2, e_name, plot_dir, save_plot = FALSE)
    
    data.frame(
      edge_idx = k,
      roi_1 = roi_names[upper_idx[k,1]],
      roi_2 = roi_names[upper_idx[k,2]],
      t(res$stats)
    )
  })
  
  edge_out_list[[ci]] <- bind_rows(chunk_res)
  setTxtProgressBar(pb, ci)
}
close(pb)

edge_out <- bind_rows(edge_out_list)
edge_out$p_Diag_FDR <- p.adjust(edge_out$p_Diag, "fdr")
edge_out$p_Age_FDR  <- p.adjust(edge_out$p_Age, "fdr")
edge_out$p_Int_FDR  <- p.adjust(edge_out$p_Int, "fdr")

cat("\nSignificant Interactions (Edge, FDR<0.05):", sum(edge_out$p_Int_FDR < 0.05, na.rm=TRUE), "\n")
write.csv(edge_out, out_edge_csv, row.names = FALSE)

output_sig_gam_main_effect <- function(res_df, feature_mat, df_demo, effect_name, feature_type, out_dir) {
  p_col <- if (effect_name == "Diagnosis") "p_Diag_FDR" else "p_Age_FDR"
  sig_data <- res_df %>% filter(.data[[p_col]] < 0.05)
  if (nrow(sig_data) == 0) {
    cat("No significant", effect_name, "results for", feature_type, "(FDR<0.05)\n")
    return(NULL)
  }

  sig_with_dir <- sig_data %>%
    rowwise() %>%
    mutate(
      idx_or_name = if (feature_type == "degree") feature else edge_idx,
      mean_TD = if (effect_name == "Diagnosis") mean(feature_mat[df_demo$Diagnosis_Fac == "TD", idx_or_name], na.rm = TRUE) else NA_real_,
      mean_DD = if (effect_name == "Diagnosis") mean(feature_mat[df_demo$Diagnosis_Fac == "DD", idx_or_name], na.rm = TRUE) else NA_real_,
      age_cor = if (effect_name == "Age") suppressWarnings(cor(df_demo$Age_Raw, feature_mat[, idx_or_name], use = "complete.obs")) else NA_real_,
      direction = if (effect_name == "Diagnosis") {
        ifelse(mean_DD > mean_TD, "DD>TD", ifelse(mean_DD < mean_TD, "DD<TD", "equal"))
      } else {
        ifelse(is.na(age_cor), NA_character_, ifelse(age_cor > 0, "Age_Pos", ifelse(age_cor < 0, "Age_Neg", "zero")))
      }
    ) %>%
    ungroup()

  base_cols <- c("feature", "edge_idx", "roi_1", "roi_2", "n", "R2", "p_Diag", "p_Age", "p_Int", "p_Diag_FDR", "p_Age_FDR", "p_Int_FDR")
  if (effect_name == "Diagnosis") {
    sig_out <- sig_with_dir %>% select(any_of(base_cols), mean_TD, mean_DD, direction)
  } else {
    sig_out <- sig_with_dir %>% select(any_of(base_cols), age_cor, direction)
  }

  out_file <- file.path(out_dir, paste0("Significant_", effect_name, "_GAM_DK68_", feature_type, "_results.csv"))
  write.csv(sig_out, out_file, row.names = FALSE)
  cat("Significant", effect_name, feature_type, "results (n =", nrow(sig_out), "):", out_file, "\n")

  return(sig_out)
}

output_sig_gam_main_effect(roi_out, degree_mat, df2, "Diagnosis", "degree", out_dir)
output_sig_gam_main_effect(edge_out, edge_mat, df2, "Diagnosis", "edge", out_dir)
output_sig_gam_main_effect(roi_out, degree_mat, df2, "Age", "degree", out_dir)
output_sig_gam_main_effect(edge_out, edge_mat, df2, "Age", "edge", out_dir)

cat("\nALL DONE.\n")