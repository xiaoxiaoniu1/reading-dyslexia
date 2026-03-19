# ============================================================
# Two-way ANOVA (Diagnosis x AgeGroup) for:
#  (A) DK-68 MIND degree (ROI-wise)
#  (B) DK-68 MIND edges (upper-triangle vector)
# Covariates: age_month, Sex, Site
# ============================================================

library(readxl)
library(dplyr)
library(purrr)
library(car)
library(stringr)
library(tools)

# ----------------------------
# 1) Paths
# ----------------------------
demo_file    <- "/data/home/tqi/data1/share/after_freesurfer/FILE/all_data_cqt.xlsx"
mind_csv_dir <- "/data/home/tqi/data1/share/after_freesurfer/fs_subjects_all/MIND_DK68_combat"

out_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/MIND_DK68_ANOVA/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_degree_csv <- file.path(out_dir, "ANOVA_DK68_degree_results.csv")
out_edge_csv   <- file.path(out_dir, "ANOVA_DK68_edge_results.csv")

# ----------------------------
# 2) Load demo + build file id rule + recode factors
# ----------------------------
df <- read_excel(demo_file, sheet = "Sheet1")

# DK68 文件命名规则：id + "_combat.csv"
# 例如：1_combat.csv, 100_combat.csv
df <- df %>%
  mutate(
    subj_id = as.character(id),
    
    # 文件名：1_combat.csv
    file_base = paste0(subj_id, "_combat"),
    
    # 完整路径
    mind_file = file.path(mind_csv_dir, paste0(file_base, ".csv")),
    
    # 组别/协变量
    Diagnosis = ifelse(group_d_or_c == 0, "TD", "DD"),
    AgeGroup  = ifelse(group_age == 1, "Adult", "Child"),
    Sex       = ifelse(sex == 1, "Male", "Female"),
    Site      = as.factor(site),
    
    Diagnosis = factor(Diagnosis, levels = c("TD", "DD")),
    AgeGroup  = factor(AgeGroup,  levels = c("Child", "Adult")),
    Sex       = factor(Sex),
    
    has_file  = file.exists(mind_file)
  ) %>%
  # 防呆：去掉任何 NA / 空 id 的行，并去掉 site=3
  filter(!is.na(id), id != "", site != 3)

# ----------------------------
# 3) Filter missing matrix subjects + write missing list
# ----------------------------
cat("Subjects in demo (after basic cleaning):", nrow(df), "\n")
cat("Subjects with existing MIND csv:", sum(df$has_file), "\n")
cat("Missing subjects (demo has but file not):", sum(!df$has_file), "\n")

missing_df <- df %>%
  filter(!has_file) %>%
  select(id, subj_id, file_base, mind_file)

write.csv(missing_df,
          file = file.path(out_dir, "missing_subjects_no_matrix.csv"),
          row.names = FALSE)

# 只保留有矩阵文件的被试
df2 <- df %>% filter(has_file)

# ----------------------------
# 4) Read matrix + helpers
# ----------------------------
read_mind_csv <- function(fp) {
  # DK68 文件格式：第一列是行名，第一行是列名
  mat <- as.matrix(read.csv(fp, row.names = 1, check.names = FALSE))
  storage.mode(mat) <- "numeric"
  mat
}

calc_degree <- function(mat) {
  diag(mat) <- NA
  rowMeans(mat, na.rm = TRUE)
}

# detect ROI count from first subject
tmp <- read_mind_csv(df2$mind_file[1])
n_roi <- nrow(tmp)
stopifnot(n_roi == ncol(tmp))
cat("Detected ROI number:", n_roi, "\n")

# upper-triangle indices mapping
upper_idx <- which(upper.tri(matrix(1, n_roi, n_roi)), arr.ind = TRUE)
n_edge <- nrow(upper_idx)
cat("Number of edges (upper triangle):", n_edge, "\n")

# ----------------------------
# 5) Build degree matrix + edge matrix
# ----------------------------
n_subj <- nrow(df2)

# 用 ROI 真名作为列名（如果你矩阵带 ROI 名）
roi_names <- rownames(tmp)
if (is.null(roi_names) || any(is.na(roi_names)) || length(roi_names) != n_roi) {
  roi_names <- paste0("ROI_", seq_len(n_roi))
}

degree_mat <- matrix(NA_real_, nrow = n_subj, ncol = n_roi,
                     dimnames = list(df2$file_base, roi_names))

edge_mat <- matrix(NA_real_, nrow = n_subj, ncol = n_edge,
                   dimnames = list(df2$file_base, paste0("E_", seq_len(n_edge))))

for (s in seq_len(n_subj)) {
  fp <- df2$mind_file[s]
  m  <- read_mind_csv(fp)
  
  # degree
  degree_mat[s, ] <- calc_degree(m)
  
  # edges (upper triangle) — 用 upper.tri(m) 得到逻辑索引
  edge_mat[s, ] <- m[upper.tri(m)]
}

cat("Degree matrix:", paste(dim(degree_mat), collapse = " x "), "\n")
cat("Edge matrix:", paste(dim(edge_mat), collapse = " x "), "\n")

# ============================================================
# 6) ROI-wise two-way ANOVA (degree)
# ============================================================
run_two_way_anova_per_feature <- function(y_mat, df_sub, feature_names,
                                          min_n = 10, var_eps = 1e-12) {
  
  res <- lapply(seq_len(ncol(y_mat)), function(k) {
    y <- y_mat[, k]
    
    # 只保留 y 非 NA 的样本
    ok <- is.finite(y)
    y2 <- y[ok]
    d2 <- df_sub[ok, , drop = FALSE]
    
    # 防呆1：有效样本太少
    if (length(y2) < min_n) {
      return(data.frame(
        feature = feature_names[k],
        n = length(y2),
        var_y = ifelse(length(y2) > 1, var(y2), NA_real_),
        p_Diagnosis = NA_real_,
        p_AgeGroup = NA_real_,
        p_Interaction = NA_real_,
        eta2_Diagnosis = NA_real_,
        eta2_AgeGroup = NA_real_,
        eta2_Interaction = NA_real_,
        cohens_f_Diagnosis = NA_real_,
        cohens_f_AgeGroup = NA_real_,
        cohens_f_Interaction = NA_real_
      ))
    }
    
    # 防呆2：y 几乎没有方差（常数）
    v <- var(y2)
    if (!is.finite(v) || v < var_eps) {
      return(data.frame(
        feature = feature_names[k],
        n = length(y2),
        var_y = v,
        p_Diagnosis = NA_real_,
        p_AgeGroup = NA_real_,
        p_Interaction = NA_real_,
        eta2_Diagnosis = NA_real_,
        eta2_AgeGroup = NA_real_,
        eta2_Interaction = NA_real_,
        cohens_f_Diagnosis = NA_real_,
        cohens_f_AgeGroup = NA_real_,
        cohens_f_Interaction = NA_real_
      ))
    }
    
    # 正常拟合
    fit <- lm(y2 ~ Diagnosis * AgeGroup + Sex, data = d2)
    
    # 防呆3：Anova 出错时不崩溃，返回 NA
    a <- tryCatch(car::Anova(fit, type = 2), error = function(e) NULL)
    if (is.null(a)) {
      return(data.frame(
        feature = feature_names[k],
        n = length(y2),
        var_y = v,
        p_Diagnosis = NA_real_,
        p_AgeGroup = NA_real_,
        p_Interaction = NA_real_,
        eta2_Diagnosis = NA_real_,
        eta2_AgeGroup = NA_real_,
        eta2_Interaction = NA_real_,
        cohens_f_Diagnosis = NA_real_,
        cohens_f_AgeGroup = NA_real_,
        cohens_f_Interaction = NA_real_
      ))
    }
    
    # 计算 Partial Eta Squared 和 Cohen's f
    # Partial η² = SS_effect / (SS_effect + SS_residual)
    # Cohen's f = sqrt(η² / (1 - η²))
    ss_residual <- sum(residuals(fit)^2)
    
    eta2_diag <- a["Diagnosis", "Sum Sq"] / (a["Diagnosis", "Sum Sq"] + ss_residual)
    eta2_age <- a["AgeGroup", "Sum Sq"] / (a["AgeGroup", "Sum Sq"] + ss_residual)
    eta2_int <- a["Diagnosis:AgeGroup", "Sum Sq"] / (a["Diagnosis:AgeGroup", "Sum Sq"] + ss_residual)
    
    cohens_f_diag <- sqrt(eta2_diag / (1 - eta2_diag))
    cohens_f_age <- sqrt(eta2_age / (1 - eta2_age))
    cohens_f_int <- sqrt(eta2_int / (1 - eta2_int))
    
    data.frame(
      feature = feature_names[k],
      n = length(y2),
      var_y = v,
      p_Diagnosis   = a["Diagnosis", "Pr(>F)"],
      p_AgeGroup    = a["AgeGroup",  "Pr(>F)"],
      p_Interaction = a["Diagnosis:AgeGroup", "Pr(>F)"],
      eta2_Diagnosis = eta2_diag,
      eta2_AgeGroup = eta2_age,
      eta2_Interaction = eta2_int,
      cohens_f_Diagnosis = cohens_f_diag,
      cohens_f_AgeGroup = cohens_f_age,
      cohens_f_Interaction = cohens_f_int
    )
  })
  
  out <- bind_rows(res)
  
  out$p_Diagnosis_FDR   <- p.adjust(out$p_Diagnosis,   method = "fdr")
  out$p_AgeGroup_FDR    <- p.adjust(out$p_AgeGroup,    method = "fdr")
  out$p_Interaction_FDR <- p.adjust(out$p_Interaction, method = "fdr")

  out$p_Interaction_FDR_diagSig <- NA_real_
  diag_sig_idx <- which(out$p_Diagnosis_FDR < 0.05 & !is.na(out$p_Interaction))
  if (length(diag_sig_idx) > 0) {
    out$p_Interaction_FDR_diagSig[diag_sig_idx] <- p.adjust(out$p_Interaction[diag_sig_idx], method = "fdr")
  }
  
  out
}

# ====== 关键：真正运行 degree 的 ANOVA，并输出结果 ======
cat("\n=== Running ROI-wise ANOVA (degree) ===\n")
degree_res <- run_two_way_anova_per_feature(degree_mat, df2, colnames(degree_mat))

cat("Degree features with NA p-values (likely constant/degenerate):",
    sum(is.na(degree_res$p_Interaction)), "\n")

cat("Significant Interaction (degree, FDR<0.05):",
    sum(degree_res$p_Interaction_FDR < 0.05, na.rm = TRUE), "\n")

write.csv(degree_res, out_degree_csv, row.names = FALSE)
cat("Saved degree results:", out_degree_csv, "\n")

# ============================================================
# 7) Edge-wise two-way ANOVA (upper-triangle edges)
#    chunked to reduce memory/time spikes
# ============================================================
run_two_way_anova_edgewise_chunked <- function(edge_mat, df_sub, upper_idx, roi_names,
                                               chunk_size = 1000,
                                               min_n = 10, var_eps = 1e-12) {
  n_edge <- ncol(edge_mat)
  chunks <- split(seq_len(n_edge), ceiling(seq_len(n_edge) / chunk_size))
  
  out_list <- vector("list", length(chunks))
  
  for (ci in seq_along(chunks)) {
    idxs <- chunks[[ci]]
    cat("Edge chunk", ci, "/", length(chunks), "- edges:", length(idxs), "\n")
    
    chunk_res <- lapply(idxs, function(k) {
      y <- edge_mat[, k]
      ok <- is.finite(y)
      y2 <- y[ok]
      d2 <- df_sub[ok, , drop = FALSE]
      
      v <- ifelse(length(y2) > 1, var(y2), NA_real_)
      
      # 少样本 / 低方差 -> NA
      if (length(y2) < min_n || !is.finite(v) || v < var_eps) {
        return(data.frame(
          edge_id = k,
          i = upper_idx[k, 1],
          j = upper_idx[k, 2],
          region_i = roi_names[upper_idx[k, 1]],
          region_j = roi_names[upper_idx[k, 2]],
          n = length(y2),
          var_y = v,
          p_Diagnosis = NA_real_,
          p_AgeGroup = NA_real_,
          p_Interaction = NA_real_,
          eta2_Diagnosis = NA_real_,
          eta2_AgeGroup = NA_real_,
          eta2_Interaction = NA_real_,
          cohens_f_Diagnosis = NA_real_,
          cohens_f_AgeGroup = NA_real_,
          cohens_f_Interaction = NA_real_
        ))
      }
      
      fit <- lm(y2 ~ Diagnosis * AgeGroup + Sex, data = d2)
      a <- tryCatch(car::Anova(fit, type = 2), error = function(e) NULL)
      
      if (is.null(a)) {
        return(data.frame(
          edge_id = k,
          i = upper_idx[k, 1],
          j = upper_idx[k, 2],
          region_i = roi_names[upper_idx[k, 1]],
          region_j = roi_names[upper_idx[k, 2]],
          n = length(y2),
          var_y = v,
          p_Diagnosis = NA_real_,
          p_AgeGroup = NA_real_,
          p_Interaction = NA_real_,
          eta2_Diagnosis = NA_real_,
          eta2_AgeGroup = NA_real_,
          eta2_Interaction = NA_real_,
          cohens_f_Diagnosis = NA_real_,
          cohens_f_AgeGroup = NA_real_,
          cohens_f_Interaction = NA_real_
        ))
      }
      
      # 计算效应量
      ss_residual <- sum(residuals(fit)^2)
      eta2_diag <- a["Diagnosis", "Sum Sq"] / (a["Diagnosis", "Sum Sq"] + ss_residual)
      eta2_age <- a["AgeGroup", "Sum Sq"] / (a["AgeGroup", "Sum Sq"] + ss_residual)
      eta2_int <- a["Diagnosis:AgeGroup", "Sum Sq"] / (a["Diagnosis:AgeGroup", "Sum Sq"] + ss_residual)
      
      cohens_f_diag <- sqrt(eta2_diag / (1 - eta2_diag))
      cohens_f_age <- sqrt(eta2_age / (1 - eta2_age))
      cohens_f_int <- sqrt(eta2_int / (1 - eta2_int))
      
      data.frame(
        edge_id = k,
        i = upper_idx[k, 1],
        j = upper_idx[k, 2],
        region_i = roi_names[upper_idx[k, 1]],
        region_j = roi_names[upper_idx[k, 2]],
        n = length(y2),
        var_y = v,
        p_Diagnosis   = a["Diagnosis", "Pr(>F)"],
        p_AgeGroup    = a["AgeGroup",  "Pr(>F)"],
        p_Interaction = a["Diagnosis:AgeGroup", "Pr(>F)"],
        eta2_Diagnosis = eta2_diag,
        eta2_AgeGroup = eta2_age,
        eta2_Interaction = eta2_int,
        cohens_f_Diagnosis = cohens_f_diag,
        cohens_f_AgeGroup = cohens_f_age,
        cohens_f_Interaction = cohens_f_int
      )
    })
    
    out_list[[ci]] <- bind_rows(chunk_res)
  }
  
  out <- bind_rows(out_list)
  
  out$p_Diagnosis_FDR   <- p.adjust(out$p_Diagnosis,   method = "fdr")
  out$p_AgeGroup_FDR    <- p.adjust(out$p_AgeGroup,    method = "fdr")
  out$p_Interaction_FDR <- p.adjust(out$p_Interaction, method = "fdr")

  out$p_Interaction_FDR_diagSig <- NA_real_
  diag_sig_idx <- which(out$p_Diagnosis_FDR < 0.05 & !is.na(out$p_Interaction))
  if (length(diag_sig_idx) > 0) {
    out$p_Interaction_FDR_diagSig[diag_sig_idx] <- p.adjust(out$p_Interaction[diag_sig_idx], method = "fdr")
  }
  
  out
}

# ====== 关键：真正运行 edge 的 ANOVA，并输出结果 ======
cat("\n=== Running Edge-wise ANOVA (upper triangle) ===\n")
edge_res <- run_two_way_anova_edgewise_chunked(edge_mat, df2, upper_idx, roi_names, chunk_size = 1000)

cat("Edge features with NA p-values (likely constant/degenerate):",
    sum(is.na(edge_res$p_Interaction)), "\n")

cat("Significant Interaction (edges, FDR<0.05):",
    sum(edge_res$p_Interaction_FDR < 0.05, na.rm = TRUE), "\n")

write.csv(edge_res, out_edge_csv, row.names = FALSE)
cat("Saved edge results:", out_edge_csv, "\n")

# ============================================================
# 8) 输出显著结果表（degree 和 edge）- 分别输出 Diagnosis、AgeGroup、Interaction
# ============================================================
cat("\n=== Generating significant results tables ===\n")

# 辅助函数：计算方向并输出
output_sig_results <- function(res_df, feature_mat, df_demo, effect_name, 
                                 feature_type = "degree", out_dir) {
  
  p_col <- paste0("p_", effect_name, "_FDR")
  sig_data <- res_df %>% filter(.data[[p_col]] < 0.05)
  
  if (nrow(sig_data) == 0) {
    cat("No significant", effect_name, "results for", feature_type, "(FDR<0.05)\n")
    return(NULL)
  }
  
  # 根据效应类型计算不同的均值
  sig_with_dir <- sig_data %>%
    rowwise() %>%
    mutate(
      # 获取当前特征的索引（degree用特征名，edge用edge_id）
      idx_or_name = if(feature_type == "degree") feature else edge_id
    )
    
  if (effect_name == "Diagnosis") {
    # === Diagnosis Effect: TD vs DD ===
    sig_with_dir <- sig_with_dir %>%
      mutate(
        mean_TD = mean(feature_mat[df_demo$Diagnosis == "TD", idx_or_name], na.rm = TRUE),
        mean_DD = mean(feature_mat[df_demo$Diagnosis == "DD", idx_or_name], na.rm = TRUE),
        direction = ifelse(mean_DD > mean_TD, "DD>TD", 
                    ifelse(mean_DD < mean_TD, "DD<TD", "equal"))
      ) %>%
      select(any_of(c("feature", "edge_id", "i", "j", "region_i", "region_j", "n")),
             starts_with("p_"), starts_with("eta2_"), starts_with("cohens_"),
             mean_TD, mean_DD, direction)
             
  } else if (effect_name == "AgeGroup") {
    # === AgeGroup Effect: Child vs Adult ===
    sig_with_dir <- sig_with_dir %>%
      mutate(
        mean_Child = mean(feature_mat[df_demo$AgeGroup == "Child", idx_or_name], na.rm = TRUE),
        mean_Adult = mean(feature_mat[df_demo$AgeGroup == "Adult", idx_or_name], na.rm = TRUE),
        direction  = ifelse(mean_Adult > mean_Child, "Adult>Child", 
                     ifelse(mean_Adult < mean_Child, "Adult<Child", "equal"))
      ) %>%
      select(any_of(c("feature", "edge_id", "i", "j", "region_i", "region_j", "n")),
             starts_with("p_"), starts_with("eta2_"), starts_with("cohens_"),
             mean_Child, mean_Adult, direction)
             
  } else {
    # === Interaction Effect: 4 Subgroups ===
    sig_with_dir <- sig_with_dir %>%
      mutate(
        mean_TD_Child = mean(feature_mat[df_demo$Diagnosis == "TD" & df_demo$AgeGroup == "Child", idx_or_name], na.rm = TRUE),
        mean_TD_Adult = mean(feature_mat[df_demo$Diagnosis == "TD" & df_demo$AgeGroup == "Adult", idx_or_name], na.rm = TRUE),
        mean_DD_Child = mean(feature_mat[df_demo$Diagnosis == "DD" & df_demo$AgeGroup == "Child", idx_or_name], na.rm = TRUE),
        mean_DD_Adult = mean(feature_mat[df_demo$Diagnosis == "DD" & df_demo$AgeGroup == "Adult", idx_or_name], na.rm = TRUE)
      ) %>%
      select(any_of(c("feature", "edge_id", "i", "j", "region_i", "region_j", "n")),
             starts_with("p_"), starts_with("eta2_"), starts_with("cohens_"),
             mean_TD_Child, mean_TD_Adult, mean_DD_Child, mean_DD_Adult)
  }

  # 解除 rowwise
  sig_with_dir <- sig_with_dir %>% ungroup() %>% arrange(.data[[p_col]])
  
  out_file <- file.path(out_dir, 
                        paste0("Significant_", effect_name, "_DK68_", 
                              feature_type, "_results.csv"))
  write.csv(sig_with_dir, out_file, row.names = FALSE)
  cat("Significant", effect_name, feature_type, "results (n =", 
      nrow(sig_with_dir), "):", out_file, "\n")
  
  return(sig_with_dir)
}

# 输出 Diagnosis 主效应显著结果
output_sig_results(degree_res, degree_mat, df2, "Diagnosis", "degree", out_dir)
output_sig_results(edge_res, edge_mat, df2, "Diagnosis", "edge", out_dir)

# 输出 AgeGroup 主效应显著结果
output_sig_results(degree_res, degree_mat, df2, "AgeGroup", "degree", out_dir)
output_sig_results(edge_res, edge_mat, df2, "AgeGroup", "edge", out_dir)

# 输出 Interaction 显著结果
output_sig_results(degree_res, degree_mat, df2, "Interaction", "degree", out_dir)
output_sig_results(edge_res, edge_mat, df2, "Interaction", "edge", out_dir)

cat("\nDONE.\n")
