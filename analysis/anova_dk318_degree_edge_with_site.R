# ============================================================
# Two-way ANOVA (Diagnosis x AgeGroup) for:
#  (A) DK-318 MIND degree (ROI-wise)
#  (B) DK-318 MIND edges (upper-triangle vector)
# Covariates: age_month, Sex, Site
# ============================================================

library(readxl)
library(dplyr)
library(purrr)
library(car)
library(stringr)
library(tools)

# ----------------------------
# 1) Paths (改这里两行就行)
# ----------------------------
demo_file    <- "/data/home/tqi/data1/share/after_freesurfer/FILE/all_data_cqt.xlsx"
mind_csv_dir <- "/data/home/tqi/data1/share/after_freesurfer/fs_subjects_all/MIND_out"

out_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/MIND_ANOVA_with_site"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_degree_csv <- file.path(out_dir, "ANOVA_DK318_degree_results.csv")
out_edge_csv   <- file.path(out_dir, "ANOVA_DK318_edge_results.csv")

# ----------------------------
# 2) Load demo + build file id rule + recode factors
# ----------------------------
df <- read_excel(demo_file, sheet = "Sheet1")

# Clean column names for easier access (e.g., "original-project" -> "original_project")
colnames(df) <- gsub("-", "_", colnames(df)) 
colnames(df) <- gsub("\\.", "_", colnames(df)) 

# Verify if original_project exists
if (!"original_project" %in% colnames(df)) {
  stop("Column 'original-project' (or 'original_project') is required for matching but not found.")
}

# DEBUG: Check Excel project names
cat("DEBUG: Unique original_project in Excel:\n")
print(unique(df$original_project))

# List all available MIND files
all_files <- list.files(mind_csv_dir, pattern = "_MIND_DK318\\.csv$", full.names = TRUE)
file_bases <- basename(all_files)

# Extract ID and Prefix from filenames
# Format: Prefix_ID_MIND_DK318.csv (e.g. ChildD_3226_MIND_DK318.csv)
# Use a more robust regex that captures everything before the last underscore+ID
extracted_ids <- sub(".*_([0-9]+)_MIND_DK318\\.csv$", "\\1", file_bases)
extracted_prefixes <- sub("^(.*)_[0-9]+_MIND_DK318\\.csv$", "\\1", file_bases)

# DEBUG: Check what we found
cat("DEBUG: Unique extracted prefixes from files:\n")
print(unique(extracted_prefixes))

# Create a mapping dataframe
file_map <- data.frame(
  id_old = as.numeric(extracted_ids),
  file_prefix = extracted_prefixes,
  mind_file = all_files,
  file_base = tools::file_path_sans_ext(file_bases),
  stringsAsFactors = FALSE
)

# Merge demo with file map using BOTH id_old and original_project
df <- df %>%
  mutate(id_old = as.numeric(id_old)) %>%
  # Match original_project from Excel with file_prefix
  inner_join(file_map, by = c("id_old", "original_project" = "file_prefix")) %>%
  mutate(
    has_file = !is.na(mind_file),

    # Groups / Covariates
    Diagnosis = ifelse(group_d_or_c == 0, "TD", "DD"),
    AgeGroup  = ifelse(group_age == 1, "Adult", "Child"),
    Sex       = ifelse(sex == 1, "Male", "Female"),
    Site      = as.factor(site),
    
    Diagnosis = factor(Diagnosis, levels = c("TD", "DD")),
    AgeGroup  = factor(AgeGroup,  levels = c("Child", "Adult")),
    Sex       = factor(Sex)
  ) %>%
  filter(!is.na(id_old)) %>%
  # Filter out site 3 as requested
  filter(site != 3)

# ----------------------------
# 3) Filter missing matrix subjects + write missing list
# ----------------------------
cat("Subjects in demo (after basic cleaning):", nrow(df), "\n")
cat("Subjects with existing MIND csv:", sum(df$has_file), "\n")
cat("Missing subjects (demo has but file not):", sum(!df$has_file), "\n")

missing_df <- df %>%
  filter(!has_file) %>%
  select(id_old, Diagnosis, AgeGroup)

write.csv(missing_df,
          file = file.path(out_dir, "missing_subjects_no_matrix.csv"),
          row.names = FALSE)

# 只保留有矩阵文件的被试
df2 <- df %>% filter(has_file)

# ----------------------------
# 4) Read matrix + helpers
# ----------------------------
read_mind_csv <- function(fp) {
  # 这里假设你存的是“labeled matrix csv”（第一列是行名，第一行是列名）
  # 用 row.names = 1 + check.names=FALSE 才能保留 ROI 名称不被改
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

    # Check factor levels to prevent "contrasts can be applied only to factors with 2 or more levels"
    n_diag <- length(unique(d2$Diagnosis))
    n_age  <- length(unique(d2$AgeGroup))
    
    if (n_diag < 2 || n_age < 2) {
       # Cannot test interaction if one of the main factors is constant
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
    
    # Build formula dynamically: only include covariates if they have > 1 level
    f_str <- "y2 ~ Diagnosis * AgeGroup"
    if (length(unique(d2$Sex)) > 1)  f_str <- paste0(f_str, " + Sex")
    if (length(unique(d2$Site)) > 1) f_str <- paste0(f_str, " + Site")
    
    # 正常拟合
    fit <- lm(as.formula(f_str), data = d2)

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

  out
}

# ====== Key: Run degree ANOVA and output results ======
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
                                               chunk_size = 2000,
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

      # Low sample size / low variance -> NA
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

      # Check factor levels
      n_diag <- length(unique(d2$Diagnosis))
      n_age  <- length(unique(d2$AgeGroup))
      
      if (n_diag < 2 || n_age < 2) {
         return(data.frame(
          edge_id = k,
          i = upper_idx[k, 1],
          j = upper_idx[k, 2],
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

      # Build formula dynamically
      f_str <- "y2 ~ Diagnosis * AgeGroup"
      if (length(unique(d2$Sex)) > 1)  f_str <- paste0(f_str, " + Sex")
      if (length(unique(d2$Site)) > 1) f_str <- paste0(f_str, " + Site")

      fit <- lm(as.formula(f_str), data = d2)
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
      
      # Calculate effect sizes
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

  out
}

# ====== Key: Run edge ANOVA (upper triangle) and output results ======
cat("\n=== Running Edge-wise ANOVA (upper triangle) ===\n")
edge_res <- run_two_way_anova_edgewise_chunked(edge_mat, df2, upper_idx, roi_names, chunk_size = 2000)

cat("Edge features with NA p-values (likely constant/degenerate):",
    sum(is.na(edge_res$p_Interaction)), "\n")

cat("Significant Interaction (edges, FDR<0.05):",
    sum(edge_res$p_Interaction_FDR < 0.05, na.rm = TRUE), "\n")

write.csv(edge_res, out_edge_csv, row.names = FALSE)
cat("Saved edge results:", out_edge_csv, "\n")

# ============================================================
# 8) Output significant results (degree and edge) - for Diagnosis, AgeGroup, Interaction
# ============================================================
cat("\n=== Generating significant results tables ===\n")

# Helper function: calculate direction and output
output_sig_results <- function(res_df, feature_mat, df_demo, effect_name, 
                                 feature_type = "degree", out_dir) {
  
  p_col <- paste0("p_", effect_name, "_FDR")
  sig_data <- res_df %>% filter(.data[[p_col]] < 0.05)
  
  if (nrow(sig_data) == 0) {
    cat("No significant", effect_name, "results for", feature_type, "(FDR<0.05)\n")
    return(NULL)
  }
  
  # Calculate means based on effect type
  sig_with_dir <- sig_data %>%
    rowwise() %>%
    mutate(
      # Get current feature index (feature name for degree, edge_id for edge)
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

  # Ungroup
  sig_with_dir <- sig_with_dir %>% ungroup() %>% arrange(.data[[p_col]])
  
  out_file <- file.path(out_dir, 
                        paste0("Significant_", effect_name, "_DK318_", 
                              feature_type, "_results.csv"))
  write.csv(sig_with_dir, out_file, row.names = FALSE)
  cat("Significant", effect_name, feature_type, "results (n =", 
      nrow(sig_with_dir), "):", out_file, "\n")
  
  return(sig_with_dir)
}

# Output Diagnosis main effect significant results
output_sig_results(degree_res, degree_mat, df2, "Diagnosis", "degree", out_dir)
output_sig_results(edge_res, edge_mat, df2, "Diagnosis", "edge", out_dir)

# Output AgeGroup main effect significant results
output_sig_results(degree_res, degree_mat, df2, "AgeGroup", "degree", out_dir)
output_sig_results(edge_res, edge_mat, df2, "AgeGroup", "edge", out_dir)

# Output Interaction significant results
output_sig_results(degree_res, degree_mat, df2, "Interaction", "degree", out_dir)
output_sig_results(edge_res, edge_mat, df2, "Interaction", "edge", out_dir)

cat("\nDONE.\n")
