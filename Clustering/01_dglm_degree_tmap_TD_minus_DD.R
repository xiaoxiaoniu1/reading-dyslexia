# ============================================================
# DK318 MIND Degree DGLM Tmap Analysis
# Purpose: Extract TD-DD t-values for each ROI in Child and Adult groups
# Output: 318 × 2 Tmap matrix for ROI clustering
# Model: y ~ Diagnosis * AgeGroup + Sex
# ============================================================

cat("Starting DK318 MIND Degree Tmap extraction...\n\n")

packages <- c("readxl", "dplyr", "dglm", "emmeans")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

library(readxl); library(dplyr); library(dglm); library(emmeans)

# Paths
demo_file <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/all_data_cqt_mean_1.5.xlsx"
mind_combat_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_combat"
out_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Clustering_5"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load demographic data
cat("Loading demographic data...\n")
df <- read_excel(demo_file, sheet = "Sheet1") %>%
  mutate(
    original_project = as.character(`original-project`),
    id_old = as.character(id_old),
    subj_prefix = paste0(original_project, "_", id_old),
    file_base = paste0(subj_prefix, "_MIND_DK318_combat"),
    degree_file = file.path(mind_combat_dir, paste0(file_base, "_degree.csv")),
    Diagnosis = factor(ifelse(group_d_or_c == 0, "TD", "DD"), levels = c("TD", "DD")),
    AgeGroup = factor(ifelse(group_age == 1, "Adult", "Child"), levels = c("Child", "Adult")),
    Sex = factor(ifelse(sex == 1, "Male", "Female"), levels = c("Female", "Male")),
    has_file = file.exists(degree_file)
  ) %>%
  filter(!is.na(original_project), !is.na(id_old), original_project != "", id_old != "", has_file)

cat("Total subjects:", nrow(df), "\n\nCell counts:\n")
print(with(df, table(Diagnosis, AgeGroup, Sex)))

# Load degree data
cat("\nLoading degree data...\n")
read_degree_csv <- function(fp) {
  d <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE)
  vals <- as.numeric(d$degree)
  names(vals) <- as.character(d$ROI)
  vals
}

tmp_deg <- read_degree_csv(df$degree_file[1])
roi_names <- names(tmp_deg)
n_roi <- length(roi_names)
n_subj <- nrow(df)

degree_mat <- matrix(NA_real_, nrow = n_subj, ncol = n_roi, dimnames = list(df$file_base, roi_names))
for (s in seq_len(n_subj)) {
  degree_mat[s, ] <- read_degree_csv(df$degree_file[s])
}

cat("Degree matrix:", paste(dim(degree_mat), collapse = " × "), "\n\n")

# Helper functions
calc_raw_stats <- function(y, data) {
  cells <- expand.grid(Diagnosis = c("TD", "DD"), AgeGroup = c("Child", "Adult"), stringsAsFactors = FALSE)
  stats <- list()
  for (i in seq_len(nrow(cells))) {
    diag <- cells$Diagnosis[i]; age <- cells$AgeGroup[i]
    idx <- data$Diagnosis == diag & data$AgeGroup == age
    yy <- y[idx]; yy <- yy[is.finite(yy)]
    prefix <- paste0("raw_", tolower(age), "_", diag)
    stats[[paste0("raw_n_", tolower(age), "_", diag)]] <- length(yy)
    stats[[paste0("raw_mean_", tolower(age), "_", diag)]] <- if (length(yy) > 0) mean(yy) else NA_real_
    stats[[paste0("raw_sd_", tolower(age), "_", diag)]] <- if (length(yy) >= 2) sd(yy) else NA_real_
  }
  stats
}

check_cells <- function(data) {
  cells <- expand.grid(Diagnosis = c("TD", "DD"), AgeGroup = c("Child", "Adult"), stringsAsFactors = FALSE)
  for (i in seq_len(nrow(cells))) {
    idx <- data$Diagnosis == cells$Diagnosis[i] & data$AgeGroup == cells$AgeGroup[i]
    if (sum(idx) == 0) return(FALSE)
  }
  TRUE
}

check_sex <- function(data) {
  sex_levels <- unique(as.character(data$Sex))
  all(c("Female", "Male") %in% sex_levels)
}

# Fit DGLM and extract contrasts for one ROI
fit_one_roi <- function(y, data, roi_name, min_n = 10, var_eps = 1e-12) {
  ok <- is.finite(y) & !is.na(data$Diagnosis) & !is.na(data$AgeGroup) & !is.na(data$Sex)
  y2 <- y[ok]; d2 <- data[ok, , drop = FALSE]
  d2$Diagnosis <- factor(as.character(d2$Diagnosis), levels = c("TD", "DD"))
  d2$AgeGroup <- factor(as.character(d2$AgeGroup), levels = c("Child", "Adult"))
  d2$Sex <- factor(as.character(d2$Sex), levels = c("Female", "Male"))
  
  result <- list(feature = roi_name, n = length(y2), fit_status = "not_run")
  result <- c(result, calc_raw_stats(y2, d2))
  
  empty_stats <- list(
    estimate_child_TD_minus_DD = NA_real_, se_child_TD_minus_DD = NA_real_,
    t_child_TD_minus_DD = NA_real_, p_child_TD_minus_DD = NA_real_,
    estimate_adult_TD_minus_DD = NA_real_, se_adult_TD_minus_DD = NA_real_,
    t_adult_TD_minus_DD = NA_real_, p_adult_TD_minus_DD = NA_real_,
    estimate_dev_adult_minus_child = NA_real_, se_dev_adult_minus_child = NA_real_,
    t_dev_adult_minus_child = NA_real_, p_dev_adult_minus_child = NA_real_
  )
  result <- c(result, empty_stats)
  
  if (length(y2) < min_n) { result$fit_status <- "low_n"; return(as.data.frame(result, stringsAsFactors = FALSE)) }
  v <- var(y2)
  if (!is.finite(v) || v < var_eps) { result$fit_status <- "near_zero_variance"; return(as.data.frame(result, stringsAsFactors = FALSE)) }
  if (!check_cells(d2)) { result$fit_status <- "empty_diagnosis_age_cell"; return(as.data.frame(result, stringsAsFactors = FALSE)) }
  if (!check_sex(d2)) { result$fit_status <- "insufficient_sex_levels"; return(as.data.frame(result, stringsAsFactors = FALSE)) }
  
  d2$y <- y2
  fit_err <- NA_character_
  fit <- tryCatch(dglm::dglm(formula = y ~ Diagnosis * AgeGroup + Sex, dformula = ~ Diagnosis * AgeGroup + Sex,
                             family = gaussian(link = "identity"), data = d2),
                  error = function(e) { fit_err <<- conditionMessage(e); NULL })
  if (is.null(fit)) { result$fit_status <- paste0("dglm_error: ", fit_err); return(as.data.frame(result, stringsAsFactors = FALSE)) }
  
  fit_mean <- fit; class(fit_mean) <- "lm"
  emm <- tryCatch(emmeans::emmeans(fit_mean, ~ Diagnosis | AgeGroup, weights = "equal"), error = function(e) NULL)
  if (is.null(emm)) { result$fit_status <- "emmeans_error"; return(as.data.frame(result, stringsAsFactors = FALSE)) }
  
  con <- tryCatch(emmeans::contrast(emm, method = list(TD_minus_DD = c(1, -1)), by = "AgeGroup", adjust = "none"), error = function(e) NULL)
  if (is.null(con)) { result$fit_status <- "contrast_error"; return(as.data.frame(result, stringsAsFactors = FALSE)) }
  
  con_df <- as.data.frame(summary(con, infer = c(FALSE, TRUE), adjust = "none"))
  child_idx <- which(con_df$AgeGroup == "Child")
  if (length(child_idx) > 0) {
    result$estimate_child_TD_minus_DD <- con_df$estimate[child_idx[1]]
    result$se_child_TD_minus_DD <- con_df$SE[child_idx[1]]
    result$t_child_TD_minus_DD <- con_df$t.ratio[child_idx[1]]
    result$p_child_TD_minus_DD <- con_df$p.value[child_idx[1]]
  }
  adult_idx <- which(con_df$AgeGroup == "Adult")
  if (length(adult_idx) > 0) {
    result$estimate_adult_TD_minus_DD <- con_df$estimate[adult_idx[1]]
    result$se_adult_TD_minus_DD <- con_df$SE[adult_idx[1]]
    result$t_adult_TD_minus_DD <- con_df$t.ratio[adult_idx[1]]
    result$p_adult_TD_minus_DD <- con_df$p.value[adult_idx[1]]
  }
  
  emm_full <- tryCatch(emmeans::emmeans(fit_mean, ~ Diagnosis * AgeGroup, weights = "equal"), error = function(e) NULL)
  if (!is.null(emm_full)) {
    con_dev <- tryCatch(emmeans::contrast(emm_full, method = list(Adult_minus_Child_change = c(-1, 1, 1, -1)), adjust = "none"),
                        error = function(e) NULL)
    if (!is.null(con_dev)) {
      con_dev_df <- as.data.frame(summary(con_dev, infer = c(FALSE, TRUE), adjust = "none"))
      if (nrow(con_dev_df) > 0) {
        result$estimate_dev_adult_minus_child <- con_dev_df$estimate[1]
        result$se_dev_adult_minus_child <- con_dev_df$SE[1]
        result$t_dev_adult_minus_child <- con_dev_df$t.ratio[1]
        result$p_dev_adult_minus_child <- con_dev_df$p.value[1]
      }
    }
  }
  
  result$fit_status <- "ok"
  as.data.frame(result, stringsAsFactors = FALSE)
}

# Run DGLM for all ROIs
cat("Fitting DGLM for each ROI...\n")
results_list <- vector("list", n_roi)
for (k in seq_len(n_roi)) {
  if (k %% 50 == 0 || k == 1 || k == n_roi) cat("  ROI", k, "/", n_roi, "\n")
  results_list[[k]] <- fit_one_roi(y = degree_mat[, k], data = df, roi_name = roi_names[k], min_n = 10, var_eps = 1e-12)
}
results_full <- dplyr::bind_rows(results_list)

cat("\nFit status summary:\n")
print(table(results_full$fit_status))

# FDR correction
cat("\nApplying FDR correction...\n")
results_full$q_child_TD_minus_DD <- p.adjust(results_full$p_child_TD_minus_DD, method = "fdr")
results_full$q_adult_TD_minus_DD <- p.adjust(results_full$p_adult_TD_minus_DD, method = "fdr")
results_full$q_dev_adult_minus_child <- p.adjust(results_full$p_dev_adult_minus_child, method = "fdr")

# Save outputs
cat("\nSaving results...\n")
full_csv <- file.path(out_dir, "DGLM_DK318_degree_Tmap_TD_minus_DD_full.csv")
write.csv(results_full, full_csv, row.names = FALSE)
cat("Saved full results:", full_csv, "\n")

tmap_318x2 <- results_full %>% select(feature, t_child_TD_minus_DD, t_adult_TD_minus_DD)
tmap_csv <- file.path(out_dir, "DGLM_DK318_degree_Tmap_TD_minus_DD_318x2.csv")
write.csv(tmap_318x2, tmap_csv, row.names = FALSE)
cat("Saved Tmap 318×2 matrix:", tmap_csv, "\n")

safe_min <- function(x) { x <- x[is.finite(x)]; if (length(x) == 0) NA_real_ else min(x) }
summary_df <- data.frame(
  contrast = c("Child_TD_minus_DD", "Adult_TD_minus_DD", "Adult_minus_Child_change"),
  n_tested = c(sum(is.finite(results_full$p_child_TD_minus_DD)), sum(is.finite(results_full$p_adult_TD_minus_DD)),
               sum(is.finite(results_full$p_dev_adult_minus_child))),
  n_sig_FDR_0.05 = c(sum(is.finite(results_full$q_child_TD_minus_DD) & results_full$q_child_TD_minus_DD < 0.05),
                     sum(is.finite(results_full$q_adult_TD_minus_DD) & results_full$q_adult_TD_minus_DD < 0.05),
                     sum(is.finite(results_full$q_dev_adult_minus_child) & results_full$q_dev_adult_minus_child < 0.05)),
  min_p = c(safe_min(results_full$p_child_TD_minus_DD), safe_min(results_full$p_adult_TD_minus_DD),
            safe_min(results_full$p_dev_adult_minus_child)),
  min_q = c(safe_min(results_full$q_child_TD_minus_DD), safe_min(results_full$q_adult_TD_minus_DD),
            safe_min(results_full$q_dev_adult_minus_child)),
  stringsAsFactors = FALSE
)

summary_csv <- file.path(out_dir, "DGLM_DK318_degree_Tmap_TD_minus_DD_summary.csv")
write.csv(summary_df, summary_csv, row.names = FALSE)
cat("Saved summary:", summary_csv, "\n\n")

cat("==============================\nSUMMARY\n==============================\n\n")
cat("Model: y ~ Diagnosis * AgeGroup + Sex\n")
cat("Contrast direction: TD - DD (positive = TD > DD)\n\n")
print(summary_df, row.names = FALSE)
cat("\n==============================\nDONE\n==============================\n")
