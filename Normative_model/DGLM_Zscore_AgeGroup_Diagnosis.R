rm(list=ls())

suppressPackageStartupMessages(library(gamlss))

result_dir <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Results/Individual_area_volume'
zscore_file <- file.path(result_dir, 'All_Zscore_summary.csv')
clinical_file <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Datasets/Datasets-diseases/Clinical_vars_control.csv'
out_file <- file.path(result_dir, 'DGLM_Zscore_AgeGroup_Diagnosis.csv')
summary_file_raw <- file.path(result_dir, 'DGLM_Zscore_AgeGroup_Diagnosis_significant_summary_raw.csv')
summary_file_fdr <- file.path(result_dir, 'DGLM_Zscore_AgeGroup_Diagnosis_significant_summary_fdr.csv')

adult_age_cutoff <- 18

max_features_env <- suppressWarnings(as.integer(Sys.getenv('DGLM_MAX_FEATURES', '0')))
max_features <- ifelse(is.na(max_features_env) || max_features_env <= 0, Inf, max_features_env)

clean_id <- function(x) {
  gsub('[()_ ]', '', tolower(as.character(x)))
}

extract_id <- function(x) {
  m <- regmatches(x, regexpr('(?<=fs_subjects_all/)[^/]+', x, perl=TRUE))
  if (length(m) == 0 || m == '') {
    m <- regmatches(x, regexpr('[^/]+(?=/stats/)', x, perl=TRUE))
  }
  if (length(m) == 0 || m == '') m <- x
  clean_id(m)
}

make_agegroup <- function(age) {
  age <- as.numeric(age)
  grp <- ifelse(is.na(age), NA_character_, ifelse(age >= adult_age_cutoff, 'Adult', 'Child'))
  factor(grp, levels = c('Child', 'Adult'))
}

normalize_sex <- function(x) {
  x <- trimws(tolower(as.character(x)))
  x[x %in% c('', 'na', 'nan')] <- NA_character_
  x[x %in% c('male', 'm')] <- 'Male'
  x[x %in% c('female', 'f')] <- 'Female'
  factor(x, levels = c('Female', 'Male'))
}

feature_anatomy_group <- function(feature_name) {
  if (feature_name %in% c('GMV', 'mean_thickness', 'sGMV', 'TCV', 'total_surface_arrea', 'Ventricles', 'WMV')) {
    return('Global')
  }
  if (startsWith(feature_name, 'Left.') || startsWith(feature_name, 'Right.')) {
    return('Subcortical')
  }
  if (grepl('_thickness$', feature_name)) return('Cortical thickness')
  if (grepl('_area$', feature_name)) return('Cortical area')
  if (grepl('_volume$', feature_name)) return('Cortical volume')
  'Other'
}

report_significant_by_anatomy <- function(res_df, summary_file_raw = NULL, summary_file_fdr = NULL) {
  tests <- list(
    mu_agegroup = c('p_mu_agegroup', 'q_mu_agegroup'),
    mu_diagnosis = c('p_mu_diagnosis', 'q_mu_diagnosis'),
    mu_interaction = c('p_mu_interaction', 'q_mu_interaction'),
    sigma_agegroup = c('p_sigma_agegroup', 'q_sigma_agegroup'),
    sigma_diagnosis = c('p_sigma_diagnosis', 'q_sigma_diagnosis'),
    sigma_interaction = c('p_sigma_interaction', 'q_sigma_interaction')
  )
  anatomy_levels <- c('Subcortical', 'Global', 'Cortical area', 'Cortical thickness', 'Cortical volume', 'Other')
  groups <- vapply(res_df$Feature, feature_anatomy_group, character(1))
  summary_rows_raw <- list()
  summary_rows_fdr <- list()

  cat('\n====================\n')
  cat('Significant results by anatomy\n')
  cat('====================\n')

  for (test_name in names(tests)) {
    p_col <- tests[[test_name]][1]
    q_col <- tests[[test_name]][2]
    cat('\n[', test_name, ']\n', sep = '')
    for (mode in c('RAW', 'FDR')) {
      sig_col <- if (mode == 'RAW') p_col else q_col
      sig_idx <- which(!is.na(res_df[[sig_col]]) & res_df[[sig_col]] < 0.05)
      cat(mode, 'significant results\n')
      for (grp in anatomy_levels) {
        feats <- res_df$Feature[sig_idx][groups[sig_idx] == grp]
        feat_text <- if (length(feats) > 0) paste(feats, collapse = '; ') else ''
        cat(' - ', grp, ' (', length(feats), '): ', if (nzchar(feat_text)) feat_text else 'None', '\n', sep = '')
        row_df <- data.frame(
          test = test_name,
          anatomy_group = grp,
          n_significant = length(feats),
          features = feat_text,
          stringsAsFactors = FALSE
        )
        if (mode == 'RAW') {
          summary_rows_raw[[length(summary_rows_raw) + 1]] <- row_df
        } else {
          summary_rows_fdr[[length(summary_rows_fdr) + 1]] <- row_df
        }
      }
    }
  }

  if (!is.null(summary_file_raw)) {
    summary_df_raw <- do.call(rbind, summary_rows_raw)
    write.csv(summary_df_raw, summary_file_raw, row.names = FALSE)
    cat('Saved RAW summary:', summary_file_raw, '\n')
  }
  if (!is.null(summary_file_fdr)) {
    summary_df_fdr <- do.call(rbind, summary_rows_fdr)
    write.csv(summary_df_fdr, summary_file_fdr, row.names = FALSE)
    cat('Saved FDR summary:', summary_file_fdr, '\n')
  }
}

lr_pvalue <- function(full, reduced) {
  if (is.null(full) || is.null(reduced)) return(NA_real_)
  p <- tryCatch(gamlss::LR.test(full, reduced)$p, error = function(e) numeric(0))
  if (length(p) == 1 && !is.na(p)) return(as.numeric(p))
  dev_full <- suppressWarnings(as.numeric(full$G.deviance))
  dev_red <- suppressWarnings(as.numeric(reduced$G.deviance))
  df_full <- suppressWarnings(as.numeric(full$df.fit))
  df_red <- suppressWarnings(as.numeric(reduced$df.fit))
  if (!is.finite(dev_full) || !is.finite(dev_red) || !is.finite(df_full) || !is.finite(df_red)) return(NA_real_)
  df_diff <- df_full - df_red
  ddev <- dev_red - dev_full
  if (!is.finite(df_diff) || df_diff <= 0) return(NA_real_)
  if (!is.finite(ddev) || ddev < 0) ddev <- 0
  stats::pchisq(ddev, df = df_diff, lower.tail = FALSE)
}

collapse_matrix_by_id <- function(df, id_vec, label = 'matrix') {
  id_vec <- as.character(id_vec)
  keep <- !is.na(id_vec) & nzchar(id_vec)
  dropped_n <- sum(!keep)
  if (dropped_n > 0) cat('[DEBUG]', label, '- dropped', dropped_n, 'rows with empty IDs\n')
  df <- df[keep, , drop = FALSE]
  id_vec <- id_vec[keep]
  if (nrow(df) == 0) {
    cat('[DEBUG]', label, '- no rows left after ID filtering\n')
    return(df)
  }

  rownames(df) <- id_vec
  dup_n <- sum(duplicated(id_vec))
  if (dup_n == 0) return(df)

  cat('[DEBUG]', label, '- collapsing', dup_n, 'duplicate rows across', length(unique(id_vec[duplicated(id_vec)])), 'IDs\n')
  out <- matrix(NA_real_, nrow = length(unique(id_vec)), ncol = ncol(df))
  rownames(out) <- unique(id_vec)
  colnames(out) <- colnames(df)

  for (j in seq_len(ncol(df))) {
    x <- suppressWarnings(as.numeric(df[[j]]))
    vals <- tapply(x, id_vec, function(v) {
      if (all(is.na(v))) NA_real_ else mean(v, na.rm = TRUE)
    })
    out[names(vals), j] <- as.numeric(vals)
  }

  as.data.frame(out, check.names = FALSE)
}

build_z_matrix_from_rds <- function(dir_path, max_features = Inf) {
  patterns <- c('_loop_our_model_individual[.]rds$', '_model_new[.]rds$')
  rds_files <- character(0)
  for (pat in patterns) {
    rds_files <- list.files(dir_path, pattern = pat, recursive = TRUE, full.names = TRUE)
    if (length(rds_files) > 0) break
  }
  if (length(rds_files) == 0) stop('No RDS files found under: ', dir_path)
  if (is.finite(max_features)) rds_files <- rds_files[seq_len(min(length(rds_files), max_features))]
  cat('Found', length(rds_files), 'rds files\n')

  Z_all <- NULL
  for (f in rds_files) {
    cat('[DEBUG] Reading:', f, '\n')
    res <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(res)) next
    feature_name <- res$i
    if (is.null(feature_name) || is.na(feature_name) || feature_name == '') next

    z_df <- NULL
    if (!is.null(res$Zscore) && length(res$Zscore) > 0) {
      if (!is.null(res$Zscore[[feature_name]])) z_df <- res$Zscore[[feature_name]]
      if (is.null(z_df) && length(res$Zscore) == 1) z_df <- res$Zscore[[1]]
    }
    if (is.null(z_df) || nrow(z_df) == 0) next
    if (ncol(z_df) != 1) z_df <- z_df[, 1, drop = FALSE]
    colnames(z_df) <- feature_name
    cat('[DEBUG]', feature_name, '- raw rows:', nrow(z_df), ' unique raw rownames:', length(unique(rownames(z_df))), '\n')
    id_clean <- vapply(rownames(z_df), extract_id, character(1))
    z_df <- collapse_matrix_by_id(z_df, id_clean, label = paste('feature', feature_name))
    cat('[DEBUG]', feature_name, '- rows after ID normalization:', nrow(z_df), '\n')

    if (is.null(Z_all)) {
      Z_all <- z_df
    } else {
      prev_n <- nrow(Z_all)
      Z_all <- merge(
        data.frame(.__id__ = rownames(Z_all), Z_all, check.names = FALSE),
        data.frame(.__id__ = rownames(z_df), z_df, check.names = FALSE),
        by = '.__id__',
        all = TRUE,
        sort = FALSE
      )
      rownames(Z_all) <- Z_all$.__id__
      Z_all$.__id__ <- NULL
      Z_all <- collapse_matrix_by_id(Z_all, rownames(Z_all), label = paste('merged after', feature_name))
      cat('[DEBUG]', feature_name, '- matrix rows:', prev_n, '->', nrow(Z_all), ' cols:', ncol(Z_all), '\n')
    }
    cat('Done:', feature_name, '\n')
  }
  if (is.null(Z_all) || nrow(Z_all) == 0 || ncol(Z_all) == 0) stop('Failed to build Z-score matrix')
  cat('[DEBUG] Final Z matrix dimensions:', nrow(Z_all), 'x', ncol(Z_all), '\n')
  Z_all
}

if (file.exists(zscore_file)) {
  cat('[DEBUG] Loading existing Z-score CSV:', zscore_file, '\n')
  Z_all <- read.csv(zscore_file, row.names = 1, check.names = FALSE)
  cat('[DEBUG] CSV dimensions before ID normalization:', nrow(Z_all), 'x', ncol(Z_all), '\n')
  id_clean <- clean_id(rownames(Z_all))
  Z_all <- collapse_matrix_by_id(Z_all, id_clean, label = 'zscore csv')
  cat('[DEBUG] CSV dimensions after ID normalization:', nrow(Z_all), 'x', ncol(Z_all), '\n')
} else {
  cat('[DEBUG] Z-score CSV not found, rebuilding from RDS files\n')
  Z_all <- build_z_matrix_from_rds(result_dir, max_features = max_features)
}

Z_all[is.infinite(as.matrix(Z_all))] <- NA

clin <- read.csv(clinical_file)
clin$ID_clean <- clean_id(clin$ID)
if (anyDuplicated(clin$ID_clean) > 0) {
  cat('[DEBUG] Clinical file has', sum(duplicated(clin$ID_clean)), 'duplicate cleaned IDs\n')
}
rownames(clin) <- clin$ID_clean

common_ids <- intersect(rownames(Z_all), rownames(clin))
cat('[DEBUG] Z matrix subjects:', nrow(Z_all), ' Clinical subjects:', nrow(clin), ' Overlap:', length(common_ids), '\n')
if (length(common_ids) == 0) stop('No overlapping subject IDs between Zscore and clinical')

Z <- Z_all[common_ids, , drop = FALSE]
clin2 <- clin[common_ids, , drop = FALSE]
cat('[DEBUG] Analysis matrix dimensions:', nrow(Z), 'x', ncol(Z), '\n')

age_col <- if ('Age_y' %in% colnames(clin2)) 'Age_y' else 'Age'
clin2$AgeGroup <- make_agegroup(clin2[[age_col]])
clin2$AgeGroup <- droplevels(as.factor(clin2$AgeGroup))
clin2$Diagnosis <- droplevels(as.factor(clin2$Diagnosis))
clin2$Diagnosis <- factor(as.character(clin2$Diagnosis), levels = c('TD', 'DD'))

features <- colnames(Z)
if (is.finite(max_features)) features <- features[seq_len(min(length(features), max_features))]

out <- vector('list', length(features))

for (k in seq_along(features)) {
  feat <- features[k]
  y <- suppressWarnings(as.numeric(Z[, feat]))
  df <- data.frame(
    Z_score = y,
    AgeGroup = clin2$AgeGroup,
    Diagnosis = clin2$Diagnosis
  )
  df <- df[is.finite(df$Z_score) & !is.na(df$AgeGroup) & !is.na(df$Diagnosis), , drop = FALSE]
  df <- df[df$Diagnosis %in% c('TD', 'DD'), , drop = FALSE]
  df$AgeGroup <- droplevels(as.factor(df$AgeGroup))
  df$Diagnosis <- droplevels(as.factor(df$Diagnosis))
  if (nrow(df) < 20 || nlevels(df$AgeGroup) < 2 || nlevels(df$Diagnosis) < 2) {
    cat('[DEBUG] Skipping', feat, '- N =', nrow(df), ' AgeGroup levels =', nlevels(df$AgeGroup), ' Diagnosis levels =', nlevels(df$Diagnosis), '\n')
    out[[k]] <- data.frame(
      Feature = feat,
      N = nrow(df),
      p_mu_agegroup = NA_real_,
      p_mu_diagnosis = NA_real_,
      p_mu_interaction = NA_real_,
      p_sigma_agegroup = NA_real_,
      p_sigma_diagnosis = NA_real_,
      p_sigma_interaction = NA_real_
    )
    next
  }

  fit_full <- tryCatch(
    gamlss(
      Z_score ~ AgeGroup * Diagnosis,
      sigma.formula = ~ AgeGroup * Diagnosis,
      family = NO,
      data = df,
      trace = FALSE
    ),
    error = function(e) NULL
  )

  if (is.null(fit_full)) {
    cat('[DEBUG] Full model failed for', feat, 'with N =', nrow(df), '\n')
    out[[k]] <- data.frame(
      Feature = feat,
      N = nrow(df),
      p_mu_agegroup = NA_real_,
      p_mu_diagnosis = NA_real_,
      p_mu_interaction = NA_real_,
      p_sigma_agegroup = NA_real_,
      p_sigma_diagnosis = NA_real_,
      p_sigma_interaction = NA_real_
    )
    next
  }

  fit_mu_no_age <- tryCatch(
    gamlss(
      Z_score ~ Diagnosis,
      sigma.formula = ~ AgeGroup * Diagnosis,
      family = NO,
      data = df,
      trace = FALSE
    ),
    error = function(e) NULL
  )
  fit_mu_no_diag <- tryCatch(
    gamlss(
      Z_score ~ AgeGroup,
      sigma.formula = ~ AgeGroup * Diagnosis,
      family = NO,
      data = df,
      trace = FALSE
    ),
    error = function(e) NULL
  )
  fit_mu_no_int <- tryCatch(
    gamlss(
      Z_score ~ AgeGroup + Diagnosis,
      sigma.formula = ~ AgeGroup * Diagnosis,
      family = NO,
      data = df,
      trace = FALSE
    ),
    error = function(e) NULL
  )
  fit_sigma_no_age <- tryCatch(
    gamlss(
      Z_score ~ AgeGroup * Diagnosis,
      sigma.formula = ~ Diagnosis,
      family = NO,
      data = df,
      trace = FALSE
    ),
    error = function(e) NULL
  )
  fit_sigma_no_diag <- tryCatch(
    gamlss(
      Z_score ~ AgeGroup * Diagnosis,
      sigma.formula = ~ AgeGroup,
      family = NO,
      data = df,
      trace = FALSE
    ),
    error = function(e) NULL
  )
  fit_sigma_no_int <- tryCatch(
    gamlss(
      Z_score ~ AgeGroup * Diagnosis,
      sigma.formula = ~ AgeGroup + Diagnosis,
      family = NO,
      data = df,
      trace = FALSE
    ),
    error = function(e) NULL
  )

  p_mu_age <- lr_pvalue(fit_full, fit_mu_no_age)
  p_mu_diag <- lr_pvalue(fit_full, fit_mu_no_diag)
  p_mu_int <- lr_pvalue(fit_full, fit_mu_no_int)
  p_sigma_age <- lr_pvalue(fit_full, fit_sigma_no_age)
  p_sigma_diag <- lr_pvalue(fit_full, fit_sigma_no_diag)
  p_sigma_int <- lr_pvalue(fit_full, fit_sigma_no_int)

  out[[k]] <- data.frame(
    Feature = feat,
    N = nrow(df),
    p_mu_agegroup = p_mu_age,
    p_mu_diagnosis = p_mu_diag,
    p_mu_interaction = p_mu_int,
    p_sigma_agegroup = p_sigma_age,
    p_sigma_diagnosis = p_sigma_diag,
    p_sigma_interaction = p_sigma_int
  )

  if (k %% 25 == 0) cat('Processed', k, '/', length(features), '\n')
}

res_df <- do.call(rbind, out)
res_df$q_mu_agegroup <- p.adjust(res_df$p_mu_agegroup, method = 'fdr')
res_df$q_mu_diagnosis <- p.adjust(res_df$p_mu_diagnosis, method = 'fdr')
res_df$q_mu_interaction <- p.adjust(res_df$p_mu_interaction, method = 'fdr')
res_df$q_sigma_agegroup <- p.adjust(res_df$p_sigma_agegroup, method = 'fdr')
res_df$q_sigma_diagnosis <- p.adjust(res_df$p_sigma_diagnosis, method = 'fdr')
res_df$q_sigma_interaction <- p.adjust(res_df$p_sigma_interaction, method = 'fdr')

write.csv(res_df, out_file, row.names = FALSE)
cat('Saved:', out_file, '\n')
report_significant_by_anatomy(
  res_df,
  summary_file_raw = summary_file_raw,
  summary_file_fdr = summary_file_fdr
)
