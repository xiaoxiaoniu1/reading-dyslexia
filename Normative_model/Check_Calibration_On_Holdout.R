rm(list = ls())
library(dplyr)

result_dir <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Results/Disease_area_volume'
application_clinical <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Datasets/Datasets-diseases/Clinical_vars_control.csv'
out_dir <- file.path(result_dir, 'calibration_check')

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

normalize_id <- function(x) {
  x <- as.character(x)
  x <- tolower(trimws(x))
  x <- gsub('[()_ ]', '', x)
  x
}

extract_subject_id <- function(df) {
  id <- rep(NA_character_, nrow(df))

  if ('ID' %in% colnames(df)) {
    id <- as.character(df$ID)
  }

  if (all(c('Freesurfer_Path2', 'Freesurfer_Path3') %in% colnames(df))) {
    tmp <- paste0(df$Freesurfer_Path2, df$Freesurfer_Path3)
    use_idx <- is.na(id) | id == ''
    id[use_idx] <- tmp[use_idx]
  }

  if (all(c('Freesufer_Path2', 'Freesufer_Path3') %in% colnames(df))) {
    tmp <- paste0(df$Freesufer_Path2, df$Freesufer_Path3)
    use_idx <- is.na(id) | id == ''
    id[use_idx] <- tmp[use_idx]
  }

  id <- ifelse(grepl('fs_subjects_all/', id),
               sub('.*fs_subjects_all/([^/]+)/stats/.*', '\\1', id),
               id)

  normalize_id(id)
}

load_subject_sets <- function() {
  app_df <- read.csv(application_clinical, header = TRUE, stringsAsFactors = FALSE)
  app_df$subject_id <- extract_subject_id(app_df)

  if (!'Diagnosis' %in% colnames(app_df)) {
    stop('Clinical file must contain a Diagnosis column.')
  }

  holdout_td_ids <- unique(app_df$subject_id[toupper(trimws(app_df$Diagnosis)) == 'TD'])

  list(
    holdout_td = holdout_td_ids,
    holdout_td_n = length(holdout_td_ids)
  )
}

read_feature_stats <- function(rds_file, holdout_ids) {
  res <- readRDS(rds_file)

  if (is.null(res$i) || is.null(res$Zscore) || is.null(res$Zscore[[res$i]])) {
    return(NULL)
  }

  feature_name <- res$i
  z_df <- res$Zscore[[feature_name]]
  if (is.null(z_df) || nrow(z_df) == 0) {
    return(NULL)
  }

  z_df <- data.frame(subject_key = rownames(z_df), Z_score = z_df[, 1], stringsAsFactors = FALSE)
  z_df$subject_id <- normalize_id(sub('.*fs_subjects_all/([^/]+)/stats/.*', '\\1', z_df$subject_key))
  z_df$subject_id[!grepl('fs_subjects_all/', z_df$subject_key)] <- normalize_id(z_df$subject_key[!grepl('fs_subjects_all/', z_df$subject_key)])

  if (!is.null(res$all_data) && nrow(res$all_data) == nrow(z_df)) {
    all_data <- res$all_data
    z_df$Diagnosis <- all_data$Diagnosis
    if ('Age' %in% colnames(all_data)) z_df$Age <- all_data$Age
    if ('Sex' %in% colnames(all_data)) z_df$Sex <- all_data$Sex
    if ('Site_ZZZ' %in% colnames(all_data)) z_df$Site_ZZZ <- all_data$Site_ZZZ
  }

  holdout_df <- z_df %>% filter(subject_id %in% holdout_ids)
  if (nrow(holdout_df) == 0) {
    return(NULL)
  }

  z_values <- holdout_df$Z_score
  z_values <- z_values[is.finite(z_values)]
  if (length(z_values) == 0) {
    return(NULL)
  }

  mean_test_p <- tryCatch(t.test(z_values, mu = 0)$p.value, error = function(e) NA_real_)

  data.frame(
    feature = feature_name,
    n = length(z_values),
    mean_z = mean(z_values),
    sd_z = sd(z_values),
    median_z = median(z_values),
    mad_z = mad(z_values),
    min_z = min(z_values),
    max_z = max(z_values),
    mean_test_p = mean_test_p,
    abs_mean = abs(mean(z_values)),
    abs_sd_minus_1 = abs(sd(z_values) - 1),
    approx_pass = abs(mean(z_values)) < 0.1 & abs(sd(z_values) - 1) < 0.1,
    stringsAsFactors = FALSE
  )
}

subject_sets <- load_subject_sets()
cat('Holdout TD subjects:', subject_sets$holdout_td_n, '\n')

rds_files <- list.files(result_dir, pattern = '_model_new[.]rds$', recursive = TRUE, full.names = TRUE)
cat('Found', length(rds_files), 'application result files\n')

feature_stats <- bind_rows(lapply(rds_files, read_feature_stats, holdout_ids = subject_sets$holdout_td))

if (nrow(feature_stats) == 0) {
  stop('No holdout TD Z-score results were found. Please confirm the disease application step has produced *_model_new.rds files.')
}

feature_stats <- feature_stats %>% arrange(abs_mean, abs_sd_minus_1)

overall_summary <- data.frame(
  holdout_td_n = subject_sets$holdout_td_n,
  feature_n = nrow(feature_stats),
  mean_of_feature_means = mean(feature_stats$mean_z, na.rm = TRUE),
  median_abs_mean = median(feature_stats$abs_mean, na.rm = TRUE),
  mean_of_feature_sds = mean(feature_stats$sd_z, na.rm = TRUE),
  median_abs_sd_minus_1 = median(feature_stats$abs_sd_minus_1, na.rm = TRUE),
  pass_rate_mean_sd_0.1 = mean(feature_stats$approx_pass, na.rm = TRUE),
  stringsAsFactors = FALSE
)

write.csv(feature_stats,
          file.path(out_dir, 'holdout_TD_feature_Zscore_summary.csv'),
          row.names = FALSE)
write.csv(overall_summary,
          file.path(out_dir, 'holdout_TD_overall_summary.csv'),
          row.names = FALSE)
write.csv(data.frame(subject_id = sort(subject_sets$holdout_td), stringsAsFactors = FALSE),
          file.path(out_dir, 'holdout_TD_subjects.csv'),
          row.names = FALSE)

cat('Saved files:\n')
cat(file.path(out_dir, 'holdout_TD_feature_Zscore_summary.csv'), '\n')
cat(file.path(out_dir, 'holdout_TD_overall_summary.csv'), '\n')
cat(file.path(out_dir, 'holdout_TD_subjects.csv'), '\n')
