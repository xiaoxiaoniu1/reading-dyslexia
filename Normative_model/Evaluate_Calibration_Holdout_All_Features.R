rm(list = ls())
library(dplyr)

base_dir <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Results/Newsite'
out_file <- file.path(base_dir, 'holdout_feature_calibration_evaluation.csv')
count_file <- file.path(base_dir, 'holdout_feature_calibration_evaluation_counts.csv')
invalid_file <- file.path(base_dir, 'holdout_feature_calibration_invalid_summaries.csv')

target_dirs <- c(
  'Global_feature',
  'aseg.vol.table',
  'lh.aparc.volume.table', 'rh.aparc.volume.table',
  'lh.aparc.thickness.table', 'rh.aparc.thickness.table',
  'lh.aparc.area.table', 'rh.aparc.area.table'
)

score_to_target <- function(z_mean, z_var) {
  abs(z_mean) + abs(z_var - 1)
}

read_one_summary <- function(csv_path, base_dir) {
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  rel_path <- sub(paste0('^', normalizePath(base_dir, winslash = '/'), '/?'), '', normalizePath(csv_path, winslash = '/'))
  feature_prefix <- sub('_holdout_zscore_summary\\.csv$', '', basename(csv_path))
  group_name <- basename(dirname(csv_path))

  if (!all(c('Model', 'N', 'Z_mean', 'Z_sd', 'Z_var') %in% colnames(df))) {
    return(list(
      valid = NULL,
      invalid = data.frame(
        group = group_name,
        feature_file_prefix = feature_prefix,
        summary_file = rel_path,
        reason = 'missing_required_columns',
        stringsAsFactors = FALSE
      )
    ))
  }

  uncal <- df %>% filter(Model == 'uncalibrated')
  cal <- df %>% filter(Model == 'calibrated')
  if (nrow(uncal) != 1 || nrow(cal) != 1) {
    return(list(
      valid = NULL,
      invalid = data.frame(
        group = group_name,
        feature_file_prefix = feature_prefix,
        summary_file = rel_path,
        reason = 'missing_uncalibrated_or_calibrated_row',
        stringsAsFactors = FALSE
      )
    ))
  }

  required_values <- c(uncal$Z_mean, uncal$Z_sd, uncal$Z_var, cal$Z_mean, cal$Z_sd, cal$Z_var)
  if (any(!is.finite(required_values))) {
    return(list(
      valid = NULL,
      invalid = data.frame(
        group = group_name,
        feature_file_prefix = feature_prefix,
        summary_file = rel_path,
        reason = 'non_finite_z_statistics',
        stringsAsFactors = FALSE
      )
    ))
  }

  uncal_dist <- score_to_target(uncal$Z_mean, uncal$Z_var)
  cal_dist <- score_to_target(cal$Z_mean, cal$Z_var)

  effect_vs_target <- if (cal_dist < uncal_dist) {
    'improved'
  } else if (cal_dist > uncal_dist) {
    'worsened'
  } else {
    'unchanged'
  }

  list(
    valid = data.frame(
      group = group_name,
      feature_file_prefix = feature_prefix,
      summary_file = rel_path,
      N_holdout = uncal$N,
      uncalibrated_Z_mean = uncal$Z_mean,
      uncalibrated_Z_sd = uncal$Z_sd,
      uncalibrated_Z_var = uncal$Z_var,
      calibrated_Z_mean = cal$Z_mean,
      calibrated_Z_sd = cal$Z_sd,
      calibrated_Z_var = cal$Z_var,
      uncalibrated_distance_to_target = uncal_dist,
      calibrated_distance_to_target = cal_dist,
      effect_vs_target = effect_vs_target,
      uncalibrated_exact_mean0_var1 = uncal$Z_mean == 0 & uncal$Z_var == 1,
      calibrated_exact_mean0_var1 = cal$Z_mean == 0 & cal$Z_var == 1,
      uncalibrated_near_mean0_var1_tol0.10 = abs(uncal$Z_mean) <= 0.10 & abs(uncal$Z_var - 1) <= 0.10,
      calibrated_near_mean0_var1_tol0.10 = abs(cal$Z_mean) <= 0.10 & abs(cal$Z_var - 1) <= 0.10,
      uncalibrated_near_mean0_var1_tol0.20 = abs(uncal$Z_mean) <= 0.20 & abs(uncal$Z_var - 1) <= 0.20,
      calibrated_near_mean0_var1_tol0.20 = abs(cal$Z_mean) <= 0.20 & abs(cal$Z_var - 1) <= 0.20,
      stringsAsFactors = FALSE
    ),
    invalid = NULL
  )
}

summary_files <- unlist(lapply(target_dirs, function(dir_name) {
  dir_path <- file.path(base_dir, dir_name)
  if (!dir.exists(dir_path)) return(character(0))
  list.files(dir_path, pattern = '_holdout_zscore_summary\\.csv$', full.names = TRUE)
}), use.names = FALSE)

cat('Found holdout summary files:', length(summary_files), '\n')

parsed_results <- lapply(summary_files, read_one_summary, base_dir = base_dir)
feature_eval <- bind_rows(lapply(parsed_results, `[[`, 'valid'))
invalid_summary <- bind_rows(lapply(parsed_results, `[[`, 'invalid'))

if (nrow(feature_eval) == 0) {
  stop('No valid *_holdout_zscore_summary.csv files were found.')
}

feature_eval <- feature_eval %>% arrange(group, calibrated_distance_to_target, feature_file_prefix)
write.csv(feature_eval, out_file, row.names = FALSE)

if (nrow(invalid_summary) > 0) {
  invalid_summary <- invalid_summary %>% arrange(group, feature_file_prefix)
  write.csv(invalid_summary, invalid_file, row.names = FALSE)
} else {
  write.csv(data.frame(), invalid_file, row.names = FALSE)
}

count_summary <- bind_rows(
  data.frame(metric = 'uncalibrated_exact_mean0_var1', count_true = sum(feature_eval$uncalibrated_exact_mean0_var1), count_false = sum(!feature_eval$uncalibrated_exact_mean0_var1)),
  data.frame(metric = 'calibrated_exact_mean0_var1', count_true = sum(feature_eval$calibrated_exact_mean0_var1), count_false = sum(!feature_eval$calibrated_exact_mean0_var1)),
  data.frame(metric = 'uncalibrated_near_mean0_var1_tol0.10', count_true = sum(feature_eval$uncalibrated_near_mean0_var1_tol0.10), count_false = sum(!feature_eval$uncalibrated_near_mean0_var1_tol0.10)),
  data.frame(metric = 'calibrated_near_mean0_var1_tol0.10', count_true = sum(feature_eval$calibrated_near_mean0_var1_tol0.10), count_false = sum(!feature_eval$calibrated_near_mean0_var1_tol0.10)),
  data.frame(metric = 'uncalibrated_near_mean0_var1_tol0.20', count_true = sum(feature_eval$uncalibrated_near_mean0_var1_tol0.20), count_false = sum(!feature_eval$uncalibrated_near_mean0_var1_tol0.20)),
  data.frame(metric = 'calibrated_near_mean0_var1_tol0.20', count_true = sum(feature_eval$calibrated_near_mean0_var1_tol0.20), count_false = sum(!feature_eval$calibrated_near_mean0_var1_tol0.20)),
  data.frame(metric = 'effect_improved', count_true = sum(feature_eval$effect_vs_target == 'improved'), count_false = NA),
  data.frame(metric = 'effect_worsened', count_true = sum(feature_eval$effect_vs_target == 'worsened'), count_false = NA),
  data.frame(metric = 'effect_unchanged', count_true = sum(feature_eval$effect_vs_target == 'unchanged'), count_false = NA),
  data.frame(metric = 'invalid_summary_files', count_true = nrow(invalid_summary), count_false = NA)
)

write.csv(count_summary, count_file, row.names = FALSE)

cat('Saved files:\n')
cat(out_file, '\n')
cat(count_file, '\n')
cat(invalid_file, '\n')
cat('Total evaluated features:', nrow(feature_eval), '\n')
cat('Total skipped invalid summaries:', nrow(invalid_summary), '\n')
