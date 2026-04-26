rm(list = ls())
library(dplyr)

base_dir <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Results/Newsite'
in_file <- file.path(base_dir, 'holdout_feature_calibration_evaluation.csv')
out_dir <- file.path(base_dir, 'holdout_feature_target_lists')

if (!file.exists(in_file)) {
  stop('Missing evaluation file: ', in_file, '\nPlease run Evaluate_Calibration_Holdout_All_Features.R first.')
}
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

df <- read.csv(in_file, stringsAsFactors = FALSE)

write.csv(df %>% filter(uncalibrated_exact_mean0_var1),
          file.path(out_dir, 'uncalibrated_exact_mean0_var1_features.csv'),
          row.names = FALSE)
write.csv(df %>% filter(calibrated_exact_mean0_var1),
          file.path(out_dir, 'calibrated_exact_mean0_var1_features.csv'),
          row.names = FALSE)

write.csv(df %>% filter(uncalibrated_near_mean0_var1_tol0.10),
          file.path(out_dir, 'uncalibrated_near_mean0_var1_tol0.10_features.csv'),
          row.names = FALSE)
write.csv(df %>% filter(calibrated_near_mean0_var1_tol0.10),
          file.path(out_dir, 'calibrated_near_mean0_var1_tol0.10_features.csv'),
          row.names = FALSE)

write.csv(df %>% filter(uncalibrated_near_mean0_var1_tol0.20),
          file.path(out_dir, 'uncalibrated_near_mean0_var1_tol0.20_features.csv'),
          row.names = FALSE)
write.csv(df %>% filter(calibrated_near_mean0_var1_tol0.20),
          file.path(out_dir, 'calibrated_near_mean0_var1_tol0.20_features.csv'),
          row.names = FALSE)

write.csv(df %>% filter(effect_vs_target == 'improved'),
          file.path(out_dir, 'features_improved_after_calibration.csv'),
          row.names = FALSE)
write.csv(df %>% filter(effect_vs_target == 'worsened'),
          file.path(out_dir, 'features_worsened_after_calibration.csv'),
          row.names = FALSE)
write.csv(df %>% filter(effect_vs_target == 'unchanged'),
          file.path(out_dir, 'features_unchanged_after_calibration.csv'),
          row.names = FALSE)

summary_df <- data.frame(
  metric = c(
    'uncalibrated_exact_mean0_var1',
    'calibrated_exact_mean0_var1',
    'uncalibrated_near_mean0_var1_tol0.10',
    'calibrated_near_mean0_var1_tol0.10',
    'uncalibrated_near_mean0_var1_tol0.20',
    'calibrated_near_mean0_var1_tol0.20',
    'effect_improved',
    'effect_worsened',
    'effect_unchanged'
  ),
  count = c(
    sum(df$uncalibrated_exact_mean0_var1),
    sum(df$calibrated_exact_mean0_var1),
    sum(df$uncalibrated_near_mean0_var1_tol0.10),
    sum(df$calibrated_near_mean0_var1_tol0.10),
    sum(df$uncalibrated_near_mean0_var1_tol0.20),
    sum(df$calibrated_near_mean0_var1_tol0.20),
    sum(df$effect_vs_target == 'improved'),
    sum(df$effect_vs_target == 'worsened'),
    sum(df$effect_vs_target == 'unchanged')
  ),
  stringsAsFactors = FALSE
)

write.csv(summary_df,
          file.path(out_dir, 'target_feature_list_summary.csv'),
          row.names = FALSE)

cat('Saved target lists to:\n')
cat(out_dir, '\n')
