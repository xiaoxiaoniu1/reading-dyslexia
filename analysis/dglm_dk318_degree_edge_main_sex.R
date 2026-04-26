# Sex-stratified wrapper for DK318 DGLM.
# It runs Male and Female separately.
# Usage:
#   Rscript dglm_dk318_degree_edge_main_sex.R [both|degree|edge]

args <- commandArgs(trailingOnly = TRUE)
analysis_type_arg <- if (length(args) >= 1) args[1] else "both"

cmd <- commandArgs(FALSE)
this_file <- sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])
script_dir <- if (!is.na(this_file) && file.exists(this_file)) {
  dirname(normalizePath(this_file))
} else {
  "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis"
}
main_script <- file.path(script_dir, "dglm_dk318_degree_edge.R")
if (!file.exists(main_script)) {
  main_script <- "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/dglm_dk318_degree_edge.R"
}
if (!file.exists(main_script)) {
  stop("Cannot find dglm_dk318_degree_edge.R. Checked: ", main_script)
}

for (sex_label in c("Male", "Female")) {
  cat("\n========================================\n")
  cat("Running sex-stratified model:", sex_label, "\n")
  cat("========================================\n")

  Sys.setenv(DK318_DGLM_OUT_ROOT = "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_DGLM_main_sex")
  Sys.setenv(DK318_DGLM_INCLUDE_AGE = "false")
  Sys.setenv(DK318_DGLM_SEX_FILTER = sex_label)

  source(main_script)
}
