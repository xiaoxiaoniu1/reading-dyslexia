# Wrapper for DK68 DGLM.
# It sources dglm_dk68_degree_edge.R from the same directory when possible.
# Usage:
#   Rscript dglm_dk68_degree_edge_sensitivity.R [both|degree|edge]

args <- commandArgs(trailingOnly = TRUE)
analysis_type_arg <- if (length(args) >= 1) args[1] else "both"

cmd <- commandArgs(FALSE)
this_file <- sub("^--file=", "", cmd[grepl("^--file=", cmd)][1])
script_dir <- if (!is.na(this_file) && file.exists(this_file)) {
  dirname(normalizePath(this_file))
} else {
  "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis"
}
main_script <- file.path(script_dir, "dglm_dk68_degree_edge.R")
if (!file.exists(main_script)) {
  main_script <- "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/dglm_dk68_degree_edge.R"
}
if (!file.exists(main_script)) {
  stop("Cannot find dglm_dk68_degree_edge.R. Checked: ", main_script)
}

Sys.setenv(DK68_DGLM_OUT_ROOT = "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK68_DGLM_sensitivity")
Sys.setenv(DK68_DGLM_INCLUDE_AGE = "true")
Sys.setenv(DK68_DGLM_SEX_FILTER = "")

source(main_script)
