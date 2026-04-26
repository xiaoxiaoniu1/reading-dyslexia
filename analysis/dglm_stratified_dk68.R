# ============================================================
# Stratified DGLM Analysis for DK-68 MIND degree and edges
#
# This is the R implementation of the stratified DK68 DGLM.
# It uses dglm() for both the mean model and dispersion model.
#
# Stratified analyses:
#   1) Within AgeGroup = Child / Adult:
#        y ~ Diagnosis * Sex + age_within_group
#   2) Within Diagnosis = DD / TD:
#        y ~ AgeGroup * Sex + age_within_group
#   3) Within Sex = Male / Female:
#        y ~ Diagnosis * AgeGroup + age_within_group
#
# Effect extraction:
#   For each 2 x 2 stratified model, marginal main effects and the
#   two-way interaction are extracted by emmeans with equal weights.
#
# Usage:
#   Rscript dglm_stratified_dk68.R [analysis_type]
#   analysis_type: both, degree, edge
#
# Optional flags:
#   --out-dir <path>
#   --demo-file <path>
#   --mind-combat-dir <path>
#   --network-dir <path>
#   --no-age-covariate
#   --run-plots true/false
# ============================================================

template_path <- "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/dglm_stratified_dk318.R"
script_txt <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

replace_once <- function(text, old, new, label) {
  if (!grepl(old, text, fixed = TRUE)) stop("Failed to locate ", label)
  sub(old, new, text, fixed = TRUE)
}

replacements <- list(
  c("DK-318", "DK-68"),
  c("dglm_stratified_dk318.R", "dglm_stratified_dk68.R"),
  c("MIND_DK318_combat", "MIND_DK68_combat"),
  c("MIND_DK318_DGLM_stratified", "MIND_DK68_DGLM_stratified"),
  c("plot_significant_DK318_stratified_brainmaps.py", "plot_significant_DK68_stratified_brainmaps.py"),
  c("violin_dk318_stratified_degree.R", "violin_dk68_stratified_degree.R"),
  c("DGLM_DK318_", "DGLM_DK68_"),
  c("DK318 combat", "DK68 combat"),
  c("DK318", "DK68")
)

for (pair in replacements) {
  script_txt <- gsub(pair[1], pair[2], script_txt, fixed = TRUE)
}

script_txt <- replace_once(
  script_txt,
  '    original_project = as.character(`original-project`),',
  '    subj_id = as.character(id),',
  'DK68 subject id line 1'
)
script_txt <- replace_once(
  script_txt,
  '    id_old = as.character(id_old),',
  '    subj_id = trimws(subj_id),',
  'DK68 subject id line 2'
)
script_txt <- replace_once(
  script_txt,
  '    subj_prefix = paste0(original_project, "_", id_old),\n',
  '',
  'DK68 subj_prefix line'
)
script_txt <- replace_once(
  script_txt,
  '    file_base = paste0(subj_prefix, "_MIND_DK68_combat"),',
  '    file_base = paste0(subj_id, "_combat"),',
  'DK68 file_base line'
)
script_txt <- replace_once(
  script_txt,
  '  filter(!is.na(original_project), !is.na(id_old), original_project != "", id_old != "")',
  '  filter(!is.na(id), subj_id != "", site != 3)',
  'DK68 demographic filter block'
)

eval(parse(text = script_txt), envir = globalenv())
