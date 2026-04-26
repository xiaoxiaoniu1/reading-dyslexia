template_path <- "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/dglm_dk318_degree_edge.R"
script_txt <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

replace_once <- function(text, old, new, label) {
  if (!grepl(old, text, fixed = TRUE)) stop("Failed to locate ", label)
  sub(old, new, text, fixed = TRUE)
}

replacements <- list(
  c("DK-318", "DK-68"),
  c("dglm_dk318_degree_edge.R", "dglm_dk68_degree_edge.R"),
  c("MIND_DK318_combat", "MIND_DK68_combat"),
  c("MIND_DK318_DGLM", "MIND_DK68_DGLM"),
  c("DGLM_DK318_", "DGLM_DK68_"),
  c("DK318_DGLM_OUT_ROOT", "DK68_DGLM_OUT_ROOT"),
  c("DK318_DGLM_INCLUDE_AGE", "DK68_DGLM_INCLUDE_AGE"),
  c("DK318_DGLM_SEX_FILTER", "DK68_DGLM_SEX_FILTER"),
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
  '  filter(\n    !is.na(original_project), !is.na(id_old),\n    original_project != "", id_old != ""\n  ) %>%',
  '  filter(!is.na(id), subj_id != "", site != 3) %>%',
  'DK68 demographic filter block'
)

eval(parse(text = script_txt), envir = globalenv())
