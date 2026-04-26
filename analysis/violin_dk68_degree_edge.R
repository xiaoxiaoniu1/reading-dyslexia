template_path <- "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/violin_dk318_degree_edge.R"
script_txt <- paste(readLines(template_path, warn = FALSE), collapse = "\n")

replace_fixed_once <- function(text, old, new, label) {
  if (!grepl(old, text, fixed = TRUE)) stop("Failed to locate ", label)
  sub(old, new, text, fixed = TRUE)
}

script_txt <- gsub("MIND_DK318/DGLM", "MIND_DK68_DGLM", script_txt, fixed = TRUE)
script_txt <- gsub("MIND_DK318_combat", "MIND_DK68_combat", script_txt, fixed = TRUE)
script_txt <- gsub("DGLM_DK318_", "DGLM_DK68_", script_txt, fixed = TRUE)
script_txt <- gsub("DK318", "DK68", script_txt, fixed = TRUE)
script_txt <- gsub("violin_dk318_degree_edge.R", "violin_dk68_degree_edge.R", script_txt, fixed = TRUE)
script_txt <- gsub("VIOLIN_Followup", "VIOLIN_Followup", script_txt, fixed = TRUE)
script_txt <- gsub("violin_followup", "violin_followup", script_txt, fixed = TRUE)

script_txt <- replace_fixed_once(
  script_txt,
  '  mutate(original_project = as.character(`original-project`),',
  '  mutate(subj_id = trimws(as.character(id)),',
  'DK68 VIOLIN demographic mutate line start'
)
script_txt <- replace_fixed_once(
  script_txt,
  '    id_old = as.character(id_old),\n    Age_Raw = as.numeric(age_month),\n    file_base = paste0(original_project,"_",id_old,"_MIND_DK68_combat"),',
  '    Age_Raw = as.numeric(age_month),\n    file_base = paste0(subj_id,"_combat"),',
  'DK68 VIOLIN file_base block'
)
script_txt <- replace_fixed_once(
  script_txt,
  '    !is.na(original_project), !is.na(id_old), original_project!="", id_old!="",\n    !is.na(Age_Raw), has_file\n  )',
  '    !is.na(id), subj_id!="", site != 3,\n    !is.na(Age_Raw), has_file\n  )',
  'DK68 VIOLIN demographic filter block'
)

eval(parse(text = script_txt), envir = globalenv())
