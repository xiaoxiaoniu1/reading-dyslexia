script_path <- "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/anova_dk318_degree_edge.R"
script_txt <- paste(readLines(script_path, warn = FALSE), collapse = "\n")

script_txt <- gsub(
  "out_dir <- \"/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_ANOVA\"",
  "out_dir <- \"/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_ANOVA_main\"",
  script_txt,
  fixed = TRUE
)
script_txt <- gsub(
  "  filter(!is.na(original_project), !is.na(id_old), original_project != \"\", id_old != \"\") %>%\n  group_by(AgeGroup) %>%\n  mutate(age_within_group = Age - mean(Age, na.rm = TRUE)) %>%\n  ungroup()",
  "  filter(!is.na(original_project), !is.na(id_old), original_project != \"\", id_old != \"\")",
  script_txt,
  fixed = TRUE
)
script_txt <- gsub(" + age_within_group", "", script_txt, fixed = TRUE)

eval(parse(text = script_txt), envir = globalenv())
