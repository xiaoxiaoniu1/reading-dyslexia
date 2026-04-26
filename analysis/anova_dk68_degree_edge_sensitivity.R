script_path <- "/data/home/tqi/data1/share/after_freesurfer/CODE/analysis/anova_dk68_degree_edge.R"
script_txt <- paste(readLines(script_path, warn = FALSE), collapse = "\n")

script_txt <- gsub(
  "out_dir <- \"/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK68_ANOVA/\"",
  "out_dir <- \"/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK68_ANOVA_sensitivity/\"",
  script_txt,
  fixed = TRUE
)

eval(parse(text = script_txt), envir = globalenv())
