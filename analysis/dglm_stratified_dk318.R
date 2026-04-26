# ============================================================
# Stratified DGLM Analysis for DK-318 MIND degree and edges
#
# This is the R implementation of the stratified DK318 DGLM.
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
#   Rscript dglm_stratified_dk318.R [analysis_type]
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

# ----------------------------
# 0) Command-line parsing
# ----------------------------
args <- commandArgs(trailingOnly = TRUE)

analysis_type <- "both"
if (length(args) >= 1 && !startsWith(args[1], "--")) {
  analysis_type <- tolower(args[1])
}
if (!(analysis_type %in% c("both", "degree", "edge"))) {
  stop("Invalid analysis_type: '", analysis_type, "'. Must be 'both', 'degree', or 'edge'.")
}

get_arg <- function(key, default = NULL) {
  idx <- match(key, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[idx + 1]
}

get_flag <- function(key) {
  any(args == key)
}

parse_bool <- function(x, default = FALSE) {
  if (is.null(x) || is.na(x)) return(default)
  tolower(as.character(x)) %in% c("1", "true", "yes", "y")
}

run_degree <- analysis_type %in% c("both", "degree")
run_edge   <- analysis_type %in% c("both", "edge")

cat("Analysis type:", analysis_type,
    "| run_degree:", run_degree,
    "| run_edge:", run_edge, "\n")

# ----------------------------
# 0b) Packages
# ----------------------------
packages <- c("readxl", "dplyr", "purrr", "dglm", "emmeans", "stringr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

library(readxl)
library(dplyr)
library(purrr)
library(dglm)
library(emmeans)
library(stringr)

# ----------------------------
# 1) Paths
# ----------------------------
demo_file_default <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/all_data_cqt_mean_1.5.xlsx"
mind_combat_dir_default <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_combat"
out_dir_default <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_DGLM_stratified"
network_brainmap_dir_default <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/brain_region_table_brainmaps"

demo_file <- get_arg("--demo-file", demo_file_default)
mind_combat_dir <- get_arg("--mind-combat-dir", mind_combat_dir_default)
out_dir <- get_arg("--out-dir", out_dir_default)
network_brainmap_dir <- get_arg("--network-dir", network_brainmap_dir_default)
include_age_covariate <- !get_flag("--no-age-covariate")
run_plots <- parse_bool(get_arg("--run-plots", "false"), default = FALSE)

threshold <- 0.05
chunk_size <- as.integer(get_arg("--chunk-size", "2000"))

cat("demo_file:", demo_file, "\n")
cat("mind_combat_dir:", mind_combat_dir, "\n")
cat("out_dir:", out_dir, "\n")
cat("include_age_covariate:", include_age_covariate, "\n")
cat("run_plots:", run_plots, "\n")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# 2) Demographics
# ----------------------------
df <- read_excel(demo_file, sheet = "Sheet1")

df <- df %>%
  mutate(
    original_project = as.character(`original-project`),
    id_old = as.character(id_old),
    Age = as.numeric(age_month),
    subj_prefix = paste0(original_project, "_", id_old),
    file_base = paste0(subj_prefix, "_MIND_DK318_combat"),
    mind_file = file.path(mind_combat_dir, paste0(file_base, ".csv")),
    degree_file = file.path(mind_combat_dir, paste0(file_base, "_degree.csv")),

    Diagnosis = ifelse(group_d_or_c == 0, "TD", "DD"),
    AgeGroup  = ifelse(group_age == 1, "Adult", "Child"),
    Sex       = ifelse(sex == 1, "Male", "Female"),

    Diagnosis = factor(Diagnosis, levels = c("TD", "DD")),
    AgeGroup  = factor(AgeGroup,  levels = c("Child", "Adult")),
    Sex       = factor(Sex,       levels = c("Female", "Male")),

    has_file = file.exists(mind_file) & file.exists(degree_file)
  ) %>%
  filter(!is.na(original_project), !is.na(id_old), original_project != "", id_old != "")

cat("Subjects in demo:", nrow(df), "\n")
cat("Subjects with matrix+degree csv:", sum(df$has_file), "\n")

df2 <- df %>% filter(has_file)
if (nrow(df2) == 0) {
  stop("No valid subjects with existing DK318 combat matrix/degree csv files were found.")
}

cat("\nCell counts among subjects with files:\n")
print(with(df2, table(Diagnosis, AgeGroup, Sex)))

# ----------------------------
# 3) Read helpers
# ----------------------------
read_mind_csv <- function(fp) {
  mat <- as.matrix(read.csv(fp, row.names = 1, check.names = FALSE))
  storage.mode(mat) <- "numeric"
  mat
}

read_degree_csv <- function(fp) {
  d <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(all(c("ROI", "degree") %in% colnames(d)))
  vals <- as.numeric(d$degree)
  names(vals) <- as.character(d$ROI)
  vals
}

tmp_mat <- read_mind_csv(df2$mind_file[1])
tmp_deg <- read_degree_csv(df2$degree_file[1])

n_roi <- nrow(tmp_mat)
stopifnot(n_roi == ncol(tmp_mat))
cat("ROI number:", n_roi, "\n")

roi_names <- colnames(tmp_mat)
if (is.null(roi_names) || any(is.na(roi_names)) || length(roi_names) != n_roi) {
  roi_names <- rownames(tmp_mat)
}
if (is.null(roi_names) || any(is.na(roi_names)) || length(roi_names) != n_roi) {
  stop("Failed to obtain ROI names from DK318 combat matrix csv.")
}

if (length(tmp_deg) != n_roi || !all(names(tmp_deg) == roi_names)) {
  stop("Degree csv does not match matrix ROI names/order for the first subject.")
}

upper_idx <- which(upper.tri(matrix(1, n_roi, n_roi)), arr.ind = TRUE)
n_edge <- nrow(upper_idx)
cat("Number of edges:", n_edge, "\n")

# ----------------------------
# 4) Build matrices
# ----------------------------
n_subj <- nrow(df2)

if (run_degree) {
  degree_mat <- matrix(
    NA_real_, nrow = n_subj, ncol = n_roi,
    dimnames = list(df2$file_base, roi_names)
  )
}
if (run_edge) {
  edge_mat <- matrix(
    NA_real_, nrow = n_subj, ncol = n_edge,
    dimnames = list(df2$file_base, paste0("E_", seq_len(n_edge)))
  )
}

for (s in seq_len(n_subj)) {
  m <- read_mind_csv(df2$mind_file[s])
  d <- read_degree_csv(df2$degree_file[s])

  if (!all(dim(m) == c(n_roi, n_roi))) {
    stop("Matrix dimension mismatch for file: ", df2$mind_file[s])
  }
  if (!identical(colnames(m), roi_names) || !identical(rownames(m), roi_names)) {
    stop("ROI names/order mismatch in matrix file: ", df2$mind_file[s])
  }
  if (length(d) != n_roi || !identical(names(d), roi_names)) {
    stop("ROI names/order mismatch in degree file: ", df2$degree_file[s])
  }

  if (run_degree) degree_mat[s, ] <- d
  if (run_edge) edge_mat[s, ] <- m[upper.tri(m)]
}

if (run_degree) cat("Degree matrix:", paste(dim(degree_mat), collapse = " x "), "\n")
if (run_edge) cat("Edge matrix:", paste(dim(edge_mat), collapse = " x "), "\n")

# ============================================================
# 5) Network ROI sets for within-network FDR
# ============================================================
sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

load_degree_network_sets <- function(base_dir) {
  if (!dir.exists(base_dir)) {
    cat("Network ROI directory not found:", base_dir, "\n")
    return(list())
  }
  network_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  out <- list()
  for (net_dir in network_dirs) {
    net_name <- basename(net_dir)
    csv_path <- file.path(net_dir, paste0(net_name, "_regions_used.csv"))
    if (!file.exists(csv_path)) next
    net_tab <- tryCatch(read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(net_tab) || !("feature" %in% colnames(net_tab))) next
    features <- unique(as.character(net_tab$feature))
    features <- features[!is.na(features) & features != ""]
    if (length(features) > 0) {
      out[[net_name]] <- features
      cat("Loaded network:", net_name, "| n_roi =", length(features), "\n")
    }
  }
  out
}

make_subset_age_covariate <- function(df_sub) {
  df_sub %>%
    group_by(AgeGroup) %>%
    mutate(age_within_group = Age - mean(Age, na.rm = TRUE)) %>%
    ungroup()
}

network_sets <- load_degree_network_sets(network_brainmap_dir)

# ============================================================
# 6) DGLM + emmeans helpers
# ============================================================
make_na_effect_values <- function(effect_names, domain) {
  vals <- list()
  for (eff in effect_names) {
    vals[[paste0("estimate_", domain, "_", eff)]] <- NA_real_
    vals[[paste0("se_", domain, "_", eff)]] <- NA_real_
    vals[[paste0("stat_", domain, "_", eff)]] <- NA_real_
    vals[[paste0("p_", domain, "_", eff)]] <- NA_real_
  }
  vals
}

extract_emm_2x2 <- function(model, factor1, factor2, effect1, effect2, interaction_effect,
                            domain, at_list = list()) {
  effect_names <- c(effect1, effect2, interaction_effect)
  vals <- make_na_effect_values(effect_names, domain)

  emm_formula <- as.formula(paste("~", factor1, "*", factor2))
  emm <- tryCatch(
    emmeans::emmeans(model, emm_formula, at = at_list, weights = "equal"),
    error = function(e) NULL
  )
  if (is.null(emm)) return(vals)

  grid <- tryCatch(as.data.frame(emm), error = function(e) NULL)
  if (is.null(grid) || nrow(grid) != 4) return(vals)

  lv1 <- if (is.factor(grid[[factor1]])) levels(grid[[factor1]]) else unique(as.character(grid[[factor1]]))
  lv2 <- if (is.factor(grid[[factor2]])) levels(grid[[factor2]]) else unique(as.character(grid[[factor2]]))
  lv1 <- lv1[lv1 %in% as.character(grid[[factor1]])]
  lv2 <- lv2[lv2 %in% as.character(grid[[factor2]])]
  # Use original model factor levels. Positive directions are the second levels:
  # Diagnosis DD, AgeGroup Adult, Sex Male.
  if (length(lv1) != 2 || length(lv2) != 2) return(vals)

  s1 <- ifelse(as.character(grid[[factor1]]) == lv1[2], 1, -1)
  s2 <- ifelse(as.character(grid[[factor2]]) == lv2[2], 1, -1)

  contrast_weights <- list()
  contrast_weights[[effect1]] <- s1 / 2
  contrast_weights[[effect2]] <- s2 / 2
  contrast_weights[[interaction_effect]] <- s1 * s2

  con <- tryCatch(
    emmeans::contrast(emm, method = contrast_weights, adjust = "none"),
    error = function(e) NULL
  )
  if (is.null(con)) return(vals)

  sm <- tryCatch(
    as.data.frame(summary(con, infer = c(FALSE, TRUE), adjust = "none")),
    error = function(e) NULL
  )
  if (is.null(sm) || !("contrast" %in% colnames(sm))) return(vals)

  stat_col <- c("t.ratio", "z.ratio", "statistic", "T.ratio", "Z.ratio")
  stat_col <- stat_col[stat_col %in% colnames(sm)][1]
  if (is.na(stat_col)) stat_col <- NULL

  for (eff in effect_names) {
    idx <- which(as.character(sm$contrast) == eff)
    if (length(idx) == 0) next
    idx <- idx[1]
    vals[[paste0("estimate_", domain, "_", eff)]] <- as.numeric(sm$estimate[idx])
    vals[[paste0("se_", domain, "_", eff)]] <- if ("SE" %in% colnames(sm)) as.numeric(sm$SE[idx]) else NA_real_
    vals[[paste0("stat_", domain, "_", eff)]] <- if (!is.null(stat_col)) as.numeric(sm[[stat_col]][idx]) else NA_real_
    vals[[paste0("p_", domain, "_", eff)]] <- if ("p.value" %in% colnames(sm)) as.numeric(sm$p.value[idx]) else NA_real_
  }
  vals
}

get_p_cols <- function(res_df) {
  p_cols <- grep("^p_(mean|disp)_", colnames(res_df), value = TRUE)
  p_cols <- p_cols[!grepl("_FDR$|_within_network", p_cols)]
  p_cols
}

apply_global_fdr <- function(res_df) {
  for (p_col in get_p_cols(res_df)) {
    res_df[[p_col]] <- as.numeric(res_df[[p_col]])
    fdr_col <- paste0(p_col, "_FDR")
    res_df[[fdr_col]] <- p.adjust(res_df[[p_col]], method = "fdr")
    res_df[[paste0("sig_", p_col, "_FDR")]] <- is.finite(res_df[[fdr_col]]) & res_df[[fdr_col]] < threshold
    res_df[[paste0("sig_", p_col, "_uncorrected")]] <- is.finite(res_df[[p_col]]) & res_df[[p_col]] < threshold
  }
  res_df
}

add_sign_columns <- function(res_df) {
  est_cols <- grep("^estimate_(mean|disp)_", colnames(res_df), value = TRUE)
  for (est_col in est_cols) {
    sign_col <- sub("^estimate_", "sign_", est_col)
    x <- as.numeric(res_df[[est_col]])
    s <- sign(x)
    s[!is.finite(s) | s == 0] <- 1
    res_df[[sign_col]] <- s
  }
  res_df
}

safe_sd <- function(x) if (sum(is.finite(x)) >= 2) sd(x, na.rm = TRUE) else NA_real_

add_raw_descriptives_2x2 <- function(y, d, factor1, factor2) {
  vals <- list()
  lv1 <- levels(d[[factor1]])
  lv2 <- levels(d[[factor2]])

  for (a in lv1) {
    yy <- y[d[[factor1]] == a]
    vals[[paste0("mean_", a)]] <- if (sum(is.finite(yy)) > 0) mean(yy, na.rm = TRUE) else NA_real_
    vals[[paste0("std_", a)]] <- safe_sd(yy)
    vals[[paste0("n_", a)]] <- sum(is.finite(yy))
  }
  for (b in lv2) {
    yy <- y[d[[factor2]] == b]
    vals[[paste0("mean_", b)]] <- if (sum(is.finite(yy)) > 0) mean(yy, na.rm = TRUE) else NA_real_
    vals[[paste0("std_", b)]] <- safe_sd(yy)
    vals[[paste0("n_", b)]] <- sum(is.finite(yy))
  }
  for (a in lv1) {
    for (b in lv2) {
      lab <- paste(a, b, sep = "_")
      yy <- y[d[[factor1]] == a & d[[factor2]] == b]
      vals[[paste0("raw_n_", lab)]] <- sum(is.finite(yy))
      vals[[paste0("raw_mean_", lab)]] <- if (sum(is.finite(yy)) > 0) mean(yy, na.rm = TRUE) else NA_real_
      vals[[paste0("raw_sd_", lab)]] <- safe_sd(yy)
    }
  }

  if (length(lv1) == 2) {
    m1 <- vals[[paste0("mean_", lv1[1])]]
    m2 <- vals[[paste0("mean_", lv1[2])]]
    s1 <- vals[[paste0("std_", lv1[1])]]
    s2 <- vals[[paste0("std_", lv1[2])]]
    vals[["mean_direction"]] <- if (!is.finite(m1) || !is.finite(m2)) NA_character_ else if (m2 > m1) paste0(lv1[2], ">", lv1[1]) else if (m2 < m1) paste0(lv1[2], "<", lv1[1]) else "equal"
    vals[["disp_direction"]] <- if (!is.finite(s1) || !is.finite(s2)) NA_character_ else if (s2 > s1) paste0(lv1[2], "_more_var") else if (s2 < s1) paste0(lv1[2], "_less_var") else "equal"
  }

  if (length(lv2) == 2) {
    m1 <- vals[[paste0("mean_", lv2[1])]]
    m2 <- vals[[paste0("mean_", lv2[2])]]
    s1 <- vals[[paste0("std_", lv2[1])]]
    s2 <- vals[[paste0("std_", lv2[2])]]
    vals[["second_factor_mean_direction"]] <- if (!is.finite(m1) || !is.finite(m2)) NA_character_ else if (m2 > m1) paste0(lv2[2], ">", lv2[1]) else if (m2 < m1) paste0(lv2[2], "<", lv2[1]) else "equal"
    vals[["second_factor_disp_direction"]] <- if (!is.finite(s1) || !is.finite(s2)) NA_character_ else if (s2 > s1) paste0(lv2[2], "_more_var") else if (s2 < s1) paste0(lv2[2], "_less_var") else "equal"
  }

  vals
}

has_min_2x2_cells <- function(d, factor1, factor2, min_cell_n = 2) {
  tab <- table(d[[factor1]], d[[factor2]])
  all(dim(tab) == c(2, 2)) && all(tab >= min_cell_n)
}

run_one_stratified_dglm <- function(y, df_sub, metadata, factor1, factor2,
                                    effect1, effect2, interaction_effect,
                                    include_age_covariate = TRUE,
                                    min_n = 10, var_eps = 1e-12,
                                    min_cell_n = 2) {
  effect_names <- c(effect1, effect2, interaction_effect)
  empty_vals <- c(make_na_effect_values(effect_names, "mean"), make_na_effect_values(effect_names, "disp"))

  ok <- is.finite(y) & !is.na(df_sub[[factor1]]) & !is.na(df_sub[[factor2]])
  if (include_age_covariate) ok <- ok & is.finite(df_sub$age_within_group)

  y2 <- y[ok]
  d2 <- df_sub[ok, , drop = FALSE]
  d2[[factor1]] <- droplevels(d2[[factor1]])
  d2[[factor2]] <- droplevels(d2[[factor2]])

  base <- c(metadata, list(n = length(y2), fit_status = "not_run"))
  raw_vals <- add_raw_descriptives_2x2(y2, d2, factor1, factor2)

  if (length(y2) < min_n) {
    base$fit_status <- "low_n"
    return(as.data.frame(c(base, raw_vals, empty_vals), stringsAsFactors = FALSE))
  }
  v <- var(y2)
  if (!is.finite(v) || v < var_eps) {
    base$fit_status <- "near_zero_variance"
    return(as.data.frame(c(base, raw_vals, empty_vals), stringsAsFactors = FALSE))
  }
  if (nlevels(d2[[factor1]]) != 2 || nlevels(d2[[factor2]]) != 2 || !has_min_2x2_cells(d2, factor1, factor2, min_cell_n)) {
    base$fit_status <- "insufficient_2x2_cells"
    return(as.data.frame(c(base, raw_vals, empty_vals), stringsAsFactors = FALSE))
  }

  d2$y <- y2

  rhs <- paste(factor1, "*", factor2)
  if (include_age_covariate) rhs <- paste(rhs, "+ age_within_group")
  mean_formula <- as.formula(paste("y ~", rhs))
  disp_formula <- as.formula(paste("~", rhs))
  at_list <- if (include_age_covariate) list(age_within_group = 0) else list()

  fit_err <- NA_character_
  fit <- tryCatch(
    dglm::dglm(
      formula = mean_formula,
      dformula = disp_formula,
      family = gaussian(link = "identity"),
      data = d2
    ),
    error = function(e) {
      fit_err <<- conditionMessage(e)
      NULL
    }
  )

  if (is.null(fit)) {
    base$fit_status <- paste0("dglm_error: ", fit_err)
    return(as.data.frame(c(base, raw_vals, empty_vals), stringsAsFactors = FALSE))
  }

  base$fit_status <- "ok"

  fit_mean <- fit
  class(fit_mean) <- "lm"

  mean_vals <- extract_emm_2x2(fit_mean, factor1, factor2, effect1, effect2, interaction_effect, "mean", at_list)
  disp_vals <- extract_emm_2x2(fit$dispersion.fit, factor1, factor2, effect1, effect2, interaction_effect, "disp", at_list)

  as.data.frame(c(base, raw_vals, mean_vals, disp_vals), stringsAsFactors = FALSE)
}

collect_matrix_results <- function(y_mat, df_sub, feature_info_df, config, progress_every = 50) {
  res <- vector("list", ncol(y_mat))
  for (k in seq_len(ncol(y_mat))) {
    if (k == 1 || k == ncol(y_mat) || k %% progress_every == 0) {
      cat("  Feature", k, "/", ncol(y_mat), "\n")
    }
    metadata <- as.list(feature_info_df[k, , drop = FALSE])
    metadata$subset_type <- config$subset_type
    metadata$subset_level <- config$subset_level
    metadata$model_label <- config$model_label
    metadata$factor1 <- config$factor1
    metadata$factor2 <- config$factor2

    res[[k]] <- run_one_stratified_dglm(
      y = y_mat[, k],
      df_sub = df_sub,
      metadata = metadata,
      factor1 = config$factor1,
      factor2 = config$factor2,
      effect1 = config$effect1,
      effect2 = config$effect2,
      interaction_effect = config$interaction_effect,
      include_age_covariate = include_age_covariate
    )
  }
  out <- bind_rows(res)
  out <- apply_global_fdr(out)
  out <- add_sign_columns(out)
  out
}

save_result_csv <- function(res_df, out_csv) {
  dir.create(dirname(out_csv), showWarnings = FALSE, recursive = TRUE)
  write.csv(res_df, out_csv, row.names = FALSE)
  cat("  Saved:", out_csv, "\n")
}

write_result_csv <- function(df_out, out_csv, keep_empty = FALSE, template_df = NULL) {
  if (is.null(df_out) || nrow(df_out) == 0) {
    if (keep_empty) {
      empty_df <- if (!is.null(template_df)) template_df[0, , drop = FALSE] else data.frame()
      write.csv(empty_df, out_csv, row.names = FALSE)
    } else if (file.exists(out_csv)) {
      file.remove(out_csv)
    }
    return(invisible(NULL))
  }
  write.csv(df_out, out_csv, row.names = FALSE)
}

save_significant_tables <- function(res_df, out_csv) {
  p_cols <- get_p_cols(res_df)
  out_base <- tools::file_path_sans_ext(basename(out_csv))
  out_parent <- dirname(out_csv)
  for (p_col in p_cols) {
    fdr_col <- paste0(p_col, "_FDR")
    sig_fdr <- res_df[is.finite(res_df[[fdr_col]]) & res_df[[fdr_col]] < threshold, , drop = FALSE]
    sig_raw <- res_df[is.finite(res_df[[p_col]]) & res_df[[p_col]] < threshold, , drop = FALSE]
    write_result_csv(sig_fdr, file.path(out_parent, paste0("Significant_", out_base, "_", p_col, ".csv")), keep_empty = TRUE, template_df = res_df)
    write_result_csv(sig_raw, file.path(out_parent, paste0("Uncorrected_Significant_", out_base, "_", p_col, ".csv")), keep_empty = FALSE, template_df = res_df)
  }
}

save_summary <- function(res_df, out_csv) {
  p_cols <- get_p_cols(res_df)
  rows <- lapply(p_cols, function(p_col) {
    fdr_col <- paste0(p_col, "_FDR")
    pvals <- as.numeric(res_df[[p_col]])
    fdrs <- as.numeric(res_df[[fdr_col]])
    data.frame(
      effect = p_col,
      n_valid = sum(is.finite(pvals)),
      min_p = if (any(is.finite(pvals))) min(pvals, na.rm = TRUE) else NA_real_,
      p01 = if (any(is.finite(pvals))) as.numeric(quantile(pvals[is.finite(pvals)], 0.01)) else NA_real_,
      p05 = if (any(is.finite(pvals))) as.numeric(quantile(pvals[is.finite(pvals)], 0.05)) else NA_real_,
      n_p_lt_0_05 = sum(is.finite(pvals) & pvals < threshold),
      min_fdr = if (any(is.finite(fdrs))) min(fdrs, na.rm = TRUE) else NA_real_,
      fdr01 = if (any(is.finite(fdrs))) as.numeric(quantile(fdrs[is.finite(fdrs)], 0.01)) else NA_real_,
      fdr05 = if (any(is.finite(fdrs))) as.numeric(quantile(fdrs[is.finite(fdrs)], 0.05)) else NA_real_,
      n_fdr_lt_0_05 = sum(is.finite(fdrs) & fdrs < threshold)
    )
  })
  summary_df <- bind_rows(rows)
  summary_csv <- file.path(dirname(out_csv), paste0("Diagnostic_Summary_", tools::file_path_sans_ext(basename(out_csv)), ".csv"))
  write.csv(summary_df, summary_csv, row.names = FALSE)
  cat("  Saved summary:", summary_csv, "\n")
  print(summary_df[, c("effect", "n_valid", "n_p_lt_0_05", "n_fdr_lt_0_05")], row.names = FALSE)
}

save_network_restricted_results <- function(res_df, network_sets, out_csv, feature_col = "feature") {
  if (length(network_sets) == 0 || !(feature_col %in% colnames(res_df))) return(invisible(NULL))

  interaction_p_cols <- get_p_cols(res_df)
  interaction_p_cols <- interaction_p_cols[grepl("Interaction", interaction_p_cols)]
  if (length(interaction_p_cols) == 0) return(invisible(NULL))

  out_base <- tools::file_path_sans_ext(basename(out_csv))
  network_out_dir <- file.path(out_dir, "network_restricted_fdr", out_base)
  dir.create(network_out_dir, showWarnings = FALSE, recursive = TRUE)

  features <- as.character(res_df[[feature_col]])
  for (network_name in names(network_sets)) {
    idx <- features %in% network_sets[[network_name]]
    if (!any(idx)) next

    net_df <- res_df[idx, , drop = FALSE]
    net_df <- cbind(network_name = network_name, net_df)

    for (p_col in interaction_p_cols) {
      fdr_col <- paste0(p_col, "_FDR_within_network")
      sig_col <- paste0("sig_", p_col, "_FDR_within_network")
      net_df[[fdr_col]] <- NA_real_
      valid <- is.finite(net_df[[p_col]])
      if (sum(valid) > 0) {
        net_df[[fdr_col]][valid] <- p.adjust(net_df[[p_col]][valid], method = "fdr")
      }
      net_df[[sig_col]] <- is.finite(net_df[[fdr_col]]) & net_df[[fdr_col]] < threshold
    }

    network_slug <- sanitize_name(network_name)
    network_csv <- file.path(network_out_dir, paste0(out_base, "_", network_slug, "_restricted_fdr.csv"))
    write.csv(net_df, network_csv, row.names = FALSE)
    cat("  Network", network_name, ": matched", nrow(net_df), "ROIs ->", network_csv, "\n")

    for (p_col in interaction_p_cols) {
      sig_col <- paste0("sig_", p_col, "_FDR_within_network")
      sig_df <- net_df[net_df[[sig_col]], , drop = FALSE]
      sig_csv <- file.path(network_out_dir, paste0("Significant_", out_base, "_", network_slug, "_", p_col, "_within_network_FDR.csv"))
      write_result_csv(sig_df, sig_csv, keep_empty = FALSE, template_df = net_df)
    }
  }
  invisible(NULL)
}

# ============================================================
# 7) Stratification configs
# ============================================================
make_config <- function(name, subset_type, subset_level, subset_var, subset_value,
                        factor1, factor2, effect1, effect2, interaction_effect,
                        output_stub, model_label) {
  list(
    name = name,
    subset_type = subset_type,
    subset_level = subset_level,
    subset_var = subset_var,
    subset_value = subset_value,
    factor1 = factor1,
    factor2 = factor2,
    effect1 = effect1,
    effect2 = effect2,
    interaction_effect = interaction_effect,
    output_stub = output_stub,
    model_label = model_label
  )
}

configs <- list(
  make_config(
    name = "Child Diagnosis x Sex",
    subset_type = "AgeGroup", subset_level = "Child", subset_var = "AgeGroup", subset_value = "Child",
    factor1 = "Diagnosis", factor2 = "Sex",
    effect1 = "Diagnosis", effect2 = "Sex", interaction_effect = "Diagnosis_Sex_Interaction",
    output_stub = "Child", model_label = "Within Child: Diagnosis * Sex"
  ),
  make_config(
    name = "Adult Diagnosis x Sex",
    subset_type = "AgeGroup", subset_level = "Adult", subset_var = "AgeGroup", subset_value = "Adult",
    factor1 = "Diagnosis", factor2 = "Sex",
    effect1 = "Diagnosis", effect2 = "Sex", interaction_effect = "Diagnosis_Sex_Interaction",
    output_stub = "Adult", model_label = "Within Adult: Diagnosis * Sex"
  ),
  make_config(
    name = "DD AgeGroup x Sex",
    subset_type = "Diagnosis", subset_level = "DD", subset_var = "Diagnosis", subset_value = "DD",
    factor1 = "AgeGroup", factor2 = "Sex",
    effect1 = "AgeGroup", effect2 = "Sex", interaction_effect = "AgeGroup_Sex_Interaction",
    output_stub = "DD_Adult_vs_Child", model_label = "Within DD: AgeGroup * Sex"
  ),
  make_config(
    name = "TD AgeGroup x Sex",
    subset_type = "Diagnosis", subset_level = "TD", subset_var = "Diagnosis", subset_value = "TD",
    factor1 = "AgeGroup", factor2 = "Sex",
    effect1 = "AgeGroup", effect2 = "Sex", interaction_effect = "AgeGroup_Sex_Interaction",
    output_stub = "TD_Adult_vs_Child", model_label = "Within TD: AgeGroup * Sex"
  ),
  make_config(
    name = "Male Diagnosis x AgeGroup",
    subset_type = "Sex", subset_level = "Male", subset_var = "Sex", subset_value = "Male",
    factor1 = "Diagnosis", factor2 = "AgeGroup",
    effect1 = "Diagnosis", effect2 = "AgeGroup", interaction_effect = "Diagnosis_AgeGroup_Interaction",
    output_stub = "Male_Diagnosis_AgeGroup", model_label = "Within Male: Diagnosis * AgeGroup"
  ),
  make_config(
    name = "Female Diagnosis x AgeGroup",
    subset_type = "Sex", subset_level = "Female", subset_var = "Sex", subset_value = "Female",
    factor1 = "Diagnosis", factor2 = "AgeGroup",
    effect1 = "Diagnosis", effect2 = "AgeGroup", interaction_effect = "Diagnosis_AgeGroup_Interaction",
    output_stub = "Female_Diagnosis_AgeGroup", model_label = "Within Female: Diagnosis * AgeGroup"
  )
)

# ============================================================
# 8) Degree analysis
# ============================================================
if (run_degree) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat("DEGREE STRATIFIED DGLM ANALYSIS\n")
  cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")

  degree_feature_info <- data.frame(feature = roi_names, stringsAsFactors = FALSE)

  for (config in configs) {
    cat("\n---", config$name, "---\n")
    mask <- as.character(df2[[config$subset_var]]) == config$subset_value
    df_g <- df2[mask, , drop = FALSE]
    y_g <- degree_mat[mask, , drop = FALSE]
    if (nrow(df_g) == 0) {
      cat("  Empty subset, skipped.\n")
      next
    }

    df_g <- make_subset_age_covariate(df_g)
    df_g[[config$factor1]] <- droplevels(df_g[[config$factor1]])
    df_g[[config$factor2]] <- droplevels(df_g[[config$factor2]])
    cat("  n =", nrow(df_g), "\n")
    print(table(df_g[[config$factor1]], df_g[[config$factor2]]))

    res <- collect_matrix_results(y_g, df_g, degree_feature_info, config, progress_every = 50)
    out_csv <- file.path(out_dir, paste0("DGLM_DK318_", config$output_stub, "_degree_results.csv"))
    save_result_csv(res, out_csv)
    save_summary(res, out_csv)
    save_significant_tables(res, out_csv)
    save_network_restricted_results(res, network_sets, out_csv, feature_col = "feature")
  }
} else {
  cat("Skipping degree analysis.\n")
}

# ============================================================
# 9) Edge analysis
# ============================================================
if (run_edge) {
  cat("\n", paste(rep("=", 60), collapse = ""), "\n", sep = "")
  cat("EDGE STRATIFIED DGLM ANALYSIS\n")
  cat(paste(rep("=", 60), collapse = ""), "\n", sep = "")

  edge_feature_info <- data.frame(
    edge_id = seq_len(n_edge),
    i = upper_idx[, 1],
    j = upper_idx[, 2],
    region_i = roi_names[upper_idx[, 1]],
    region_j = roi_names[upper_idx[, 2]],
    stringsAsFactors = FALSE
  )

  for (config in configs) {
    cat("\n---", config$name, "---\n")
    mask <- as.character(df2[[config$subset_var]]) == config$subset_value
    df_g <- df2[mask, , drop = FALSE]
    y_g <- edge_mat[mask, , drop = FALSE]
    if (nrow(df_g) == 0) {
      cat("  Empty subset, skipped.\n")
      next
    }

    df_g <- make_subset_age_covariate(df_g)
    df_g[[config$factor1]] <- droplevels(df_g[[config$factor1]])
    df_g[[config$factor2]] <- droplevels(df_g[[config$factor2]])
    cat("  n =", nrow(df_g), "\n")
    print(table(df_g[[config$factor1]], df_g[[config$factor2]]))

    chunks <- split(seq_len(n_edge), ceiling(seq_len(n_edge) / chunk_size))
    chunk_res <- vector("list", length(chunks))
    for (ci in seq_along(chunks)) {
      idxs <- chunks[[ci]]
      cat("  Edge chunk", ci, "/", length(chunks), "| edges", min(idxs), "-", max(idxs), "\n")
      chunk_info <- edge_feature_info[idxs, , drop = FALSE]
      chunk_res[[ci]] <- collect_matrix_results(y_g[, idxs, drop = FALSE], df_g, chunk_info, config, progress_every = 1000000)
    }

    res <- bind_rows(chunk_res)
    res <- apply_global_fdr(res)
    res <- add_sign_columns(res)
    out_csv <- file.path(out_dir, paste0("DGLM_DK318_", config$output_stub, "_edge_results.csv"))
    save_result_csv(res, out_csv)
    save_summary(res, out_csv)
    save_significant_tables(res, out_csv)
  }
} else {
  cat("Skipping edge analysis.\n")
}

# ============================================================
# 10) Optional plotting
# ============================================================
cat("\nALL DONE.\n")

if (run_degree && run_plots) {
  brainmap_script <- "/data/home/tqi/data1/share/after_freesurfer/CODE/degree/plot_significant_DK318_stratified_brainmaps.py"
  violin_script <- "/data/home/tqi/data1/share/after_freesurfer/CODE/degree/violin_dk318_stratified_degree.R"

  if (file.exists(brainmap_script)) {
    cat("\nRunning stratified degree brainmap plotting...\n")
    cmd <- c("python3", brainmap_script, "--result-dir", out_dir, "--demo-file", demo_file, "--mind-dir", mind_combat_dir)
    print(cmd)
    try(system2(cmd[1], cmd[-1], wait = TRUE))
  }

  if (file.exists(violin_script)) {
    cat("\nRunning stratified interaction violin plotting...\n")
    cmd <- c("Rscript", violin_script, "--result-dir", out_dir, "--demo-file", demo_file, "--mind-combat-dir", mind_combat_dir)
    print(cmd)
    try(system2(cmd[1], cmd[-1], wait = TRUE))
  }
}
