# ============================================================
# DGLM Analysis for DK-318 MIND degree and edges
#
# Default main analysis:
#   y ~ Diagnosis * AgeGroup * Sex + age_within_group
#   dformula = same as mean model
#
# Effects are estimated from the fitted model by equal-weight
# marginal contrasts over the model-predicted cells.
#
# Usage:
#   Rscript dglm_dk318_degree_edge.R [analysis_type]
#   analysis_type: both (default), degree, edge
#
# Optional environment variables for wrapper scripts:
#   DK318_DGLM_OUT_ROOT
#   DK318_DGLM_INCLUDE_AGE        true/false
#   DK318_DGLM_SEX_FILTER         empty, Male, or Female
# ============================================================

args <- commandArgs(trailingOnly = TRUE)
analysis_type <- if (length(args) >= 1) tolower(args[1]) else "both"

if (!(analysis_type %in% c("both", "degree", "edge"))) {
  stop("Invalid analysis_type: '", analysis_type, "'. Must be 'both', 'degree', or 'edge'.")
}

run_degree <- analysis_type %in% c("both", "degree")
run_edge   <- analysis_type %in% c("both", "edge")

cat("Analysis type:", analysis_type,
    "| run_degree:", run_degree,
    "| run_edge:", run_edge, "\n")

packages <- c("readxl", "dplyr", "purrr", "dglm", "stringr", "emmeans")
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
library(stringr)
library(emmeans)

# ----------------------------
# 1) Paths
# Input paths are unchanged.
# ----------------------------
demo_file       <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/all_data_cqt_mean_1.5.xlsx"
mind_combat_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_combat"
network_brainmap_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/brain_region_table_brainmaps"

default_out_root <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_DGLM"
out_root <- Sys.getenv("DK318_DGLM_OUT_ROOT", unset = default_out_root)

include_age <- tolower(Sys.getenv("DK318_DGLM_INCLUDE_AGE", unset = "true")) %in% c("true", "1", "yes", "y")
sex_filter <- Sys.getenv("DK318_DGLM_SEX_FILTER", unset = "")
sex_filter <- trimws(sex_filter)
if (!(sex_filter %in% c("", "Male", "Female"))) {
  stop("DK318_DGLM_SEX_FILTER must be empty, Male, or Female. Got: ", sex_filter)
}

if (sex_filter != "") {
  out_root <- file.path(out_root, paste0(tolower(sex_filter), "_sex"))
}

out_degree_dir  <- file.path(out_root, "degree")
out_edge_dir    <- file.path(out_root, "edge")
out_summary_dir <- file.path(out_root, "summary")

dir.create(out_root, showWarnings = FALSE, recursive = TRUE)
dir.create(out_degree_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_edge_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_summary_dir, showWarnings = FALSE, recursive = TRUE)

out_degree_csv <- file.path(out_degree_dir, "DGLM_DK318_degree_results.csv")
out_edge_csv   <- file.path(out_edge_dir,   "DGLM_DK318_edge_results.csv")

threshold <- 0.05

target_factors <- if (sex_filter == "") {
  c("Diagnosis", "AgeGroup", "Sex")
} else {
  c("Diagnosis", "AgeGroup")
}

effect_specs <- list(
  Diagnosis = c("Diagnosis"),
  AgeGroup = c("AgeGroup"),
  Sex = c("Sex"),
  Diagnosis_AgeGroup = c("Diagnosis", "AgeGroup"),
  Diagnosis_Sex = c("Diagnosis", "Sex"),
  AgeGroup_Sex = c("AgeGroup", "Sex"),
  Diagnosis_AgeGroup_Sex = c("Diagnosis", "AgeGroup", "Sex")
)
effect_specs <- effect_specs[vapply(effect_specs, function(x) all(x %in% target_factors), logical(1))]
effect_names <- names(effect_specs)

model_factor_rhs <- paste(target_factors, collapse = " * ")
model_rhs <- model_factor_rhs
if (include_age) model_rhs <- paste(model_rhs, "+ age_within_group")

mean_formula <- as.formula(paste("y ~", model_rhs))
disp_formula <- as.formula(paste("~", model_rhs))
emm_formula <- as.formula(paste("~", paste(target_factors, collapse = " * ")))

cat("Output root:", out_root, "\n")
cat("Sex filter:", ifelse(sex_filter == "", "none", sex_filter), "\n")
cat("Include age_within_group:", include_age, "\n")
cat("Mean formula:", deparse(mean_formula), "\n")
cat("Dispersion formula:", deparse(disp_formula), "\n")
cat("Effects:", paste(effect_names, collapse = ", "), "\n")

# ----------------------------
# 2) Load demo + recode factors
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

    has_file  = file.exists(mind_file) & file.exists(degree_file)
  ) %>%
  filter(
    !is.na(original_project), !is.na(id_old),
    original_project != "", id_old != ""
  ) %>%
  group_by(AgeGroup) %>%
  mutate(age_within_group = Age - mean(Age, na.rm = TRUE)) %>%
  ungroup()

cat("Subjects in demo:", nrow(df), "\n")
cat("Subjects with matrix+degree csv:", sum(df$has_file), "\n")

df2 <- df %>% filter(has_file)
if (sex_filter != "") {
  df2 <- df2 %>% filter(Sex == sex_filter)
}

if (nrow(df2) == 0) {
  stop("No valid subjects with existing DK318 combat matrix/degree csv files were found.")
}

cat("\nCell counts among subjects used:\n")
if (sex_filter == "") {
  print(with(df2, table(Diagnosis, AgeGroup, Sex)))
} else {
  print(with(df2, table(Diagnosis, AgeGroup)))
}

# ----------------------------
# 3) Read matrix/degree helpers
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

if (length(tmp_deg) != n_roi) {
  stop("Degree length does not match matrix ROI count for first DK318 subject.")
}
if (!all(names(tmp_deg) == roi_names)) {
  stop("ROI names in degree csv do not match ROI names in matrix csv for first DK318 subject.")
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
    stop("Matrix dimension mismatch for file: ", df2$mind_file[s],
         " | expected ", n_roi, "x", n_roi,
         ", got ", paste(dim(m), collapse = "x"))
  }

  if (!identical(colnames(m), roi_names) || !identical(rownames(m), roi_names)) {
    stop("ROI names/order mismatch in matrix file: ", df2$mind_file[s])
  }

  if (length(d) != n_roi || !identical(names(d), roi_names)) {
    stop("ROI names/order mismatch in degree file: ", df2$degree_file[s])
  }

  if (run_degree) degree_mat[s, ] <- d
  if (run_edge)   edge_mat[s, ] <- m[upper.tri(m)]
}

if (run_degree) cat("Degree matrix:", paste(dim(degree_mat), collapse = " x "), "\n")
if (run_edge)   cat("Edge matrix:", paste(dim(edge_mat), collapse = " x "), "\n")

# ============================================================
# 5) Degree network helpers
# ============================================================
load_degree_network_sets <- function(base_dir) {
  if (!dir.exists(base_dir)) {
    cat("Network ROI directory not found:", base_dir, "\n")
    return(list())
  }

  network_dirs <- list.dirs(base_dir, recursive = FALSE, full.names = TRUE)
  network_dirs <- network_dirs[basename(network_dirs) != ""]

  if (length(network_dirs) == 0) {
    cat("No network directories found in:", base_dir, "\n")
    return(list())
  }

  out <- list()
  for (net_dir in network_dirs) {
    net_name <- basename(net_dir)
    csv_path <- file.path(net_dir, paste0(net_name, "_regions_used.csv"))
    if (!file.exists(csv_path)) {
      cat("Skipping network without regions_used csv:", net_name, "\n")
      next
    }

    net_tab <- tryCatch(
      read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (is.null(net_tab) || !("feature" %in% colnames(net_tab))) {
      cat("Skipping network with invalid csv format:", net_name, "\n")
      next
    }

    features <- unique(as.character(net_tab$feature))
    features <- features[!is.na(features) & features != ""]
    if (length(features) == 0) {
      cat("Skipping empty network:", net_name, "\n")
      next
    }

    out[[net_name]] <- features
    cat("Loaded network:", net_name, "| n_roi =", length(features), "\n")
  }

  out
}

sanitize_network_name <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  tolower(x)
}

get_raw_p_cols <- function(res_df) {
  p_cols <- grep("^p_(mean|disp)_", colnames(res_df), value = TRUE)
  p_cols <- p_cols[!grepl("_FDR$|_networkFDR_", p_cols)]
  p_cols
}

add_global_fdr <- function(res_df) {
  p_cols <- get_raw_p_cols(res_df)
  for (p_col in p_cols) {
    res_df[[p_col]] <- as.numeric(res_df[[p_col]])
    res_df[[paste0(p_col, "_FDR")]] <- p.adjust(res_df[[p_col]], method = "fdr")
  }
  res_df
}

add_sign_columns <- function(res_df) {
  for (domain in c("mean", "disp")) {
    for (eff in effect_names) {
      est_col <- paste0("estimate_", domain, "_", eff)
      if (est_col %in% colnames(res_df)) {
        sign_col <- paste0("sign_", domain, "_", eff)
        x <- as.numeric(res_df[[est_col]])
        s <- sign(x)
        s[!is.finite(s) | s == 0] <- 1
        res_df[[sign_col]] <- s
      }
    }
  }
  res_df
}

add_network_fdr <- function(res_df, p_cols, network_sets, feature_col = "feature") {
  if (length(network_sets) == 0 || nrow(res_df) == 0) return(res_df)
  if (!(feature_col %in% colnames(res_df))) return(res_df)

  features <- as.character(res_df[[feature_col]])

  for (p_col in p_cols) {
    if (!(p_col %in% colnames(res_df))) next

    for (net_name in names(network_sets)) {
      members <- intersect(features, network_sets[[net_name]])
      if (length(members) == 0) next

      idx <- which(features %in% members & is.finite(res_df[[p_col]]))
      if (length(idx) == 0) next

      out_col <- paste0(p_col, "_networkFDR_", sanitize_network_name(net_name))
      res_df[[out_col]] <- NA_real_
      res_df[[out_col]][idx] <- p.adjust(res_df[[p_col]][idx], method = "fdr")
    }
  }

  res_df
}

export_network_degree_results <- function(out_degree_dir, degree_res, network_sets) {
  if (length(network_sets) == 0 || nrow(degree_res) == 0) return(invisible(NULL))

  network_root <- file.path(out_degree_dir, "network_FDR")
  dir.create(network_root, showWarnings = FALSE, recursive = TRUE)

  p_cols <- get_raw_p_cols(degree_res)

  for (net_name in names(network_sets)) {
    net_features <- unique(network_sets[[net_name]])
    net_features <- net_features[!is.na(net_features) & net_features != ""]
    if (length(net_features) == 0) next

    net_tag <- sanitize_network_name(net_name)
    net_dir <- file.path(network_root, net_name)
    dir.create(net_dir, showWarnings = FALSE, recursive = TRUE)

    idx <- degree_res$feature %in% net_features
    net_res <- degree_res[idx, , drop = FALSE]
    if (nrow(net_res) == 0) next

    keep_cols <- c("feature", "n", "fit_status")
    for (p_col in p_cols) {
      keep_cols <- c(
        keep_cols,
        sub("^p_", "estimate_", p_col),
        sub("^p_", "se_", p_col),
        sub("^p_", "stat_", p_col),
        p_col,
        paste0(p_col, "_FDR"),
        paste0(p_col, "_networkFDR_", net_tag)
      )
    }
    keep_cols <- unique(keep_cols[keep_cols %in% colnames(net_res)])
    net_res_out <- net_res[, keep_cols, drop = FALSE]

    write.csv(
      net_res_out,
      file.path(net_dir, "DGLM_DK318_degree_network.csv"),
      row.names = FALSE
    )

    manifest <- data.frame(
      network = net_name,
      network_tag = net_tag,
      n_roi_in_definition = length(net_features),
      n_roi_in_results = nrow(net_res)
    )
    write.csv(manifest, file.path(net_dir, "network_manifest.csv"), row.names = FALSE)

    cat("Exported network degree results:", net_name, "| n_roi in results =", nrow(net_res), "\n")
  }

  invisible(NULL)
}

network_sets <- load_degree_network_sets(network_brainmap_dir)

# ============================================================
# 6) Equal-weight marginal contrasts
# ============================================================
factor_levels <- list(
  Diagnosis = c("TD", "DD"),
  AgeGroup = c("Child", "Adult"),
  Sex = c("Female", "Male")
)

factor_positive <- list(
  Diagnosis = "DD",
  AgeGroup = "Adult",
  Sex = "Male"
)

factor_negative <- list(
  Diagnosis = "TD",
  AgeGroup = "Child",
  Sex = "Female"
)

cell_specs <- do.call(
  expand.grid,
  c(
    lapply(target_factors, function(f) factor_levels[[f]]),
    list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  )
)
names(cell_specs) <- target_factors

cell_label <- function(row) {
  paste(vapply(target_factors, function(f) as.character(row[[f]]), character(1)), collapse = "_")
}

safe_sd <- function(x) {
  if (sum(is.finite(x)) >= 2) sd(x, na.rm = TRUE) else NA_real_
}

raw_cell_descriptives <- function(y, d) {
  vals <- list()

  for (r in seq_len(nrow(cell_specs))) {
    row <- cell_specs[r, , drop = FALSE]
    lab <- cell_label(row)

    idx <- rep(TRUE, nrow(d))
    for (f in target_factors) {
      idx <- idx & as.character(d[[f]]) == as.character(row[[f]][1])
    }

    yy <- y[idx]
    vals[[paste0("raw_n_", lab)]] <- sum(is.finite(yy))
    vals[[paste0("raw_mean_", lab)]] <- if (sum(is.finite(yy)) >= 1) mean(yy, na.rm = TRUE) else NA_real_
    vals[[paste0("raw_sd_", lab)]] <- safe_sd(yy)
  }

  vals
}

has_all_target_cells <- function(d) {
  for (r in seq_len(nrow(cell_specs))) {
    row <- cell_specs[r, , drop = FALSE]
    idx <- rep(TRUE, nrow(d))
    for (f in target_factors) {
      idx <- idx & as.character(d[[f]]) == as.character(row[[f]][1])
    }
    if (sum(idx) == 0) return(FALSE)
  }
  TRUE
}

make_empty_domain_values <- function(domain) {
  vals <- list()
  for (eff in effect_names) {
    vals[[paste0("estimate_", domain, "_", eff)]] <- NA_real_
    vals[[paste0("se_",       domain, "_", eff)]] <- NA_real_
    vals[[paste0("stat_",     domain, "_", eff)]] <- NA_real_
    vals[[paste0("p_",        domain, "_", eff)]] <- NA_real_
  }
  vals
}

make_empty_effect_values <- function() {
  c(make_empty_domain_values("mean"), make_empty_domain_values("disp"))
}

extract_equal_contrasts <- function(model, domain) {
  vals <- make_empty_domain_values(domain)

  at_arg <- if (include_age) list(age_within_group = 0) else NULL

  emm <- tryCatch(
    emmeans::emmeans(
      model,
      specs = emm_formula,
      at = at_arg,
      weights = "equal"
    ),
    error = function(e) NULL
  )
  if (is.null(emm)) return(vals)

  grid <- tryCatch(as.data.frame(emm), error = function(e) NULL)
  if (is.null(grid)) return(vals)
  if (!all(target_factors %in% colnames(grid))) return(vals)

  n_cells <- nrow(grid)
  if (n_cells != 2 ^ length(target_factors)) return(vals)

  sign_by_factor <- list()
  for (f in target_factors) {
    sign_by_factor[[f]] <- ifelse(
      as.character(grid[[f]]) == factor_positive[[f]], 1,
      ifelse(as.character(grid[[f]]) == factor_negative[[f]], -1, NA_real_)
    )
    if (any(!is.finite(sign_by_factor[[f]]))) return(vals)
  }

  contrast_weights <- list()
  for (eff in effect_names) {
    eff_factors <- effect_specs[[eff]]
    w <- rep(1, n_cells)
    for (f in eff_factors) {
      w <- w * sign_by_factor[[f]]
    }
    w <- w / (2 ^ (length(target_factors) - length(eff_factors)))
    contrast_weights[[eff]] <- w
  }

  con <- tryCatch(
    emmeans::contrast(emm, method = contrast_weights, adjust = "none"),
    error = function(e) NULL
  )
  if (is.null(con)) return(vals)

  sm <- tryCatch(
    as.data.frame(summary(con, infer = c(FALSE, TRUE), adjust = "none")),
    error = function(e) NULL
  )
  if (is.null(sm)) return(vals)

  stat_col <- c("t.ratio", "z.ratio", "statistic", "T.ratio", "Z.ratio")
  stat_col <- stat_col[stat_col %in% colnames(sm)][1]
  if (is.na(stat_col)) stat_col <- NULL

  for (eff in effect_names) {
    idx <- which(as.character(sm$contrast) == eff)
    if (length(idx) == 0) next
    idx <- idx[1]

    vals[[paste0("estimate_", domain, "_", eff)]] <- as.numeric(sm$estimate[idx])
    vals[[paste0("se_",       domain, "_", eff)]] <- if ("SE" %in% colnames(sm)) as.numeric(sm$SE[idx]) else NA_real_
    vals[[paste0("stat_",     domain, "_", eff)]] <- if (!is.null(stat_col)) as.numeric(sm[[stat_col]][idx]) else NA_real_
    vals[[paste0("p_",        domain, "_", eff)]] <- if ("p.value" %in% colnames(sm)) as.numeric(sm$p.value[idx]) else NA_real_
  }

  vals
}

# ============================================================
# 7) One-feature model
# ============================================================
run_one_dglm <- function(y, df_sub, metadata, min_n = 10, var_eps = 1e-12) {
  ok <- is.finite(y) &
    !is.na(df_sub$Diagnosis) &
    !is.na(df_sub$AgeGroup) &
    !is.na(df_sub$Sex)

  if (include_age) ok <- ok & is.finite(df_sub$age_within_group)

  y2 <- y[ok]
  d2 <- df_sub[ok, , drop = FALSE]

  d2$Diagnosis <- factor(as.character(d2$Diagnosis), levels = c("TD", "DD"))
  d2$AgeGroup  <- factor(as.character(d2$AgeGroup),  levels = c("Child", "Adult"))
  d2$Sex       <- factor(as.character(d2$Sex),       levels = c("Female", "Male"))

  base <- c(metadata, list(n = length(y2), fit_status = "not_run"))
  raw_vals <- raw_cell_descriptives(y2, d2)
  empty_vals <- make_empty_effect_values()

  if (length(y2) < min_n) {
    base$fit_status <- "low_n"
    return(as.data.frame(c(base, raw_vals, empty_vals), stringsAsFactors = FALSE))
  }

  v <- var(y2)
  if (!is.finite(v) || v < var_eps) {
    base$fit_status <- "near_zero_variance"
    return(as.data.frame(c(base, raw_vals, empty_vals), stringsAsFactors = FALSE))
  }

  if (!has_all_target_cells(d2)) {
    base$fit_status <- "empty_target_cell"
    return(as.data.frame(c(base, raw_vals, empty_vals), stringsAsFactors = FALSE))
  }

  d2$y <- y2

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

  mean_vals <- extract_equal_contrasts(fit_mean, domain = "mean")
  disp_vals <- extract_equal_contrasts(fit$dispersion.fit, domain = "disp")

  as.data.frame(c(base, raw_vals, mean_vals, disp_vals), stringsAsFactors = FALSE)
}

run_dglm_per_feature <- function(y_mat, df_sub, feature_names,
                                 min_n = 10, var_eps = 1e-12,
                                 progress_every = 25) {
  res <- vector("list", ncol(y_mat))

  for (k in seq_len(ncol(y_mat))) {
    if (k %% progress_every == 0 || k == 1 || k == ncol(y_mat)) {
      cat("Feature", k, "/", ncol(y_mat), "\n")
    }
    metadata <- list(feature = feature_names[k])
    res[[k]] <- run_one_dglm(
      y = y_mat[, k],
      df_sub = df_sub,
      metadata = metadata,
      min_n = min_n,
      var_eps = var_eps
    )
  }

  out <- dplyr::bind_rows(res)
  out <- add_global_fdr(out)
  out <- add_sign_columns(out)
  out
}

run_dglm_edgewise_chunked <- function(edge_mat, df_sub, upper_idx, roi_names,
                                      chunk_size = 2000, min_n = 10, var_eps = 1e-12) {
  n_edge <- ncol(edge_mat)
  chunks <- split(seq_len(n_edge), ceiling(seq_len(n_edge) / chunk_size))
  out_list <- vector("list", length(chunks))

  for (ci in seq_along(chunks)) {
    idxs <- chunks[[ci]]
    cat("Edge chunk", ci, "/", length(chunks),
        "| edges", min(idxs), "-", max(idxs), "\n")

    chunk_res <- lapply(idxs, function(k) {
      region_i <- roi_names[upper_idx[k, 1]]
      region_j <- roi_names[upper_idx[k, 2]]
      metadata <- list(
        edge_id = k,
        i = upper_idx[k, 1],
        j = upper_idx[k, 2],
        region_i = region_i,
        region_j = region_j,
        feature = paste(region_i, region_j, sep = "__")
      )
      run_one_dglm(
        y = edge_mat[, k],
        df_sub = df_sub,
        metadata = metadata,
        min_n = min_n,
        var_eps = var_eps
      )
    })

    out_list[[ci]] <- dplyr::bind_rows(chunk_res)
  }

  out <- dplyr::bind_rows(out_list)
  out <- add_global_fdr(out)
  out <- add_sign_columns(out)
  out
}

# ============================================================
# 8) Summaries
# ============================================================
safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else min(x)
}

write_significant_summaries <- function(res_df, result_type, id_cols, out_summary_dir,
                                        threshold = 0.05) {
  p_cols <- get_raw_p_cols(res_df)

  summary_df <- dplyr::bind_rows(lapply(p_cols, function(p_col) {
    fdr_col <- paste0(p_col, "_FDR")
    data.frame(
      result_type = result_type,
      effect = p_col,
      n_tested = sum(is.finite(res_df[[p_col]])),
      n_sig_FDR = if (fdr_col %in% colnames(res_df)) {
        sum(is.finite(res_df[[fdr_col]]) & res_df[[fdr_col]] < threshold)
      } else {
        NA_integer_
      },
      min_p_raw = safe_min(res_df[[p_col]]),
      min_p_FDR = if (fdr_col %in% colnames(res_df)) safe_min(res_df[[fdr_col]]) else NA_real_
    )
  }))

  summary_csv <- file.path(
    out_summary_dir,
    paste0("DGLM_DK318_", result_type, "_FDR_summary.csv")
  )
  write.csv(summary_df, summary_csv, row.names = FALSE)

  sig_rows <- list()
  for (p_col in p_cols) {
    fdr_col <- paste0(p_col, "_FDR")
    if (!(fdr_col %in% colnames(res_df))) next

    idx <- which(is.finite(res_df[[fdr_col]]) & res_df[[fdr_col]] < threshold)
    if (length(idx) == 0) next

    est_col <- sub("^p_", "estimate_", p_col)
    se_col  <- sub("^p_", "se_", p_col)
    st_col  <- sub("^p_", "stat_", p_col)

    keep_cols <- c(id_cols, "n", est_col, se_col, st_col, p_col, fdr_col)
    keep_cols <- keep_cols[keep_cols %in% colnames(res_df)]

    tmp <- res_df[idx, keep_cols, drop = FALSE]
    tmp$effect <- p_col
    tmp$p_raw <- res_df[[p_col]][idx]
    tmp$p_FDR <- res_df[[fdr_col]][idx]

    sig_rows[[p_col]] <- tmp
  }

  sig_df <- if (length(sig_rows) > 0) dplyr::bind_rows(sig_rows) else data.frame()

  sig_csv <- file.path(
    out_summary_dir,
    paste0("DGLM_DK318_", result_type, "_significant_FDR",
           gsub("\\.", "p", as.character(threshold)), ".csv")
  )
  write.csv(sig_df, sig_csv, row.names = FALSE)

  cat("\n==============================\n")
  cat("Significant", result_type, "results, FDR <", threshold, "\n")
  cat("==============================\n")
  print(summary_df, row.names = FALSE)
  cat("Saved summary:", summary_csv, "\n")
  cat("Saved significant table:", sig_csv, "\n")
  invisible(list(summary = summary_df, significant = sig_df))
}

# ============================================================
# 9) Diagnostic on first ROI
# ============================================================
if (run_degree) {
  cat("\n=== Diagnostic: first ROI ===\n")
  diag_y <- degree_mat[, 1]
  diag_res <- run_one_dglm(
    y = diag_y,
    df_sub = df2,
    metadata = list(feature = colnames(degree_mat)[1]),
    min_n = 10,
    var_eps = 1e-12
  )
  diag_p_cols <- c(paste0("p_mean_", effect_names), paste0("p_disp_", effect_names))
  diag_keep <- c("feature", "n", "fit_status", diag_p_cols)
  diag_keep <- diag_keep[diag_keep %in% colnames(diag_res)]
  print(diag_res[, diag_keep, drop = FALSE], row.names = FALSE)
}

# ============================================================
# 10) Run degree
# ============================================================
if (run_degree) {
  cat("\n=== Running DGLM degree ===\n")

  degree_res <- run_dglm_per_feature(
    y_mat = degree_mat,
    df_sub = df2,
    feature_names = colnames(degree_mat),
    min_n = 10,
    var_eps = 1e-12,
    progress_every = 25
  )

  p_cols_degree <- get_raw_p_cols(degree_res)
  degree_res <- add_network_fdr(
    res_df = degree_res,
    p_cols = p_cols_degree,
    network_sets = network_sets,
    feature_col = "feature"
  )

  write.csv(degree_res, out_degree_csv, row.names = FALSE)
  cat("Saved degree results:", out_degree_csv, "\n")

  export_network_degree_results(
    out_degree_dir = out_degree_dir,
    degree_res = degree_res,
    network_sets = network_sets
  )

  write_significant_summaries(
    res_df = degree_res,
    result_type = "degree",
    id_cols = c("feature"),
    out_summary_dir = out_summary_dir,
    threshold = threshold
  )
}

# ============================================================
# 11) Run edge
# ============================================================
if (run_edge) {
  cat("\n=== Running DGLM edges ===\n")

  edge_res <- run_dglm_edgewise_chunked(
    edge_mat = edge_mat,
    df_sub = df2,
    upper_idx = upper_idx,
    roi_names = roi_names,
    chunk_size = 2000,
    min_n = 10,
    var_eps = 1e-12
  )

  write.csv(edge_res, out_edge_csv, row.names = FALSE)
  cat("Saved edge results:", out_edge_csv, "\n")

  write_significant_summaries(
    res_df = edge_res,
    result_type = "edge",
    id_cols = c("edge_id", "i", "j", "region_i", "region_j"),
    out_summary_dir = out_summary_dir,
    threshold = threshold
  )
}

cat("\nDONE.\n")
