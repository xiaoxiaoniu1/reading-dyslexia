# DGLM within clusters from DK318 degree_Tmap_cluster_K*.csv
# Usage:
#   Rscript 05_dglm_by_cluster_network.R
#   Rscript 05_dglm_by_cluster_network.R degree
#   Rscript 05_dglm_by_cluster_network.R /path/to/degree_Tmap_cluster_K3.csv degree

args <- commandArgs(trailingOnly = TRUE)
default_cluster_csv <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Clustering_3/clustering/degree_Tmap_cluster_K3.csv"
valid_analysis_types <- c("both", "degree", "edge")

if (length(args) == 0) {
  cluster_csv <- default_cluster_csv
  analysis_type <- "both"
} else if (tolower(args[1]) %in% valid_analysis_types) {
  cluster_csv <- default_cluster_csv
  analysis_type <- tolower(args[1])
} else {
  cluster_csv <- args[1]
  analysis_type <- if (length(args) >= 2) tolower(args[2]) else "both"
}

if (!(analysis_type %in% valid_analysis_types)) stop("analysis_type must be both, degree, or edge")
if (!file.exists(cluster_csv)) stop("Cannot find cluster_csv: ", cluster_csv)

pkgs <- c("readxl", "dplyr", "dglm", "stringr", "emmeans")
for (p in pkgs) if (!require(p, character.only = TRUE)) install.packages(p, repos = "https://cran.r-project.org")
invisible(lapply(pkgs, library, character.only = TRUE))

base_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5"
demo_file <- file.path(base_dir, "all_data_cqt_mean_1.5.xlsx")
mind_dir <- file.path(base_dir, "MIND_DK318_combat")
out_root <- Sys.getenv("DK318_CLUSTER_DGLM_OUT_ROOT", unset = file.path(base_dir, "MIND_DK318_DGLM_by_cluster", tools::file_path_sans_ext(basename(cluster_csv))))
include_age <- tolower(Sys.getenv("DK318_DGLM_INCLUDE_AGE", unset = "false")) %in% c("true", "1", "yes", "y")
include_iq <- tolower(Sys.getenv("DK318_DGLM_INCLUDE_IQ", unset = "false")) %in% c("true", "1", "yes", "y")
sex_filter <- trimws(Sys.getenv("DK318_DGLM_SEX_FILTER", unset = ""))
if (!(sex_filter %in% c("", "Male", "Female"))) stop("DK318_DGLM_SEX_FILTER must be empty, Male, or Female")
if (sex_filter != "") out_root <- file.path(out_root, paste0(tolower(sex_filter), "_sex"))
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

cluster_tab <- read.csv(cluster_csv, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("feature", "cluster") %in% names(cluster_tab))) stop("cluster csv must contain feature and cluster columns")
cluster_tab <- cluster_tab %>% filter(!is.na(feature), feature != "", !is.na(cluster)) %>% mutate(feature = as.character(feature), cluster = as.character(cluster))
cluster_sets <- split(cluster_tab$feature, cluster_tab$cluster)
cluster_sets <- lapply(cluster_sets, unique)

cat("Cluster csv:", cluster_csv, "\n")
cat("Output root:", out_root, "\n")
cat("Clusters:", paste(names(cluster_sets), lengths(cluster_sets), sep = "=", collapse = "; "), "\n")

factors_used <- if (sex_filter == "") c("Diagnosis", "AgeGroup", "Sex") else c("Diagnosis", "AgeGroup")
effect_specs <- list(Diagnosis = "Diagnosis", AgeGroup = "AgeGroup", Sex = "Sex", Diagnosis_AgeGroup = c("Diagnosis", "AgeGroup"), Diagnosis_Sex = c("Diagnosis", "Sex"), AgeGroup_Sex = c("AgeGroup", "Sex"), Diagnosis_AgeGroup_Sex = c("Diagnosis", "AgeGroup", "Sex"))
effect_specs <- effect_specs[vapply(effect_specs, function(x) all(x %in% factors_used), logical(1))]
effect_names <- names(effect_specs)
rhs <- paste(factors_used, collapse = " * ")
if (include_age) rhs <- paste(rhs, "+ age_within_group")
if (include_iq) rhs <- paste(rhs, "+ IQ")
mean_formula <- as.formula(paste("y ~", rhs))
disp_formula <- as.formula(paste("~", rhs))
emm_formula <- as.formula(paste("~", paste(factors_used, collapse = " * ")))
cat("Mean formula:", deparse(mean_formula), "\n")

read_mind <- function(fp) { m <- as.matrix(read.csv(fp, row.names = 1, check.names = FALSE)); storage.mode(m) <- "numeric"; m }
read_degree <- function(fp) { d <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE); vals <- as.numeric(d$degree); names(vals) <- as.character(d$ROI); vals }

subj <- readxl::read_excel(demo_file, sheet = "Sheet1")
iq_col <- names(subj)[stringr::str_detect(names(subj), stringr::regex("iq", ignore_case = TRUE))][1]
subj <- subj %>% mutate(
  original_project = as.character(`original-project`), id_old = as.character(id_old), Age = as.numeric(age_month),
  IQ = if (is.na(iq_col)) NA_real_ else as.numeric(.data[[iq_col]]), subj_prefix = paste0(original_project, "_", id_old),
  file_base = paste0(subj_prefix, "_MIND_DK318_combat"), mind_file = file.path(mind_dir, paste0(file_base, ".csv")), degree_file = file.path(mind_dir, paste0(file_base, "_degree.csv")),
  Diagnosis = factor(ifelse(group_d_or_c == 0, "TD", "DD"), levels = c("TD", "DD")),
  AgeGroup = factor(ifelse(group_age == 1, "Adult", "Child"), levels = c("Child", "Adult")),
  Sex = factor(ifelse(sex == 1, "Male", "Female"), levels = c("Female", "Male")),
  has_file = file.exists(mind_file) & file.exists(degree_file)
) %>% filter(!is.na(original_project), !is.na(id_old), original_project != "", id_old != "") %>% group_by(AgeGroup) %>% mutate(age_within_group = Age - mean(Age, na.rm = TRUE)) %>% ungroup() %>% filter(has_file)
if (sex_filter != "") subj <- subj %>% filter(Sex == sex_filter)
if (include_iq) subj <- subj %>% filter(is.finite(IQ))
if (nrow(subj) == 0) stop("No valid subjects")
print(if (sex_filter == "") with(subj, table(Diagnosis, AgeGroup, Sex)) else with(subj, table(Diagnosis, AgeGroup)))

m0 <- read_mind(subj$mind_file[1]); d0 <- read_degree(subj$degree_file[1])
roi_names <- colnames(m0); if (is.null(roi_names)) roi_names <- rownames(m0)
if (!identical(names(d0), roi_names)) stop("ROI names differ between matrix and degree csv")

empty_domain <- function(domain) setNames(as.list(rep(NA_real_, 4 * length(effect_names))), as.vector(sapply(effect_names, function(e) paste0(c("estimate_", "se_", "stat_", "p_"), domain, "_", e))))
empty_effects <- function() c(empty_domain("mean"), empty_domain("disp"))
cell_grid <- do.call(expand.grid, c(lapply(factors_used, function(f) list(Diagnosis = c("TD", "DD"), AgeGroup = c("Child", "Adult"), Sex = c("Female", "Male"))[[f]]), list(KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)))
names(cell_grid) <- factors_used
has_all_cells <- function(d) all(apply(cell_grid, 1, function(r) sum(Reduce(`&`, Map(function(f) as.character(d[[f]]) == r[[f]], factors_used))) > 0))
raw_cells <- function(y, d) {
  out <- list()
  for (i in seq_len(nrow(cell_grid))) {
    r <- cell_grid[i, , drop = FALSE]; lab <- paste(as.character(r[1, factors_used]), collapse = "_")
    idx <- Reduce(`&`, Map(function(f) as.character(d[[f]]) == as.character(r[[f]][1]), factors_used)); yy <- y[idx]
    out[[paste0("raw_n_", lab)]] <- sum(is.finite(yy)); out[[paste0("raw_mean_", lab)]] <- if (sum(is.finite(yy)) > 0) mean(yy, na.rm = TRUE) else NA_real_; out[[paste0("raw_sd_", lab)]] <- if (sum(is.finite(yy)) > 1) sd(yy, na.rm = TRUE) else NA_real_
  }
  out
}

extract_contrasts <- function(model, domain) {
  vals <- empty_domain(domain); at_arg <- if (include_age) list(age_within_group = 0) else NULL
  emm <- tryCatch(emmeans::emmeans(model, specs = emm_formula, at = at_arg, weights = "equal"), error = function(e) NULL)
  if (is.null(emm)) return(vals)
  grid <- as.data.frame(emm); pos <- list(Diagnosis = "DD", AgeGroup = "Adult", Sex = "Male"); neg <- list(Diagnosis = "TD", AgeGroup = "Child", Sex = "Female")
  signs <- lapply(factors_used, function(f) ifelse(as.character(grid[[f]]) == pos[[f]], 1, ifelse(as.character(grid[[f]]) == neg[[f]], -1, NA_real_))); names(signs) <- factors_used
  if (any(!is.finite(unlist(signs)))) return(vals)
  weights <- lapply(effect_specs, function(eff) Reduce(`*`, signs[eff]) / (2 ^ (length(factors_used) - length(eff))))
  con <- tryCatch(emmeans::contrast(emm, method = weights, adjust = "none"), error = function(e) NULL); if (is.null(con)) return(vals)
  sm <- as.data.frame(summary(con, infer = c(FALSE, TRUE), adjust = "none")); stat_col <- intersect(c("t.ratio", "z.ratio", "statistic"), names(sm))[1]
  for (eff in effect_names) { j <- match(eff, as.character(sm$contrast)); if (!is.na(j)) { vals[[paste0("estimate_", domain, "_", eff)]] <- sm$estimate[j]; vals[[paste0("se_", domain, "_", eff)]] <- if ("SE" %in% names(sm)) sm$SE[j] else NA_real_; vals[[paste0("stat_", domain, "_", eff)]] <- if (!is.na(stat_col)) sm[[stat_col]][j] else NA_real_; vals[[paste0("p_", domain, "_", eff)]] <- if ("p.value" %in% names(sm)) sm$p.value[j] else NA_real_ } }
  vals
}

run_one <- function(y, meta, min_n = 10) {
  ok <- is.finite(y) & !is.na(subj$Diagnosis) & !is.na(subj$AgeGroup) & !is.na(subj$Sex); if (include_age) ok <- ok & is.finite(subj$age_within_group); if (include_iq) ok <- ok & is.finite(subj$IQ)
  y2 <- y[ok]; d2 <- subj[ok, , drop = FALSE]; base <- c(meta, list(n = length(y2), fit_status = "not_run")); raw <- raw_cells(y2, d2); emp <- empty_effects()
  if (length(y2) < min_n) { base$fit_status <- "low_n"; return(as.data.frame(c(base, raw, emp), stringsAsFactors = FALSE)) }
  if (!is.finite(var(y2)) || var(y2) < 1e-12) { base$fit_status <- "near_zero_variance"; return(as.data.frame(c(base, raw, emp), stringsAsFactors = FALSE)) }
  if (!has_all_cells(d2)) { base$fit_status <- "empty_target_cell"; return(as.data.frame(c(base, raw, emp), stringsAsFactors = FALSE)) }
  d2$y <- y2; err <- NA_character_; fit <- tryCatch(dglm::dglm(mean_formula, dformula = disp_formula, family = gaussian(link = "identity"), data = d2), error = function(e) { err <<- conditionMessage(e); NULL })
  if (is.null(fit)) { base$fit_status <- paste0("dglm_error: ", err); return(as.data.frame(c(base, raw, emp), stringsAsFactors = FALSE)) }
  mf <- fit; class(mf) <- "lm"; as.data.frame(c(meta, list(n = length(y2), fit_status = "ok"), raw, extract_contrasts(mf, "mean"), extract_contrasts(fit$dispersion.fit, "disp")), stringsAsFactors = FALSE)
}

add_fdr <- function(x) { for (pc in grep("^p_(mean|disp)_", names(x), value = TRUE)) x[[paste0(pc, "_FDR_cluster")]] <- p.adjust(as.numeric(x[[pc]]), "fdr"); x }

for (cl in names(cluster_sets)) {
  features <- intersect(cluster_sets[[cl]], roi_names); tag <- paste0("cluster_", gsub("[^A-Za-z0-9]+", "_", cl)); cl_dir <- file.path(out_root, tag); dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(data.frame(cluster = cl, feature = features), file.path(cl_dir, "cluster_regions_used.csv"), row.names = FALSE)
  cat("\n===", tag, "n_roi =", length(features), "===\n")
  if (length(features) < 1) next
  idx <- match(features, roi_names)
  if (analysis_type %in% c("both", "degree")) {
    res <- lapply(seq_along(idx), function(k) { if (k %% 25 == 0 || k == 1 || k == length(idx)) cat("degree", k, "/", length(idx), "\n"); y <- vapply(subj$degree_file, function(fp) read_degree(fp)[idx[k]], numeric(1)); run_one(y, list(cluster = cl, feature = features[k], result_scope = "within_cluster_degree")) })
    write.csv(add_fdr(dplyr::bind_rows(res)), file.path(cl_dir, "DGLM_DK318_degree_within_cluster.csv"), row.names = FALSE)
  }
  if (analysis_type %in% c("both", "edge") && length(features) >= 2) {
    pairs <- t(combn(idx, 2)); res <- vector("list", nrow(pairs))
    for (k in seq_len(nrow(pairs))) { if (k %% 500 == 0 || k == 1 || k == nrow(pairs)) cat("edge", k, "/", nrow(pairs), "\n"); ii <- pairs[k, 1]; jj <- pairs[k, 2]; y <- vapply(subj$mind_file, function(fp) read_mind(fp)[ii, jj], numeric(1)); res[[k]] <- run_one(y, list(cluster = cl, edge_id = k, i = ii, j = jj, region_i = roi_names[ii], region_j = roi_names[jj], feature = paste(roi_names[ii], roi_names[jj], sep = "__"), result_scope = "within_cluster_edge")) }
    write.csv(add_fdr(dplyr::bind_rows(res)), file.path(cl_dir, "DGLM_DK318_edge_within_cluster.csv"), row.names = FALSE)
  }
}
cat("\nDONE. Results saved to:", out_root, "\n")
