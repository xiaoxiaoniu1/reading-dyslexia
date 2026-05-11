# ============================================================
# DK318 MIND Degree NMF - DD Only
# Purpose: Run NMF only in DD participants to explore subtypes.
# Outputs: K metrics, W/H matrices, dominant subtype assignments, top ROIs, QC.
# ============================================================

cat("Starting DK318 MIND Degree NMF - DD only...\n\n")

packages <- c("readxl", "dplyr", "NMF", "car", "ggplot2", "reshape2")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}
library(readxl); library(dplyr); library(NMF); library(car); library(ggplot2); library(reshape2)
options(contrasts = c("contr.sum", "contr.poly"))

demo_file <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/all_data_cqt_mean_1.5.xlsx"
mind_combat_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_combat"
out_dir <- Sys.getenv("DD_NMF_OUT_DIR", "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NMF_DD_only")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

K_range <- 2:10
K_final <- as.integer(Sys.getenv("NMF_K_FINAL", "4"))
nrun <- as.integer(Sys.getenv("NMF_NRUN", "100"))
nmf_method <- Sys.getenv("NMF_METHOD", "brunet")
nmf_seed <- as.integer(Sys.getenv("NMF_SEED", "2026"))
if (!K_final %in% K_range) stop("NMF_K_FINAL must be in ", paste(K_range, collapse = ", "))

cat("Output directory:", out_dir, "\n")
cat("K range:", paste(K_range, collapse = ", "), "| K final:", K_final, "| nrun:", nrun, "\n\n")

cat("Loading demographic data...\n")
demo_raw <- read_excel(demo_file, sheet = "Sheet1")
df_all <- demo_raw %>%
  mutate(
    original_project = as.character(`original-project`),
    id_old = as.character(id_old),
    subj_prefix = paste0(original_project, "_", id_old),
    file_base = paste0(subj_prefix, "_MIND_DK318_combat"),
    degree_file = file.path(mind_combat_dir, paste0(file_base, "_degree.csv")),
    Diagnosis = factor(ifelse(group_d_or_c == 0, "TD", "DD"), levels = c("TD", "DD")),
    AgeGroup = factor(ifelse(group_age == 1, "Adult", "Child"), levels = c("Child", "Adult")),
    Sex = factor(ifelse(sex == 1, "Male", "Female"), levels = c("Female", "Male")),
    has_file = file.exists(degree_file)
  )

df <- df_all %>%
  filter(!is.na(original_project), !is.na(id_old), original_project != "", id_old != "",
         has_file, Diagnosis == "DD", !is.na(AgeGroup), !is.na(Sex)) %>%
  droplevels()

cat("Total subjects:", nrow(demo_raw), "\n")
cat("Subjects with degree files:", sum(df_all$has_file, na.rm = TRUE), "\n")
cat("DD subjects used:", nrow(df), "\n")
print(with(df, table(AgeGroup, Sex)))
cat("\n")
if (nrow(df) < max(K_range)) stop("Too few DD subjects for requested K range.")

read_degree_csv <- function(fp) {
  d <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE)
  vals <- as.numeric(d$degree)
  names(vals) <- as.character(d$ROI)
  vals
}

cat("Loading degree data...\n")
roi_names <- names(read_degree_csv(df$degree_file[1]))
X <- matrix(NA_real_, nrow = nrow(df), ncol = length(roi_names), dimnames = list(df$file_base, roi_names))
for (s in seq_len(nrow(df))) {
  if (s %% 50 == 0 || s == 1 || s == nrow(df)) cat("  Subject", s, "/", nrow(df), "\n")
  X[s, ] <- read_degree_csv(df$degree_file[s])
}
if (any(is.na(X))) stop("Degree matrix contains NA values.")

min_raw <- min(X); max_raw <- max(X); n_negative <- sum(X < 0)
cat("Raw degree range:", min_raw, "to", max_raw, "| negative values:", n_negative, "\n")

scaling_info <- data.frame(ROI = roi_names, min_value = NA_real_, max_value = NA_real_, range = NA_real_, was_constant = FALSE)
X_scaled <- X
for (j in seq_along(roi_names)) {
  mn <- min(X[, j]); mx <- max(X[, j]); rg <- mx - mn
  scaling_info[j, c("min_value", "max_value", "range")] <- c(mn, mx, rg)
  if (rg < 1e-12) {
    scaling_info$was_constant[j] <- TRUE
    X_scaled[, j] <- 0
  } else {
    X_scaled[, j] <- (X[, j] - mn) / rg
  }
}
write.csv(scaling_info, file.path(out_dir, "NMF_DD_only_degree_scaling_info.csv"), row.names = FALSE)
cat("Scaled range:", min(X_scaled), "to", max(X_scaled), "\n\n")

cat("Evaluating K...\n")
k_metrics <- data.frame()
nmf_fits <- list()
for (k in K_range) {
  cat("  Fitting K =", k, "\n")
  fit <- nmf(X_scaled, rank = k, method = nmf_method, nrun = nrun, seed = nmf_seed)
  nmf_fits[[as.character(k)]] <- fit
  k_metrics <- rbind(k_metrics, data.frame(
    K = k,
    reconstruction_error = fit@residuals,
    rss = sum((X_scaled - fitted(fit))^2),
    cophenetic = tryCatch(cophcor(fit), error = function(e) NA_real_),
    dispersion = tryCatch(dispersion(fit), error = function(e) NA_real_)
  ))
  print(tail(k_metrics, 1))
}
write.csv(k_metrics, file.path(out_dir, "NMF_DD_only_K_selection_metrics.csv"), row.names = FALSE)

plot_metric <- function(y, ylab, file, col) {
  if (all(is.na(y))) return(invisible(NULL))
  png(file.path(out_dir, file), width = 800, height = 600, res = 120)
  plot(k_metrics$K, y, type = "b", pch = 19, xlab = "K", ylab = ylab,
       main = paste("DD-only NMF K Selection:", ylab), col = col, lwd = 2)
  grid(); dev.off()
}
plot_metric(k_metrics$reconstruction_error, "Reconstruction Error", "NMF_DD_only_K_selection_reconstruction_error.png", "steelblue")
plot_metric(k_metrics$cophenetic, "Cophenetic", "NMF_DD_only_K_selection_cophenetic.png", "darkgreen")
plot_metric(k_metrics$dispersion, "Dispersion", "NMF_DD_only_K_selection_dispersion.png", "purple")

cat("Extracting final model K =", K_final, "\n")
fit_final <- nmf_fits[[as.character(K_final)]]
W <- basis(fit_final); H <- coef(fit_final)
component_names <- paste0("SubtypeComponent_", seq_len(K_final))
colnames(W) <- component_names; rownames(H) <- component_names

cat("Testing AgeGroup/Sex effects on component scores...\n")
anova_results <- data.frame()
for (k in seq_len(K_final)) {
  lm_fit <- lm(score ~ AgeGroup + Sex, data = df %>% mutate(score = W[, k]))
  anova_k <- car::Anova(lm_fit, type = 3)
  for (term in rownames(anova_k)[-1]) {
    anova_results <- rbind(anova_results, data.frame(
      component = component_names[k], term = term, df = anova_k[term, "Df"],
      statistic = anova_k[term, "F value"], p_value = anova_k[term, "Pr(>F)"]
    ))
  }
}
anova_results <- anova_results %>% group_by(term) %>% mutate(q_value_by_term = p.adjust(p_value, "fdr")) %>% ungroup()
anova_results$q_value_all_terms <- p.adjust(anova_results$p_value, "fdr")
write.csv(anova_results, file.path(out_dir, "NMF_DD_only_component_lm_anova_results.csv"), row.names = FALSE)

W_prop <- W / rowSums(W)
colnames(W_prop) <- paste0(component_names, "_proportion")
dom_idx <- max.col(W_prop, ties.method = "first")
W_df <- df %>% select(file_base, original_project, id_old, Diagnosis, AgeGroup, Sex) %>%
  bind_cols(as.data.frame(W)) %>% bind_cols(as.data.frame(W_prop)) %>%
  mutate(dominant_subtype = component_names[dom_idx],
         dominant_loading = W[cbind(seq_len(nrow(W)), dom_idx)],
         dominant_proportion = W_prop[cbind(seq_len(nrow(W_prop)), dom_idx)],
         second_proportion = apply(W_prop, 1, function(x) sort(x, decreasing = TRUE)[2]),
         subtype_margin = dominant_proportion - second_proportion)
write.csv(W_df, file.path(out_dir, "NMF_DD_only_W_subject_loadings_and_subtypes.csv"), row.names = FALSE)
write.csv(W_df, file.path(out_dir, "NMF_W_subject_loadings.csv"), row.names = FALSE)

subtype_summary <- W_df %>% group_by(dominant_subtype) %>% summarise(
  n_subjects = n(), percent = 100 * n() / nrow(W_df),
  n_child = sum(AgeGroup == "Child"), n_adult = sum(AgeGroup == "Adult"),
  n_female = sum(Sex == "Female"), n_male = sum(Sex == "Male"),
  mean_dominant_proportion = mean(dominant_proportion), mean_subtype_margin = mean(subtype_margin), .groups = "drop")
write.csv(subtype_summary, file.path(out_dir, "NMF_DD_only_subtype_summary.csv"), row.names = FALSE)

H_wide <- as.data.frame(H); H_wide$component <- component_names; H_wide <- H_wide[, c("component", roi_names)]
H_long <- reshape2::melt(H, varnames = c("component", "ROI"), value.name = "weight")
H_long$component <- component_names[H_long$component]
H_long <- H_long %>% group_by(component) %>% arrange(desc(weight)) %>% mutate(rank_within_component = row_number()) %>% ungroup()
write.csv(H_wide, file.path(out_dir, "NMF_DD_only_H_component_roi_weights.csv"), row.names = FALSE)
write.csv(H_long, file.path(out_dir, "NMF_DD_only_H_component_roi_weights_long.csv"), row.names = FALSE)
write.csv(H_wide, file.path(out_dir, "NMF_H_component_roi_weights.csv"), row.names = FALSE)
write.csv(H_long, file.path(out_dir, "NMF_H_component_roi_weights_long.csv"), row.names = FALSE)

top_rois <- data.frame()
for (comp in component_names) {
  comp_weights <- H_long %>% filter(component == comp)
  top_rois <- rbind(top_rois,
                    comp_weights %>% slice(1:min(10, n())) %>% mutate(top_type = "top10"),
                    comp_weights %>% slice(1:min(20, n())) %>% mutate(top_type = "top20"),
                    comp_weights %>% slice(1:ceiling(0.1 * length(roi_names))) %>% mutate(top_type = "top10percent"))
}
write.csv(top_rois, file.path(out_dir, "NMF_DD_only_component_top_ROIs.csv"), row.names = FALSE)
write.csv(top_rois, file.path(out_dir, "NMF_component_top_ROIs.csv"), row.names = FALSE)

roi_assignment <- data.frame(ROI = roi_names, assigned_component = NA_character_, max_weight = NA_real_, second_max_weight = NA_real_, assignment_margin = NA_real_)
for (j in seq_along(roi_names)) {
  sw <- sort(H[, j], decreasing = TRUE); idx <- which.max(H[, j])
  roi_assignment$assigned_component[j] <- component_names[idx]
  roi_assignment$max_weight[j] <- sw[1]; roi_assignment$second_max_weight[j] <- sw[2]
  roi_assignment$assignment_margin[j] <- sw[1] - sw[2]
}
write.csv(roi_assignment, file.path(out_dir, "NMF_DD_only_ROI_hard_assignment.csv"), row.names = FALSE)
write.csv(roi_assignment, file.path(out_dir, "NMF_ROI_hard_assignment.csv"), row.names = FALSE)

W_long <- W_df %>% select(file_base, AgeGroup, Sex, dominant_subtype, all_of(component_names)) %>%
  reshape2::melt(id.vars = c("file_base", "AgeGroup", "Sex", "dominant_subtype"), variable.name = "component", value.name = "score")
ggsave(file.path(out_dir, "NMF_DD_only_component_scores_boxplot.png"),
       ggplot(W_long, aes(AgeGroup, score, fill = Sex)) + geom_boxplot(outlier.size = 0.5) +
         geom_point(position = position_jitterdodge(jitter.width = 0.15), alpha = 0.35, size = 0.8) +
         facet_wrap(~ component, scales = "free_y", ncol = 2) + theme_bw() + labs(title = "DD-only NMF Component Scores"),
       width = 10, height = 8, dpi = 150)
ggsave(file.path(out_dir, "NMF_DD_only_dominant_subtype_counts.png"),
       ggplot(W_df, aes(dominant_subtype, fill = AgeGroup)) + geom_bar(position = "dodge") + theme_bw() +
         labs(title = "DD Dominant NMF Subtypes", x = "Dominant subtype", y = "Number of DD subjects"),
       width = 8, height = 6, dpi = 150)

qc_summary <- data.frame(
  metric = c("population", "n_subjects_total_in_demo", "n_subjects_with_degree_file", "n_DD_subjects_used", "n_roi", "min_raw_degree", "max_raw_degree", "min_scaled_degree", "max_scaled_degree", "n_negative_raw_values", "n_constant_rois", "K_range_tested", "K_final", "NMF_method", "nrun", "seed"),
  value = as.character(c("DD_only", nrow(demo_raw), sum(df_all$has_file, na.rm = TRUE), nrow(df), length(roi_names), min_raw, max_raw, min(X_scaled), max(X_scaled), n_negative, sum(scaling_info$was_constant), paste(K_range, collapse = "-"), K_final, nmf_method, nrun, nmf_seed))
)
write.csv(qc_summary, file.path(out_dir, "NMF_DD_only_degree_QC_summary.csv"), row.names = FALSE)

cat("\n==============================\nDD-ONLY NMF COMPLETE\n==============================\n")
cat("Output directory:", out_dir, "\n")
cat("DD subjects:", nrow(df), "| ROIs:", length(roi_names), "| K final:", K_final, "\n\n")
print(subtype_summary, row.names = FALSE)
cat("\nDONE\n")
