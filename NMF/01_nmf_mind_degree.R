# ============================================================
# DK318 MIND Degree NMF Analysis
# Purpose: Data-driven NMF decomposition of MIND degree patterns
# Output: NMF components and statistical tests for AgeGroup × Diagnosis
# Model: score_k ~ AgeGroup * Diagnosis + Sex
# ============================================================

cat("Starting DK318 MIND Degree NMF Analysis...\n\n")

# ============================================================
# 1. Load and Install Packages
# ============================================================

packages <- c("readxl", "dplyr", "stringr", "NMF", "broom", "car", "ggplot2", "reshape2")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

library(readxl); library(dplyr); library(stringr); library(NMF)
library(broom); library(car); library(ggplot2); library(reshape2)

# Set contrasts for Type III ANOVA
options(contrasts = c("contr.sum", "contr.poly"))

# ============================================================
# 2. Define Paths and Options
# ============================================================

demo_file <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/all_data_cqt_mean_1.5.xlsx"
mind_combat_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_combat"
out_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NMF_3"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# NMF parameters
K_range <- 2:10
nrun <- 100
nmf_method <- "brunet"
nmf_seed <- 2026
K_final <- 3 # Default, can be adjusted after K selection

cat("Output directory:", out_dir, "\n")
cat("K range:", paste(K_range, collapse = ", "), "\n")
cat("NMF runs per K:", nrun, "\n\n")

# ============================================================
# 3. Load Demographic Data and Recode Variables
# ============================================================

cat("Loading demographic data...\n")
df <- read_excel(demo_file, sheet = "Sheet1") %>%
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
  ) %>%
  filter(!is.na(original_project), !is.na(id_old), 
         original_project != "", id_old != "", has_file,
         !is.na(Diagnosis), !is.na(AgeGroup), !is.na(Sex))

n_subjects_total <- nrow(read_excel(demo_file, sheet = "Sheet1"))
n_subjects_used <- nrow(df)

cat("Total subjects in demo file:", n_subjects_total, "\n")
cat("Subjects with valid degree files and covariates:", n_subjects_used, "\n\n")

cat("Cell counts:\n")
print(with(df, table(Diagnosis, AgeGroup, Sex)))
cat("\n")

# ============================================================
# 4. Read Degree Files and Build X Matrix
# ============================================================

cat("Loading degree data...\n")

read_degree_csv <- function(fp) {
  d <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE)
  vals <- as.numeric(d$degree)
  names(vals) <- as.character(d$ROI)
  vals
}

# Read first file to get ROI names
tmp_deg <- read_degree_csv(df$degree_file[1])
roi_names <- names(tmp_deg)
n_roi <- length(roi_names)
n_subj <- nrow(df)

cat("Building", n_subj, "×", n_roi, "degree matrix...\n")

X <- matrix(NA_real_, nrow = n_subj, ncol = n_roi, 
            dimnames = list(df$file_base, roi_names))

for (s in seq_len(n_subj)) {
  if (s %% 50 == 0 || s == 1 || s == n_subj) {
    cat("  Reading subject", s, "/", n_subj, "\n")
  }
  X[s, ] <- read_degree_csv(df$degree_file[s])
}

cat("\nDegree matrix dimensions:", paste(dim(X), collapse = " × "), "\n")
cat("Missing values:", sum(is.na(X)), "\n")
cat("Raw degree range: [", min(X, na.rm = TRUE), ",", max(X, na.rm = TRUE), "]\n\n")

# Check for NA values
if (any(is.na(X))) {
  stop("ERROR: Degree matrix contains NA values. Please check input files.")
}

# ============================================================
# 5. Check Non-negativity and Scale X
# ============================================================

cat("Checking non-negativity...\n")
min_raw <- min(X, na.rm = TRUE)
max_raw <- max(X, na.rm = TRUE)
n_negative <- sum(X < 0, na.rm = TRUE)

cat("Min raw degree:", min_raw, "\n")
cat("Max raw degree:", max_raw, "\n")
cat("Number of negative values:", n_negative, "\n\n")

# Scale each ROI to [0, 1] using min-max scaling
cat("Applying min-max scaling to [0, 1] for each ROI...\n")

scaling_info <- data.frame(
  ROI = roi_names,
  min_value = NA_real_,
  max_value = NA_real_,
  range = NA_real_,
  was_constant = FALSE,
  stringsAsFactors = FALSE
)

X_scaled <- X
for (j in seq_len(n_roi)) {
  col_vals <- X[, j]
  col_min <- min(col_vals, na.rm = TRUE)
  col_max <- max(col_vals, na.rm = TRUE)
  col_range <- col_max - col_min
  
  scaling_info$min_value[j] <- col_min
  scaling_info$max_value[j] <- col_max
  scaling_info$range[j] <- col_range
  
  if (col_range < 1e-12) {
    scaling_info$was_constant[j] <- TRUE
    X_scaled[, j] <- 0
  } else {
    X_scaled[, j] <- (col_vals - col_min) / col_range
  }
}

n_constant <- sum(scaling_info$was_constant)
cat("Number of constant ROIs:", n_constant, "\n")
cat("Scaled degree range: [", min(X_scaled, na.rm = TRUE), ",", 
    max(X_scaled, na.rm = TRUE), "]\n\n")

# Save scaling info
scaling_file <- file.path(out_dir, "NMF_degree_scaling_info.csv")
write.csv(scaling_info, scaling_file, row.names = FALSE)
cat("Saved scaling info:", scaling_file, "\n\n")

# ============================================================
# 6. Evaluate K from 2 to 10
# ============================================================

cat("Evaluating NMF for K =", paste(K_range, collapse = ", "), "...\n")
cat("This may take several minutes...\n\n")

k_metrics <- data.frame(
  K = integer(),
  reconstruction_error = numeric(),
  rss = numeric(),
  cophenetic = numeric(),
  stringsAsFactors = FALSE
)

nmf_fits <- list()

for (k in K_range) {
  cat("Fitting NMF with K =", k, "(", nrun, "runs)...\n")
  
  fit <- nmf(X_scaled, rank = k, method = nmf_method, nrun = nrun, seed = nmf_seed)
  nmf_fits[[as.character(k)]] <- fit
  
  # Extract metrics
  recon_error <- fit@residuals
  rss_val <- sum((X_scaled - fitted(fit))^2)
  
  # Cophenetic coefficient (if available)
  coph <- tryCatch({
    cophcor(fit)
  }, error = function(e) NA_real_)
  
  k_metrics <- rbind(k_metrics, data.frame(
    K = k,
    reconstruction_error = recon_error,
    rss = rss_val,
    cophenetic = coph,
    stringsAsFactors = FALSE
  ))
  
  cat("  K =", k, "| Reconstruction error:", round(recon_error, 4), 
      "| RSS:", round(rss_val, 2), "| Cophenetic:", round(coph, 3), "\n")
}

cat("\n")

# Save K selection metrics
k_metrics_file <- file.path(out_dir, "NMF_K_selection_metrics.csv")
write.csv(k_metrics, k_metrics_file, row.names = FALSE)
cat("Saved K selection metrics:", k_metrics_file, "\n\n")

# Plot reconstruction error
png(file.path(out_dir, "NMF_K_selection_reconstruction_error.png"), 
    width = 800, height = 600, res = 120)
plot(k_metrics$K, k_metrics$reconstruction_error, type = "b", pch = 19,
     xlab = "Number of Components (K)", ylab = "Reconstruction Error",
     main = "NMF K Selection: Reconstruction Error", col = "steelblue", lwd = 2)
grid()
dev.off()

# Plot cophenetic coefficient if available
if (any(!is.na(k_metrics$cophenetic))) {
  png(file.path(out_dir, "NMF_K_selection_cophenetic.png"), 
      width = 800, height = 600, res = 120)
  plot(k_metrics$K, k_metrics$cophenetic, type = "b", pch = 19,
       xlab = "Number of Components (K)", ylab = "Cophenetic Coefficient",
       main = "NMF K Selection: Cophenetic Coefficient", col = "darkgreen", lwd = 2)
  grid()
  dev.off()
  cat("Saved K selection plots\n\n")
}

# ============================================================
# 7. Fit Final NMF Model with K_final
# ============================================================

cat("Fitting final NMF model with K =", K_final, "...\n")

if (as.character(K_final) %in% names(nmf_fits)) {
  fit_final <- nmf_fits[[as.character(K_final)]]
  cat("Using previously fitted model for K =", K_final, "\n\n")
} else {
  fit_final <- nmf(X_scaled, rank = K_final, method = nmf_method, 
                   nrun = nrun, seed = nmf_seed)
  cat("Fitted new model for K =", K_final, "\n\n")
}

# ============================================================
# 8. Extract W and H
# ============================================================

cat("Extracting W (subject loadings) and H (component ROI weights)...\n")

W <- basis(fit_final)      # n × K
H <- coef(fit_final)       # K × 318

cat("W dimensions:", paste(dim(W), collapse = " × "), "\n")
cat("H dimensions:", paste(dim(H), collapse = " × "), "\n\n")

# Verify dimensions
if (nrow(W) != n_subj || ncol(W) != K_final) {
  stop("ERROR: W matrix dimensions incorrect")
}
if (nrow(H) != K_final || ncol(H) != n_roi) {
  stop("ERROR: H matrix dimensions incorrect")
}

# Name components
component_names <- paste0("Component_", 1:K_final)
colnames(W) <- component_names
rownames(H) <- component_names

# ============================================================
# 9. Run lm or ANOVA for Each Component Score
# ============================================================

cat("Running statistical tests for each component...\n")
cat("Model: score_k ~ AgeGroup * Diagnosis + Sex\n\n")

anova_results <- data.frame()

for (k in 1:K_final) {
  score_k <- W[, k]
  
  # Build data frame
  test_df <- df %>%
    mutate(score = score_k)
  
  # Fit linear model
  lm_fit <- lm(score ~ AgeGroup * Diagnosis + Sex, data = test_df)
  
  # Type III ANOVA
  anova_k <- car::Anova(lm_fit, type = 3)
  
  # Extract results
  terms <- rownames(anova_k)[-1]  # Exclude intercept
  for (term in terms) {
    anova_results <- rbind(anova_results, data.frame(
      component = component_names[k],
      term = term,
      df = anova_k[term, "Df"],
      statistic = anova_k[term, "F value"],
      p_value = anova_k[term, "Pr(>F)"],
      stringsAsFactors = FALSE
    ))
  }
  
  cat("  Component", k, "done\n")
}

cat("\n")

# ============================================================
# 10. Apply FDR Correction
# ============================================================

cat("Applying FDR correction...\n")

# FDR by term
anova_results <- anova_results %>%
  group_by(term) %>%
  mutate(q_value_by_term = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

# FDR across all tests
anova_results$q_value_all_terms <- p.adjust(anova_results$p_value, method = "fdr")

cat("FDR correction complete\n\n")

# Save ANOVA results
anova_file <- file.path(out_dir, "NMF_component_lm_anova_results.csv")
write.csv(anova_results, anova_file, row.names = FALSE)
cat("Saved ANOVA results:", anova_file, "\n\n")

# Print summary
cat("Significant results (q_value_by_term < 0.05):\n")
sig_results <- anova_results %>% filter(q_value_by_term < 0.05)
if (nrow(sig_results) > 0) {
  print(sig_results, row.names = FALSE)
} else {
  cat("  No significant results at FDR < 0.05\n")
}
cat("\n")

# ============================================================
# 11. Save W, H, Top ROI, ROI Assignment
# ============================================================

cat("Saving W matrix (subject loadings)...\n")

W_df <- df %>%
  select(file_base, original_project, id_old, Diagnosis, AgeGroup, Sex) %>%
  bind_cols(as.data.frame(W))

w_file <- file.path(out_dir, "NMF_W_subject_loadings.csv")
write.csv(W_df, w_file, row.names = FALSE)
cat("Saved:", w_file, "\n\n")

cat("Saving H matrix (component ROI weights)...\n")

# Wide format
H_wide <- as.data.frame(H)
H_wide$component <- component_names
H_wide <- H_wide[, c("component", roi_names)]

h_wide_file <- file.path(out_dir, "NMF_H_component_roi_weights.csv")
write.csv(H_wide, h_wide_file, row.names = FALSE)
cat("Saved:", h_wide_file, "\n")

# Long format
H_long <- reshape2::melt(H, varnames = c("component", "ROI"), value.name = "weight")
H_long$component <- component_names[H_long$component]
H_long <- H_long %>%
  group_by(component) %>%
  arrange(desc(weight)) %>%
  mutate(rank_within_component = row_number()) %>%
  ungroup()

h_long_file <- file.path(out_dir, "NMF_H_component_roi_weights_long.csv")
write.csv(H_long, h_long_file, row.names = FALSE)
cat("Saved:", h_long_file, "\n\n")

# Top ROIs
cat("Extracting top ROIs for each component...\n")

top_rois <- data.frame()

for (k in 1:K_final) {
  comp_weights <- H_long %>% filter(component == component_names[k])
  
  # Top 10
  top10 <- comp_weights %>% slice(1:min(10, n())) %>% mutate(top_type = "top10")
  
  # Top 20
  top20 <- comp_weights %>% slice(1:min(20, n())) %>% mutate(top_type = "top20")
  
  # Top 10%
  n_top10pct <- ceiling(0.1 * n_roi)
  top10pct <- comp_weights %>% slice(1:n_top10pct) %>% mutate(top_type = "top10percent")
  
  top_rois <- rbind(top_rois, top10, top20, top10pct)
}

top_rois_file <- file.path(out_dir, "NMF_component_top_ROIs.csv")
write.csv(top_rois, top_rois_file, row.names = FALSE)
cat("Saved:", top_rois_file, "\n\n")

# ROI hard assignment
cat("Assigning each ROI to primary component...\n")

roi_assignment <- data.frame(
  ROI = roi_names,
  assigned_component = NA_character_,
  max_weight = NA_real_,
  second_max_weight = NA_real_,
  assignment_margin = NA_real_,
  stringsAsFactors = FALSE
)

for (j in 1:n_roi) {
  weights_j <- H[, j]
  sorted_weights <- sort(weights_j, decreasing = TRUE)
  max_idx <- which.max(weights_j)
  
  roi_assignment$assigned_component[j] <- component_names[max_idx]
  roi_assignment$max_weight[j] <- sorted_weights[1]
  roi_assignment$second_max_weight[j] <- if (length(sorted_weights) > 1) sorted_weights[2] else NA_real_
  roi_assignment$assignment_margin[j] <- sorted_weights[1] - 
    (if (length(sorted_weights) > 1) sorted_weights[2] else 0)
}

assignment_file <- file.path(out_dir, "NMF_ROI_hard_assignment.csv")
write.csv(roi_assignment, assignment_file, row.names = FALSE)
cat("Saved:", assignment_file, "\n\n")

# ============================================================
# 12. Generate Plots
# ============================================================

cat("Generating plots...\n")

# Subject loading boxplots
cat("  Creating component scores boxplot...\n")

W_long <- W_df %>%
  select(file_base, Diagnosis, AgeGroup, Sex, all_of(component_names)) %>%
  reshape2::melt(id.vars = c("file_base", "Diagnosis", "AgeGroup", "Sex"),
                 variable.name = "component", value.name = "score")

p_boxplot <- ggplot(W_long, aes(x = AgeGroup, y = score, fill = Diagnosis)) +
  geom_boxplot(outlier.size = 0.5) +
  facet_wrap(~ component, scales = "free_y", ncol = 2) +
  labs(title = "NMF Component Scores by AgeGroup and Diagnosis",
       x = "Age Group", y = "Component Score") +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "NMF_component_scores_boxplot.png"), 
       p_boxplot, width = 10, height = 8, dpi = 150)

# Component ROI weights
cat("  Creating component ROI weights plot...\n")

p_weights <- ggplot(H_long, aes(x = rank_within_component, y = weight, color = component)) +
  geom_line(size = 1) +
  facet_wrap(~ component, scales = "free_y", ncol = 2) +
  labs(title = "NMF Component ROI Weights (Ranked)",
       x = "ROI Rank", y = "Weight") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "NMF_component_roi_weights.png"), 
       p_weights, width = 10, height = 8, dpi = 150)

cat("Plots saved\n\n")

# ============================================================
# 13. Save QC Summary
# ============================================================

cat("Saving QC summary...\n")

qc_summary <- data.frame(
  metric = c("n_subjects_total_in_demo", "n_subjects_with_degree_file", 
             "n_subjects_used", "n_roi", "min_raw_degree", "max_raw_degree",
             "min_scaled_degree", "max_scaled_degree", "n_missing_values",
             "n_negative_raw_values", "n_constant_rois", "K_range_tested",
             "K_final", "NMF_method", "nrun", "seed"),
  value = as.character(c(n_subjects_total, sum(read_excel(demo_file, sheet = "Sheet1")$has_file, na.rm = TRUE),
            n_subjects_used, n_roi, min_raw, max_raw,
            min(X_scaled, na.rm = TRUE), max(X_scaled, na.rm = TRUE),
            sum(is.na(X)), n_negative, n_constant,
            paste(K_range, collapse = "-"), K_final, nmf_method, nrun, nmf_seed)),
  stringsAsFactors = FALSE
)

qc_file <- file.path(out_dir, "NMF_degree_QC_summary.csv")
write.csv(qc_summary, qc_file, row.names = FALSE)
cat("Saved:", qc_file, "\n\n")

# ============================================================
# Final Summary
# ============================================================

cat("==============================\n")
cat("NMF ANALYSIS COMPLETE\n")
cat("==============================\n\n")

cat("Output directory:", out_dir, "\n\n")

cat("Key results:\n")
cat("  - K final:", K_final, "\n")
cat("  - Subjects analyzed:", n_subjects_used, "\n")
cat("  - ROIs:", n_roi, "\n\n")

cat("Significant components (AgeGroup:Diagnosis, q < 0.05):\n")
sig_interaction <- anova_results %>% 
  filter(term == "AgeGroup:Diagnosis", q_value_by_term < 0.05)
if (nrow(sig_interaction) > 0) {
  print(sig_interaction %>% select(component, statistic, p_value, q_value_by_term), 
        row.names = FALSE)
} else {
  cat("  No significant interactions at FDR < 0.05\n")
  cat("  Lowest p-value:\n")
  lowest_p <- anova_results %>% 
    filter(term == "AgeGroup:Diagnosis") %>%
    arrange(p_value) %>%
    slice(1)
  print(lowest_p %>% select(component, statistic, p_value, q_value_by_term), 
        row.names = FALSE)
}

cat("\n==============================\n")
cat("DONE\n")
cat("==============================\n")
