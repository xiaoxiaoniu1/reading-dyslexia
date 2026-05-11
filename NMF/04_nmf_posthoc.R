# ============================================================
# DK318 MIND Degree NMF - Post-hoc Contrasts
# Purpose: Detailed post-hoc tests for significant components
# Input: NMF W matrix and demographic data
# Output: Pairwise comparisons and effect sizes
# ============================================================

cat("Starting NMF Post-hoc Analysis...\n\n")

# ============================================================
# 1. Load Packages
# ============================================================

packages <- c("dplyr", "emmeans", "effectsize", "broom")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

library(dplyr); library(emmeans); library(effectsize); library(broom)

# ============================================================
# 2. Define Paths
# ============================================================

nmf_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NMF_3"
out_dir <- nmf_dir

cat("NMF directory:", nmf_dir, "\n\n")

# ============================================================
# 3. Load Data
# ============================================================

cat("Loading data...\n")

# Load W matrix (subject loadings)
W_df <- read.csv(file.path(nmf_dir, "NMF_W_subject_loadings.csv"), stringsAsFactors = FALSE)

# Load ANOVA results
anova_results <- read.csv(file.path(nmf_dir, "NMF_component_lm_anova_results.csv"), stringsAsFactors = FALSE)

# Get component names
component_cols <- grep("^Component_", colnames(W_df), value = TRUE)
n_components <- length(component_cols)

# Ensure factors are properly coded
W_df$Diagnosis <- factor(W_df$Diagnosis, levels = c("TD", "DD"))
W_df$AgeGroup <- factor(W_df$AgeGroup, levels = c("Child", "Adult"))
W_df$Sex <- factor(W_df$Sex, levels = c("Female", "Male"))

cat("Loaded data for", n_components, "components\n")
cat("Number of subjects:", nrow(W_df), "\n\n")

# ============================================================
# 4. Identify Components for Post-hoc Testing
# ============================================================

cat("Identifying components for post-hoc testing...\n")

# Components with significant interaction (q < 0.05)
sig_interaction <- anova_results %>%
  filter(term == "AgeGroup:Diagnosis", q_value_by_term < 0.05) %>%
  pull(component)

# Components with significant main effects (q < 0.05)
sig_diagnosis <- anova_results %>%
  filter(term == "Diagnosis", q_value_by_term < 0.05) %>%
  pull(component)

sig_agegroup <- anova_results %>%
  filter(term == "AgeGroup", q_value_by_term < 0.05) %>%
  pull(component)

# All components with any significant effect
sig_components <- unique(c(sig_interaction, sig_diagnosis, sig_agegroup))

cat("Components with significant interaction:", length(sig_interaction), "\n")
cat("Components with significant Diagnosis effect:", length(sig_diagnosis), "\n")
cat("Components with significant AgeGroup effect:", length(sig_agegroup), "\n")
cat("Total components for post-hoc testing:", length(sig_components), "\n\n")

if (length(sig_components) == 0) {
  cat("No significant components found. Testing all components...\n\n")
  sig_components <- component_cols
}

# ============================================================
# 5. Post-hoc Contrasts
# ============================================================

cat("Running post-hoc contrasts...\n\n")

posthoc_results <- data.frame()

for (comp in sig_components) {
  cat("  Processing", comp, "...\n")
  
  # Fit linear model
  formula_str <- paste(comp, "~ AgeGroup * Diagnosis + Sex")
  lm_fit <- lm(as.formula(formula_str), data = W_df)
  
  # Estimated marginal means
  emm <- emmeans(lm_fit, ~ Diagnosis | AgeGroup)
  
  # Pairwise contrasts: TD vs DD within each age group
  contrasts_by_age <- pairs(emm, adjust = "none")
  contrasts_df <- as.data.frame(summary(contrasts_by_age, infer = c(TRUE, TRUE)))
  
  # Add component name
  contrasts_df$component <- comp
  
  # Calculate Cohen's d for each contrast
  for (i in 1:nrow(contrasts_df)) {
    age_group <- contrasts_df$AgeGroup[i]
    
    # Get data for this age group
    td_data <- W_df[[comp]][W_df$AgeGroup == age_group & W_df$Diagnosis == "TD"]
    dd_data <- W_df[[comp]][W_df$AgeGroup == age_group & W_df$Diagnosis == "DD"]
    
    # Calculate Cohen's d
    pooled_sd <- sqrt(((length(td_data) - 1) * var(td_data) + 
                       (length(dd_data) - 1) * var(dd_data)) / 
                      (length(td_data) + length(dd_data) - 2))
    
    cohens_d <- (mean(td_data) - mean(dd_data)) / pooled_sd
    
    contrasts_df$cohens_d[i] <- cohens_d
  }
  
  posthoc_results <- rbind(posthoc_results, contrasts_df)
  
  # Also test the interaction contrast: (Adult_TD - Adult_DD) - (Child_TD - Child_DD)
  emm_full <- emmeans(lm_fit, ~ Diagnosis * AgeGroup)
  interaction_contrast <- contrast(emm_full, 
                                   method = list("Adult_minus_Child_change" = c(-1, 1, 1, -1)),
                                   adjust = "none")
  
  int_df <- as.data.frame(summary(interaction_contrast, infer = c(TRUE, TRUE)))
  int_df$component <- comp
  int_df$contrast_type <- "Interaction"
  
  # Store interaction contrast separately
  if (!exists("interaction_results")) {
    interaction_results <- int_df
  } else {
    interaction_results <- rbind(interaction_results, int_df)
  }
}

cat("\n")

# ============================================================
# 6. Apply FDR Correction
# ============================================================

cat("Applying FDR correction to post-hoc tests...\n")

posthoc_results$q_value <- p.adjust(posthoc_results$p.value, method = "fdr")

if (exists("interaction_results")) {
  interaction_results$q_value <- p.adjust(interaction_results$p.value, method = "fdr")
}

# ============================================================
# 7. Save Results
# ============================================================

cat("Saving post-hoc results...\n")

# Save pairwise contrasts
posthoc_file <- file.path(out_dir, "NMF_posthoc_pairwise_contrasts.csv")
write.csv(posthoc_results, posthoc_file, row.names = FALSE)
cat("Saved:", posthoc_file, "\n")

# Save interaction contrasts
if (exists("interaction_results")) {
  interaction_file <- file.path(out_dir, "NMF_posthoc_interaction_contrasts.csv")
  write.csv(interaction_results, interaction_file, row.names = FALSE)
  cat("Saved:", interaction_file, "\n")
}

cat("\n")

# ============================================================
# 8. Descriptive Statistics
# ============================================================

cat("Calculating descriptive statistics...\n")

descriptives <- data.frame()

for (comp in component_cols) {
  comp_desc <- W_df %>%
    group_by(AgeGroup, Diagnosis) %>%
    summarise(
      n = n(),
      mean = mean(.data[[comp]], na.rm = TRUE),
      sd = sd(.data[[comp]], na.rm = TRUE),
      se = sd / sqrt(n),
      median = median(.data[[comp]], na.rm = TRUE),
      min = min(.data[[comp]], na.rm = TRUE),
      max = max(.data[[comp]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(component = comp)
  
  descriptives <- rbind(descriptives, comp_desc)
}

# Save descriptives
desc_file <- file.path(out_dir, "NMF_descriptive_statistics.csv")
write.csv(descriptives, desc_file, row.names = FALSE)
cat("Saved:", desc_file, "\n\n")

# ============================================================
# 9. Effect Size Summary
# ============================================================

cat("Creating effect size summary...\n")

effect_size_summary <- posthoc_results %>%
  select(component, AgeGroup, contrast, estimate, cohens_d, p.value, q_value) %>%
  mutate(
    effect_size_interpretation = case_when(
      abs(cohens_d) < 0.2 ~ "negligible",
      abs(cohens_d) < 0.5 ~ "small",
      abs(cohens_d) < 0.8 ~ "medium",
      TRUE ~ "large"
    )
  )

effect_file <- file.path(out_dir, "NMF_effect_sizes.csv")
write.csv(effect_size_summary, effect_file, row.names = FALSE)
cat("Saved:", effect_file, "\n\n")

# ============================================================
# 10. Summary Report
# ============================================================

cat("==============================\n")
cat("POST-HOC ANALYSIS SUMMARY\n")
cat("==============================\n\n")

cat("Significant pairwise contrasts (q < 0.05):\n")
sig_posthoc <- posthoc_results %>%
  filter(q_value < 0.05) %>%
  arrange(q_value)

if (nrow(sig_posthoc) > 0) {
  print(sig_posthoc %>% 
          select(component, AgeGroup, contrast, estimate, cohens_d, p.value, q_value), 
        row.names = FALSE)
} else {
  cat("  No significant pairwise contrasts at FDR < 0.05\n")
}

cat("\n")

if (exists("interaction_results")) {
  cat("Interaction contrasts (Adult-Child developmental change):\n")
  sig_int <- interaction_results %>%
    filter(q_value < 0.05) %>%
    arrange(q_value)
  
  if (nrow(sig_int) > 0) {
    print(sig_int %>% 
            select(component, estimate, SE, t.ratio, p.value, q_value), 
          row.names = FALSE)
  } else {
    cat("  No significant interaction contrasts at FDR < 0.05\n")
  }
}

cat("\n")

cat("Effect sizes (Cohen's d) for significant contrasts:\n")
if (nrow(sig_posthoc) > 0) {
  print(sig_posthoc %>% 
          select(component, AgeGroup, cohens_d) %>%
          mutate(cohens_d = round(cohens_d, 3)), 
        row.names = FALSE)
}

cat("\n==============================\n")
cat("DONE\n")
cat("==============================\n")
