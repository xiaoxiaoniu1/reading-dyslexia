# ============================================================
# DK318 MIND Degree NMF - DD-only Post-hoc / Descriptives
# Purpose: Post-hoc tests for DD-only NMF component scores.
# Model context: score_k ~ AgeGroup + Sex
# Outputs: AgeGroup/Sex pairwise contrasts, Cohen's d, descriptives.
# ============================================================

cat("Starting DD-only NMF Post-hoc Analysis...\n\n")

packages <- c("dplyr", "emmeans", "effectsize")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}
library(dplyr); library(emmeans); library(effectsize)

nmf_dir <- Sys.getenv("DD_NMF_OUT_DIR", "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NMF_DD_only")
out_dir <- nmf_dir
cat("DD-only NMF directory:", nmf_dir, "\n\n")

W_df <- read.csv(file.path(nmf_dir, "NMF_DD_only_W_subject_loadings_and_subtypes.csv"), stringsAsFactors = FALSE)
anova_results <- read.csv(file.path(nmf_dir, "NMF_DD_only_component_lm_anova_results.csv"), stringsAsFactors = FALSE)

component_cols <- grep("^SubtypeComponent_[0-9]+$", colnames(W_df), value = TRUE)
if (length(component_cols) == 0) component_cols <- grep("^Component_[0-9]+$", colnames(W_df), value = TRUE)
if (length(component_cols) == 0) stop("No component columns found in W file.")

W_df$AgeGroup <- factor(W_df$AgeGroup, levels = c("Child", "Adult"))
W_df$Sex <- factor(W_df$Sex, levels = c("Female", "Male"))
W_df$dominant_subtype <- factor(W_df$dominant_subtype, levels = component_cols)

cat("Loaded", nrow(W_df), "DD subjects and", length(component_cols), "components.\n\n")

sig_components <- anova_results %>%
  filter(term %in% c("AgeGroup", "Sex"), q_value_by_term < 0.05) %>%
  pull(component) %>% unique()
if (length(sig_components) == 0) {
  cat("No significant AgeGroup/Sex effects found. Testing all components for descriptive completeness.\n\n")
  sig_components <- component_cols
}

calc_cohens_d_two_group <- function(values, group) {
  keep <- !is.na(values) & !is.na(group)
  values <- values[keep]
  group <- droplevels(factor(group[keep]))
  if (nlevels(group) != 2) return(NA_real_)
  g1 <- values[group == levels(group)[1]]
  g2 <- values[group == levels(group)[2]]
  if (length(g1) < 2 || length(g2) < 2) return(NA_real_)
  pooled_sd <- sqrt(((length(g1) - 1) * var(g1) + (length(g2) - 1) * var(g2)) / (length(g1) + length(g2) - 2))
  if (is.na(pooled_sd) || pooled_sd == 0) return(NA_real_)
  (mean(g2) - mean(g1)) / pooled_sd
}

cat("Running post-hoc contrasts...\n")
posthoc_results <- data.frame()

for (comp in sig_components) {
  cat("  Processing", comp, "\n")
  lm_fit <- lm(as.formula(paste(comp, "~ AgeGroup + Sex")), data = W_df)

  age_df <- as.data.frame(summary(pairs(emmeans(lm_fit, ~ AgeGroup), adjust = "none"), infer = c(TRUE, TRUE)))
  age_df$component <- comp
  age_df$term <- "AgeGroup"
  age_df$cohens_d <- calc_cohens_d_two_group(W_df[[comp]], W_df$AgeGroup)

  sex_df <- as.data.frame(summary(pairs(emmeans(lm_fit, ~ Sex), adjust = "none"), infer = c(TRUE, TRUE)))
  sex_df$component <- comp
  sex_df$term <- "Sex"
  sex_df$cohens_d <- calc_cohens_d_two_group(W_df[[comp]], W_df$Sex)

  posthoc_results <- rbind(posthoc_results, age_df, sex_df)
}

posthoc_results$q_value <- p.adjust(posthoc_results$p.value, method = "fdr")
posthoc_results <- posthoc_results %>%
  mutate(effect_size_interpretation = case_when(
    is.na(cohens_d) ~ NA_character_,
    abs(cohens_d) < 0.2 ~ "negligible",
    abs(cohens_d) < 0.5 ~ "small",
    abs(cohens_d) < 0.8 ~ "medium",
    TRUE ~ "large"
  ))

posthoc_file <- file.path(out_dir, "NMF_DD_only_posthoc_age_sex_contrasts.csv")
write.csv(posthoc_results, posthoc_file, row.names = FALSE)
write.csv(posthoc_results, file.path(out_dir, "NMF_posthoc_pairwise_contrasts.csv"), row.names = FALSE)
cat("Saved:", posthoc_file, "\n")

cat("Calculating descriptive statistics...\n")
descriptives_age_sex <- data.frame()
for (comp in component_cols) {
  comp_desc <- W_df %>%
    group_by(AgeGroup, Sex) %>%
    summarise(
      n = n(),
      mean = mean(.data[[comp]], na.rm = TRUE),
      sd = sd(.data[[comp]], na.rm = TRUE),
      se = sd / sqrt(n),
      median = median(.data[[comp]], na.rm = TRUE),
      min = min(.data[[comp]], na.rm = TRUE),
      max = max(.data[[comp]], na.rm = TRUE),
      .groups = "drop"
    ) %>% mutate(component = comp)
  descriptives_age_sex <- rbind(descriptives_age_sex, comp_desc)
}

descriptives_by_subtype <- W_df %>%
  group_by(dominant_subtype, AgeGroup, Sex) %>%
  summarise(
    n = n(),
    mean_dominant_proportion = mean(dominant_proportion, na.rm = TRUE),
    mean_subtype_margin = mean(subtype_margin, na.rm = TRUE),
    .groups = "drop"
  )

age_sex_desc_file <- file.path(out_dir, "NMF_DD_only_descriptive_statistics_age_sex.csv")
subtype_desc_file <- file.path(out_dir, "NMF_DD_only_descriptive_statistics_by_subtype.csv")
write.csv(descriptives_age_sex, age_sex_desc_file, row.names = FALSE)
write.csv(descriptives_by_subtype, subtype_desc_file, row.names = FALSE)
write.csv(descriptives_age_sex, file.path(out_dir, "NMF_descriptive_statistics.csv"), row.names = FALSE)
cat("Saved:", age_sex_desc_file, "\n")
cat("Saved:", subtype_desc_file, "\n")

sig_posthoc <- posthoc_results %>% filter(q_value < 0.05) %>% arrange(q_value)
cat("\n==============================\n")
cat("DD-ONLY POST-HOC SUMMARY\n")
cat("==============================\n\n")
cat("Significant AgeGroup/Sex contrasts, q < 0.05:\n")
if (nrow(sig_posthoc) > 0) {
  print(sig_posthoc %>% select(component, term, contrast, estimate, cohens_d, p.value, q_value, effect_size_interpretation), row.names = FALSE)
} else {
  cat("  No significant post-hoc contrasts at FDR < 0.05\n")
}
cat("\nDONE\n")
