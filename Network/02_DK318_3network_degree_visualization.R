# ============================================================
# DK318 MIND Degree 3-Network Visualization
# Purpose: Generate figures from analysis results
# ============================================================

cat("Starting DK318 3-network degree visualization...\n\n")

# ============================================================
# 1. Load Packages
# ============================================================

cat("Loading packages...\n")
pkgs <- c("dplyr", "readr", "ggplot2", "tidyr", "magrittr")
for (p in pkgs) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) {
    install.packages(p, repos = "https://cran.r-project.org")
    library(p, character.only = TRUE)
  }
}
library(dplyr); library(readr); library(ggplot2); library(tidyr); library(magrittr)

# ============================================================
# 2. Define Paths
# ============================================================

out_root <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_3network_analysis"
analysis_dir <- file.path(out_root, "analysis")
fig_dir <- file.path(out_root, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 3. Helper Functions
# ============================================================

save2 <- function(p, n, w = 8, h = 6) {
  ggsave(file.path(fig_dir, paste0(n, ".png")), p, width = w, height = h, dpi = 300)
}

star <- function(q) {
  case_when(
    is.na(q) ~ "",
    q < 0.001 ~ "***",
    q < 0.01 ~ "**",
    q < 0.05 ~ "*",
    TRUE ~ ""
  )
}

check_file <- function(fp, desc) {
  if (!file.exists(fp)) {
    warning("Required file not found: ", desc, " (", fp, ")")
    return(FALSE)
  }
  return(TRUE)
}

# ============================================================
# 4. Check Required Input Files
# ============================================================

cat("Checking required input files...\n")
required_files <- list(
  membership = file.path(analysis_dir, "network_membership_3net.csv"),
  mixed_emm = file.path(analysis_dir, "mixed_model_3network_emmeans_adjusted_means.csv"),
  long_data = file.path(analysis_dir, "DK318_3network_degree_long.csv"),
  mixed_td_dd = file.path(analysis_dir, "mixed_model_3network_TD_minus_DD_by_network_age.csv"),
  mixed_dev = file.path(analysis_dir, "mixed_model_3network_developmental_contrast.csv"),
  per_td_dd = file.path(analysis_dir, "per_network_TD_minus_DD_by_age.csv"),
  per_dev = file.path(analysis_dir, "per_network_developmental_contrast.csv")
)

all_exist <- TRUE
for (name in names(required_files)) {
  if (!check_file(required_files[[name]], name)) {
    all_exist <- FALSE
  }
}

if (!all_exist) {
  stop("Some required input files are missing. Please run analysis script first.")
}

cat("All required files found.\n\n")

# ============================================================
# 5. Figure 1: Network Overlap Membership Bar
# ============================================================

cat("Generating Figure 1: Network overlap membership bar...\n")
mem <- read_csv(required_files$membership, show_col_types = FALSE)
p1 <- mem %>%
  count(membership_pattern, name = "n") %>%
  ggplot(aes(reorder(membership_pattern, n), n, fill = membership_pattern)) +
  geom_col() +
  coord_flip() +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Membership pattern", y = "ROI count", title = "3-network overlap membership")
save2(p1, "fig_01_network_overlap_membership_bar", 8, 5)

# ============================================================
# 6. Figure 2: Mixed Model Adjusted Means
# ============================================================

cat("Generating Figure 2: Mixed model adjusted means...\n")
emm <- read_csv(required_files$mixed_emm, show_col_types = FALSE) %>%
  mutate(
    Network = factor(Network, levels = c("DM", "MD", "Reading")),
    AgeGroup = factor(AgeGroup, levels = c("Child", "Adult")),
    Diagnosis = factor(Diagnosis, levels = c("TD", "DD"))
  )

p2 <- ggplot(emm, aes(AgeGroup, emmean, color = Diagnosis, group = Diagnosis)) +
  geom_point(position = position_dodge(0.2), size = 2) +
  geom_line(position = position_dodge(0.2)) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, position = position_dodge(0.2)) +
  facet_wrap(~Network, nrow = 1) +
  theme_bw() +
  labs(y = "Adjusted mean degree", title = "Mixed-model adjusted means by network")
save2(p2, "fig_02_mixed_model_adjusted_means_by_network", 10, 4)

# ============================================================
# 7. Figure 3: Raw Network Degree Distributions
# ============================================================

cat("Generating Figure 3: Raw network degree distributions...\n")
raw <- read_csv(required_files$long_data, show_col_types = FALSE) %>%
  mutate(
    Network = factor(Network, levels = c("DM", "MD", "Reading")),
    AgeGroup = factor(AgeGroup, levels = c("Child", "Adult")),
    Diagnosis = factor(Diagnosis, levels = c("TD", "DD"))
  )

p3 <- ggplot(raw, aes(AgeGroup, Degree, fill = Diagnosis, color = Diagnosis)) +
  geom_violin(alpha = 0.25, position = position_dodge(0.8), trim = FALSE) +
  geom_boxplot(width = 0.18, position = position_dodge(0.8), outlier.shape = NA, alpha = 0.5) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.8), alpha = 0.5, size = 1) +
  facet_wrap(~Network, nrow = 1, scales = "free_y") +
  theme_bw() +
  labs(title = "Raw network degree distributions")
save2(p3, "fig_03_raw_network_degree_distribution", 11, 4.5)

# ============================================================
# 8. Figure 4: Mixed Model TD-DD t Heatmap
# ============================================================

cat("Generating Figure 4: Mixed model TD-DD t heatmap...\n")
c_age <- read_csv(required_files$mixed_td_dd, show_col_types = FALSE) %>%
  mutate(column = paste(AgeGroup, "TD-DD"), label = paste0(round(t.ratio, 2), star(q_all)))

c_dev <- read_csv(required_files$mixed_dev, show_col_types = FALSE) %>%
  mutate(column = "Adult(TD-DD)-Child(TD-DD)", label = paste0(round(t.ratio, 2), star(q_value)))

ht <- bind_rows(
  c_age %>% transmute(Network, column, t.ratio, label),
  c_dev %>% transmute(Network, column, t.ratio, label)
) %>%
  mutate(
    Network = factor(Network, levels = c("DM", "MD", "Reading")),
    column = factor(column, levels = c("Child TD-DD", "Adult TD-DD", "Adult(TD-DD)-Child(TD-DD)"))
  )

p4 <- ggplot(ht, aes(column, Network, fill = t.ratio)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient2(low = "#3b4cc0", mid = "white", high = "#b40426", midpoint = 0) +
  theme_bw() +
  labs(x = NULL, y = NULL, title = "Mixed-model TD-DD t heatmap")
save2(p4, "fig_04_mixed_model_TD_minus_DD_t_heatmap", 9, 4)

# ============================================================
# 9. Figure 5: Mixed Model TD-DD Contrast Profile
# ============================================================

cat("Generating Figure 5: Mixed model TD-DD contrast profile...\n")
p5 <- c_age %>%
  mutate(
    AgeGroup = factor(AgeGroup, levels = c("Child", "Adult")),
    Network = factor(Network, levels = c("DM", "MD", "Reading"))
  ) %>%
  ggplot(aes(AgeGroup, estimate, color = Network, group = Network)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE), width = 0.08) +
  theme_bw() +
  labs(y = "TD - DD estimate", title = "Mixed-model TD-DD contrast profile")
save2(p5, "fig_05_mixed_model_TD_minus_DD_contrast_profile", 7, 5)

# ============================================================
# 10. Figure 6: Per-Network TD-DD t Heatmap
# ============================================================

cat("Generating Figure 6: Per-network TD-DD t heatmap...\n")
p_age <- read_csv(required_files$per_td_dd, show_col_types = FALSE) %>%
  mutate(column = paste(AgeGroup, "TD-DD"), label = paste0(round(t.ratio, 2), star(q_all)))

p_dev <- read_csv(required_files$per_dev, show_col_types = FALSE) %>%
  mutate(column = "Adult(TD-DD)-Child(TD-DD)", label = paste0(round(t.ratio, 2), star(q_value)))

ht2 <- bind_rows(
  p_age %>% transmute(Network, column, t.ratio, label),
  p_dev %>% transmute(Network, column, t.ratio, label)
) %>%
  mutate(
    Network = factor(Network, levels = c("DM", "MD", "Reading")),
    column = factor(column, levels = c("Child TD-DD", "Adult TD-DD", "Adult(TD-DD)-Child(TD-DD)"))
  )

p6 <- ggplot(ht2, aes(column, Network, fill = t.ratio)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient2(low = "#3b4cc0", mid = "white", high = "#b40426", midpoint = 0) +
  theme_bw() +
  labs(x = NULL, y = NULL, title = "Per-network TD-DD t heatmap")
save2(p6, "fig_06_per_network_TD_minus_DD_t_heatmap", 9, 4)

# ============================================================
# PART 2: PER-NETWORK MODEL RESULTS
# ============================================================

# ============================================================
# 7. Figure 7: Per-Network Model Adjusted Means
# ============================================================

cat("Generating Figure 7: Per-network model adjusted means...\n")
per_emm <- read_csv(file.path(analysis_dir, "per_network_emmeans_adjusted_means.csv"), show_col_types = FALSE) %>%
  mutate(
    Network = factor(Network, levels = c("DM", "MD", "Reading")),
    AgeGroup = factor(AgeGroup, levels = c("Child", "Adult")),
    Diagnosis = factor(Diagnosis, levels = c("TD", "DD"))
  )

p7 <- ggplot(per_emm, aes(AgeGroup, emmean, color = Diagnosis, group = Diagnosis)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  geom_line(position = position_dodge(0.2), size = 1) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, position = position_dodge(0.2)) +
  facet_wrap(~Network, nrow = 1, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "top") +
  labs(
    y = "Adjusted mean degree",
    title = "Per-network model: Adjusted means",
    subtitle = "Independent linear model for each network"
  )
save2(p7, "fig_07_per_network_adjusted_means", 11, 5)

# ============================================================
# 8. Figure 8: Per-Network TD-DD t Heatmap (moved from Fig 6)
# ============================================================

cat("Generating Figure 8: Per-network TD-DD t heatmap...\n")
p_age <- read_csv(required_files$per_td_dd, show_col_types = FALSE) %>%
  mutate(column = paste(AgeGroup, "TD-DD"), label = paste0(round(t.ratio, 2), star(q_all)))

p_dev <- read_csv(required_files$per_dev, show_col_types = FALSE) %>%
  mutate(column = "Adult(TD-DD)-Child(TD-DD)", label = paste0(round(t.ratio, 2), star(q_value)))

ht2 <- bind_rows(
  p_age %>% transmute(Network, column, t.ratio, label),
  p_dev %>% transmute(Network, column, t.ratio, label)
) %>%
  mutate(
    Network = factor(Network, levels = c("DM", "MD", "Reading")),
    column = factor(column, levels = c("Child TD-DD", "Adult TD-DD", "Adult(TD-DD)-Child(TD-DD)"))
  )

p8 <- ggplot(ht2, aes(column, Network, fill = t.ratio)) +
  geom_tile(color = "white") +
  geom_text(aes(label = label), size = 4) +
  scale_fill_gradient2(low = "#3b4cc0", mid = "white", high = "#b40426", midpoint = 0) +
  theme_bw() +
  labs(x = NULL, y = NULL, title = "Per-network model: TD-DD t heatmap")
save2(p8, "fig_08_per_network_TD_minus_DD_t_heatmap", 9, 4)

# ============================================================
# 9. Figure 9: Per-Network TD-DD Contrast Profile
# ============================================================

cat("Generating Figure 9: Per-network TD-DD contrast profile...\n")
p9 <- p_age %>%
  mutate(
    AgeGroup = factor(AgeGroup, levels = c("Child", "Adult")),
    Network = factor(Network, levels = c("DM", "MD", "Reading"))
  ) %>%
  ggplot(aes(AgeGroup, estimate, color = Network, group = Network)) +
  geom_hline(yintercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 2) +
  geom_line() +
  geom_errorbar(aes(ymin = estimate - SE, ymax = estimate + SE), width = 0.08) +
  theme_bw() +
  labs(y = "TD - DD estimate", title = "Per-network model: TD-DD contrast profile")
save2(p9, "fig_09_per_network_TD_minus_DD_profile", 7, 5)

# ============================================================
# 10. Figure 10: Per-Network Diagnosis×AgeGroup Interaction
# ============================================================

cat("Generating Figure 10: Per-network Diagnosis×AgeGroup interaction...\n")
p10 <- ggplot(per_emm, aes(AgeGroup, emmean, color = Diagnosis, group = Diagnosis)) +
  geom_point(size = 3, position = position_dodge(0.2)) +
  geom_line(position = position_dodge(0.2), size = 1) +
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, position = position_dodge(0.2)) +
  facet_wrap(~Network, nrow = 1, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "top") +
  labs(
    y = "Adjusted mean degree",
    title = "Per-network model: Diagnosis × AgeGroup interaction",
    subtitle = "Parallel lines indicate no interaction"
  )
save2(p10, "fig_10_per_network_diagnosis_agegroup_interaction", 11, 5)

# ============================================================
# PART 3: MODEL COMPARISON
# ============================================================

# ============================================================
# 11. Figure 11: Mixed vs Per-Network Model Comparison
# ============================================================

cat("Generating Figure 11: Model comparison scatter plot...\n")
mixed_comp <- read_csv(required_files$mixed_td_dd, show_col_types = FALSE) %>%
  transmute(Network, AgeGroup, mixed_t = t.ratio, mixed_p = p.value, mixed_est = estimate)

per_comp <- read_csv(required_files$per_td_dd, show_col_types = FALSE) %>%
  transmute(Network, AgeGroup, per_t = t.ratio, per_p = p.value, per_est = estimate)

comp <- left_join(mixed_comp, per_comp, by = c("Network", "AgeGroup")) %>%
  mutate(
    Network = factor(Network, levels = c("DM", "MD", "Reading")),
    AgeGroup = factor(AgeGroup, levels = c("Child", "Adult")),
    label = paste0(Network, "\n", AgeGroup)
  )

p7 <- ggplot(comp, aes(mixed_t, per_t, color = Network, shape = AgeGroup)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey50") +
  geom_point(size = 4, alpha = 0.8) +
  geom_text(aes(label = label), vjust = -0.8, size = 3, show.legend = FALSE) +
  theme_bw() +
  labs(
    x = "Mixed-effects model t-value",
    y = "Per-network model t-value",
    title = "Model comparison: TD-DD contrasts",
    subtitle = "Points near diagonal indicate high consistency"
  ) +
  coord_fixed()
save2(p11, "fig_11_model_comparison_scatter", 7, 7)

# ============================================================
# 12. Figure 12: Forest Plot for All TD-DD Contrasts (Both Models)
# ============================================================

cat("Generating Figure 12: Forest plot for TD-DD contrasts (both models)...\n")
forest_mixed <- read_csv(required_files$mixed_td_dd, show_col_types = FALSE) %>%
  mutate(
    Network = factor(Network, levels = c("Reading", "MD", "DM")),
    AgeGroup = factor(AgeGroup, levels = c("Child", "Adult")),
    group_label = paste(Network, AgeGroup, sep = " - "),
    ci_lower = estimate - 1.96 * SE,
    ci_upper = estimate + 1.96 * SE,
    sig = ifelse(q_all < 0.05, "Significant", "Non-significant"),
    model = "Mixed-effects"
  )

forest_per <- read_csv(required_files$per_td_dd, show_col_types = FALSE) %>%
  mutate(
    Network = factor(Network, levels = c("Reading", "MD", "DM")),
    AgeGroup = factor(AgeGroup, levels = c("Child", "Adult")),
    group_label = paste(Network, AgeGroup, sep = " - "),
    ci_lower = estimate - 1.96 * SE,
    ci_upper = estimate + 1.96 * SE,
    sig = ifelse(q_all < 0.05, "Significant", "Non-significant"),
    model = "Per-network"
  )

forest_combined <- bind_rows(forest_mixed, forest_per) %>%
  mutate(model = factor(model, levels = c("Mixed-effects", "Per-network")))

p12 <- ggplot(forest_combined, aes(estimate, group_label, color = model, shape = model)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
  geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.2, size = 0.8, 
                 position = position_dodge(height = 0.5)) +
  geom_point(size = 3, position = position_dodge(height = 0.5)) +
  scale_color_manual(values = c("Mixed-effects" = "#2166ac", "Per-network" = "#b2182b")) +
  theme_bw() +
  theme(legend.position = "top") +
  labs(
    x = "Effect size (TD - DD)",
    y = NULL,
    color = "Model",
    shape = "Model",
    title = "Model comparison: Forest plot of TD-DD contrasts",
    subtitle = "Error bars show 95% confidence intervals"
  )
save2(p12, "fig_12_forest_plot_both_models", 9, 6)

# ============================================================
# 13. Figure 13: Effect Size Comparison
# ============================================================

cat("Generating Figure 13: Effect size comparison...\n")
cat("Generating Figure 13: Effect size comparison...\n")
comp <- left_join(mixed_comp, per_comp, by = c("Network", "AgeGroup")) %>%
  mutate(
    Network = factor(Network, levels = c("DM", "MD", "Reading")),
    AgeGroup = factor(AgeGroup, levels = c("Child", "Adult")),
    label = paste0(Network, "\n", AgeGroup)
  )

effect_comp <- comp %>%
  select(Network, AgeGroup, mixed_est, per_est) %>%
  pivot_longer(c(mixed_est, per_est), names_to = "model", values_to = "effect_size") %>%
  mutate(
    model = recode(model, mixed_est = "Mixed-effects", per_est = "Per-network"),
    model = factor(model, levels = c("Mixed-effects", "Per-network")),
    group = paste(Network, AgeGroup, sep = "\n")
  )

p13 <- ggplot(effect_comp, aes(group, effect_size, fill = model)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_fill_manual(values = c("Mixed-effects" = "#2166ac", "Per-network" = "#b2182b")) +
  theme_bw() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 9)
  ) +
  labs(
    x = NULL,
    y = "Effect size (TD - DD)",
    fill = "Model type",
    title = "Model comparison: Effect size comparison"
  )
save2(p13, "fig_13_effect_size_comparison", 9, 5)

# ============================================================
# 14. Figure 14: Combined Heatmap (moved from Fig 11)
# ============================================================

cat("Generating Figure 14: Combined model heatmap...\n")
ht_mixed <- ht %>%
  filter(column != "Adult(TD-DD)-Child(TD-DD)") %>%
  mutate(model = "Mixed-effects")

ht_per <- ht2 %>%
  filter(column != "Adult(TD-DD)-Child(TD-DD)") %>%
  mutate(model = "Per-network")

ht_combined <- bind_rows(ht_mixed, ht_per) %>%
  mutate(
    model = factor(model, levels = c("Mixed-effects", "Per-network")),
    Network = factor(Network, levels = c("DM", "MD", "Reading")),
    column = factor(column, levels = c("Child TD-DD", "Adult TD-DD"))
  )

p14 <- ggplot(ht_combined, aes(column, Network, fill = t.ratio)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = label), size = 3.5) +
  facet_wrap(~model, nrow = 1) +
  scale_fill_gradient2(low = "#3b4cc0", mid = "white", high = "#b40426", midpoint = 0, 
                       name = "t-value") +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  ) +
  labs(
    x = NULL, 
    y = NULL, 
    title = "Model comparison: TD-DD t-values",
    subtitle = "High consistency between models"
  )
save2(p14, "fig_14_combined_model_heatmap", 10, 4.5)

# ============================================================
# 15. Complete
# ============================================================

cat("\n==============================\n")
cat("VISUALIZATION COMPLETE\n")
cat("==============================\n")
cat("Total figures generated: 14\n")
cat("  - Mixed-effects model: Fig 01-06\n")
cat("  - Per-network model: Fig 07-10\n")
cat("  - Model comparison: Fig 11-14\n")
cat("All figures saved to:", fig_dir, "\n")
