# ============================================================
# DK318 MIND Degree NMF - Visualization Script
# Purpose: Create detailed visualizations for NMF results
# Input: NMF results from 01_nmf_mind_degree.R
# Output: Enhanced plots and heatmaps
# ============================================================

cat("Starting NMF Visualization...\n\n")

# ============================================================
# 1. Load Packages
# ============================================================

packages <- c("dplyr", "ggplot2", "reshape2", "RColorBrewer", "pheatmap", "gridExtra")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

library(dplyr); library(ggplot2); library(reshape2)
library(RColorBrewer); library(pheatmap); library(gridExtra)

# ============================================================
# 2. Define Paths
# ============================================================

nmf_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NMF_3"
out_dir <- file.path(nmf_dir, "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("NMF directory:", nmf_dir, "\n")
cat("Output directory:", out_dir, "\n\n")

# ============================================================
# 3. Load Data
# ============================================================

cat("Loading NMF results...\n")

# Load W matrix (subject loadings)
W_df <- read.csv(file.path(nmf_dir, "NMF_W_subject_loadings.csv"), stringsAsFactors = FALSE)

# Load H matrix (component ROI weights)
H_long <- read.csv(file.path(nmf_dir, "NMF_H_component_roi_weights_long.csv"), stringsAsFactors = FALSE)

# Load ANOVA results
anova_results <- read.csv(file.path(nmf_dir, "NMF_component_lm_anova_results.csv"), stringsAsFactors = FALSE)

# Load top ROIs
top_rois <- read.csv(file.path(nmf_dir, "NMF_component_top_ROIs.csv"), stringsAsFactors = FALSE)

# Get component names
component_cols <- grep("^Component_", colnames(W_df), value = TRUE)
n_components <- length(component_cols)

cat("Loaded data for", n_components, "components\n")
cat("Number of subjects:", nrow(W_df), "\n\n")

# ============================================================
# 4. Component Scores - Individual Plots
# ============================================================

cat("Creating individual component score plots...\n")

for (comp in component_cols) {
  # Get ANOVA results for this component
  comp_anova <- anova_results %>% filter(component == comp)
  
  # Extract p-values for annotation
  p_diag <- comp_anova$p_value[comp_anova$term == "Diagnosis"]
  p_age <- comp_anova$p_value[comp_anova$term == "AgeGroup"]
  p_int <- comp_anova$p_value[comp_anova$term == "AgeGroup:Diagnosis"]
  
  q_diag <- comp_anova$q_value_by_term[comp_anova$term == "Diagnosis"]
  q_age <- comp_anova$q_value_by_term[comp_anova$term == "AgeGroup"]
  q_int <- comp_anova$q_value_by_term[comp_anova$term == "AgeGroup:Diagnosis"]
  
  # Create annotation text
  anno_text <- sprintf(
    "Diagnosis: p=%.3f, q=%.3f\nAgeGroup: p=%.3f, q=%.3f\nInteraction: p=%.3f, q=%.3f",
    p_diag, q_diag, p_age, q_age, p_int, q_int
  )
  
  # Create plot
  p <- ggplot(W_df, aes(x = AgeGroup, y = .data[[comp]], fill = Diagnosis)) +
    geom_boxplot(outlier.size = 1, alpha = 0.7) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2), 
               alpha = 0.3, size = 0.8) +
    scale_fill_manual(values = c("TD" = "#4DAF4A", "DD" = "#E41A1C")) +
    labs(title = paste(comp, "Scores by Age Group and Diagnosis"),
         subtitle = anno_text,
         x = "Age Group", y = "Component Score") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom",
          plot.subtitle = element_text(size = 9, face = "italic"))
  
  ggsave(file.path(out_dir, paste0(comp, "_boxplot.png")), 
         p, width = 8, height = 6, dpi = 300)
}

cat("Saved", n_components, "individual component plots\n\n")

# ============================================================
# 5. Component Scores Heatmap
# ============================================================

cat("Creating component scores heatmap...\n")

# Prepare data for heatmap
W_matrix <- as.matrix(W_df[, component_cols])
rownames(W_matrix) <- W_df$file_base

# Create annotation for rows
row_annotation <- data.frame(
  Diagnosis = W_df$Diagnosis,
  AgeGroup = W_df$AgeGroup,
  Sex = W_df$Sex,
  row.names = W_df$file_base
)

# Define colors
anno_colors <- list(
  Diagnosis = c("TD" = "#4DAF4A", "DD" = "#E41A1C"),
  AgeGroup = c("Child" = "#377EB8", "Adult" = "#FF7F00"),
  Sex = c("Female" = "#984EA3", "Male" = "#FFFF33")
)

# Create heatmap
png(file.path(out_dir, "component_scores_heatmap.png"), 
    width = 10, height = 12, units = "in", res = 300)

pheatmap(t(W_matrix),
         annotation_col = row_annotation,
         annotation_colors = anno_colors,
         scale = "row",
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = FALSE,
         main = "NMF Component Scores Across Subjects",
         fontsize = 10,
         color = colorRampPalette(c("blue", "white", "red"))(100))

dev.off()

cat("Saved component scores heatmap\n\n")

# ============================================================
# 6. Top ROI Heatmap for Each Component
# ============================================================

cat("Creating top ROI heatmaps...\n")

# Get top 20 ROIs for each component
top20_rois <- top_rois %>%
  filter(top_type == "top20") %>%
  select(component, ROI, weight, rank_within_component)

for (comp in unique(top20_rois$component)) {
  comp_top <- top20_rois %>% 
    filter(component == comp) %>%
    arrange(rank_within_component)
  
  # Create barplot
  p <- ggplot(comp_top, aes(x = reorder(ROI, weight), y = weight)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = paste(comp, "- Top 20 ROIs"),
         x = "ROI", y = "Weight") +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 8))
  
  ggsave(file.path(out_dir, paste0(comp, "_top20_rois.png")), 
         p, width = 8, height = 10, dpi = 300)
}

cat("Saved", length(unique(top20_rois$component)), "top ROI plots\n\n")

# ============================================================
# 7. ANOVA Results Heatmap
# ============================================================

cat("Creating ANOVA results heatmap...\n")

# Prepare data: -log10(p-value) for visualization
anova_wide <- anova_results %>%
  select(component, term, p_value) %>%
  mutate(neg_log10_p = -log10(p_value)) %>%
  select(-p_value) %>%
  tidyr::pivot_wider(names_from = term, values_from = neg_log10_p)

anova_matrix <- as.matrix(anova_wide[, -1])
rownames(anova_matrix) <- anova_wide$component

# Create heatmap
png(file.path(out_dir, "anova_results_heatmap.png"), 
    width = 8, height = 6, units = "in", res = 300)

pheatmap(anova_matrix,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         display_numbers = TRUE,
         number_format = "%.2f",
         main = "ANOVA Results: -log10(p-value)",
         fontsize = 12,
         color = colorRampPalette(c("white", "yellow", "red"))(100),
         breaks = seq(0, max(anova_matrix, na.rm = TRUE), length.out = 101))

dev.off()

cat("Saved ANOVA results heatmap\n\n")

# ============================================================
# 8. Component Correlation Plot
# ============================================================

cat("Creating component correlation plot...\n")

# Calculate correlation between components
W_cor <- cor(W_matrix)

# Create correlation heatmap
png(file.path(out_dir, "component_correlation.png"), 
    width = 8, height = 7, units = "in", res = 300)

pheatmap(W_cor,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         display_numbers = TRUE,
         number_format = "%.2f",
         main = "Correlation Between NMF Components",
         fontsize = 12,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = seq(-1, 1, length.out = 101))

dev.off()

cat("Saved component correlation plot\n\n")

# ============================================================
# 9. Violin Plots for Significant Components
# ============================================================

cat("Creating violin plots for significant components...\n")

# Find significant interactions
sig_interactions <- anova_results %>%
  filter(term == "AgeGroup:Diagnosis", q_value_by_term < 0.05) %>%
  pull(component)

if (length(sig_interactions) > 0) {
  for (comp in sig_interactions) {
    p <- ggplot(W_df, aes(x = AgeGroup, y = .data[[comp]], fill = Diagnosis)) +
      geom_violin(alpha = 0.5, trim = FALSE) +
      geom_boxplot(width = 0.2, position = position_dodge(0.9), outlier.size = 0.5) +
      scale_fill_manual(values = c("TD" = "#4DAF4A", "DD" = "#E41A1C")) +
      labs(title = paste(comp, "- Significant AgeGroup × Diagnosis Interaction"),
           x = "Age Group", y = "Component Score") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom")
    
    ggsave(file.path(out_dir, paste0(comp, "_violin_significant.png")), 
           p, width = 8, height = 6, dpi = 300)
  }
  
  cat("Saved", length(sig_interactions), "violin plots for significant components\n\n")
} else {
  cat("No significant interactions found (q < 0.05)\n\n")
}

# ============================================================
# 10. Summary Plot: All Components Overview
# ============================================================

cat("Creating summary overview plot...\n")

# Prepare data for faceted plot
W_long <- W_df %>%
  select(file_base, Diagnosis, AgeGroup, all_of(component_cols)) %>%
  reshape2::melt(id.vars = c("file_base", "Diagnosis", "AgeGroup"),
                 variable.name = "component", value.name = "score")

# Add significance stars
sig_labels <- anova_results %>%
  filter(term == "AgeGroup:Diagnosis") %>%
  mutate(sig_label = case_when(
    q_value_by_term < 0.001 ~ "***",
    q_value_by_term < 0.01 ~ "**",
    q_value_by_term < 0.05 ~ "*",
    TRUE ~ ""
  )) %>%
  select(component, sig_label)

# Create faceted plot
p_summary <- ggplot(W_long, aes(x = AgeGroup, y = score, fill = Diagnosis)) +
  geom_boxplot(outlier.size = 0.3) +
  facet_wrap(~ component, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("TD" = "#4DAF4A", "DD" = "#E41A1C")) +
  labs(title = "NMF Component Scores: Overview",
       subtitle = "Significance: *** q<0.001, ** q<0.01, * q<0.05 (AgeGroup:Diagnosis)",
       x = "Age Group", y = "Component Score") +
  theme_bw(base_size = 11) +
  theme(legend.position = "bottom",
        strip.background = element_rect(fill = "lightblue"))

ggsave(file.path(out_dir, "all_components_overview.png"), 
       p_summary, width = 12, height = 10, dpi = 300)

cat("Saved summary overview plot\n\n")

# ============================================================
# Summary
# ============================================================

cat("==============================\n")
cat("VISUALIZATION COMPLETE\n")
cat("==============================\n\n")

cat("Output directory:", out_dir, "\n\n")

cat("Generated plots:\n")
cat("  - Individual component boxplots:", n_components, "files\n")
cat("  - Component scores heatmap: 1 file\n")
cat("  - Top ROI plots:", n_components, "files\n")
cat("  - ANOVA results heatmap: 1 file\n")
cat("  - Component correlation plot: 1 file\n")
cat("  - Violin plots for significant components:", length(sig_interactions), "files\n")
cat("  - Summary overview plot: 1 file\n\n")

cat("==============================\n")
cat("DONE\n")
cat("==============================\n")
