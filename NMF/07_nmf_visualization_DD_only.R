# ============================================================
# DK318 MIND Degree NMF - DD-only Visualization
# Purpose: Create detailed plots for DD-only NMF subtype results.
# ============================================================

cat("Starting DD-only NMF Visualization...\n\n")

packages <- c("dplyr", "ggplot2", "reshape2", "pheatmap", "tidyr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}
library(dplyr); library(ggplot2); library(reshape2); library(pheatmap); library(tidyr)

nmf_dir <- Sys.getenv("DD_NMF_OUT_DIR", "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NMF_DD_only")
out_dir <- file.path(nmf_dir, "plots")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("DD-only NMF directory:", nmf_dir, "\n")
cat("Plot directory:", out_dir, "\n\n")

W_df <- read.csv(file.path(nmf_dir, "NMF_DD_only_W_subject_loadings_and_subtypes.csv"), stringsAsFactors = FALSE)
H_long <- read.csv(file.path(nmf_dir, "NMF_DD_only_H_component_roi_weights_long.csv"), stringsAsFactors = FALSE)
anova_results <- read.csv(file.path(nmf_dir, "NMF_DD_only_component_lm_anova_results.csv"), stringsAsFactors = FALSE)
top_rois <- read.csv(file.path(nmf_dir, "NMF_DD_only_component_top_ROIs.csv"), stringsAsFactors = FALSE)
subtype_summary_file <- file.path(nmf_dir, "NMF_DD_only_subtype_summary.csv")
subtype_summary <- if (file.exists(subtype_summary_file)) read.csv(subtype_summary_file, stringsAsFactors = FALSE) else NULL

component_cols <- grep("^SubtypeComponent_[0-9]+$", colnames(W_df), value = TRUE)
if (length(component_cols) == 0) component_cols <- grep("^Component_[0-9]+$", colnames(W_df), value = TRUE)
if (length(component_cols) == 0) stop("No component columns found in W file.")

W_df$AgeGroup <- factor(W_df$AgeGroup, levels = c("Child", "Adult"))
W_df$Sex <- factor(W_df$Sex, levels = c("Female", "Male"))
W_df$dominant_subtype <- factor(W_df$dominant_subtype, levels = component_cols)

cat("Loaded", nrow(W_df), "DD subjects and", length(component_cols), "components.\n\n")

get_pq <- function(comp, term) {
  x <- anova_results %>% filter(component == comp, term == !!term)
  if (nrow(x) == 0) return(c(NA_real_, NA_real_))
  c(x$p_value[1], x$q_value_by_term[1])
}

cat("Creating individual component score plots...\n")
for (comp in component_cols) {
  age_pq <- get_pq(comp, "AgeGroup")
  sex_pq <- get_pq(comp, "Sex")
  anno_text <- sprintf("AgeGroup: p=%.3f, q=%.3f\nSex: p=%.3f, q=%.3f", age_pq[1], age_pq[2], sex_pq[1], sex_pq[2])
  p <- ggplot(W_df, aes(x = AgeGroup, y = .data[[comp]], fill = Sex)) +
    geom_boxplot(outlier.size = 1, alpha = 0.75) +
    geom_point(position = position_jitterdodge(jitter.width = 0.2), alpha = 0.45, size = 1.0) +
    scale_fill_manual(values = c("Female" = "#984EA3", "Male" = "#4DAF4A")) +
    labs(title = paste(comp, "scores within DD"), subtitle = anno_text, x = "Age group", y = "Component score") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 9, face = "italic"))
  ggsave(file.path(out_dir, paste0(comp, "_DD_only_boxplot.png")), p, width = 8, height = 6, dpi = 300)
}

cat("Creating component scores heatmap...\n")
W_matrix <- as.matrix(W_df[, component_cols])
rownames(W_matrix) <- W_df$file_base
row_annotation <- data.frame(
  AgeGroup = W_df$AgeGroup,
  Sex = W_df$Sex,
  DominantSubtype = W_df$dominant_subtype,
  row.names = W_df$file_base
)
anno_colors <- list(
  AgeGroup = c("Child" = "#377EB8", "Adult" = "#FF7F00"),
  Sex = c("Female" = "#984EA3", "Male" = "#4DAF4A")
)
png(file.path(out_dir, "DD_only_component_scores_heatmap.png"), width = 10, height = 12, units = "in", res = 300)
pheatmap(t(W_matrix), annotation_col = row_annotation, annotation_colors = anno_colors,
         scale = "row", cluster_rows = TRUE, cluster_cols = TRUE, show_colnames = FALSE,
         main = "DD-only NMF Component Scores", fontsize = 10,
         color = colorRampPalette(c("blue", "white", "red"))(100))
dev.off()

cat("Creating top ROI plots...\n")
top20_rois <- top_rois %>% filter(top_type == "top20") %>% select(component, ROI, weight, rank_within_component)
for (comp in unique(top20_rois$component)) {
  comp_top <- top20_rois %>% filter(component == comp) %>% arrange(rank_within_component)
  p <- ggplot(comp_top, aes(x = reorder(ROI, weight), y = weight)) +
    geom_bar(stat = "identity", fill = "steelblue") + coord_flip() +
    labs(title = paste(comp, "- DD-only Top 20 ROIs"), x = "ROI", y = "Weight") +
    theme_bw(base_size = 10) + theme(axis.text.y = element_text(size = 8))
  ggsave(file.path(out_dir, paste0(comp, "_DD_only_top20_rois.png")), p, width = 8, height = 10, dpi = 300)
}

cat("Creating ANOVA heatmap...\n")
anova_wide <- anova_results %>%
  select(component, term, p_value) %>% mutate(neg_log10_p = -log10(pmax(p_value, .Machine$double.xmin))) %>%
  select(-p_value) %>% tidyr::pivot_wider(names_from = term, values_from = neg_log10_p)
anova_matrix <- as.matrix(anova_wide[, -1, drop = FALSE])
rownames(anova_matrix) <- anova_wide$component
png(file.path(out_dir, "DD_only_anova_results_heatmap.png"), width = 8, height = 6, units = "in", res = 300)
pheatmap(anova_matrix, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE,
         number_format = "%.2f", main = "DD-only ANOVA: -log10(p-value)", fontsize = 12,
         color = colorRampPalette(c("white", "yellow", "red"))(100))
dev.off()

cat("Creating component correlation plot...\n")
png(file.path(out_dir, "DD_only_component_correlation.png"), width = 8, height = 7, units = "in", res = 300)
pheatmap(cor(W_matrix), cluster_rows = TRUE, cluster_cols = TRUE, display_numbers = TRUE,
         number_format = "%.2f", main = "DD-only Component Correlations", fontsize = 12,
         color = colorRampPalette(c("blue", "white", "red"))(100), breaks = seq(-1, 1, length.out = 101))
dev.off()

cat("Creating subtype distribution plots...\n")
p_counts <- ggplot(W_df, aes(x = dominant_subtype, fill = AgeGroup)) +
  geom_bar(position = "dodge") + scale_fill_manual(values = c("Child" = "#377EB8", "Adult" = "#FF7F00")) +
  labs(title = "DD dominant NMF subtype counts", x = "Dominant subtype", y = "Number of DD subjects") +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "bottom")
ggsave(file.path(out_dir, "DD_only_dominant_subtype_counts_by_age.png"), p_counts, width = 9, height = 6, dpi = 300)

p_prop <- ggplot(W_df, aes(x = dominant_subtype, y = dominant_proportion, fill = dominant_subtype)) +
  geom_boxplot(outlier.size = 0.8, alpha = 0.75) + geom_jitter(width = 0.15, alpha = 0.45, size = 1.0) +
  labs(title = "DD dominant subtype proportion", x = "Dominant subtype", y = "Dominant proportion") +
  theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "none")
ggsave(file.path(out_dir, "DD_only_dominant_subtype_proportion.png"), p_prop, width = 9, height = 6, dpi = 300)

W_long <- W_df %>% select(file_base, AgeGroup, Sex, dominant_subtype, all_of(component_cols)) %>%
  reshape2::melt(id.vars = c("file_base", "AgeGroup", "Sex", "dominant_subtype"), variable.name = "component", value.name = "score")
p_summary <- ggplot(W_long, aes(x = AgeGroup, y = score, fill = Sex)) +
  geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
  geom_point(position = position_jitterdodge(jitter.width = 0.15), alpha = 0.35, size = 0.7) +
  facet_wrap(~ component, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Female" = "#984EA3", "Male" = "#4DAF4A")) +
  labs(title = "DD-only NMF component scores overview", subtitle = "Model tested in main script: score ~ AgeGroup + Sex", x = "Age group", y = "Component score") +
  theme_bw(base_size = 11) + theme(legend.position = "bottom", strip.background = element_rect(fill = "lightblue"))
ggsave(file.path(out_dir, "DD_only_all_components_overview.png"), p_summary, width = 12, height = 10, dpi = 300)

if (!is.null(subtype_summary)) {
  write.csv(subtype_summary, file.path(out_dir, "DD_only_subtype_summary_copy.csv"), row.names = FALSE)
}

cat("\n==============================\n")
cat("DD-ONLY VISUALIZATION COMPLETE\n")
cat("==============================\n")
cat("Output directory:", out_dir, "\n")
cat("Generated individual boxplots:", length(component_cols), "\n")
cat("Generated top ROI plots:", length(unique(top20_rois$component)), "\n")
cat("DONE\n")
