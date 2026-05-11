# ============================================================
# DK318 MIND Degree T-map Clustering - Network Enrichment Analysis & Visualization
# Purpose: Combined script for enrichment testing and comprehensive visualization
# Input: Cluster assignments and network definitions
# Output: Enrichment test results, heatmaps, bar plots, pie charts, and detailed reports
# ============================================================

cat("Starting Cluster Network Enrichment Analysis & Visualization...\n\n")

# ============================================================
# PART 1: ENRICHMENT ANALYSIS
# ============================================================
cat("==============================\n")
cat("PART 1: ENRICHMENT ANALYSIS\n")
cat("==============================\n\n")

# Load Packages
packages <- c("dplyr", "readr", "tidyr", "ggplot2", "gridExtra", "scales")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

library(dplyr); library(readr); library(tidyr); library(ggplot2); library(gridExtra); library(scales)

# Define Paths
cluster_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Clustering_5/clustering"
network_base_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/brain_region_table_brainmaps"
yeo_base_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/DK-318"
out_dir <- cluster_dir

cluster_file <- file.path(cluster_dir, "degree_Tmap_cluster_K5.csv")

network_files <- list(
  DMN = file.path(network_base_dir, "default/default_regions_used.csv"),
  MD = file.path(network_base_dir, "multiple_demand/multiple_demand_regions_used.csv"),
  Reading = file.path(network_base_dir, "mean_gt1p0_union_reading_v4_v5/mean_gt1p0_union_reading_v4_v5_regions_used.csv")
)

# Yeo network mapping files
yeo_files <- list(
  Yeo7_lh = file.path(yeo_base_dir, "lh_DK318_to_Yeo7_mapping.csv"),
  Yeo7_rh = file.path(yeo_base_dir, "rh_DK318_to_Yeo7_mapping.csv"),
  Yeo17_lh = file.path(yeo_base_dir, "lh_DK318_to_Yeo17_mapping.csv"),
  Yeo17_rh = file.path(yeo_base_dir, "rh_DK318_to_Yeo17_mapping.csv")
)

cat("Cluster directory:", cluster_dir, "\n")
cat("Network base directory:", network_base_dir, "\n")
cat("Yeo network directory:", yeo_base_dir, "\n\n")

# Load Cluster Results
cat("Loading cluster results...\n")
cluster_data <- read.csv(cluster_file, stringsAsFactors = FALSE)
cat("Loaded cluster data:", nrow(cluster_data), "ROIs\n")

clusters <- sort(unique(cluster_data$cluster))
n_clusters <- length(clusters)
cat("Number of clusters:", n_clusters, "\n\n")

# Load Network Definitions
cat("Loading network definitions...\n")
networks <- list()

# Load original networks (DMN, MD, Reading)
for (net_name in names(network_files)) {
  net_file <- network_files[[net_name]]
  if (file.exists(net_file)) {
    net_data <- read.csv(net_file, stringsAsFactors = FALSE)
    if ("feature" %in% colnames(net_data)) {
      networks[[net_name]] <- unique(as.character(net_data[["feature"]]))
      cat("  ", net_name, ":", length(networks[[net_name]]), "ROIs\n")
    }
  }
}

# Load Yeo7 networks
cat("\nLoading Yeo7 networks...\n")
yeo7_networks <- list()
for (yeo_key in c("Yeo7_lh", "Yeo7_rh")) {
  yeo_file <- yeo_files[[yeo_key]]
  hemi_prefix <- ifelse(grepl("_lh$", yeo_key), "lh_", "rh_")
  if (file.exists(yeo_file)) {
    yeo_data <- read.csv(yeo_file, stringsAsFactors = FALSE)
    if ("DK318_region" %in% colnames(yeo_data) && "primary_Yeo7_network" %in% colnames(yeo_data)) {
      yeo_data$feature <- paste0(hemi_prefix, yeo_data$DK318_region)
      for (net in unique(yeo_data$primary_Yeo7_network)) {
        net_rois <- yeo_data$feature[yeo_data$primary_Yeo7_network == net]
        net_name <- paste0("Yeo7_", net)
        if (net_name %in% names(yeo7_networks)) {
          yeo7_networks[[net_name]] <- unique(c(yeo7_networks[[net_name]], net_rois))
        } else {
          yeo7_networks[[net_name]] <- unique(net_rois)
        }
      }
    }
  }
}
for (net_name in names(yeo7_networks)) {
  cat("  ", net_name, ":", length(yeo7_networks[[net_name]]), "ROIs\n")
}

# Load Yeo17 networks
cat("\nLoading Yeo17 networks...\n")
yeo17_networks <- list()
for (yeo_key in c("Yeo17_lh", "Yeo17_rh")) {
  yeo_file <- yeo_files[[yeo_key]]
  hemi_prefix <- ifelse(grepl("_lh$", yeo_key), "lh_", "rh_")
  if (file.exists(yeo_file)) {
    yeo_data <- read.csv(yeo_file, stringsAsFactors = FALSE)
    if ("DK318_region" %in% colnames(yeo_data) && "primary_Yeo17_network" %in% colnames(yeo_data)) {
      yeo_data$feature <- paste0(hemi_prefix, yeo_data$DK318_region)
      for (net in unique(yeo_data$primary_Yeo17_network)) {
        net_rois <- yeo_data$feature[yeo_data$primary_Yeo17_network == net]
        net_name <- paste0("Yeo17_", net)
        if (net_name %in% names(yeo17_networks)) {
          yeo17_networks[[net_name]] <- unique(c(yeo17_networks[[net_name]], net_rois))
        } else {
          yeo17_networks[[net_name]] <- unique(net_rois)
        }
      }
    }
  }
}
for (net_name in names(yeo17_networks)) {
  cat("  ", net_name, ":", length(yeo17_networks[[net_name]]), "ROIs\n")
}

# Combine all networks
networks <- c(networks, yeo7_networks, yeo17_networks)
cat("\nTotal networks loaded:", length(networks), "\n")

network_match_summary <- data.frame(
  network = names(networks),
  n_defined_rois = sapply(networks, length),
  n_matched_to_cluster_rois = sapply(networks, function(x) sum(x %in% cluster_data$feature)),
  stringsAsFactors = FALSE
)
cat("\nNetwork ROI matching summary:\n")
print(network_match_summary, row.names = FALSE)
write.csv(network_match_summary, file.path(out_dir, "network_roi_matching_summary.csv"), row.names = FALSE)
cat("Saved: network_roi_matching_summary.csv\n\n")

# Fisher's Exact Test Function
test_enrichment <- function(cluster_rois, network_rois, all_rois) {
  cluster_rois <- intersect(cluster_rois, all_rois)
  network_rois <- intersect(network_rois, all_rois)
  
  a <- sum(cluster_rois %in% network_rois)
  b <- sum(!(cluster_rois %in% network_rois))
  c <- sum((all_rois %in% network_rois) & !(all_rois %in% cluster_rois))
  d <- sum(!(all_rois %in% network_rois) & !(all_rois %in% cluster_rois))
  
  contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  test_result <- fisher.test(contingency_table, alternative = "greater")
  
  n_cluster <- a + b
  n_network <- a + c
  n_total <- a + b + c + d
  expected <- (n_cluster * n_network) / n_total
  enrichment_ratio <- if (expected > 0) a / expected else NA_real_
  
  list(
    n_cluster_rois = n_cluster, n_network_rois = n_network, n_overlap = a,
    expected_overlap = expected, enrichment_ratio = enrichment_ratio,
    p_value = test_result$p.value, odds_ratio = test_result$estimate
  )
}

# Test Enrichment
cat("Testing enrichment for all clusters...\n")
enrichment_results <- data.frame()
all_rois <- cluster_data$feature

for (clust in clusters) {
  clust_rois <- cluster_data$feature[cluster_data$cluster == clust]
  for (net_name in names(networks)) {
    result <- test_enrichment(clust_rois, networks[[net_name]], all_rois)
    enrichment_results <- rbind(enrichment_results, data.frame(
      cluster = clust, network = net_name,
      n_cluster_rois = result$n_cluster_rois, n_network_rois = result$n_network_rois,
      n_overlap = result$n_overlap, expected_overlap = result$expected_overlap,
      enrichment_ratio = result$enrichment_ratio, odds_ratio = result$odds_ratio,
      p_value = result$p_value, stringsAsFactors = FALSE
    ))
  }
}

enrichment_results <- enrichment_results %>%
  mutate(
    network_category = case_when(
      grepl("^Yeo7_", network) ~ "Yeo7",
      grepl("^Yeo17_", network) ~ "Yeo17",
      TRUE ~ "Original"
    )
  ) %>%
  group_by(network_category) %>%
  mutate(q_value = p.adjust(p_value, method = "fdr")) %>%
  ungroup()

write.csv(enrichment_results, file.path(out_dir, "cluster_network_enrichment.csv"), row.names = FALSE)
cat("Saved: cluster_network_enrichment.csv\n")
cat("FDR correction: applied separately within Original, Yeo7, and Yeo17 network systems.\n\n")

# Add Cluster Characteristics
cluster_chars <- cluster_data %>%
  group_by(cluster) %>%
  summarise(mean_t_child = mean(t_child_TD_minus_DD),
            mean_t_adult = mean(t_adult_TD_minus_DD), .groups = "drop") %>%
  mutate(pattern = case_when(
    mean_t_child > 0 & mean_t_adult > 0 ~ "Persistent TD>DD",
    mean_t_child < 0 & mean_t_adult < 0 ~ "Persistent DD>TD",
    mean_t_child < 0 & mean_t_adult > 0 ~ "Developmental Catch-up",
    mean_t_child > 0 & mean_t_adult < 0 ~ "Developmental Reversal",
    TRUE ~ "Mixed"
  ))

enrichment_with_pattern <- enrichment_results %>% left_join(cluster_chars, by = "cluster")
write.csv(enrichment_with_pattern, file.path(out_dir, "cluster_network_enrichment_with_pattern.csv"), row.names = FALSE)
cat("Saved: cluster_network_enrichment_with_pattern.csv\n\n")

# Detailed Overlap
cat("Extracting detailed overlap information...\n")
overlap_details <- data.frame()
for (clust in clusters) {
  clust_rois <- cluster_data$feature[cluster_data$cluster == clust]
  for (net_name in names(networks)) {
    overlap_rois <- intersect(clust_rois, networks[[net_name]])
    if (length(overlap_rois) > 0) {
      for (roi in overlap_rois) {
        overlap_details <- rbind(overlap_details, data.frame(
          cluster = clust, network = net_name, ROI = roi, stringsAsFactors = FALSE))
      }
    }
  }
}

overlap_details <- overlap_details %>%
  left_join(cluster_data %>% select(feature, t_child_TD_minus_DD, t_adult_TD_minus_DD),
            by = c("ROI" = "feature"))

write.csv(overlap_details, file.path(out_dir, "cluster_network_overlap_details.csv"), row.names = FALSE)
cat("Saved: cluster_network_overlap_details.csv\n\n")

# ============================================================
# PART 2: ENHANCED VISUALIZATIONS
# ============================================================
cat("==============================\n")
cat("PART 2: ENHANCED VISUALIZATIONS\n")
cat("==============================\n\n")

# 1. Enhanced Heatmap with Numbers
cat("Creating enhanced heatmap...\n")

enrichment_plot <- enrichment_with_pattern %>%
  mutate(
    cluster_label = paste0("Cluster ", cluster, "\n", pattern),
    sig_label = case_when(
      q_value < 0.001 ~ "***",
      q_value < 0.01 ~ "**", 
      q_value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    display_label = paste0(n_overlap, "\n(", round(enrichment_ratio, 2), ")", sig_label),
    network_category = case_when(
      grepl("^Yeo7_", network) ~ "Yeo7",
      grepl("^Yeo17_", network) ~ "Yeo17",
      TRUE ~ "Original"
    )
  )

# Create separate heatmaps for each network category
# Original networks (DMN, MD, Reading)
enrichment_original <- enrichment_plot %>% filter(network_category == "Original")

cat("Network categories summary:\n")
cat("  Original networks:", nrow(enrichment_original), "rows\n")
cat("  Yeo7 networks:", nrow(enrichment_plot %>% filter(network_category == "Yeo7")), "rows\n")
cat("  Yeo17 networks:", nrow(enrichment_plot %>% filter(network_category == "Yeo17")), "rows\n")
cat("  Unique networks in data:", paste(unique(enrichment_plot$network), collapse = ", "), "\n\n")

if (nrow(enrichment_original) > 0) {
  p1 <- ggplot(enrichment_original, aes(x = network, y = cluster_label)) +
    geom_tile(aes(fill = enrichment_ratio), color = "white", linewidth = 1.5) +
    geom_text(aes(label = display_label), size = 4.5, fontface = "bold") +
    scale_fill_gradient2(
      low = "#2166AC", mid = "white", high = "#B2182B",
      midpoint = 1,
      limits = c(0.5, max(enrichment_original$enrichment_ratio, na.rm = TRUE)),
      name = "Enrichment\nRatio",
      breaks = c(0.5, 1.0, 1.5, 2.0)
    ) +
    labs(
      title = "Cluster Network Enrichment Analysis - Original Networks",
      subtitle = "Numbers show: Overlap Count\n(Enrichment Ratio) Significance\n* network-system FDR q<0.05, ** q<0.01, *** q<0.001",
      x = "Functional Network",
      y = "Cluster (Developmental Pattern)"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(size = 13, face = "bold", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 11),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      panel.grid = element_blank(),
      legend.position = "right"
    )
  
  ggsave(file.path(out_dir, "enrichment_heatmap_enhanced_original.png"), 
         p1, width = 10, height = 7, dpi = 300)
  cat("Saved: enrichment_heatmap_enhanced_original.png\n")
}

# Yeo7 networks
enrichment_yeo7 <- enrichment_plot %>% filter(network_category == "Yeo7")

if (nrow(enrichment_yeo7) > 0 && any(!is.na(enrichment_yeo7$enrichment_ratio))) {
  max_ratio_yeo7 <- max(enrichment_yeo7$enrichment_ratio, na.rm = TRUE)
  if (is.finite(max_ratio_yeo7) && max_ratio_yeo7 > 0.5) {
    p1_yeo7 <- ggplot(enrichment_yeo7, aes(x = network, y = cluster_label)) +
      geom_tile(aes(fill = enrichment_ratio), color = "white", linewidth = 1.5) +
      geom_text(aes(label = display_label), size = 3.5, fontface = "bold") +
      scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 1,
        limits = c(0.5, max_ratio_yeo7),
        name = "Enrichment\nRatio"
      ) +
      labs(
        title = "Cluster Network Enrichment Analysis - Yeo7 Networks",
        subtitle = "Numbers show: Overlap Count\n(Enrichment Ratio) Significance\n* network-system FDR q<0.05, ** q<0.01, *** q<0.001",
        x = "Yeo7 Network",
        y = "Cluster (Developmental Pattern)"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(size = 10, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 9),
        panel.grid = element_blank(),
        legend.position = "right"
      )
    
    ggsave(file.path(out_dir, "enrichment_heatmap_enhanced_yeo7.png"), 
           p1_yeo7, width = 12, height = 7, dpi = 300)
    cat("Saved: enrichment_heatmap_enhanced_yeo7.png\n")
  } else {
    cat("Skipped Yeo7 heatmap: invalid enrichment ratios\n")
  }
} else {
  cat("Skipped Yeo7 heatmap: no data available\n")
}

# Yeo17 networks
enrichment_yeo17 <- enrichment_plot %>% filter(network_category == "Yeo17")

if (nrow(enrichment_yeo17) > 0 && any(!is.na(enrichment_yeo17$enrichment_ratio))) {
  max_ratio_yeo17 <- max(enrichment_yeo17$enrichment_ratio, na.rm = TRUE)
  if (is.finite(max_ratio_yeo17) && max_ratio_yeo17 > 0.5) {
    p1_yeo17 <- ggplot(enrichment_yeo17, aes(x = network, y = cluster_label)) +
      geom_tile(aes(fill = enrichment_ratio), color = "white", linewidth = 1.5) +
      geom_text(aes(label = display_label), size = 2.8, fontface = "bold") +
      scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 1,
        limits = c(0.5, max_ratio_yeo17),
        name = "Enrichment\nRatio"
      ) +
      labs(
        title = "Cluster Network Enrichment Analysis - Yeo17 Networks",
        subtitle = "Numbers show: Overlap Count\n(Enrichment Ratio) Significance\n* network-system FDR q<0.05, ** q<0.01, *** q<0.001",
        x = "Yeo17 Network",
        y = "Cluster (Developmental Pattern)"
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.x = element_text(size = 9, face = "bold", angle = 45, hjust = 1),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        plot.subtitle = element_text(hjust = 0.5, size = 8),
        panel.grid = element_blank(),
        legend.position = "right"
      )
    
    ggsave(file.path(out_dir, "enrichment_heatmap_enhanced_yeo17.png"), 
           p1_yeo17, width = 16, height = 7, dpi = 300)
    cat("Saved: enrichment_heatmap_enhanced_yeo17.png\n")
  } else {
    cat("Skipped Yeo17 heatmap: invalid enrichment ratios\n")
  }
} else {
  cat("Skipped Yeo17 heatmap: no data available\n")
}

cat("\n")

# 2. Bar Plot: Enrichment Ratio by Cluster
cat("Creating enrichment ratio bar plots by network category...\n")

enrichment_plot2 <- enrichment_with_pattern %>%
  mutate(
    cluster_label = paste0("Cluster ", cluster, "\n", pattern),
    sig = ifelse(q_value < 0.05, "Significant", "Not Significant"),
    network_category = case_when(
      grepl("^Yeo7_", network) ~ "Yeo7",
      grepl("^Yeo17_", network) ~ "Yeo17",
      TRUE ~ "Original"
    ),
    network_label = case_when(
      network_category == "Yeo7" ~ sub("^Yeo7_", "", network),
      network_category == "Yeo17" ~ sub("^Yeo17_", "", network),
      TRUE ~ network
    ),
    overlap_label = ifelse(is.na(enrichment_ratio), "", as.character(n_overlap))
  )

create_enrichment_barplot <- function(plot_data, filename, title, x_title, width, height, x_text_size = 10, label_size = 3.5) {
  plot_data <- plot_data %>% filter(!is.na(enrichment_ratio), is.finite(enrichment_ratio))
  
  if (nrow(plot_data) == 0) {
    cat("Skipped:", filename, "- no valid enrichment ratios\n")
    return(invisible(NULL))
  }
  
  p <- ggplot(plot_data, aes(x = network_label, y = enrichment_ratio, fill = sig)) +
    geom_col(color = "black", linewidth = 0.4, width = 0.75) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red", linewidth = 0.8) +
    geom_text(aes(label = overlap_label), vjust = -0.35, size = label_size, fontface = "bold") +
    facet_wrap(~cluster_label, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = c("Significant" = "#E41A1C", "Not Significant" = "gray70")) +
    labs(
      title = title,
      subtitle = "Numbers above bars show overlap ROI count; red dashed line = no enrichment (ratio = 1)",
      x = x_title,
      y = "Enrichment Ratio",
      fill = "Network-system FDR q < 0.05"
    ) +
    theme_bw(base_size = 12) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = x_text_size),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 9),
      legend.position = "bottom",
      panel.grid.major.x = element_blank()
    )
  
  ggsave(file.path(out_dir, filename), p, width = width, height = height, dpi = 300)
  cat("Saved:", filename, "\n")
}

create_enrichment_barplot(
  enrichment_plot2 %>% filter(network_category == "Original"),
  "enrichment_barplot_by_cluster_original.png",
  "Network Enrichment Ratio by Cluster - Original Networks",
  "Original Network",
  width = 8,
  height = 10,
  x_text_size = 11,
  label_size = 4
)

create_enrichment_barplot(
  enrichment_plot2 %>% filter(network_category == "Yeo7"),
  "enrichment_barplot_by_cluster_yeo7.png",
  "Network Enrichment Ratio by Cluster - Yeo7 Networks",
  "Yeo7 Network",
  width = 10,
  height = 10,
  x_text_size = 10,
  label_size = 3.5
)

create_enrichment_barplot(
  enrichment_plot2 %>% filter(network_category == "Yeo17"),
  "enrichment_barplot_by_cluster_yeo17.png",
  "Network Enrichment Ratio by Cluster - Yeo17 Networks",
  "Yeo17 Network",
  width = 16,
  height = 10,
  x_text_size = 8,
  label_size = 3
)

cat("\n")

# 3. Pie Charts: Network Distribution per Cluster
cat("Creating network distribution pie charts...\n")

network_dist <- overlap_details %>%
  group_by(cluster, network) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(cluster) %>%
  mutate(
    total = sum(n),
    percentage = n / total * 100,
    label = paste0(network, "\n", n, " (", round(percentage, 1), "%)")
  ) %>%
  left_join(enrichment_with_pattern %>% select(cluster, pattern) %>% distinct(), by = "cluster")

pie_plots <- lapply(sort(unique(network_dist$cluster)), function(clust) {
  data <- network_dist %>% filter(cluster == clust)
  pattern <- unique(data$pattern)
  
  ggplot(data, aes(x = "", y = n, fill = network)) +
    geom_bar(stat = "identity", width = 1, color = "white", size = 1.5) +
    coord_polar("y", start = 0) +
    geom_text(aes(label = paste0(n, "\n(", round(percentage, 1), "%)")),
              position = position_stack(vjust = 0.5), size = 4.5, fontface = "bold") +
    scale_fill_manual(values = c("DMN" = "#1B9E77", "MD" = "#D95F02", "Reading" = "#7570B3")) +
    labs(title = paste0("Cluster ", clust, "\n", pattern)) +
    theme_void(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 11),
      legend.position = "bottom",
      legend.title = element_blank()
    )
})

p3 <- grid.arrange(grobs = pie_plots, ncol = 2,
                   top = grid::textGrob("Network Distribution by Cluster\n(ROI counts in overlapping regions)",
                                       gp = grid::gpar(fontsize = 16, fontface = "bold")))

ggsave(file.path(out_dir, "network_distribution_pie_charts.png"),
       p3, width = 12, height = 10, dpi = 300)
cat("Saved: network_distribution_pie_charts.png\n\n")

# 4. Basic Heatmap (log2 enrichment)
cat("Creating basic enrichment heatmap...\n")

heatmap_data <- enrichment_with_pattern %>%
  mutate(log_enrichment = log2(enrichment_ratio),
         significance = case_when(q_value < 0.001 ~ "***", q_value < 0.01 ~ "**",
                                  q_value < 0.05 ~ "*", TRUE ~ ""),
         cluster_label = paste0("Cluster ", cluster, "\n", pattern))

p_heatmap <- ggplot(heatmap_data, aes(x = network, y = cluster_label, fill = log_enrichment)) +
  geom_tile(color = "white", size = 1) +
  geom_text(aes(label = significance), size = 8, vjust = 0.75) +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0,
                       limits = c(-max(abs(heatmap_data$log_enrichment), na.rm = TRUE),
                                  max(abs(heatmap_data$log_enrichment), na.rm = TRUE)),
                       name = "Log2\nEnrichment\nRatio") +
  labs(title = "Cluster Network Enrichment Analysis",
       subtitle = "Fisher's Exact Test\n* network-system FDR q<0.05, ** q<0.01, *** q<0.001",
       x = "Network", y = "Cluster") +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5), panel.grid = element_blank())

ggsave(file.path(out_dir, "cluster_network_enrichment_heatmap.png"), p_heatmap, width = 8, height = 6, dpi = 300)
cat("Saved: cluster_network_enrichment_heatmap.png\n\n")

# 5. Bar Plot for Significant Enrichments Only
if (any(enrichment_with_pattern$q_value < 0.05)) {
  cat("Creating bar plot for significant enrichments...\n")
  sig_data <- enrichment_with_pattern %>% filter(q_value < 0.05) %>%
    mutate(cluster_label = paste0("Cluster ", cluster, "\n(", pattern, ")"))
  
  p_bar <- ggplot(sig_data, aes(x = network, y = enrichment_ratio, fill = network)) +
    geom_bar(stat = "identity", color = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", color = "red") +
    facet_wrap(~cluster_label, ncol = 2) +
    labs(title = "Significant Network Enrichments by Cluster", subtitle = "Network-system FDR q < 0.05",
         x = "Network", y = "Enrichment Ratio") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5), legend.position = "none")
  
  ggsave(file.path(out_dir, "cluster_network_enrichment_barplot.png"), p_bar, width = 10, height = 8, dpi = 300)
  cat("Saved: cluster_network_enrichment_barplot.png\n\n")
}

# ============================================================
# PART 3: DETAILED SUMMARY TABLES AND REPORTS
# ============================================================
cat("==============================\n")
cat("PART 3: SUMMARY TABLES & REPORTS\n")
cat("==============================\n\n")

# Detailed Summary Table
cat("Creating detailed summary table...\n")

summary_table <- enrichment_with_pattern %>%
  arrange(cluster, p_value) %>%
  mutate(
    sig = ifelse(q_value < 0.05, "***", ""),
    enrichment_ratio = round(enrichment_ratio, 3),
    p_value_fmt = format.pval(p_value, digits = 3),
    q_value_fmt = format.pval(q_value, digits = 3)
  ) %>%
  select(cluster, pattern, network_category, network, n_cluster_rois, n_network_rois, 
         n_overlap, enrichment_ratio, p_value_fmt, q_value_fmt, sig)

write.csv(summary_table, file.path(out_dir, "enrichment_summary_table.csv"), row.names = FALSE)
cat("Saved: enrichment_summary_table.csv\n\n")

# Text Summary Report
cat("Generating text summary report...\n")

sink(file.path(out_dir, "enrichment_summary_report.txt"))

cat(strrep("=", 70), "\n")
cat("CLUSTER NETWORK ENRICHMENT ANALYSIS - SUMMARY REPORT\n")
cat(strrep("=", 70), "\n\n")

cat("Analysis Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Overall statistics
cat("OVERALL STATISTICS\n")
cat(strrep("-", 70), "\n")
cat("Total clusters:", length(unique(enrichment_with_pattern$cluster)), "\n")
cat("Total networks tested:", length(unique(enrichment_with_pattern$network)), "\n")
cat("Total tests performed:", nrow(enrichment_with_pattern), "\n")
cat("FDR correction: separately within Original, Yeo7, and Yeo17 network systems\n")
cat("Significant enrichments (network-system FDR q < 0.05):", sum(enrichment_with_pattern$q_value < 0.05), "\n\n")

cat("Tests by network system:\n")
print(enrichment_with_pattern %>% count(network_category, name = "n_tests"), row.names = FALSE)
cat("\n")

# By cluster
for (clust in sort(unique(enrichment_with_pattern$cluster))) {
  clust_data <- enrichment_with_pattern %>% filter(cluster == clust)
  pattern <- unique(clust_data$pattern)
  
  cat("\n")
  cat(strrep("=", 70), "\n")
  cat("CLUSTER", clust, ":", pattern, "\n")
  cat(strrep("=", 70), "\n\n")
  
  cat("Cluster size:", unique(clust_data$n_cluster_rois), "ROIs\n")
  cat("Mean t-value (child):", round(unique(clust_data$mean_t_child), 3), "\n")
  cat("Mean t-value (adult):", round(unique(clust_data$mean_t_adult), 3), "\n\n")
  
  cat("Network Enrichment Results:\n")
  cat(strrep("-", 70), "\n")
  
  for (i in 1:nrow(clust_data)) {
    row <- clust_data[i, ]
    cat(sprintf("  %s Network:\n", row$network))
    cat(sprintf("    - Overlap: %d ROIs (expected: %.1f)\n", 
                row$n_overlap, row$expected_overlap))
    cat(sprintf("    - Enrichment Ratio: %.3f\n", row$enrichment_ratio))
    cat(sprintf("    - p-value: %.4f, q-value: %.4f", row$p_value, row$q_value))
    if (row$q_value < 0.05) {
      cat(" ***SIGNIFICANT***")
    }
    cat("\n\n")
  }
  
  # Top ROIs in each network
  clust_overlap <- overlap_details %>% filter(cluster == clust)
  if (nrow(clust_overlap) > 0) {
    cat("Top ROIs by Network:\n")
    cat(strrep("-", 70), "\n")
    for (net in unique(clust_overlap$network)) {
      net_rois <- clust_overlap %>% filter(network == net) %>% arrange(desc(abs(t_adult_TD_minus_DD)))
      cat(sprintf("  %s (%d ROIs):\n", net, nrow(net_rois)))
      if (nrow(net_rois) > 0) {
        top_n <- min(5, nrow(net_rois))
        for (j in 1:top_n) {
          roi <- net_rois[j, ]
          cat(sprintf("    %d. %s (t_child=%.2f, t_adult=%.2f)\n",
                      j, roi$ROI, roi$t_child_TD_minus_DD, roi$t_adult_TD_minus_DD))
        }
      }
      cat("\n")
    }
  }
}

cat("\n")
cat(strrep("=", 70), "\n")
cat("END OF REPORT\n")
cat(strrep("=", 70), "\n")

sink()

cat("Saved: enrichment_summary_report.txt\n\n")

# ============================================================
# PART 4: CONSOLE SUMMARY
# ============================================================
cat("==============================\n")
cat("ENRICHMENT ANALYSIS SUMMARY\n")
cat("==============================\n\n")

cat("Cluster patterns:\n")
print(cluster_chars, row.names = FALSE)
cat("\n")

sig_enrichment <- enrichment_with_pattern %>% 
  filter(q_value < 0.05) %>% arrange(q_value) %>%
  select(cluster, pattern, network_category, network, n_overlap, enrichment_ratio, p_value, q_value)

cat("Significant enrichments (network-system FDR q < 0.05):\n")
if (nrow(sig_enrichment) > 0) {
  print(sig_enrichment, row.names = FALSE)
  cat("\n")
  for (i in 1:nrow(sig_enrichment)) {
    row <- sig_enrichment[i, ]
    cat(sprintf("Cluster %d (%s) -> %s Network\n", 
                row$cluster, row$pattern, row$network))
    cat(sprintf("  Enrichment Ratio: %.3f, q-value: %.4f\n", 
                row$enrichment_ratio, row$q_value))
    cat(sprintf("  %d ROIs overlap (expected: %.1f)\n\n", 
                row$n_overlap, 
                enrichment_with_pattern$expected_overlap[enrichment_with_pattern$cluster == row$cluster & 
                                                          enrichment_with_pattern$network == row$network]))
  }
} else {
  cat("  No significant enrichments at network-system FDR < 0.05\n\n")
}

cat("All enrichments ranked by p-value:\n")
ranked <- enrichment_with_pattern %>% 
  arrange(p_value) %>%
  select(cluster, pattern, network_category, network, n_overlap, enrichment_ratio, p_value, q_value)
print(ranked, row.names = FALSE)

cat("\n")
cat(strrep("=", 70), "\n")
cat("ANALYSIS COMPLETE\n")
cat(strrep("=", 70), "\n\n")

cat("Generated files:\n")
cat("  Data files:\n")
cat("    - cluster_network_enrichment.csv\n")
cat("    - cluster_network_enrichment_with_pattern.csv\n")
cat("    - cluster_network_overlap_details.csv\n")
cat("    - enrichment_summary_table.csv\n")
cat("    - enrichment_summary_report.txt\n\n")
cat("  Visualization files:\n")
cat("    - enrichment_heatmap_enhanced_original.png (original networks: DMN, MD, Reading)\n")
cat("    - enrichment_heatmap_enhanced_yeo7.png (Yeo7 networks)\n")
cat("    - enrichment_heatmap_enhanced_yeo17.png (Yeo17 networks)\n")
cat("    - cluster_network_enrichment_heatmap.png (log2 enrichment)\n")
cat("    - enrichment_barplot_by_cluster_original.png (original networks: DMN, MD, Reading)\n")
cat("    - enrichment_barplot_by_cluster_yeo7.png (Yeo7 networks)\n")
cat("    - enrichment_barplot_by_cluster_yeo17.png (Yeo17 networks)\n")
if (any(enrichment_with_pattern$q_value < 0.05)) {
  cat("    - cluster_network_enrichment_barplot.png (significant only)\n")
}
cat("    - network_distribution_pie_charts.png\n\n")

cat("Key findings:\n")
cat("  - Check enrichment_heatmap_enhanced_original.png for DMN/MD/Reading patterns\n")
cat("  - Check enrichment_heatmap_enhanced_yeo7.png for Yeo7 network patterns\n")
cat("  - Check enrichment_heatmap_enhanced_yeo17.png for Yeo17 network patterns\n")
cat("  - Numbers in heatmap show: overlap count (enrichment ratio) significance\n")
cat("  - Enrichment ratio > 1 means more overlap than expected by chance\n")
cat("  - Enrichment ratio < 1 means less overlap than expected\n")
cat("  - Read enrichment_summary_report.txt for detailed results\n\n")

cat("==============================\n")
cat("DONE\n")
cat("==============================\n")
