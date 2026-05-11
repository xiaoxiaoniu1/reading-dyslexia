# ============================================================
# DK318 MIND Degree NMF - Network Enrichment Analysis
# Purpose: Test if NMF components are enriched in predefined networks
# Input: NMF component top ROIs and network definitions
# Output: Enrichment test results (Fisher's exact test)
# ============================================================

cat("Starting NMF Component Network Enrichment Analysis...\n\n")

# ============================================================
# 1. Load Packages
# ============================================================

packages <- c("dplyr", "readr", "tidyr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

library(dplyr); library(readr); library(tidyr)

# ============================================================
# 2. Define Paths
# ============================================================

nmf_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NMF_3"
network_base_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/brain_region_table_brainmaps"
out_dir <- nmf_dir

# Input files
roi_assignment_file <- file.path(nmf_dir, "NMF_ROI_hard_assignment.csv")
top_rois_file <- file.path(nmf_dir, "NMF_component_top_ROIs.csv")

# Network definition files - updated paths
network_files <- list(
  DMN = file.path(network_base_dir, "default/default_regions_used.csv"),
  MD = file.path(network_base_dir, "multiple_demand/multiple_demand_regions_used.csv"),
  Reading = file.path(network_base_dir, "mean_gt1p0_union_reading_v4_v5/mean_gt1p0_union_reading_v4_v5_regions_used.csv")
)

cat("NMF directory:", nmf_dir, "\n")
cat("Network base directory:", network_base_dir, "\n\n")

# ============================================================
# 3. Load NMF Results
# ============================================================

cat("Loading NMF results...\n")

# Load ROI hard assignment
roi_assignment <- read.csv(roi_assignment_file, stringsAsFactors = FALSE)
cat("Loaded ROI assignment:", nrow(roi_assignment), "ROIs\n")

# Load top ROIs
top_rois <- read.csv(top_rois_file, stringsAsFactors = FALSE)
cat("Loaded top ROIs:", nrow(top_rois), "entries\n\n")

# Get unique components
components <- unique(roi_assignment$assigned_component)
n_components <- length(components)
cat("Number of components:", n_components, "\n\n")

# ============================================================
# 4. Load Network Definitions
# ============================================================

cat("Loading network definitions...\n")

networks <- list()
network_names <- names(network_files)

for (net_name in network_names) {
  net_file <- network_files[[net_name]]
  
  if (file.exists(net_file)) {
    net_data <- read.csv(net_file, stringsAsFactors = FALSE)
    
    # The column name is "feature" in these network files
    roi_col <- "feature"
    
    if (roi_col %in% colnames(net_data)) {
      networks[[net_name]] <- unique(as.character(net_data[[roi_col]]))
      cat("  ", net_name, ":", length(networks[[net_name]]), "ROIs\n")
    } else {
      cat("  WARNING:", net_name, "- could not find 'feature' column\n")
      cat("    Available columns:", paste(colnames(net_data), collapse = ", "), "\n")
    }
  } else {
    cat("  WARNING:", net_name, "file not found:", net_file, "\n")
  }
}

cat("\n")

if (length(networks) == 0) {
  stop("ERROR: No network definitions loaded. Please check network file paths.")
}

# ============================================================
# 5. Fisher's Exact Test for Enrichment
# ============================================================

cat("Running Fisher's exact test for enrichment...\n\n")

# Function to perform Fisher's exact test
test_enrichment <- function(component_rois, network_rois, all_rois) {
  # 2x2 contingency table:
  #                In Network    Not in Network
  # In Component        a              b
  # Not in Component    c              d
  
  a <- sum(component_rois %in% network_rois)
  b <- sum(!(component_rois %in% network_rois))
  c <- sum((all_rois %in% network_rois) & !(all_rois %in% component_rois))
  d <- sum(!(all_rois %in% network_rois) & !(all_rois %in% component_rois))
  
  contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  
  # Fisher's exact test
  test_result <- fisher.test(contingency_table, alternative = "greater")
  
  # Calculate enrichment metrics
  n_component <- a + b
  n_network <- a + c
  n_total <- a + b + c + d
  
  expected <- (n_component * n_network) / n_total
  enrichment_ratio <- if (expected > 0) a / expected else NA_real_
  
  list(
    n_component_rois = n_component,
    n_network_rois = n_network,
    n_overlap = a,
    expected_overlap = expected,
    enrichment_ratio = enrichment_ratio,
    p_value = test_result$p.value,
    odds_ratio = test_result$estimate
  )
}

# ============================================================
# 6. Test Enrichment for Hard Assignment
# ============================================================

cat("Testing enrichment for hard-assigned ROIs...\n")

enrichment_hard <- data.frame()

all_rois <- roi_assignment$ROI

for (comp in components) {
  comp_rois <- roi_assignment$ROI[roi_assignment$assigned_component == comp]
  
  for (net_name in names(networks)) {
    net_rois <- networks[[net_name]]
    
    result <- test_enrichment(comp_rois, net_rois, all_rois)
    
    enrichment_hard <- rbind(enrichment_hard, data.frame(
      component = comp,
      network = net_name,
      n_component_rois = result$n_component_rois,
      n_network_rois = result$n_network_rois,
      n_overlap = result$n_overlap,
      expected_overlap = result$expected_overlap,
      enrichment_ratio = result$enrichment_ratio,
      odds_ratio = result$odds_ratio,
      p_value = result$p_value,
      stringsAsFactors = FALSE
    ))
  }
}

# FDR correction
enrichment_hard$q_value <- p.adjust(enrichment_hard$p_value, method = "fdr")

# Save results
hard_file <- file.path(out_dir, "NMF_network_enrichment_hard_assignment.csv")
write.csv(enrichment_hard, hard_file, row.names = FALSE)
cat("Saved:", hard_file, "\n\n")

# ============================================================
# 7. Test Enrichment for Top ROIs
# ============================================================

cat("Testing enrichment for top ROIs...\n")

enrichment_top <- data.frame()

top_types <- unique(top_rois$top_type)

for (comp in components) {
  for (top_type in top_types) {
    comp_top_rois <- top_rois$ROI[top_rois$component == comp & top_rois$top_type == top_type]
    
    for (net_name in names(networks)) {
      net_rois <- networks[[net_name]]
      
      result <- test_enrichment(comp_top_rois, net_rois, all_rois)
      
      enrichment_top <- rbind(enrichment_top, data.frame(
        component = comp,
        top_type = top_type,
        network = net_name,
        n_component_rois = result$n_component_rois,
        n_network_rois = result$n_network_rois,
        n_overlap = result$n_overlap,
        expected_overlap = result$expected_overlap,
        enrichment_ratio = result$enrichment_ratio,
        odds_ratio = result$odds_ratio,
        p_value = result$p_value,
        stringsAsFactors = FALSE
      ))
    }
  }
}

# FDR correction
enrichment_top$q_value <- p.adjust(enrichment_top$p_value, method = "fdr")

# Save results
top_file <- file.path(out_dir, "NMF_network_enrichment_top_rois.csv")
write.csv(enrichment_top, top_file, row.names = FALSE)
cat("Saved:", top_file, "\n\n")

# ============================================================
# 8. Summary
# ============================================================

cat("==============================\n")
cat("ENRICHMENT ANALYSIS SUMMARY\n")
cat("==============================\n\n")

cat("Significant enrichments (hard assignment, q < 0.05):\n")
sig_hard <- enrichment_hard %>% 
  filter(q_value < 0.05) %>%
  arrange(q_value)

if (nrow(sig_hard) > 0) {
  print(sig_hard %>% select(component, network, n_overlap, enrichment_ratio, p_value, q_value), 
        row.names = FALSE)
} else {
  cat("  No significant enrichments at FDR < 0.05\n")
}

cat("\n")

cat("Significant enrichments (top ROIs, q < 0.05):\n")
sig_top <- enrichment_top %>% 
  filter(q_value < 0.05) %>%
  arrange(q_value)

if (nrow(sig_top) > 0) {
  print(sig_top %>% select(component, top_type, network, n_overlap, enrichment_ratio, p_value, q_value), 
        row.names = FALSE)
} else {
  cat("  No significant enrichments at FDR < 0.05\n")
}

cat("\n==============================\n")
cat("DONE\n")
cat("==============================\n")
