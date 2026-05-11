# ============================================================
# DK318 MIND Degree NMF - DD-only Network Enrichment
# Purpose: Test whether DD-only NMF subtype components are enriched
#          in predefined DMN, MD, and Reading networks.
# ============================================================

cat("Starting DD-only NMF Network Enrichment Analysis...\n\n")

packages <- c("dplyr", "readr", "tidyr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}
library(dplyr); library(readr); library(tidyr)

nmf_dir <- Sys.getenv("DD_NMF_OUT_DIR", "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NMF_DD_only")
network_base_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/brain_region_table_brainmaps"
out_dir <- nmf_dir

roi_assignment_file <- file.path(nmf_dir, "NMF_DD_only_ROI_hard_assignment.csv")
top_rois_file <- file.path(nmf_dir, "NMF_DD_only_component_top_ROIs.csv")

network_files <- list(
  DMN = file.path(network_base_dir, "default/default_regions_used.csv"),
  MD = file.path(network_base_dir, "multiple_demand/multiple_demand_regions_used.csv"),
  Reading = file.path(network_base_dir, "mean_gt1p0_union_reading_v4_v5/mean_gt1p0_union_reading_v4_v5_regions_used.csv")
)

cat("DD-only NMF directory:", nmf_dir, "\n")
cat("Network base directory:", network_base_dir, "\n\n")

if (!file.exists(roi_assignment_file)) stop("Missing ROI assignment file: ", roi_assignment_file)
if (!file.exists(top_rois_file)) stop("Missing top ROIs file: ", top_rois_file)

roi_assignment <- read.csv(roi_assignment_file, stringsAsFactors = FALSE)
top_rois <- read.csv(top_rois_file, stringsAsFactors = FALSE)
components <- unique(roi_assignment$assigned_component)
all_rois <- roi_assignment$ROI

cat("Loaded", length(all_rois), "ROIs and", length(components), "components.\n\n")

networks <- list()
for (net_name in names(network_files)) {
  net_file <- network_files[[net_name]]
  if (!file.exists(net_file)) {
    cat("WARNING:", net_name, "file not found:", net_file, "\n")
    next
  }
  net_data <- read.csv(net_file, stringsAsFactors = FALSE)
  if (!"feature" %in% colnames(net_data)) {
    cat("WARNING:", net_name, "does not contain feature column. Available:", paste(colnames(net_data), collapse = ", "), "\n")
    next
  }
  networks[[net_name]] <- unique(as.character(net_data$feature))
  cat("  ", net_name, ":", length(networks[[net_name]]), "ROIs\n")
}
cat("\n")
if (length(networks) == 0) stop("No network definitions loaded.")

test_enrichment <- function(component_rois, network_rois, all_rois) {
  a <- sum(component_rois %in% network_rois)
  b <- sum(!(component_rois %in% network_rois))
  c <- sum((all_rois %in% network_rois) & !(all_rois %in% component_rois))
  d <- sum(!(all_rois %in% network_rois) & !(all_rois %in% component_rois))
  contingency_table <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
  test_result <- fisher.test(contingency_table, alternative = "greater")
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

cat("Testing enrichment for hard-assigned ROIs...\n")
enrichment_hard <- data.frame()
for (comp in components) {
  comp_rois <- roi_assignment$ROI[roi_assignment$assigned_component == comp]
  for (net_name in names(networks)) {
    result <- test_enrichment(comp_rois, networks[[net_name]], all_rois)
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
enrichment_hard$q_value <- p.adjust(enrichment_hard$p_value, method = "fdr")
write.csv(enrichment_hard, file.path(out_dir, "NMF_DD_only_network_enrichment_hard_assignment.csv"), row.names = FALSE)
write.csv(enrichment_hard, file.path(out_dir, "NMF_network_enrichment_hard_assignment.csv"), row.names = FALSE)

cat("Testing enrichment for top ROIs...\n")
enrichment_top <- data.frame()
for (comp in components) {
  for (top_type in unique(top_rois$top_type)) {
    comp_top_rois <- top_rois$ROI[top_rois$component == comp & top_rois$top_type == top_type]
    for (net_name in names(networks)) {
      result <- test_enrichment(comp_top_rois, networks[[net_name]], all_rois)
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
enrichment_top$q_value <- p.adjust(enrichment_top$p_value, method = "fdr")
write.csv(enrichment_top, file.path(out_dir, "NMF_DD_only_network_enrichment_top_rois.csv"), row.names = FALSE)
write.csv(enrichment_top, file.path(out_dir, "NMF_network_enrichment_top_rois.csv"), row.names = FALSE)

cat("\nSignificant DD-only enrichments, hard assignment q < 0.05:\n")
sig_hard <- enrichment_hard %>% filter(q_value < 0.05) %>% arrange(q_value)
if (nrow(sig_hard) > 0) {
  print(sig_hard %>% select(component, network, n_overlap, enrichment_ratio, p_value, q_value), row.names = FALSE)
} else {
  cat("  No significant enrichments at FDR < 0.05\n")
}

cat("\nSignificant DD-only enrichments, top ROIs q < 0.05:\n")
sig_top <- enrichment_top %>% filter(q_value < 0.05) %>% arrange(q_value)
if (nrow(sig_top) > 0) {
  print(sig_top %>% select(component, top_type, network, n_overlap, enrichment_ratio, p_value, q_value), row.names = FALSE)
} else {
  cat("  No significant enrichments at FDR < 0.05\n")
}

cat("\nDD-only network enrichment complete.\n")
