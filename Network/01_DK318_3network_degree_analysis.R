# ============================================================
# DK318 MIND Degree 3-Network Aggregation Analysis
# Purpose: Aggregate ROI-level degree into DM/MD/Reading networks
#          Test Diagnosis × AgeGroup effects
# ============================================================

cat("Starting DK318 3-network degree analysis...\n\n")

# ============================================================
# 1. Load Packages
# ============================================================

cat("Loading packages...\n")
pkgs <- c("readxl", "dplyr", "tidyr", "tibble", "purrr", "car", "lme4", "lmerTest", "emmeans", "magrittr")
for (p in pkgs) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) {
    install.packages(p, repos = "https://cran.r-project.org")
    library(p, character.only = TRUE)
  }
}
library(readxl); library(dplyr); library(tidyr); library(tibble); library(purrr)
library(car); library(lme4); library(lmerTest); library(emmeans); library(magrittr)
options(contrasts = c("contr.sum", "contr.poly"))

# ============================================================
# 2. Define Paths
# ============================================================

demo_file <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/all_data_cqt_mean_1.5.xlsx"
mind_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_combat"
reading_file <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/brain_region_table_brainmaps/mean_gt1p0_union_reading_v4_v5/mean_gt1p0_union_reading_v4_v5_regions_used.csv"
md_file <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/brain_region_table_brainmaps/multiple_demand/multiple_demand_regions_used.csv"
dm_file <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network/brain_region_table_brainmaps/default/default_regions_used.csv"
out_root <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_3network_analysis"
analysis_dir <- file.path(out_root, "analysis")
fig_dir <- file.path(out_root, "figures")
dir.create(analysis_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
standardize_roi_before_network_mean <- FALSE

cat("Output directory:", analysis_dir, "\n\n")

# ============================================================
# 3. Helper Functions
# ============================================================

read_network <- function(fp, net) {
  x <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"feature" %in% names(x)) stop("Missing feature column: ", fp)
  tibble(feature = trimws(as.character(x$feature)), Network = net) %>%
    filter(!is.na(feature), feature != "") %>% distinct()
}

read_degree <- function(fp) {
  x <- read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE)
  if (!all(c("ROI", "degree") %in% names(x))) stop("Bad degree file: ", fp)
  v <- as.numeric(x$degree)
  names(v) <- as.character(x$ROI)
  return(v)
}

add_fdr_by_age <- function(x) {
  x %>% group_by(AgeGroup) %>% mutate(q_by_age = p.adjust(p.value, "fdr")) %>% ungroup()
}

model_status <- function(name, status, msg) {
  tibble(model_name = name, status = status, message = msg)
}

# Section 1 complete

# ============================================================
# 4. Load Demographics
# ============================================================

cat("Loading demographic data...\n")
raw_demo <- read_excel(demo_file, sheet = "Sheet1")
n_total <- nrow(raw_demo)

df <- raw_demo %>%
  mutate(
    original_project = as.character(`original-project`),
    id_old = as.character(id_old),
    subj_prefix = paste0(original_project, "_", id_old),
    file_base = paste0(subj_prefix, "_MIND_DK318_combat"),
    degree_file = file.path(mind_dir, paste0(file_base, "_degree.csv")),
    SubjectID = subj_prefix,
    Diagnosis = factor(ifelse(group_d_or_c == 0, "TD", "DD"), levels = c("TD", "DD")),
    AgeGroup = factor(ifelse(group_age == 1, "Adult", "Child"), levels = c("Child", "Adult")),
    Sex = factor(ifelse(sex == 1, "Male", "Female"), levels = c("Female", "Male")),
    has_file = file.exists(degree_file)
  ) %>%
  filter(!is.na(original_project), original_project != "", !is.na(id_old), id_old != "",
         !is.na(Diagnosis), !is.na(AgeGroup), !is.na(Sex), has_file) %>%
  select(SubjectID, subj_prefix, file_base, Diagnosis, AgeGroup, Sex, degree_file)

n_used <- nrow(df)
if (n_used == 0) stop("No subjects entered final analysis.")

all_files <- file.path(mind_dir, paste0(as.character(raw_demo$`original-project`), "_", 
                                        as.character(raw_demo$id_old), "_MIND_DK318_combat_degree.csv"))
n_with_file <- sum(file.exists(all_files), na.rm = TRUE)

cat("Total subjects:", n_total, "\n")
cat("With degree files:", n_with_file, "\n")
cat("Final analysis:", n_used, "\n\n")
print(with(df, table(Diagnosis, AgeGroup, Sex)))

write.csv(df %>% count(Diagnosis, AgeGroup, Sex, name = "n"), 
          file.path(analysis_dir, "QC_subject_cell_counts.csv"), row.names = FALSE)
write.csv(tibble(setting = c("standardize_roi_before_network_mean", "n_subjects_total", 
                              "n_subjects_with_degree_file", "n_subjects_used"),
                 value = c(as.character(standardize_roi_before_network_mean), n_total, n_with_file, n_used)),
          file.path(analysis_dir, "QC_analysis_settings.csv"), row.names = FALSE)

# ============================================================
# 5. Load Degree Matrix
# ============================================================

cat("\nLoading degree matrix...\n")
v1 <- read_degree(df$degree_file[1])
roi_names <- names(v1)
n_roi <- length(roi_names)
cat("Number of ROIs:", n_roi, "\n")

X <- matrix(NA_real_, nrow(df), n_roi, dimnames = list(df$SubjectID, roi_names))
X[1, ] <- v1[roi_names]

if (nrow(df) > 1) {
  for (i in 2:nrow(df)) {
    if (i %% 50 == 0 || i == nrow(df)) cat("  Reading subject", i, "/", nrow(df), "\n")
    v <- read_degree(df$degree_file[i])
    miss <- setdiff(roi_names, names(v))
    if (length(miss)) stop("ROI missing in ", df$degree_file[i], ": ", paste(head(miss, 10), collapse = ", "))
    X[i, ] <- v[roi_names]
  }
}

if (any(is.na(X))) stop("Degree matrix contains NA values.")
cat("Degree range: [", min(X), ", ", max(X), "]\n")

# Section 2 complete

# ============================================================
# 6. Load Networks and Compute Overlap
# ============================================================

cat("\nLoading network definitions...\n")
nets <- bind_rows(read_network(dm_file, "DM"), read_network(md_file, "MD"), read_network(reading_file, "Reading"))
all_rois <- sort(unique(nets$feature))

mem <- tibble(feature = all_rois) %>%
  left_join(nets %>% mutate(flag = TRUE) %>% pivot_wider(names_from = Network, values_from = flag, values_fill = FALSE), 
            by = "feature") %>%
  mutate(DM = ifelse(is.na(DM), FALSE, DM), MD = ifelse(is.na(MD), FALSE, MD), 
         Reading = ifelse(is.na(Reading), FALSE, Reading),
         in_DM = DM, in_MD = MD, in_Reading = Reading,
         n_networks = as.integer(in_DM) + as.integer(in_MD) + as.integer(in_Reading),
         membership_pattern = case_when(
           in_DM & !in_MD & !in_Reading ~ "DM_only",
           !in_DM & in_MD & !in_Reading ~ "MD_only",
           !in_DM & !in_MD & in_Reading ~ "Reading_only",
           in_DM & in_MD & !in_Reading ~ "DM+MD",
           in_DM & !in_MD & in_Reading ~ "Reading+DM",
           !in_DM & in_MD & in_Reading ~ "Reading+MD",
           in_DM & in_MD & in_Reading ~ "DM+MD+Reading",
           TRUE ~ "None")) %>%
  select(feature, in_DM, in_MD, in_Reading, n_networks, membership_pattern)

write.csv(mem, file.path(analysis_dir, "network_membership_3net.csv"), row.names = FALSE)

sets <- list(DM = mem$feature[mem$in_DM], MD = mem$feature[mem$in_MD], Reading = mem$feature[mem$in_Reading])
cat("Network sizes: DM =", length(sets$DM), ", MD =", length(sets$MD), ", Reading =", length(sets$Reading), "\n")

qc_roi <- map_dfr(names(sets), function(nm) {
  fd <- intersect(sets[[nm]], colnames(X))
  ms <- setdiff(sets[[nm]], colnames(X))
  if (!length(fd)) stop("No ROI found for network: ", nm)
  tibble(Network = nm, n_roi_in_network_file = length(sets[[nm]]), 
         n_roi_found_in_degree_matrix = length(fd), n_roi_missing_from_degree_matrix = length(ms),
         missing_roi_list = if (length(ms)) paste(ms, collapse = "; ") else "")
})
write.csv(qc_roi, file.path(analysis_dir, "QC_network_roi_coverage.csv"), row.names = FALSE)

ov <- function(a, b) length(intersect(sets[[a]], sets[[b]]))
overlap_summary <- bind_rows(
  tibble(comparison = "Reading_vs_MD", n_overlap = ov("Reading", "MD"), n_network_1 = length(sets$Reading), 
         n_network_2 = length(sets$MD), percent_of_network_1 = 100 * ov("Reading", "MD") / length(sets$Reading), 
         percent_of_network_2 = 100 * ov("Reading", "MD") / length(sets$MD)),
  tibble(comparison = "Reading_vs_DM", n_overlap = ov("Reading", "DM"), n_network_1 = length(sets$Reading), 
         n_network_2 = length(sets$DM), percent_of_network_1 = 100 * ov("Reading", "DM") / length(sets$Reading), 
         percent_of_network_2 = 100 * ov("Reading", "DM") / length(sets$DM)),
  tibble(comparison = "MD_vs_DM", n_overlap = ov("MD", "DM"), n_network_1 = length(sets$MD), 
         n_network_2 = length(sets$DM), percent_of_network_1 = 100 * ov("MD", "DM") / length(sets$MD), 
         percent_of_network_2 = 100 * ov("MD", "DM") / length(sets$DM)),
  tibble(comparison = "All_three", n_overlap = length(Reduce(intersect, sets)), n_network_1 = length(sets$DM), 
         n_network_2 = NA, percent_of_network_1 = 100 * length(Reduce(intersect, sets)) / length(sets$DM), 
         percent_of_network_2 = NA)
)
write.csv(overlap_summary, file.path(analysis_dir, "network_overlap_summary_3net.csv"), row.names = FALSE)

# ============================================================
# 7. Compute Network Degree
# ============================================================

cat("\nComputing network-level degree...\n")
Xu <- if (standardize_roi_before_network_mean) {
  z <- as.matrix(scale(X))
  z[is.na(z)] <- 0
  z
} else X

keep <- lapply(sets, intersect, y = colnames(Xu))
wide <- df %>%
  mutate(DM_degree = rowMeans(Xu[, keep$DM, drop = FALSE]),
         MD_degree = rowMeans(Xu[, keep$MD, drop = FALSE]),
         Reading_degree = rowMeans(Xu[, keep$Reading, drop = FALSE])) %>%
  select(SubjectID, subj_prefix, file_base, Diagnosis, AgeGroup, Sex, DM_degree, MD_degree, Reading_degree)

long <- wide %>%
  pivot_longer(c(DM_degree, MD_degree, Reading_degree), names_to = "Network", values_to = "Degree") %>%
  mutate(Network = dplyr::recode(Network, DM_degree = "DM", MD_degree = "MD", Reading_degree = "Reading"),
         Network = factor(Network, levels = c("DM", "MD", "Reading"))) %>%
  select(SubjectID, subj_prefix, file_base, Diagnosis, AgeGroup, Sex, Network, Degree)

write.csv(wide, file.path(analysis_dir, "DK318_3network_degree_wide.csv"), row.names = FALSE)
write.csv(long, file.path(analysis_dir, "DK318_3network_degree_long.csv"), row.names = FALSE)

# Section 3 complete

# ============================================================
# 8. Mixed-Effects Model
# ============================================================

cat("\nFitting mixed-effects model...\n")
status_list <- list()

fit_mixed <- tryCatch(lmer(Degree ~ Network * Diagnosis * AgeGroup + Sex + (1 | SubjectID), 
                           data = long, REML = FALSE), error = function(e) e)

if (inherits(fit_mixed, "error")) {
  cat("ERROR: Mixed model failed.\n")
  status_list[[1]] <- model_status("mixed_model_3network", "failed", conditionMessage(fit_mixed))
} else {
  cat("Mixed model fitted successfully.\n")
  status_list[[1]] <- model_status("mixed_model_3network", "success", "ok")
  
  a <- as.data.frame(anova(fit_mixed, type = 3)) %>% rownames_to_column("term")
  write.csv(a, file.path(analysis_dir, "mixed_model_3network_type3_anova.csv"), row.names = FALSE)
  
  writeLines(capture.output(summary(fit_mixed)), file.path(analysis_dir, "mixed_model_3network_summary.txt"))
  
  emm <- as.data.frame(emmeans(fit_mixed, ~ Diagnosis * AgeGroup | Network))
  write.csv(emm, file.path(analysis_dir, "mixed_model_3network_emmeans_adjusted_means.csv"), row.names = FALSE)
  
  c1 <- as.data.frame(contrast(emmeans(fit_mixed, ~ Diagnosis | Network * AgeGroup), 
                                method = list(TD_minus_DD = c(1, -1)), by = c("Network", "AgeGroup"))) %>%
    mutate(q_all = p.adjust(p.value, "fdr")) %>% add_fdr_by_age() %>%
    select(Network, AgeGroup, contrast, estimate, SE, df, t.ratio, p.value, q_by_age, q_all)
  write.csv(c1, file.path(analysis_dir, "mixed_model_3network_TD_minus_DD_by_network_age.csv"), row.names = FALSE)
  
  c2 <- as.data.frame(contrast(emmeans(fit_mixed, ~ Diagnosis * AgeGroup | Network), 
                                method = list("Adult(TD-DD) - Child(TD-DD)" = c(-1, 1, 1, -1)), by = "Network")) %>%
    mutate(q_value = p.adjust(p.value, "fdr")) %>%
    select(Network, contrast, estimate, SE, df, t.ratio, p.value, q_value)
  write.csv(c2, file.path(analysis_dir, "mixed_model_3network_developmental_contrast.csv"), row.names = FALSE)
}

# Section 4 complete

# ============================================================
# 9. Per-Network Models
# ============================================================

cat("\nFitting per-network linear models...\n")
pa <- list(); pe <- list(); pt <- list(); pd <- list(); idx <- 2

for (net in c("DM", "MD", "Reading")) {
  cat("  Network:", net, "\n")
  fit <- tryCatch(lm(Degree ~ Diagnosis * AgeGroup + Sex, data = filter(long, Network == net)), error = function(e) e)
  
  if (inherits(fit, "error")) {
    status_list[[idx]] <- model_status(paste0("per_network_", net), "failed", conditionMessage(fit))
    idx <- idx + 1
    next
  }
  
  status_list[[idx]] <- model_status(paste0("per_network_", net), "success", "ok")
  idx <- idx + 1
  
  pa[[net]] <- as.data.frame(car::Anova(fit, type = 3)) %>% rownames_to_column("term") %>% 
    filter(term != "(Intercept)") %>% transmute(Network = net, term, statistic = `F value`, df = Df, p.value = `Pr(>F)`)
  pe[[net]] <- as.data.frame(emmeans(fit, ~ Diagnosis * AgeGroup)) %>% mutate(Network = net) %>% 
    select(Network, Diagnosis, AgeGroup, emmean, SE, df, lower.CL, upper.CL)
  pt[[net]] <- as.data.frame(contrast(emmeans(fit, ~ Diagnosis | AgeGroup), method = list(TD_minus_DD = c(1, -1)), 
                                      by = "AgeGroup")) %>% mutate(Network = net) %>% 
    select(Network, AgeGroup, contrast, estimate, SE, df, t.ratio, p.value)
  pd[[net]] <- as.data.frame(contrast(emmeans(fit, ~ Diagnosis * AgeGroup), 
                                      method = list("Adult(TD-DD) - Child(TD-DD)" = c(-1, 1, 1, -1)))) %>% 
    mutate(Network = net) %>% select(Network, contrast, estimate, SE, df, t.ratio, p.value)
}

pa <- bind_rows(pa)
if (nrow(pa)) {
  pa$q_value <- NA_real_
  ix <- pa$term == "Diagnosis:AgeGroup"
  if (any(ix)) pa$q_value[ix] <- p.adjust(pa$p.value[ix], "fdr")
  pa$q_all <- p.adjust(pa$p.value, "fdr")
  write.csv(pa, file.path(analysis_dir, "per_network_lm_type3_anova.csv"), row.names = FALSE)
}

pe <- bind_rows(pe)
if (nrow(pe)) write.csv(pe, file.path(analysis_dir, "per_network_emmeans_adjusted_means.csv"), row.names = FALSE)

pt <- bind_rows(pt)
if (nrow(pt)) {
  pt <- pt %>% mutate(q_all = p.adjust(p.value, "fdr")) %>% add_fdr_by_age() %>% 
    select(Network, AgeGroup, contrast, estimate, SE, df, t.ratio, p.value, q_by_age, q_all)
  write.csv(pt, file.path(analysis_dir, "per_network_TD_minus_DD_by_age.csv"), row.names = FALSE)
}

pd <- bind_rows(pd)
if (nrow(pd)) {
  pd <- pd %>% mutate(q_value = p.adjust(p.value, "fdr")) %>% 
    select(Network, contrast, estimate, SE, df, t.ratio, p.value, q_value)
  write.csv(pd, file.path(analysis_dir, "per_network_developmental_contrast.csv"), row.names = FALSE)
}

# ============================================================
# 10. Save Status and Complete
# ============================================================

write.csv(bind_rows(status_list), file.path(analysis_dir, "QC_model_status.csv"), row.names = FALSE)

cat("\n==============================\n")
cat("ANALYSIS COMPLETE\n")
cat("==============================\n")
cat("Outputs saved to:", analysis_dir, "\n")

