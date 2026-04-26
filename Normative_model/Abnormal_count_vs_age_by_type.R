rm(list=ls())

library(ggplot2)

results_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Results/Individual_area_volume"
clinical_path <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Datasets/Datasets-diseases/Clinical_vars_control.csv"
out_dir <- file.path(results_dir, "_summary_abnormality")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
csv_path <- file.path(out_dir, "abnormal_count_by_age_by_type.csv")

z_thr_thickness <- 2
z_thr_area <- 2
z_thr_volume <- 2

subj_summary <- NULL

gather_type <- function(type) {
  dirs <- character()
  if (type == "thickness") {
    dirs <- c("lh.aparc.thickness.table", "rh.aparc.thickness.table")
  } else if (type == "area") {
    dirs <- c("lh.aparc.area.table", "rh.aparc.area.table")
  } else if (type == "volume") {
    dirs <- c("lh.aparc.volume.table", "rh.aparc.volume.table", "aseg.vol.table")
  } else {
    stop("unknown type")
  }

  files <- unlist(lapply(dirs, function(d) {
    list.files(file.path(results_dir, d), pattern = "_loop_our_model_individual\\.rds$", full.names = TRUE)
  }))
  files <- files[file.exists(files)]
  if (length(files) == 0) return(NULL)
  files
}

if (file.exists(csv_path)) {
  subj_summary <- read.csv(csv_path, header = TRUE, stringsAsFactors = FALSE)
} else {
  clinical <- read.csv(clinical_path, header = TRUE, stringsAsFactors = FALSE)
  clinical$ID <- as.character(clinical$ID)
  clinical$path_key <- paste0(clinical$Freesurfer_Path2, clinical$Freesurfer_Path3)
  clinical$Age <- as.numeric(clinical$Age)

  ids <- clinical$ID
  id_set <- unique(ids)
  path_to_id <- setNames(clinical$ID, clinical$path_key)
  age_by_id <- setNames(clinical$Age, clinical$ID)

  map_rowname_to_id <- function(rn) {
    rn <- as.character(rn)
    out <- rn
    not_id <- !(out %in% id_set)
    if (any(not_id)) {
      out[not_id] <- unname(path_to_id[rn[not_id]])
    }
    out
  }

  accumulate_counts <- function(files, z_thr, type) {
    if (length(files) == 0) return(NULL)
    counts_total <- setNames(integer(length(id_set)), id_set)
    counts_abn <- setNames(integer(length(id_set)), id_set)

    for (k in seq_along(files)) {
      f <- files[k]
      x <- readRDS(f)
      if (is.null(x$Zscore) || length(x$Zscore) == 0) next
      zdf <- x$Zscore[[1]]
      if (is.null(zdf) || !("Z_score" %in% colnames(zdf))) next

      rn <- rownames(zdf)
      subj_id <- map_rowname_to_id(rn)
      ok <- !is.na(subj_id) & (subj_id %in% id_set)
      if (!any(ok)) next

      subj_id <- subj_id[ok]
      z <- as.numeric(zdf$Z_score[ok])
      if (length(z) != length(subj_id)) next

      counts_total[subj_id] <- counts_total[subj_id] + 1L
      counts_abn[subj_id] <- counts_abn[subj_id] + as.integer(abs(z) >= z_thr)

      if (k %% 50 == 0) {
        cat(type, "processed", k, "of", length(files), "\n")
      }
    }

    out <- data.frame(
      type = type,
      ID = names(counts_total),
      Age = unname(age_by_id[names(counts_total)]),
      n_features = as.integer(counts_total),
      n_abnormal = as.integer(counts_abn),
      stringsAsFactors = FALSE
    )
    out <- out[!is.na(out$Age) & out$n_features > 0, ]
    out$prop_abnormal <- out$n_abnormal / out$n_features
    out
  }

  files_thk <- gather_type("thickness")
  files_area <- gather_type("area")
  files_vol <- gather_type("volume")

  subj_summary <- rbind(
    accumulate_counts(files_thk, z_thr_thickness, "thickness"),
    accumulate_counts(files_area, z_thr_area, "area"),
    accumulate_counts(files_vol, z_thr_volume, "volume")
  )

  if (is.null(subj_summary) || nrow(subj_summary) == 0) stop("No subject summaries computed")
  write.csv(subj_summary, csv_path, row.names = FALSE)
}

subj_summary$type <- factor(subj_summary$type, levels = c("thickness", "area", "volume"))
type_colors <- c(thickness = "#1b9e77", area = "#7570b3", volume = "#d95f02")

subj_summary$Age_bin <- round(subj_summary$Age)

bin_mean_n <- aggregate(n_abnormal ~ type + Age_bin, data = subj_summary, FUN = mean)
bin_sd_n <- aggregate(n_abnormal ~ type + Age_bin, data = subj_summary, FUN = sd)
bin_n_n <- aggregate(n_abnormal ~ type + Age_bin, data = subj_summary, FUN = length)
bin_stats_n <- merge(bin_mean_n, bin_sd_n, by = c("type", "Age_bin"), suffixes = c("_mean", "_sd"))
bin_stats_n <- merge(bin_stats_n, bin_n_n, by = c("type", "Age_bin"))
colnames(bin_stats_n)[colnames(bin_stats_n) == "n_abnormal"] <- "n"
bin_stats_n$se <- bin_stats_n$n_abnormal_sd / sqrt(bin_stats_n$n)
bin_stats_n$ci <- 1.96 * bin_stats_n$se
bin_stats_n$ymin <- pmax(bin_stats_n$n_abnormal_mean - bin_stats_n$ci, 0)
bin_stats_n$ymax <- bin_stats_n$n_abnormal_mean + bin_stats_n$ci

p1 <- ggplot(subj_summary, aes(x = Age, y = n_abnormal)) +
  geom_point(aes(color = type), alpha = 0.18, size = 1, position = position_jitter(width = 0.12, height = 0.15)) +
  geom_ribbon(
    data = bin_stats_n,
    aes(x = Age_bin, ymin = ymin, ymax = ymax, fill = type),
    alpha = 0.15,
    color = NA,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = bin_stats_n,
    aes(x = Age_bin, y = n_abnormal_mean, color = type),
    linewidth = 1.1,
    inherit.aes = FALSE
  ) +
  facet_wrap(~type, scales = "free_y") +
  scale_color_manual(values = type_colors) +
  scale_fill_manual(values = type_colors) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(x = "Age", y = "Abnormal region count (|Z| ≥ 2)", title = "Abnormal region count vs Age")

ggsave(file.path(out_dir, "abnormal_count_vs_age_by_type.png"), p1, width = 10, height = 4.5, dpi = 200)

bin_mean_p <- aggregate(prop_abnormal ~ type + Age_bin, data = subj_summary, FUN = mean)
bin_sd_p <- aggregate(prop_abnormal ~ type + Age_bin, data = subj_summary, FUN = sd)
bin_n_p <- aggregate(prop_abnormal ~ type + Age_bin, data = subj_summary, FUN = length)
bin_stats_p <- merge(bin_mean_p, bin_sd_p, by = c("type", "Age_bin"), suffixes = c("_mean", "_sd"))
bin_stats_p <- merge(bin_stats_p, bin_n_p, by = c("type", "Age_bin"))
colnames(bin_stats_p)[colnames(bin_stats_p) == "prop_abnormal"] <- "n"
bin_stats_p$se <- bin_stats_p$prop_abnormal_sd / sqrt(bin_stats_p$n)
bin_stats_p$ci <- 1.96 * bin_stats_p$se
bin_stats_p$ymin <- pmax(bin_stats_p$prop_abnormal_mean - bin_stats_p$ci, 0)
bin_stats_p$ymax <- pmin(bin_stats_p$prop_abnormal_mean + bin_stats_p$ci, 1)

p2 <- ggplot(subj_summary, aes(x = Age, y = prop_abnormal)) +
  geom_point(aes(color = type), alpha = 0.18, size = 1, position = position_jitter(width = 0.12, height = 0.005)) +
  geom_ribbon(
    data = bin_stats_p,
    aes(x = Age_bin, ymin = ymin, ymax = ymax, fill = type),
    alpha = 0.15,
    color = NA,
    inherit.aes = FALSE
  ) +
  geom_line(
    data = bin_stats_p,
    aes(x = Age_bin, y = prop_abnormal_mean, color = type),
    linewidth = 1.1,
    inherit.aes = FALSE
  ) +
  facet_wrap(~type, scales = "free_y") +
  scale_color_manual(values = type_colors) +
  scale_fill_manual(values = type_colors) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "none"
  ) +
  labs(x = "Age", y = "Abnormal proportion (|Z| ≥ 2)", title = "Abnormal proportion vs Age")

ggsave(file.path(out_dir, "abnormal_prop_vs_age_by_type.png"), p2, width = 10, height = 4.5, dpi = 200)

cat(out_dir, "\n")
