# ============================================================
# 基于 318x2 T-map 矩阵对 ROI 进行聚类
# 识别不同的发育差异模式，并评估不同 K 值的聚类质量
# ============================================================

cat("T-map 聚类分析\n")
cat("==============\n\n")

packages <- c("cluster", "factoextra", "ggplot2", "dplyr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

library(cluster)
library(factoextra)
library(ggplot2)
library(dplyr)

K_main <- 5
K_range <- 2:8

tmap_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Clustering_5"
output_dir <- file.path(tmap_dir, "clustering")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 1. 读取数据
# ============================================================
cat("读取 T-map 数据...\n")
tmap <- read.csv(file.path(tmap_dir, "DGLM_DK318_degree_Tmap_TD_minus_DD_318x2.csv"))
full_results <- read.csv(file.path(tmap_dir, "DGLM_DK318_degree_Tmap_TD_minus_DD_full.csv"))

tmap_ok <- tmap %>%
  left_join(full_results %>% select(feature, fit_status), by = "feature") %>%
  filter(fit_status == "ok") %>%
  filter(!is.na(t_child_TD_minus_DD), !is.na(t_adult_TD_minus_DD))

cat("总 ROI 数:", nrow(tmap), "\n")
cat("用于聚类的 ROI 数:", nrow(tmap_ok), "\n\n")

X <- as.matrix(tmap_ok[, c("t_child_TD_minus_DD", "t_adult_TD_minus_DD")])
rownames(X) <- tmap_ok$feature

cat("聚类矩阵维度:", paste(dim(X), collapse = " x "), "\n\n")

if (!(K_main %in% K_range)) {
  stop("K_main 必须在 K_range 中。")
}

# ============================================================
# 2. K 选择分析
# ============================================================
cat("K 选择分析 (K = ", min(K_range), " 到 ", max(K_range), ")...\n", sep = "")
set.seed(123)

# 执行 K-means 聚类
all_kmeans <- list()
for (k in K_range) {
  cat("  K =", k, "\n")
  all_kmeans[[paste0("K", k)]] <- kmeans(X, centers = k, nstart = 50, iter.max = 100)
}

# 计算距离矩阵
d <- dist(X)

# 提取聚类质量指标
wss <- numeric(length(K_range))
bss <- numeric(length(K_range))
tss <- numeric(length(K_range))
sil_scores <- numeric(length(K_range))

for (i in seq_along(K_range)) {
  k <- K_range[i]
  km <- all_kmeans[[paste0("K", k)]]
  wss[i] <- km$tot.withinss
  bss[i] <- km$betweenss
  tss[i] <- km$totss
  sil_scores[i] <- mean(silhouette(km$cluster, d)[, 3])
}

bss_tss_ratio <- bss / tss

ch_index <- numeric(length(K_range))
for (i in seq_along(K_range)) {
  k <- K_range[i]
  ch_index[i] <- (bss[i] / (k - 1)) / (wss[i] / (nrow(X) - k))
}

calc_db <- function(km, X) {
  centers <- km$centers
  k <- nrow(centers)
  scatters <- sapply(seq_len(k), function(i) {
    Xi <- X[km$cluster == i, , drop = FALSE]
    if (nrow(Xi) == 0) return(0)
    center_mat <- matrix(centers[i, ], nrow = nrow(Xi), ncol = ncol(X), byrow = TRUE)
    mean(sqrt(rowSums((Xi - center_mat)^2)))
  })
  db_vals <- sapply(seq_len(k), function(i) {
    max(sapply(setdiff(seq_len(k), i), function(j) {
      center_dist <- sqrt(sum((centers[i, ] - centers[j, ])^2))
      if (center_dist == 0) return(Inf)
      (scatters[i] + scatters[j]) / center_dist
    }))
  })
  mean(db_vals)
}

db_index <- numeric(length(K_range))
size_cv <- numeric(length(K_range))
for (i in seq_along(K_range)) {
  k <- K_range[i]
  km <- all_kmeans[[paste0("K", k)]]
  db_index[i] <- calc_db(km, X)
  sizes <- as.numeric(table(km$cluster))
  size_cv[i] <- sd(sizes) / mean(sizes)
}

k_selection <- data.frame(
  K = K_range,
  WSS = wss,
  BSS = bss,
  TSS = tss,
  BSS_TSS_ratio = bss_tss_ratio,
  Silhouette = sil_scores,
  Calinski_Harabasz = ch_index,
  Davies_Bouldin = db_index,
  Size_CV = size_cv
)

write.csv(k_selection, file.path(output_dir, "k_selection_metrics.csv"), row.names = FALSE)
cat("\n聚类质量指标:\n")
print(k_selection)
cat("\n保存: k_selection_metrics.csv\n\n")

# ============================================================
# 3. K 选择可视化
# ============================================================
cat("生成 K 选择图表...\n")

png(file.path(output_dir, "k_selection_comprehensive.png"), width = 1600, height = 1000)
par(mfrow = c(2, 3), mar = c(4, 4, 3, 2))
plot(K_range, wss, type = "b", pch = 19, col = "blue", lwd = 2,
     xlab = "K", ylab = "WSS", main = "Elbow Method")
grid()
plot(K_range, sil_scores, type = "b", pch = 19, col = "red", lwd = 2,
     xlab = "K", ylab = "Average Silhouette", main = "Silhouette")
abline(h = 0.5, lty = 3, col = "gray50")
grid()
plot(K_range, bss_tss_ratio, type = "b", pch = 19, col = "green3", lwd = 2,
     xlab = "K", ylab = "BSS/TSS", main = "Explained Variance")
grid()
plot(K_range, ch_index, type = "b", pch = 19, col = "purple", lwd = 2,
     xlab = "K", ylab = "Calinski-Harabasz", main = "C-H Index")
grid()
plot(K_range, db_index, type = "b", pch = 19, col = "orange", lwd = 2,
     xlab = "K", ylab = "Davies-Bouldin", main = "D-B Index")
grid()
plot(K_range, size_cv, type = "b", pch = 19, col = "cyan4", lwd = 2,
     xlab = "K", ylab = "Size CV", main = "Cluster Size Balance")
grid()
dev.off()
cat("保存: k_selection_comprehensive.png\n")

png(file.path(output_dir, "elbow_plot.png"), width = 800, height = 600)
plot(K_range, wss, type = "b", pch = 19, col = "blue", lwd = 2,
     xlab = "K", ylab = "Within-cluster Sum of Squares", main = "Elbow Method")
grid()
dev.off()

png(file.path(output_dir, "silhouette_plot.png"), width = 800, height = 600)
plot(K_range, sil_scores, type = "b", pch = 19, col = "red", lwd = 2,
     xlab = "K", ylab = "Average Silhouette", main = "Silhouette Method")
abline(h = max(sil_scores), lty = 2, col = "gray")
grid()
dev.off()

png(file.path(output_dir, "bss_tss_ratio_plot.png"), width = 800, height = 600)
plot(K_range, bss_tss_ratio, type = "b", pch = 19, col = "green3", lwd = 2,
     xlab = "K", ylab = "BSS/TSS Ratio", main = "BSS/TSS Ratio")
grid()
dev.off()

# ============================================================
# 4. 综合评分
# ============================================================
normalize <- function(x) {
  if (max(x) == min(x)) return(rep(1, length(x)))
  (x - min(x)) / (max(x) - min(x))
}

k_selection_norm <- k_selection %>%
  mutate(
    norm_Silhouette = normalize(Silhouette),
    norm_BSS_TSS = normalize(BSS_TSS_ratio),
    norm_CH = normalize(Calinski_Harabasz),
    norm_DB = 1 - normalize(Davies_Bouldin),
    norm_Size = 1 - normalize(Size_CV),
    Composite_Score = (norm_Silhouette + norm_BSS_TSS + norm_CH + norm_DB + norm_Size) / 5
  )

write.csv(k_selection_norm, file.path(output_dir, "k_selection_normalized.csv"), row.names = FALSE)
best_k_composite <- k_selection_norm$K[which.max(k_selection_norm$Composite_Score)]
cat("综合评分推荐 K =", best_k_composite, "\n\n")

png(file.path(output_dir, "composite_score.png"), width = 1200, height = 600)
par(mfrow = c(1, 2))
plot(k_selection_norm$K, k_selection_norm$Composite_Score, type = "b", pch = 19,
     col = "darkblue", lwd = 2, xlab = "K", ylab = "Composite Score",
     main = "Composite Score")
abline(v = best_k_composite, lty = 2, col = "red", lwd = 2)
abline(v = K_main, lty = 2, col = "green3", lwd = 2)
legend("topright",
       legend = c(paste0("Best K=", best_k_composite), paste0("Current K=", K_main)),
       col = c("red", "green3"), lty = 2, lwd = 2)
grid()
barplot(t(as.matrix(k_selection_norm[, c("norm_Silhouette", "norm_BSS_TSS", "norm_CH", "norm_DB", "norm_Size")])),
        beside = TRUE, names.arg = k_selection_norm$K,
        col = c("red", "green3", "purple", "orange", "cyan4"),
        xlab = "K", ylab = "Normalized Score", main = "Normalized Metrics",
        legend.text = c("Sil", "BSS/TSS", "C-H", "D-B", "Size"),
        args.legend = list(x = "topright", cex = 0.75))
dev.off()
cat("保存: composite_score.png\n\n")

# ============================================================
# 5. 主分析
# ============================================================
cat(paste0("K=", K_main, " 主分析...\n"))
km_main <- all_kmeans[[paste0("K", K_main)]]

tmap_clustered <- tmap_ok
tmap_clustered$cluster <- km_main$cluster

write.csv(tmap_clustered, file.path(output_dir, paste0("degree_Tmap_cluster_K", K_main, ".csv")), row.names = FALSE)
cat("保存:", paste0("degree_Tmap_cluster_K", K_main, ".csv"), "\n")

centroids <- as.data.frame(km_main$centers)
centroids$cluster <- seq_len(K_main)
centroids$size <- km_main$size
write.csv(centroids, file.path(output_dir, paste0("degree_Tmap_cluster_K", K_main, "_centroids.csv")), row.names = FALSE)
cat("保存:", paste0("degree_Tmap_cluster_K", K_main, "_centroids.csv"), "\n\n")

cluster_summary <- tmap_clustered %>%
  group_by(cluster) %>%
  summarise(
    n = n(),
    mean_t_child = mean(t_child_TD_minus_DD),
    mean_t_adult = mean(t_adult_TD_minus_DD),
    median_t_child = median(t_child_TD_minus_DD),
    median_t_adult = median(t_adult_TD_minus_DD),
    sd_t_child = sd(t_child_TD_minus_DD),
    sd_t_adult = sd(t_adult_TD_minus_DD),
    mean_overall_TD_minus_DD = (mean_t_child + mean_t_adult) / 2,
    mean_developmental_change = mean_t_adult - mean_t_child,
    .groups = "drop"
  )

print(cluster_summary)
write.csv(cluster_summary, file.path(output_dir, paste0("degree_Tmap_cluster_K", K_main, "_summary.csv")), row.names = FALSE)
cat("保存:", paste0("degree_Tmap_cluster_K", K_main, "_summary.csv"), "\n\n")

cat("聚类解释:\n")
for (i in seq_len(K_main)) {
  t_child <- cluster_summary$mean_t_child[i]
  t_adult <- cluster_summary$mean_t_adult[i]
  pattern <- if (t_child > 0 & t_adult > 0) {
    "持续 TD > DD"
  } else if (t_child < 0 & t_adult < 0) {
    "持续 DD > TD"
  } else if (t_child < 0 & t_adult > 0) {
    "发育追赶"
  } else if (t_child > 0 & t_adult < 0) {
    "发育逆转"
  } else {
    "混合模式"
  }
  cat("\n聚类", i, "(n =", cluster_summary$n[i], "):\n")
  cat("  儿童组 t 均值:", round(t_child, 2), "\n")
  cat("  成人组 t 均值:", round(t_adult, 2), "\n")
  cat("  发育变化:", round(cluster_summary$mean_developmental_change[i], 2), "\n")
  cat("  模式:", pattern, "\n")
}
cat("\n")

# ============================================================
# 6. 聚类可视化
# ============================================================
base_colors <- c("red", "blue", "green3", "orange", "purple", "brown", "cyan4", "magenta")
cluster_colors <- setNames(base_colors[seq_len(K_main)], as.character(seq_len(K_main)))

png(file.path(output_dir, paste0("degree_Tmap_cluster_K", K_main, "_scatter.png")), width = 1200, height = 1000)
plot(tmap_clustered$t_child_TD_minus_DD, tmap_clustered$t_adult_TD_minus_DD,
     col = cluster_colors[as.character(tmap_clustered$cluster)], pch = 16, cex = 1.5,
     xlab = "Child TD-DD t", ylab = "Adult TD-DD t",
     main = paste0("DK318 MIND Degree T-map Clustering (K=", K_main, ")"))
abline(h = 0, v = 0, col = "gray30", lwd = 2, lty = 2)
abline(0, 1, col = "gray50", lwd = 1, lty = 3)
points(centroids$t_child_TD_minus_DD, centroids$t_adult_TD_minus_DD,
       col = "black", pch = 4, cex = 3, lwd = 3)
legend("topleft", legend = paste0("Cluster ", seq_len(K_main), " (n=", centroids$size, ")"),
       col = cluster_colors, pch = 16, cex = 1.2)
grid()
dev.off()
cat("保存:", paste0("degree_Tmap_cluster_K", K_main, "_scatter.png"), "\n")

sil_main <- silhouette(km_main$cluster, d)
png(file.path(output_dir, paste0("degree_Tmap_silhouette_K", K_main, ".png")), width = 900, height = 800)
plot(sil_main, col = cluster_colors,
     main = paste0("Silhouette Plot (K=", K_main, ", avg=", round(mean(sil_main[, 3]), 3), ")"))
dev.off()
cat("保存:", paste0("degree_Tmap_silhouette_K", K_main, ".png"), "\n\n")

sil_by_cluster <- data.frame(
  cluster = seq_len(K_main),
  mean_silhouette = as.numeric(tapply(sil_main[, 3], sil_main[, 1], mean)),
  min_silhouette = as.numeric(tapply(sil_main[, 3], sil_main[, 1], min)),
  n_negative = as.numeric(tapply(sil_main[, 3], sil_main[, 1], function(x) sum(x < 0))),
  n_total = as.numeric(table(sil_main[, 1]))
)
sil_by_cluster$pct_negative <- sil_by_cluster$n_negative / sil_by_cluster$n_total * 100
write.csv(sil_by_cluster, file.path(output_dir, paste0("degree_Tmap_silhouette_by_cluster_K", K_main, ".csv")), row.names = FALSE)

# ============================================================
# 7. 总结
# ============================================================
main_idx <- match(K_main, K_range)
cat("================\n")
cat("聚类分析完成\n")
cat("================\n\n")
cat("输出目录:", output_dir, "\n")
cat("当前 K =", K_main, "\n")
cat("综合评分推荐 K =", best_k_composite, "\n")
cat("平均 Silhouette:", round(sil_scores[main_idx], 3), "\n")
cat("BSS/TSS:", round(bss_tss_ratio[main_idx], 3), "\n")
cat("Calinski-Harabasz:", round(ch_index[main_idx], 1), "\n")
cat("Davies-Bouldin:", round(db_index[main_idx], 3), "\n\n")
cat("重点查看图表:\n")
cat("  - k_selection_comprehensive.png\n")
cat("  - composite_score.png\n")
cat("  -", paste0("degree_Tmap_cluster_K", K_main, "_scatter.png"), "\n")
cat("  -", paste0("degree_Tmap_silhouette_K", K_main, ".png"), "\n")
