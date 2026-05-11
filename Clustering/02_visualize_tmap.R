# ============================================================
# T-map 结果可视化和初步检查
# 用于快速查看 T-map 结果的分布和模式
# ============================================================

cat("T-map 结果可视化\n")
cat("================\n\n")

# 加载包
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# 设置路径
tmap_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Clustering_5"
output_dir <- file.path(tmap_dir, "figures")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# 读取结果
cat("读取 T-map 结果...\n")
full_results <- read.csv(file.path(tmap_dir, "DGLM_DK318_degree_Tmap_TD_minus_DD_full.csv"))
tmap_318x2 <- read.csv(file.path(tmap_dir, "DGLM_DK318_degree_Tmap_TD_minus_DD_318x2.csv"))
summary_stats <- read.csv(file.path(tmap_dir, "DGLM_DK318_degree_Tmap_TD_minus_DD_summary.csv"))

cat("数据维度:\n")
cat("  完整结果:", nrow(full_results), "行 ×", ncol(full_results), "列\n")
cat("  318×2 矩阵:", nrow(tmap_318x2), "行 ×", ncol(tmap_318x2), "列\n\n")

# ============================================================
# 1. 拟合状态统计
# ============================================================
cat("1. 拟合状态统计\n")
cat("----------------\n")
fit_status_table <- table(full_results$fit_status)
print(fit_status_table)
cat("\n")

fit_status_pct <- prop.table(fit_status_table) * 100
cat("百分比:\n")
print(round(fit_status_pct, 1))
cat("\n")

# 保存拟合状态图
png(file.path(output_dir, "01_fit_status.png"), width = 800, height = 600)
barplot(fit_status_table, 
        main = "DGLM 拟合状态分布",
        xlab = "状态", ylab = "ROI 数量",
        col = c("ok" = "green3", "low_n" = "orange", 
                "near_zero_variance" = "red", 
                "empty_diagnosis_age_cell" = "purple")[names(fit_status_table)],
        las = 2)
dev.off()
cat("保存: 01_fit_status.png\n\n")

# ============================================================
# 2. T 值分布
# ============================================================
cat("2. T 值分布统计\n")
cat("----------------\n")

# 只使用成功拟合的 ROI
ok_results <- full_results %>% filter(fit_status == "ok")
cat("成功拟合的 ROI 数量:", nrow(ok_results), "\n\n")

# 儿童组 t 值
cat("儿童组 TD-DD t 值:\n")
cat("  范围:", round(range(ok_results$t_child_TD_minus_DD, na.rm = TRUE), 2), "\n")
cat("  均值:", round(mean(ok_results$t_child_TD_minus_DD, na.rm = TRUE), 2), "\n")
cat("  中位数:", round(median(ok_results$t_child_TD_minus_DD, na.rm = TRUE), 2), "\n")
cat("  正值数量:", sum(ok_results$t_child_TD_minus_DD > 0, na.rm = TRUE), "\n")
cat("  负值数量:", sum(ok_results$t_child_TD_minus_DD < 0, na.rm = TRUE), "\n\n")

# 成人组 t 值
cat("成人组 TD-DD t 值:\n")
cat("  范围:", round(range(ok_results$t_adult_TD_minus_DD, na.rm = TRUE), 2), "\n")
cat("  均值:", round(mean(ok_results$t_adult_TD_minus_DD, na.rm = TRUE), 2), "\n")
cat("  中位数:", round(median(ok_results$t_adult_TD_minus_DD, na.rm = TRUE), 2), "\n")
cat("  正值数量:", sum(ok_results$t_adult_TD_minus_DD > 0, na.rm = TRUE), "\n")
cat("  负值数量:", sum(ok_results$t_adult_TD_minus_DD < 0, na.rm = TRUE), "\n\n")

# T 值分布直方图
png(file.path(output_dir, "02_t_value_distributions.png"), width = 1200, height = 600)
par(mfrow = c(1, 2))

hist(ok_results$t_child_TD_minus_DD, breaks = 30,
     main = "儿童组 TD-DD t 值分布",
     xlab = "t 值", ylab = "频数",
     col = "skyblue", border = "white")
abline(v = 0, col = "red", lwd = 2, lty = 2)

hist(ok_results$t_adult_TD_minus_DD, breaks = 30,
     main = "成人组 TD-DD t 值分布",
     xlab = "t 值", ylab = "频数",
     col = "lightcoral", border = "white")
abline(v = 0, col = "red", lwd = 2, lty = 2)

dev.off()
cat("保存: 02_t_value_distributions.png\n\n")

# ============================================================
# 3. T-map 散点图（核心可视化）
# ============================================================
cat("3. T-map 散点图\n")
cat("----------------\n")

# 合并数据
tmap_plot_data <- tmap_318x2 %>%
  left_join(full_results %>% select(feature, fit_status, 
                                     q_child_TD_minus_DD, q_adult_TD_minus_DD),
            by = "feature") %>%
  filter(fit_status == "ok")

# 确定象限
tmap_plot_data <- tmap_plot_data %>%
  mutate(
    quadrant = case_when(
      t_child_TD_minus_DD > 0 & t_adult_TD_minus_DD > 0 ~ "I: 持续 TD>DD",
      t_child_TD_minus_DD < 0 & t_adult_TD_minus_DD > 0 ~ "II: 发育追赶",
      t_child_TD_minus_DD < 0 & t_adult_TD_minus_DD < 0 ~ "III: 持续 DD>TD",
      t_child_TD_minus_DD > 0 & t_adult_TD_minus_DD < 0 ~ "IV: 发育逆转",
      TRUE ~ "原点"
    ),
    sig_both = q_child_TD_minus_DD < 0.05 & q_adult_TD_minus_DD < 0.05,
    sig_child_only = q_child_TD_minus_DD < 0.05 & q_adult_TD_minus_DD >= 0.05,
    sig_adult_only = q_child_TD_minus_DD >= 0.05 & q_adult_TD_minus_DD < 0.05,
    sig_neither = q_child_TD_minus_DD >= 0.05 & q_adult_TD_minus_DD >= 0.05
  )

# 象限统计
cat("象限分布:\n")
quadrant_table <- table(tmap_plot_data$quadrant)
print(quadrant_table)
cat("\n")

cat("显著性统计:\n")
cat("  两组都显著:", sum(tmap_plot_data$sig_both, na.rm = TRUE), "\n")
cat("  仅儿童组显著:", sum(tmap_plot_data$sig_child_only, na.rm = TRUE), "\n")
cat("  仅成人组显著:", sum(tmap_plot_data$sig_adult_only, na.rm = TRUE), "\n")
cat("  都不显著:", sum(tmap_plot_data$sig_neither, na.rm = TRUE), "\n\n")

# 基础散点图
png(file.path(output_dir, "03_tmap_scatter_basic.png"), width = 1000, height = 1000)
plot(tmap_plot_data$t_child_TD_minus_DD, tmap_plot_data$t_adult_TD_minus_DD,
     xlab = "儿童组 TD-DD t 值", ylab = "成人组 TD-DD t 值",
     main = "DK318 MIND Degree T-map (318 ROIs)",
     pch = 16, col = rgb(0, 0, 0, 0.5), cex = 1.2)
abline(h = 0, v = 0, col = "red", lwd = 2, lty = 2)
abline(0, 1, col = "blue", lwd = 2, lty = 2)
grid()
legend("topleft", 
       legend = c("x = 0 (儿童无差异)", "y = 0 (成人无差异)", "y = x (一致性)"),
       col = c("red", "red", "blue"), lty = 2, lwd = 2)
dev.off()
cat("保存: 03_tmap_scatter_basic.png\n\n")

# 按象限着色的散点图
png(file.path(output_dir, "04_tmap_scatter_quadrants.png"), width = 1000, height = 1000)
quadrant_colors <- c("I: 持续 TD>DD" = "green3", 
                     "II: 发育追赶" = "blue",
                     "III: 持续 DD>TD" = "red",
                     "IV: 发育逆转" = "orange")
plot(tmap_plot_data$t_child_TD_minus_DD, tmap_plot_data$t_adult_TD_minus_DD,
     xlab = "儿童组 TD-DD t 值", ylab = "成人组 TD-DD t 值",
     main = "T-map 按象限着色",
     pch = 16, col = quadrant_colors[tmap_plot_data$quadrant], cex = 1.2)
abline(h = 0, v = 0, col = "gray30", lwd = 2, lty = 2)
abline(0, 1, col = "gray50", lwd = 1, lty = 3)
grid()
legend("topleft", legend = names(quadrant_colors), 
       col = quadrant_colors, pch = 16, cex = 0.9)
dev.off()
cat("保存: 04_tmap_scatter_quadrants.png\n\n")

# 按显著性着色的散点图
png(file.path(output_dir, "05_tmap_scatter_significance.png"), width = 1000, height = 1000)
sig_colors <- ifelse(tmap_plot_data$sig_both, "red",
                     ifelse(tmap_plot_data$sig_child_only, "blue",
                            ifelse(tmap_plot_data$sig_adult_only, "green",
                                   rgb(0.7, 0.7, 0.7, 0.5))))
plot(tmap_plot_data$t_child_TD_minus_DD, tmap_plot_data$t_adult_TD_minus_DD,
     xlab = "儿童组 TD-DD t 值", ylab = "成人组 TD-DD t 值",
     main = "T-map 按显著性着色 (FDR < 0.05)",
     pch = 16, col = sig_colors, cex = 1.2)
abline(h = 0, v = 0, col = "black", lwd = 2, lty = 2)
abline(0, 1, col = "gray50", lwd = 1, lty = 3)
grid()
legend("topleft", 
       legend = c("两组都显著", "仅儿童组显著", "仅成人组显著", "都不显著"),
       col = c("red", "blue", "green", "gray"), pch = 16, cex = 0.9)
dev.off()
cat("保存: 05_tmap_scatter_significance.png\n\n")

# ============================================================
# 4. 相关性分析
# ============================================================
cat("4. 儿童组与成人组 t 值相关性\n")
cat("------------------------------\n")
cor_result <- cor.test(tmap_plot_data$t_child_TD_minus_DD, 
                       tmap_plot_data$t_adult_TD_minus_DD)
cat("Pearson 相关系数:", round(cor_result$estimate, 3), "\n")
cat("p 值:", format.pval(cor_result$p.value, digits = 3), "\n\n")

# ============================================================
# 5. 显著 ROI 列表
# ============================================================
cat("5. 显著 ROI 列表\n")
cat("----------------\n")

# 两组都显著的 ROI
sig_both_rois <- tmap_plot_data %>%
  filter(sig_both) %>%
  arrange(desc(abs(t_child_TD_minus_DD) + abs(t_adult_TD_minus_DD))) %>%
  select(feature, t_child_TD_minus_DD, t_adult_TD_minus_DD, 
         q_child_TD_minus_DD, q_adult_TD_minus_DD)

if (nrow(sig_both_rois) > 0) {
  cat("\n两组都显著的 ROI (前10个):\n")
  print(head(sig_both_rois, 10))
  write.csv(sig_both_rois, 
            file.path(output_dir, "significant_ROIs_both_groups.csv"),
            row.names = FALSE)
  cat("\n保存: significant_ROIs_both_groups.csv\n")
}

# ============================================================
# 6. 汇总报告
# ============================================================
cat("\n")
cat("================\n")
cat("可视化完成\n")
cat("================\n\n")

cat("生成的图表:\n")
cat("  1. 01_fit_status.png - 拟合状态分布\n")
cat("  2. 02_t_value_distributions.png - t 值分布直方图\n")
cat("  3. 03_tmap_scatter_basic.png - 基础散点图\n")
cat("  4. 04_tmap_scatter_quadrants.png - 按象限着色\n")
cat("  5. 05_tmap_scatter_significance.png - 按显著性着色\n\n")

cat("输出目录:", output_dir, "\n\n")

cat("关键发现:\n")
cat("  - 成功拟合:", nrow(ok_results), "/", nrow(full_results), "ROI\n")
cat("  - 儿童组 t 值范围:", round(range(ok_results$t_child_TD_minus_DD, na.rm = TRUE), 2), "\n")
cat("  - 成人组 t 值范围:", round(range(ok_results$t_adult_TD_minus_DD, na.rm = TRUE), 2), "\n")
cat("  - 相关系数:", round(cor_result$estimate, 3), "\n")
cat("  - 两组都显著的 ROI:", sum(tmap_plot_data$sig_both, na.rm = TRUE), "\n\n")

cat("下一步: 使用 318×2 矩阵进行 K-means 聚类分析\n")
