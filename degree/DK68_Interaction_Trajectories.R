# ==============================================================================
# GAMLSS 组间交互作用显著脑区的发育轨迹分析 (DK68 版)
# 功能：针对 ANOVA 中 Interaction 显著的脑区，分别拟合 TD 和 DD 的发育曲线
# ==============================================================================

library(readxl)
library(dplyr)
library(gamlss) # 使用 gamlss
library(ggplot2)
library(tidyr)

# 1. 环境配置与路径定义
base_dir     <- "/data/home/tqi/data1/share/after_freesurfer"
demo_file    <- file.path(base_dir, "FILE/all_data_cqt.xlsx")
mind_csv_dir <- file.path(base_dir, "fs_subjects_all/MIND_DK68_combat")
anova_file   <- file.path(base_dir, "FILE/MIND_DK68_ANOVA/ANOVA_DK68_degree_results.csv")
plot_out_dir <- file.path(base_dir, "FILE/MIND_DK68_ANOVA/Growth_Trajectories_Interaction")

dir.create(plot_out_dir, showWarnings = FALSE, recursive = TRUE)

# 2. 识别 Interaction 显著的脑区 (p < 0.05)
anova_res <- read.csv(anova_file)
sig_rois <- anova_res$feature[which(anova_res$p_Interaction < 0.05)]
cat("发现 Interaction 显著的脑区数量:", length(sig_rois), "\n")

if (length(sig_rois) == 0) {
  stop("未发现显著的 Interaction 效应脑区。")
}

# 3. 加载被试信息 (包含年龄、性别、分组等)
# 自动匹配 MIND_DK68_combat 目录下的文件 (格式: ID_combat.csv)
matrix_files <- list.files(mind_csv_dir, pattern = "_combat\\.csv$")
file_ids <- sub("_combat\\.csv$", "", matrix_files)
file_map <- data.frame(id = file_ids, file_name = matrix_files, stringsAsFactors = FALSE)

df_demo <- read_excel(demo_file, sheet = "Sheet1") %>%
  mutate(id = as.character(id)) %>%
  inner_join(file_map, by = "id") %>%
  mutate(
    group_d_or_c = as.numeric(group_d_or_c), # 0=TD, 1=DD
    sex = as.numeric(sex),
    age_month = as.numeric(age_month),
    file_path = file.path(mind_csv_dir, file_name)
  )

cat("匹配到脑网络数据的被试数量:", nrow(df_demo), "\n")

# 定义计算 Degree 的辅助函数
calc_degree <- function(mat) {
  diag(mat) <- NA
  rowMeans(mat, na.rm = TRUE)
}

# 4. 加载所有被试显著脑区的 Degree 数据
cat("正在读取脑网络数据并提取显著脑区指标...\n")
all_data_list <- list()

for (i in 1:nrow(df_demo)) {
  subj_file <- df_demo$file_path[i]
  tmp_deg <- tryCatch({
    mat <- as.matrix(read.csv(subj_file, row.names = 1, check.names = FALSE))
    deg <- calc_degree(mat)
    data.frame(region = names(deg), degree = as.numeric(deg), stringsAsFactors = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(tmp_deg)) next
  
  # 提取显著脑区
  tmp_sig <- tmp_deg %>% filter(region %in% sig_rois)
  if (nrow(tmp_sig) > 0) {
    # 合并人口学信息
    tmp_sig$id <- df_demo$id[i]
    tmp_sig$group_d_or_c <- df_demo$group_d_or_c[i]
    tmp_sig$age_month <- df_demo$age_month[i]
    tmp_sig$sex <- df_demo$sex[i]
    all_data_list[[i]] <- tmp_sig
  }
}

df_long <- bind_rows(all_data_list)
colnames(df_long)[colnames(df_long) == "region"] <- "ROI"

# 5. 循环处理每个显著脑区并绘图
for (roi_name in sig_rois) {
  cat("正在处理脑区:", roi_name, "\n")
  
  # 筛选当前脑区数据
  roi_data <- df_long %>%
    filter(ROI == roi_name) %>%
    drop_na(degree, age_month, group_d_or_c)
  
  if (nrow(roi_data) < 30) {
    cat("  数据量过少，跳过:", roi_name, "\n")
    next
  }
  
  # 分组
  roi_td <- roi_data %>% filter(group_d_or_c == 0)
  roi_dd <- roi_data %>% filter(group_d_or_c == 1)
  
  # 分别拟合线性模型
  m_td <- tryCatch({
    lm(degree ~ age_month, data = roi_td)
  }, error = function(e) NULL)
  
  m_dd <- tryCatch({
    lm(degree ~ age_month, data = roi_dd)
  }, error = function(e) NULL)
  
  if (is.null(m_td) || is.null(m_dd)) {
    cat("  模型拟合失败，跳过:", roi_name, "\n")
    next
  }
  
  # 生成预测区间
  age_range <- seq(min(roi_data$age_month), max(roi_data$age_month), length.out = 100)
  
  # TD 预测
  pred_td <- predict(m_td, newdata = data.frame(age_month = age_range), interval = "prediction", level = 0.90)
  plot_td <- data.frame(
    age_month = age_range,
    mu = pred_td[, "fit"],
    lower = pred_td[, "lwr"],
    upper = pred_td[, "upr"],
    Group = "TD"
  )
  
  # DD 预测
  pred_dd <- predict(m_dd, newdata = data.frame(age_month = age_range), interval = "prediction", level = 0.90)
  plot_dd <- data.frame(
    age_month = age_range,
    mu = pred_dd[, "fit"],
    lower = pred_dd[, "lwr"],
    upper = pred_dd[, "upr"],
    Group = "DD"
  )
  
  plot_combined <- bind_rows(plot_td, plot_dd)
  
  # 原始散点分组标签
  roi_data_plot <- roi_data %>%
    mutate(Group = ifelse(group_d_or_c == 0, "TD", "DD"))
  
  # 绘图
  p <- ggplot() +
    # 原始散点
    geom_point(data = roi_data_plot, aes(x = age_month, y = degree, color = Group), alpha = 0.3, size = 1.5) +
    # 均值曲线
    geom_line(data = plot_combined, aes(x = age_month, y = mu, color = Group), size = 1.2) +
    # 5%-95% 发育参考区间
    geom_ribbon(data = plot_combined, aes(x = age_month, ymin = lower, ymax = upper, fill = Group), alpha = 0.15) +
    # 样式美化
    scale_color_manual(values = c("TD" = "firebrick", "DD" = "steelblue")) +
    scale_fill_manual(values = c("TD" = "firebrick", "DD" = "steelblue")) +
    labs(x = "Age (Years)", y = "MSN Degree", title = roi_name) +
    theme_minimal(base_size = 14) +
    theme(
      legend.position = "right",
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      axis.text = element_text(color = "black"),
      panel.grid.major = element_line(color = "grey90"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  # 保存图片
  safe_roi_name <- gsub("/", "_", roi_name)
  ggsave(file.path(plot_out_dir, paste0(safe_roi_name, "_trajectory.png")), plot = p, width = 6, height = 5, dpi = 300)
}

cat("分析完成！所有显著脑区的轨迹图已保存至:", plot_out_dir, "\n")
