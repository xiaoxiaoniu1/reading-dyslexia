#/Volumes/Ting7/2025_7_16_331_subjects/analysis.R
#阅读障碍vs.健康对照（成人、儿童）

# 331个被试的人口统计学信息
################################################################################
library(readxl)
library(dplyr)

df <- read_excel("/Volumes/Ting7/RD/0526all_data.xlsx", sheet = "Sheet1")       #读取人口统计学数据

# 将诊断组别、年龄组和性别转换为标签（可选，方便阅读）
df <- df %>%
  mutate(
    group_d_or_c = ifelse(group_d_or_c == 0, "TD", "DD"),
    group_age = ifelse(group_age == 1, "Adult", "Child"),
    sex = ifelse(sex == 1, "Male", "Female")
  )

# 统计每个组（诊断组、年龄组）在每个站点中的男女数量
group_summary <- df %>%
  group_by(site, group_d_or_c, group_age, sex) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(site, group_d_or_c, group_age, sex)
print(group_summary)# 查看结果

# 统计每个组（诊断组、年龄组）在每个站点中的数量、平均年龄
group_summary <- df %>%
  group_by(group_d_or_c, group_age) %>%
  summarise(
    n = n(),
    mean_age = mean(age_month, na.rm = TRUE),
    sd_age = sd(age_month, na.rm = TRUE),  # 可选
    .groups = "drop"
  ) %>%
  arrange(group_d_or_c, group_age)
print(group_summary)
################################################################################


#Averaged MIND patterns in DD and TD 
################################################################################
library(readxl)
library(tidyverse)

# 1.读取人口学数据
df <- read_excel("/Volumes/Ting7/RD/0526all_data.xlsx", sheet = "Sheet1") %>%
  mutate(
    group_d_or_c = ifelse(group_d_or_c == 0, "TD", "DD"),
    group_age = ifelse(group_age == 1, "Adult", "Child"),
    sex = ifelse(sex == 1, "Male", "Female")
  )

# 2.定义 MIND 文件夹路径
mind_dir <- "/Volumes/Ting7/2025_7_16_331_subjects/MIND"

# 3.写一个函数读取单个被试的 MIND 矩阵（去掉第一行第一列的标签）
read_mind_matrix <- function(id) {
  file <- file.path(mind_dir, paste0(id, ".csv"))
  mat <- as.matrix(read.csv(file, row.names = 1, check.names = FALSE))
  return(mat)
}

# 4.计算每个组的平均 MIND
group_avg_mind <- df %>%
  group_by(group_d_or_c, group_age) %>%
  group_map(~ {
    ids <- .x$id
    mats <- lapply(ids, read_mind_matrix)
    avg_mat <- Reduce("+", mats) / length(mats)
    list(group = paste(.y$group_age, .y$group_d_or_c, sep = "_"), mat = avg_mat)
  })

# 5.保存结果（可选：存成 csv）
for (res in group_avg_mind) {
  write.csv(res$mat, file = paste0("/Volumes/Ting7/2025_7_16_331_subjects/MIND_avg_", res$group, ".csv"))
}

library(pheatmap)

# 读取平均 MIND 矩阵
adult_dd <- as.matrix(read.csv("/Volumes/Ting7/2025_7_16_331_subjects/MIND_avg_Adult_DD.csv", row.names = 1))
adult_td <- as.matrix(read.csv("/Volumes/Ting7/2025_7_16_331_subjects/MIND_avg_Adult_TD.csv", row.names = 1))
child_dd <- as.matrix(read.csv("/Volumes/Ting7/2025_7_16_331_subjects/MIND_avg_Child_DD.csv", row.names = 1))
child_td <- as.matrix(read.csv("/Volumes/Ting7/2025_7_16_331_subjects/MIND_avg_Child_TD.csv", row.names = 1))

# 可视化 Adult DD, Adult TD, Children DD, Children TD 4个组的平均MIND
pheatmap(adult_dd, main = "Adult DD", cluster_rows = FALSE, cluster_cols = FALSE, 
         color = colorRampPalette(c("white","red"))(100),breaks = seq(0, 0.5, length.out = 101), legend_breaks = seq(0, 0.5, 0.1))
pheatmap(adult_td, main = "Adult TD", cluster_rows = FALSE, cluster_cols = FALSE, 
         color = colorRampPalette(c("white","red"))(100),breaks = seq(0, 0.5, length.out = 101), legend_breaks = seq(0, 0.5, 0.1))
pheatmap(child_dd, main = "Child DD", cluster_rows = FALSE, cluster_cols = FALSE, 
         color = colorRampPalette(c("white","red"))(100),breaks = seq(0, 0.5, length.out = 101), legend_breaks = seq(0, 0.5, 0.1))
pheatmap(child_td, main = "Child TD", cluster_rows = FALSE, cluster_cols = FALSE, 
         color = colorRampPalette(c("white","red"))(100),breaks = seq(0, 0.5, length.out = 101), legend_breaks = seq(0, 0.5, 0.1))

#可视化差异图
diff_mat1 <- child_dd - child_td  # children DD vs. TD组间差异矩阵
diff_mat2 <- adult_dd - adult_td  # adult DD vs. TD组间差异矩阵
diff_mat3 <- adult_dd - child_dd  # DD  adult vs. children组间差异矩阵
diff_mat4 <- adult_td - child_td  # TD  adult vs. children组间差异矩阵

pheatmap(
  diff_mat4,
  # main = "Child DD - Child TD",
  # main = "Adult DD - Adult TD",
  # main = "DD Adult - DD Child",
  main = "TD Adult - TD Child",
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  color = colorRampPalette(c("blue", "white", "red"))(100),
  breaks = seq(-0.05, 0.05, length.out = 101),  # 根据差值范围设置
  legend_breaks = seq(-0.05, 0.05, 0.01)
)

#进行T test统计差异
################################################################################
library(readxl)
library(dplyr)
library(pheatmap)
library(purrr)
library(stringr)

# ============ 1) 读取人口统计学信息 ============
demo_file <- "/Volumes/Ting7/RD/0526all_data.xlsx"
df <- read_excel(demo_file, sheet = "Sheet1") %>%
  mutate(
    group_d_or_c = ifelse(group_d_or_c == 0, "TD", "DD"),
    group_age    = ifelse(group_age == 1, "Adult", "Child"),
    sex          = ifelse(sex == 1, "Male", "Female")
  )
  # ) %>%
  # filter(site %in% c(1, 2))   #筛选指定site的被试

# ============ 2) 读取 MIND 矩阵 ============
mind_dir <- "/Volumes/Ting7/2025_7_16_331_subjects/MIND"

# 读单个被试矩阵（第一列是行名，保留脑区标签；自动转为矩阵）
read_mind_matrix <- function(id) {
  fp <- file.path(mind_dir, paste0(id, ".csv"))
  if (!file.exists(fp)) return(NULL)
  mat <- as.matrix(read.csv(fp, row.names = 1, check.names = FALSE))
  storage.mode(mat) <- "numeric"
  return(mat)
}

# 取有文件的被试（避免因缺文件报错）
df$has_file <- vapply(df$id, function(x) file.exists(file.path(mind_dir, paste0(x, ".csv"))), logical(1))
df_avail <- df %>% filter(has_file)

# ============ 3) 分组并读入矩阵列表 ============
get_group_mats <- function(demo, diag_label, age_label) {
  ids <- demo %>%
    filter(group_d_or_c == diag_label, group_age == age_label) %>%
    pull(id)
  mats <- lapply(ids, read_mind_matrix)
  mats <- mats[!vapply(mats, is.null, logical(1))]  # 去掉读取失败
  if (length(mats) == 0) stop(paste("该组没有可用矩阵：", age_label, diag_label))
  # 确保都按同一行列名顺序（以第一个为模板）
  rn <- rownames(mats[[1]]); cn <- colnames(mats[[1]])
  mats <- lapply(mats, function(m) m[rn, cn, drop = FALSE])
  return(mats)
}

child_dd_list <- get_group_mats(df_avail, "DD", "Child")
child_td_list <- get_group_mats(df_avail, "TD", "Child")
adult_dd_list <- get_group_mats(df_avail, "DD", "Adult")
adult_td_list <- get_group_mats(df_avail, "TD", "Adult")

# ============ 4) 逐连接 t 检验（返回 t、p、FDR、显著 t-map） ============
ttest_edgewise <- function(g1_list, g2_list) {
  n_region <- nrow(g1_list[[1]])
  rn <- rownames(g1_list[[1]]); cn <- colnames(g1_list[[1]])
  
  t_mat  <- matrix(0, n_region, n_region, dimnames = list(rn, cn))
  p_mat  <- matrix(1, n_region, n_region, dimnames = list(rn, cn))
  
  # 逐连接 Welch t 检验（独立样本）
  for (i in 1:n_region) {
    for (j in 1:n_region) {
      vals1 <- vapply(g1_list, function(x) x[i, j], numeric(1))
      vals2 <- vapply(g2_list, function(x) x[i, j], numeric(1))
      tt <- t.test(vals1, vals2)   # Welch (var.equal = FALSE)
      t_mat[i, j] <- unname(tt$statistic)                # t 值（带正负号）
      p_mat[i, j] <- tt$p.value                          # p 值（用于校正）
    }
  }
  
  # FDR 多重比较校正（对所有边整体校正）
  p_adj <- matrix(p.adjust(as.vector(p_mat), method = "fdr"),
                  nrow = n_region, ncol = n_region, dimnames = list(rn, cn))
  
  # 显著性掩膜：非显著处设为 NA，但保留 t 值作为显示量
  t_sig <- t_mat
  t_sig[p_adj > 0.05] <- NA
  
  # 可选：对角线不显示（一般对角线=1或无意义）
  diag(t_sig) <- NA
  diag(t_mat) <- NA
  
  list(t_mat = t_mat, p_mat = p_mat, p_adj = p_adj, t_sig = t_sig)
}

# 运行（DD 相对 TD；正值 = DD > TD）
res_child <- ttest_edgewise(child_dd_list, child_td_list)
res_adult <- ttest_edgewise(adult_dd_list, adult_td_list)

# ============ 5) 可视化：显示 t 值（显著处保留） ============
plot_tsig <- function(t_sig, title = "") {
  # 颜色对称：以 |t| 的最大值定尺，红=DD>TD，蓝=TD>DD
  tmax <- max(abs(t_sig), na.rm = TRUE)
  tmax <- ifelse(is.finite(tmax), tmax, 1)
  breaks <- seq(-tmax, tmax, length.out = 101)
  
  pheatmap(
    t_sig,
    main = title,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = colorRampPalette(c("blue", "white", "red"))(100),
    breaks = breaks,
    na_col = "white",
    border_color = "grey",
    legend_breaks = pretty(c(-tmax, tmax), n = 5)
  )
}

plot_tsig(res_child$t_sig, "Child: DD - TD (t-map, FDR<0.05)")
plot_tsig(res_adult$t_sig, "Adult: DD - TD (t-map, FDR<0.05)")








#MIND degree
################################################################################
# 将每位被试的 MSN 特征提取为向量
library(readxl)
library(tidyverse)

# 读取元数据
df <- read_excel("/Volumes/Ting7/RD/0526all_data.xlsx", sheet = "Sheet1")

# 准备一个空矩阵存储所有被试的 MSN 向量化结果
subject_ids <- df$id
n_subjects <- length(subject_ids)
n_features <- 68 * (68 - 1) / 2
msn_data <- matrix(NA, nrow = n_subjects, ncol = n_features)

# 提取上三角（不含对角线）索引
get_upper_triangle <- function(mat) {
  mat[upper.tri(mat)]
}

# 读取每位被试的 MSN
for (i in seq_along(subject_ids)) {
  id <- subject_ids[i]
  file_path <- file.path("/Volumes/Ting7/2025_7_16_331_subjects/MIND", paste0(id, ".csv"))
  mat <- as.matrix(read.csv(file_path, row.names = 1))
  vec <- get_upper_triangle(mat)
  msn_data[i, ] <- vec
}

#####################
n_roi <- 68
n_subj <- nrow(msn_data)

# 辅助函数：把向量还原为 68x68 对称矩阵
vec_to_matrix <- function(vec, n = n_roi) {
  mat <- matrix(0, n, n)
  mat[upper.tri(mat)] <- vec
  mat <- mat + t(mat)
  diag(mat) <- 0  # 对角线设 NA
  return(mat)
}

# 存储度值矩阵（331 × 68）
msn_degree <- matrix(NA, nrow = n_subj, ncol = n_roi)

for (i in 1:n_subj) {
  mat <- vec_to_matrix(msn_data[i, ])
  msn_degree[i, ] <- rowMeans(mat, na.rm = TRUE)  # 按行均值（加权度）
}

# msn_degree 就是你要的 331 × 68
dim(msn_degree)

#构建t test模型
run_ttest_degree <- function(data, df, group_age_value, group_label, p_thresh = 0.05, use_fdr = FALSE) {
  idx <- df$group_age == group_age_value
  data_sub <- data[idx, ]
  group_sub <- df$group_d_or_c[idx]
  
  results <- apply(data_sub, 2, function(feature) {
    test <- t.test(feature ~ group_sub)
    c(t = test$statistic,
      p = test$p.value,
      mean_DD = mean(feature[group_sub == 1]),
      mean_TD = mean(feature[group_sub == 0]))
  })
  
  res_df <- as.data.frame(t(results))
  res_df$fdr <- p.adjust(res_df$p, method = "fdr")
  
  sig_df <- res_df[if (use_fdr) res_df$fdr < p_thresh else res_df$p < p_thresh, ]
  
  message("Found ", nrow(sig_df), " significant ROIs in ", group_label)
  
  return(list(
    all_results = res_df,
    sig_results = sig_df
  ))
}
#构建线性回归模型
run_lm_degree <- function(data, df, group_age_value, group_label, p_thresh = 0.05, use_fdr = FALSE) {
  # 选择对应的年龄组
  idx <- df$group_age == group_age_value
  data_sub <- data[idx, ]
  df_sub <- df[idx, ]
  
  results <- apply(data_sub, 2, function(feature) {
    # 构建线性模型: group + age + sex
    model <- lm(feature ~ group_d_or_c + age_month + sex, data = df_sub)
    coef_summary <- summary(model)$coefficients
    
    # 提取 group 的统计量
    c(
      t = coef_summary["group_d_or_c", "t value"],
      p = coef_summary["group_d_or_c", "Pr(>|t|)"],
      mean_DD = mean(feature[df_sub$group_d_or_c == 1]),
      mean_TD = mean(feature[df_sub$group_d_or_c == 0])
    )
  })
  
  res_df <- as.data.frame(t(results))
  res_df$fdr <- p.adjust(res_df$p, method = "fdr")
  
  sig_df <- res_df[if (use_fdr) res_df$fdr < p_thresh else res_df$p < p_thresh, ]
  
  message("Found ", nrow(sig_df), " significant ROIs in ", group_label, " (controlling age & sex)")
  
  return(list(
    all_results = res_df,
    sig_results = sig_df
  ))
}
child_combat_result <- run_lm_degree(harmonized_data, df, group_age_value = 2, group_label = "Children (ComBat)", p_thresh = 0.05, use_fdr = TRUE)
adult_combat_result <- run_lm_degree(harmonized_data, df, group_age_value = 1, group_label = "Adults (ComBat)", p_thresh = 0.05, use_fdr = TRUE)

#用残差
run_ttest_resid <- function(data, df, group_age_value, group_label,
                            p_thresh = 0.05, use_fdr = FALSE) {
  # 选择对应年龄组
  idx <- df$group_age == group_age_value
  data_sub <- data[idx, ]
  df_sub <- df[idx, ]
  group_sub <- df_sub$group_d_or_c
  
  results <- apply(data_sub, 2, function(feature) {
    # 先控制 age, sex
    model <- lm(feature ~ age_month + sex, data = df_sub)
    resid_feature <- residuals(model)
    
    # 在残差上做 t test
    test <- t.test(resid_feature ~ group_sub)
    c(
      t = test$statistic,
      p = test$p.value,
      mean_DD = mean(resid_feature[group_sub == 1]),
      mean_TD = mean(resid_feature[group_sub == 0])
    )
  })
  
  res_df <- as.data.frame(t(results))
  res_df$fdr <- p.adjust(res_df$p, method = "fdr")
  
  sig_df <- res_df[if (use_fdr) res_df$fdr < p_thresh else res_df$p < p_thresh, ]
  
  message("Found ", nrow(sig_df), " significant ROIs in ",
          group_label, " (residualized for age & sex)")
  
  return(list(
    all_results = res_df,
    sig_results = sig_df
  ))
}
child_combat_result <- run_ttest_resid(harmonized_data, df, group_age_value = 2, group_label = "Children (ComBat)", p_thresh = 0.05, use_fdr = TRUE)
adult_combat_result <- run_ttest_resid(harmonized_data, df, group_age_value = 1, group_label = "Adults (ComBat)", p_thresh = 0.05, use_fdr = TRUE)



#######

child_deg_result <- run_ttest_degree(msn_degree, df, 2, "Children", p_thresh = 0.05, use_fdr = TRUE)
adult_deg_result <- run_ttest_degree(msn_degree, df, 1, "Adults", p_thresh = 0.05, use_fdr = TRUE)


path <- "/Users/qiantingcheng/testproject/aparc__t_value/age/"
t_data_area <- read_freesurfer_table(paste0(path, "aparc_area_t_value.table"), measure = "area")
roi_labels <- t_data_area$label 
# 创建绘图数据
# 选children or adult
# age_group <- child_deg_result
# age_group <- adult_deg_result
# 
age_group <- child_combat_result
age_group <- adult_combat_result
# 从 t 检验结果中取显著 ROI
sig_mask <- age_group$all_results$p < 0.05
t_vals_sig <- age_group$all_results$t
t_vals_sig[!sig_mask] <- NA  # 非显著置为 NA

plot_data <- data.frame(
  label = roi_labels,  # 这里必须和 atlas 中的 label 名称匹配
  value = t_vals_sig  # 对应每个 ROI 的 t 值
)

# 画图
p <- ggseg(plot_data,                   
           mapping = aes(fill = value),color='black',size=0.1,
           show.legend = TRUE) +
  scale_fill_viridis_c(
    option = "viridis",
    limits = c(-1, 5),
    oob = scales::squish
  ) +
  labs(title = " ", fill = "t-value")


print(p)










#####################


#仅展示组之间的情况
################################################################################
library(dplyr)
library(ggseg)

# ---- 1. 定义函数：向量 -> 矩阵
vec_to_matrix <- function(vec, n_nodes = 68) {
  mat <- matrix(0, n_nodes, n_nodes)
  mat[upper.tri(mat)] <- vec
  mat[lower.tri(mat)] <- t(mat)[lower.tri(mat)]
  return(mat)
}

# ---- 2. 分组计算平均 MSN
group_labels <- paste(df$group_age, df$group_d_or_c, sep = "_")  # e.g., "Child_DD"
unique_groups <- unique(group_labels)

group_avg_degree <- list()

for (grp in unique_groups) {
  idx <- which(group_labels == grp)
  avg_vec <- colMeans(msn_data[idx, ], na.rm = TRUE)   # 平均向量
  avg_mat <- vec_to_matrix(avg_vec)                    # 转回 68×68 矩阵
  degree_vals <- rowMeans(avg_mat)                     # 每个脑区的平均连接强度（加权度）
  group_avg_degree[[grp]] <- degree_vals
}

# ---- 3. 转换为 ggseg 可用数据框
# 你需要有一个包含脑区名称的向量（和68个ROI顺序对应）

path <- "/Users/qiantingcheng/testproject/aparc__t_value/age/"
t_data_area <- read_freesurfer_table(paste0(path, "aparc_area_t_value.table"), measure = "area")
roi_labels <- t_data_area$label  # 例如 Desikan-Killiany 的 label

plot_data <- do.call(rbind, lapply(names(group_avg_degree), function(grp) {
  data.frame(
    label = roi_labels,
    value = group_avg_degree[[grp]],
    group = grp
  )
}))

# ---- 4. ggseg 绘图
ggseg(data = plot_data,
      atlas = dk,                 # 如果是 Desikan-Killiany
      mapping = aes(fill = value),
      show.legend = TRUE) +
  facet_wrap(~group) +
  scale_fill_viridis_c(option = "plasma") +
  theme_void()



library(dplyr)
library(ggseg)
library(viridis)
library(scales)


# atlas 改列名成 label 以匹配 ROI 名称
atlas_data <- dk %>% rename(label = region)

plots <- list()  # 用来存储所有组的图
# 循环每个组，分别处理 & 绘图
for (grp in names(group_avg_degree)) {
  
  # --- 1. 构建该组的 68×3 数据 ---
  plot_data <- data.frame(
    label = roi_labels,                 # 68 个 ROI 名称
    value = group_avg_degree[[grp]],    # 对应的度值
    group = grp                         # 当前组名
  )
  

  
  # --- 3. 绘制单独一张 ggseg 图 ---
  # p <- ggseg(plot_data,
  #            mapping = aes(fill = value),
  #            show.legend = TRUE) +
  #   # scale_fill_viridis_c(option = "turbo",
  #   #                      limits = c(0.05, 0.25), oob = scales::squish, n =5) +
  #   scale_fill_manual(values = viridis(5, option = "turbo")) +
  #   labs(title = grp) +
  #   theme_void() +
  #   theme(
  #     plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
  #   )
  library(dplyr)
  library(viridis)
  library(scales)
  
  # 将连续值分成 5 个区间
  plot_data <- plot_data %>%
    mutate(value_cat = cut(value,
                           breaks = seq(0.0, 0.25, length.out = 11),
                           include.lowest = TRUE))
  
  # 绘图
  p <- ggseg(plot_data,
             mapping = aes(fill = value_cat),
             show.legend = TRUE) +
    scale_fill_manual(values = viridis(10, option = "turbo")) +
    labs(title = grp) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold")
    )

  
  
  # 存入结构体（list）
  plots[[grp]] <- p
  # # 显示
  # print(p)
  # 
  # # 如果要保存
  # ggsave(paste0("ggseg_", grp, ".png"), p, width = 6, height = 4, dpi = 300)
}
# 之后可以单独访问：
# plots[["Child_DD"]]
# 或者一次性显示：
library(patchwork)
wrap_plots(plots, ncol = 1)





################################################################################




#不进行combat 直接组间比较
######################################
library(pheatmap)

run_ttest_raw <- function(data, df, group_age_value, group_label, p_thresh = 0.05, use_fdr = TRUE) {
  idx <- df$group_age == group_age_value
  data_sub <- data[idx, ]
  group_sub <- df$group_d_or_c[idx]
  
  results <- apply(data_sub, 2, function(feature) {
    test <- t.test(feature ~ group_sub)
    c(t = test$statistic,
      p = test$p.value,
      mean_DD = mean(feature[group_sub == 1]),
      mean_TD = mean(feature[group_sub == 0]))
  })
  
  res_df <- as.data.frame(t(results))
  res_df$fdr <- p.adjust(res_df$p, method = "fdr")
  
  n_roi <- 68
  t_mat <- matrix(0, n_roi, n_roi)
  p_mat <- matrix(1, n_roi, n_roi)
  tri_idx <- which(upper.tri(t_mat))
  
  t_mat[tri_idx] <- res_df$t
  p_mat[tri_idx] <- res_df$p
  t_mat <- t_mat + t(t_mat)
  p_mat <- p_mat + t(p_mat)
  
  ### ✅ 绘制所有连接的热图（不做筛选）：
  # range_vals <- range(t_mat, na.rm = TRUE)
  # if (range_vals[1] == range_vals[2]) range_vals <- c(-1, 1)
  t_max <- max(abs(t_mat), na.rm = TRUE)
  colormap <- colorRampPalette(c("blue", "white", "red"))(100)
  breaks_seq <- seq(-t_max, t_max, length.out = 101)
  # pheatmap(t_mat,
  #          cluster_rows = FALSE,
  #          cluster_cols = FALSE,
  #          color = colormap,
  #          breaks = breaks_seq,
  #          main = paste0("All Connections T-values (", group_label, ")"))
  
  ### ✅ 仅显著连接热图（FDR 或 p）
  sig_stat <- if (use_fdr) res_df$fdr else res_df$p
  sig_mask <- sig_stat < p_thresh
  sig_df <- res_df[sig_mask, , drop = FALSE]
  message("Found ", nrow(sig_df), " significant connections in ", group_label)
  
  if (nrow(sig_df) > 0) {
    t_sig <- matrix(0, n_roi, n_roi)
    t_sig[tri_idx[sig_mask]] <- res_df$t[sig_mask]
    t_sig <- t_sig + t(t_sig)
    
    range_vals2 <- range(t_sig[t_sig != 0], na.rm = TRUE)
    if (length(range_vals2) == 0 || range_vals2[1] == range_vals2[2]) {
      range_vals2 <- c(-1, 1)
    }
    
    pheatmap(t_sig,
             cluster_rows = FALSE,
             cluster_cols = FALSE,
             color = colormap,
             breaks = breaks_seq,
             main = paste0("Significant Connections T-values (", group_label, ", ",
                           ifelse(use_fdr, "FDR", "p"), " <", p_thresh, ")"))
  } else {
    message("No significant connections to plot (", group_label, ")")
  }
  
  return(list(
    all_results = res_df,
    sig_results = sig_df,
    t_matrix = t_mat,
    p_matrix = p_mat
  ))
}
# 儿童组（group_age == 2）
child_raw_result <- run_ttest_raw(msn_data, df, group_age_value = 2, group_label = "Children", p_thresh = 0.05, use_fdr = TRUE)
adult_raw_result <- run_ttest_raw(msn_data, df, group_age_value = 1, group_label = "Adults", p_thresh = 0.05, use_fdr = TRUE)

#msn_degree
child_raw_result <- run_ttest_raw(msn_degree, df, group_age_value = 2, group_label = "Children", p_thresh = 0.05, use_fdr = TRUE)
adult_raw_result <- run_ttest_raw(msn_degree, df, group_age_value = 1, group_label = "Adults", p_thresh = 0.05, use_fdr = TRUE)






















################################################################################
#MIND degree
################################################################################
# 将每位被试的 MSN 特征提取为向量
library(readxl)
library(tidyverse)

# 读取元数据
df <- read_excel("/Volumes/Ting7/RD/0526all_data.xlsx", sheet = "Sheet1")

# 准备一个空矩阵存储所有被试的 MSN 向量化结果
subject_ids <- df$id
n_subjects <- length(subject_ids)
n_features <- 68 * (68 - 1) / 2
msn_data <- matrix(NA, nrow = n_subjects, ncol = n_features)

# 提取上三角（不含对角线）索引
get_upper_triangle <- function(mat) {
  mat[upper.tri(mat)]
}

# 读取每位被试的 MSN
for (i in seq_along(subject_ids)) {
  id <- subject_ids[i]
  file_path <- file.path("/Volumes/Ting7/2025_7_16_331_subjects/MIND", paste0(id, ".csv"))
  mat <- as.matrix(read.csv(file_path, row.names = 1))
  vec <- get_upper_triangle(mat)
  msn_data[i, ] <- vec
}

#####################
n_roi <- 68
n_subj <- nrow(msn_data)

# 辅助函数：把向量还原为 68x68 对称矩阵
vec_to_matrix <- function(vec, n = n_roi) {
  mat <- matrix(0, n, n)
  mat[upper.tri(mat)] <- vec
  mat <- mat + t(mat)
  diag(mat) <- 0  # 对角线设 NA
  return(mat)
}

# 存储度值矩阵（331 × 68）
msn_degree <- matrix(NA, nrow = n_subj, ncol = n_roi)

for (i in 1:n_subj) {
  mat <- vec_to_matrix(msn_data[i, ])
  msn_degree[i, ] <- rowMeans(mat, na.rm = TRUE)  # 按行均值（加权度）
}

# msn_degree 就是你要的 331 × 68
dim(msn_degree)


library(tidyverse)
library(mgcv)
library(gamlss)

# 把 msn_degree 转成数据框，并加上被试信息
colnames(msn_degree) <- paste0("ROI", 1:ncol(msn_degree))
df_degree <- bind_cols(df, as.data.frame(msn_degree))

# 转成长格式，便于批量建模
df_long <- df_degree %>%
  pivot_longer(cols = starts_with("ROI"),
               names_to = "ROI",
               values_to = "degree")

# 示例：对某个 ROI 拟合 GAM 曲线
m1 <- gam(degree ~ s(age_month, k=5) + sex, data = filter(df_long, ROI == "ROI1"), method = "REML")
summary(m1)
plot(m1, shade = TRUE, main = "ROI1 growth curve")

#批量拟合所有 ROI，并做多重比较校正
rois <- unique(df_long$ROI)
results <- tibble(ROI = rois, edf = NA_real_, p_val = NA_real_)

for (i in seq_along(rois)) {
  sub <- filter(df_long, ROI == rois[i])
  m <- gam(degree ~ s(age_month, k=5) + sex, data = sub, method = "REML")
  s <- summary(m)
  results$edf[i] <- s$edf[1]          # 平滑项的自由度
  results$p_val[i] <- s$s.pv[1]       # 平滑项的 p 值
}
results <- results %>% mutate(p_fdr = p.adjust(p_val, "fdr"))

#如果想要 “百分位生长曲线”（类似身高体重曲线）
# 以 ROI1 为例，gamlss 可以建均值和方差随年龄变化
roi1 <- filter(df_long, ROI == "ROI1")
m2 <- gamlss(degree ~ pb(age_month), sigma.formula = ~ pb(age_month),
             data = roi1, family = NO)
roi1 <- df_long %>%
  filter(ROI == "ROI1") %>%
  select(age_month, degree, sex, group_d_or_c)   # 只保留需要的列
roi1 <- drop_na(roi1)        # 再确保没有 NA

m2 <- gamlss(
  degree ~ pb(age_month),
  sigma.formula = ~ pb(age_month),
  data = roi1,
  family = NO
)

# 预测百分位曲线
ages <- seq(min(roi1$age_month), max(roi1$age_month), by = 0.5)
new_data <- data.frame(age_month = ages)

pred_mu <- predict(m2, newdata = new_data, what = "mu", type = "response")
pred_sigma <- predict(m2, newdata = new_data, what = "sigma", type = "response")

# 百分位曲线
p5  <- pred_mu + qnorm(0.05) * pred_sigma
p50 <- pred_mu
p95 <- pred_mu + qnorm(0.95) * pred_sigma
plot(roi1$age_month, roi1$degree, main = "ROI1 Percentile growth curve",
     xlab = "Age (months)", ylab = "Degree")
lines(ages, p50, col="blue", lwd=2)
lines(ages, p5,  col="red", lty=2)
lines(ages, p95, col="red", lty=2)

plot_df <- data.frame(
  age_month = ages,
  p5 = p5,
  p50 = p50,
  p95 = p95
)

# 绘图
ggplot(roi1, aes(x = age_month, y = degree, color = factor(group_d_or_c))) +
  geom_point(alpha = 0.6) +
  geom_line(data = plot_df, aes(x = age_month, y = p50), color = "black", size = 1) +
  geom_line(data = plot_df, aes(x = age_month, y = p5), color = "black", linetype = "dashed") +
  geom_line(data = plot_df, aes(x = age_month, y = p95), color = "black", linetype = "dashed") +
  scale_color_manual(values = c("0" = "steelblue", "1" = "firebrick"),
                     labels = c("TD", "DD"),
                     name = "Group") +
  labs(x = "Age (months)", y = "Degree", 
       title = "Normative Growth Model of MSN Degree") +
  theme_minimal(base_size = 14)




#可视化
library(ggplot2)
ggplot(filter(df_long, ROI=="ROI1"), aes(x=age_month, y=degree)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="gam", formula = y ~ s(x, k=5), color="blue") +
  theme_minimal() +
  labs(title="ROI1 growth curve", y="Degree")

ggplot(df_long %>% filter(ROI %in% paste0("ROI", 1:12)),  # 先画前12个ROI
       aes(x=age_month, y=degree)) +
  geom_point(alpha=0.4) +
  geom_smooth(method="gam", formula = y ~ s(x, k=5), color="blue") +
  facet_wrap(~ROI, scales="free_y") +
  theme_bw()









#DD TD分别拟合
############################
library(gamlss)
library(dplyr)
library(ggplot2)

# 假设 df_long 有 age_month, degree, group_d_or_c
# group_d_or_c: 0 = TD, 1 = DD

# =============================
# 1) 拆分数据
df_td <- df_long %>% filter(group_d_or_c == 0)
df_dd <- df_long %>% filter(group_d_or_c == 1)

# ROI1  ROI3 ROI8  ROI10      ROI26  ROI28  ROI29   ROI35
# 针对某一个 ROI (比如 ROI1)
roi1 <- df_long %>%
  filter(ROI == "ROI1") %>%                  # 只保留 ROI1
  select(age_month, degree, sex, group_d_or_c) %>%  # 只保留必要列
  drop_na()                                  # 删除 NA

# 检查是否还存在 NA
colSums(is.na(roi1))

# 分组数据
roi1_td <- roi1 %>% filter(group_d_or_c == 0)
roi1_dd <- roi1 %>% filter(group_d_or_c == 1)

# 2) 分别拟合 GAMLSS 模型
m_td <- gamlss(degree ~ pb(age_month),
               sigma.formula = ~pb(age_month),
               data = roi1_td, family = NO)

m_dd <- gamlss(degree ~ pb(age_month),
               sigma.formula = ~pb(age_month),
               data = roi1_dd, family = NO)

# =============================
# 3) 预测轨迹
ages <- seq(min(df_long$age_month), max(df_long$age_month), by = 1)

# TD
pred_mu_td <- predict(m_td, newdata = data.frame(age_month = ages),
                      what = "mu", type = "response")
pred_sigma_td <- predict(m_td, newdata = data.frame(age_month = ages),
                         what = "sigma", type = "response")
plot_td <- data.frame(
  age_month = ages,
  mu = pred_mu_td,
  lower = pred_mu_td + qnorm(0.05) * pred_sigma_td,
  upper = pred_mu_td + qnorm(0.95) * pred_sigma_td,
  group = "TD"
)

# DD
pred_mu_dd <- predict(m_dd, newdata = data.frame(age_month = ages),
                      what = "mu", type = "response")
pred_sigma_dd <- predict(m_dd, newdata = data.frame(age_month = ages),
                         what = "sigma", type = "response")
plot_dd <- data.frame(
  age_month = ages,
  mu = pred_mu_dd,
  lower = pred_mu_dd + qnorm(0.05) * pred_sigma_dd,
  upper = pred_mu_dd + qnorm(0.95) * pred_sigma_dd,
  group = "DD"
)

# 合并
plot_group <- bind_rows(plot_td, plot_dd)

# =============================
# 4) 可视化
# ggplot() +
#   # 散点：原始每个被试的数据
#   geom_point(data = roi1,
#              aes(x = age_month, y = degree, color = factor(group_d_or_c)),
#              alpha = 0.6, size = 2) +
#   
#   # TD / DD 均值轨迹
#   geom_line(data = plot_group,
#             aes(x = age_month, y = mu, color = group),
#             size = 1.2) +
#   
#   # TD / DD 区间带
#   geom_ribbon(data = plot_group,
#               aes(x = age_month, ymin = lower, ymax = upper, fill = group),
#               alpha = 0.2) +
#   
#   scale_color_manual(values = c("0" = "steelblue", "1" = "firebrick",
#                                 "TD" = "steelblue", "DD" = "firebrick"),
#                      labels = c("TD", "DD"),
#                      name = "Group") +
# 
#   scale_fill_manual(values = c("TD" = "steelblue", "DD" = "firebrick"),
#                     guide = "none") +   # fill 图例隐藏
#   
#   labs(x = "Age", y = "MSN Degree",
#        title = "Group-specific Growth Trajectories (TD vs DD)") +
#   theme_minimal(base_size = 14) +
#   theme(legend.position = "top")

# 数据预处理，统一分组变量
roi1 <- roi1 %>%
  mutate(Group = ifelse(group_d_or_c == 0, "TD", "DD"))

plot_group <- plot_group %>%
  mutate(Group = ifelse(group %in% c("0","TD"), "TD", "DD"))

# 绘图
p <- ggplot() +
  # 散点：原始每个被试
  geom_point(data = roi1,
             aes(x = age_month, y = degree, color = Group),
             alpha = 0.6, size = 2) +
  
  # TD / DD 均值轨迹
  geom_line(data = plot_group,
            aes(x = age_month, y = mu, color = Group),
            size = 1.2) +
  
  # TD / DD 置信区间
  geom_ribbon(data = plot_group,
              aes(x = age_month, ymin = lower, ymax = upper, fill = Group),
              alpha = 0.2) +
  
  # 只保留一份图例
  scale_color_manual(values = c("TD" = "steelblue", "DD" = "firebrick"),
                     name = "Group") +
  scale_fill_manual(values = c("TD" = "steelblue", "DD" = "firebrick"),
                    guide = "none") +
  
  labs(x = "Age", y = "MSN Degree",
       title = " ") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.8), # 画出坐标轴边框
    axis.line = element_line(color = "black", linewidth = 0.6),                 # 加上坐标轴线
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(color = "black"),
    panel.grid.major = element_blank(),  # 去掉主网格线
    panel.grid.minor = element_blank(),  # 去掉次网格线
    plot.title = element_text(hjust = 0.5)
  )
  # theme(legend.position = "right")
print(p)
# 保存画布大小（单位：英寸）
ggsave("/Volumes/Ting7/2025_7_16_331_subjects/2025.9.9_growth_trajectory/zhongqi/ROI1.pdf", plot = p, width = 4, height = 4, dpi = 300, bg = "white")












































