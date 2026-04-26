rm(list=ls())

result_dir <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Results/Disease_area_volume'
output_file <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Results/Disease_area_volume/All_Zscore_summary.csv'
output_quant <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Results/Disease_area_volume/All_Quant_summary.csv'

# 找到所有 rds 文件
rds_files <- list.files(result_dir, pattern = '_model_new[.]rds$', recursive = TRUE, full.names = TRUE)
cat('Found', length(rds_files), 'rds files\n')

Z_all <- NULL
Q_all <- NULL

for (f in rds_files) {
  tryCatch({
    res <- readRDS(f)
    feature_name <- res$i  # 特征名
    
    # 取最后一个特征的 Z-score（每个rds只包含当前特征）
    if (!is.null(res$Zscore) && length(res$Zscore) > 0) {
      z_df <- res$Zscore[[feature_name]]
      if (!is.null(z_df) && nrow(z_df) > 0) {
        colnames(z_df) <- feature_name
        if (is.null(Z_all)) {
          Z_all <- z_df
        } else {
          # 按行名合并
          Z_all <- merge(Z_all, z_df, by = 'row.names', all = TRUE)
          rownames(Z_all) <- Z_all$Row.names
          Z_all$Row.names <- NULL
        }
      }
    }
    
    # 取 Quant_score
    if (!is.null(res$Quant_data) && length(res$Quant_data) > 0) {
      q_df <- res$Quant_data[[feature_name]]
      if (!is.null(q_df) && nrow(q_df) > 0) {
        colnames(q_df) <- feature_name
        if (is.null(Q_all)) {
          Q_all <- q_df
        } else {
          Q_all <- merge(Q_all, q_df, by = 'row.names', all = TRUE)
          rownames(Q_all) <- Q_all$Row.names
          Q_all$Row.names <- NULL
        }
      }
    }
    
    cat('Done:', feature_name, '\n')
  }, error = function(e) {
    cat('Error in', f, ':', e$message, '\n')
  })
}

cat('\nZ-score matrix:', nrow(Z_all), 'subjects x', ncol(Z_all), 'features\n')

# 修复行名：从路径提取 subject ID，统一去掉括号、空格、下划线
extract_id <- function(x) {
  m <- regmatches(x, regexpr('(?<=fs_subjects_all/)[^/]+', x, perl=TRUE))
  if(length(m)==0 || m=='') {
    m <- regmatches(x, regexpr('[^/]+(?=/stats/)', x, perl=TRUE))
  }
  if(length(m)==0 || m=='') return(tolower(x))
  gsub('[()_ ]', '', tolower(m))
}
rownames(Z_all) <- sapply(rownames(Z_all), extract_id)
rownames(Q_all) <- sapply(rownames(Q_all), extract_id)
cat('Row name example:', head(rownames(Z_all), 3), '\n')

# 保存
write.csv(Z_all, output_file)
write.csv(Q_all, output_quant)
cat('Saved to:\n', output_file, '\n', output_quant, '\n')
