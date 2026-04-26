rm(list=ls())
library(ggplot2)
library(dplyr)
library(tidyr)

zscore_file <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Results/Disease_area_volume/All_Zscore_summary.csv'
clinical_file <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Datasets/Datasets-diseases/Clinical_vars_control.csv'
out_dir <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Results/Disease_area_volume'

Z <- read.csv(zscore_file, row.names=1, check.names=FALSE)
Z[is.infinite(as.matrix(Z))] <- NA
clin <- read.csv(clinical_file)

extract_id <- function(x) {
  m <- regmatches(x, regexpr('(?<=fs_subjects_all/)[^/]+', x, perl=TRUE))
  if(length(m)==0 || m=='') {
    m <- regmatches(x, regexpr('[^/]+(?=/stats/)', x, perl=TRUE))
  }
  if(length(m)==0 || m=='') return(tolower(x))
  gsub('[()_ ]', '', tolower(m))
}
rownames(Z) <- sapply(rownames(Z), extract_id)
rownames(clin) <- gsub('[()_ ]', '', tolower(as.character(clin$ID)))

common_ids <- intersect(rownames(Z), rownames(clin))
Z <- Z[common_ids, ]
diag_group <- clin[common_ids, 'Diagnosis']
cat('DD:', sum(diag_group=='DD'), ' TD:', sum(diag_group=='TD'), '\n')

Z_DD <- Z[diag_group=='DD', ]
Z_TD <- Z[diag_group=='TD', ]

# 分类型做 t-test 和 FDR
get_feature_type <- function(f) {
  # Thickness
  if(grepl('_thickness$', f) || f == 'mean_thickness') return('thickness')
  # Area
  if(grepl('_area$', f) || grepl('_arrea$', f) || f == 'total_surface_area' || f == 'total_surface_arrea') return('area')
  # Volume
  if(grepl('_volume$', f)) return('volume')
  
  return('global')
}
feature_types <- sapply(colnames(Z), get_feature_type)

pval <- sapply(colnames(Z), function(f) {
  tryCatch(t.test(Z_DD[,f], Z_TD[,f])$p.value, error=function(e) NA)
})

# 按类型分别做 FDR 校正
pval_fdr <- rep(NA, length(pval))
for(type in c('volume','thickness','area','global')) {
  idx <- which(feature_types == type)
  if(length(idx) > 0) {
    pval_fdr[idx] <- p.adjust(pval[idx], method='fdr')
    sig <- sum(pval_fdr[idx] < 0.05, na.rm=TRUE)
    cat(type, ':', length(idx), 'features,', sig, 'significant (FDR<0.05)\n')
  }
}

result_df <- data.frame(
  Feature=colnames(Z),
  Type=feature_types,
  Mean_DD=colMeans(Z_DD, na.rm=TRUE),
  Mean_TD=colMeans(Z_TD, na.rm=TRUE),
  Mean_diff=colMeans(Z_DD, na.rm=TRUE) - colMeans(Z_TD, na.rm=TRUE),
  p_value=pval,
  p_fdr=pval_fdr,
  stringsAsFactors=FALSE
)
result_df <- result_df[order(result_df$p_value), ]
write.csv(result_df, file.path(out_dir, 'DD_vs_TD_Zscore_comparison_by_type.csv'), row.names=FALSE)

# 按四种类型分别输出：完整结果 + FDR显著结果
for(type in c('volume','thickness','area','global')) {
  sub_df <- result_df[result_df$Type == type, ]
  sub_df <- sub_df[order(sub_df$p_value), ]
  write.csv(sub_df, file.path(out_dir, paste0('DD_vs_TD_', type, '_all_results.csv')), row.names=FALSE)

  sig_df <- sub_df[!is.na(sub_df$p_fdr) & sub_df$p_fdr < 0.05, ]
  write.csv(sig_df, file.path(out_dir, paste0('DD_vs_TD_', type, '_FDRsig_results.csv')), row.names=FALSE)

  cat(type, 'all:', nrow(sub_df), 'features; FDR sig:', nrow(sig_df), '\n')
}

# Top20 箱线图（全类型混合）
top20 <- head(result_df$Feature, 20)
Z_long <- Z[, top20, drop=FALSE]
Z_long$Diagnosis <- diag_group
Z_long <- pivot_longer(Z_long, cols=all_of(top20), names_to='Feature', values_to='Z_score')
Z_long$Feature <- factor(Z_long$Feature, levels=top20)

p <- ggplot(Z_long, aes(x=Feature, y=Z_score, fill=Diagnosis)) +
  geom_boxplot(outlier.size=0.5, alpha=0.8) +
  geom_hline(yintercept=c(-1.96, 1.96), linetype='dashed', color='red', linewidth=0.4) +
  geom_hline(yintercept=0, color='gray40', linewidth=0.4) +
  scale_fill_manual(values=c('DD'='#E74C3C', 'TD'='#3498DB')) +
  labs(title='Top 20 Brain Regions: DD vs TD Z-scores (FDR by type)', x='', y='Z-score', fill='Group') +
  theme_bw(base_size=11) +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=8),
        plot.title=element_text(hjust=0.5))

ggsave(file.path(out_dir, 'DD_vs_TD_Top20_Zscore_boxplot_by_type.png'), p, width=14, height=6, dpi=150)
cat('Done!\n')
