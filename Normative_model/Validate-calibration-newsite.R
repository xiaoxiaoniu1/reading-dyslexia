rm(list=ls())
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(gamlss)

# ============================================================
# Paths - adjust as needed
# ============================================================
datapath      <- '/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/Source-codes'
clinical_path <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Datasets/Datasets-new/Clinical_vars.csv'
MR_path       <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Datasets/Datasets-new/MR_measures.xlsx'
modelpath_raw <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Models/GAMLSS/DK'  # original uncalibrated models
modelpath_cal <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Results/Newsite'   # calibrated models
savepath      <- '/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Results/Validation'

if(!dir.exists(savepath)) dir.create(savepath, recursive=TRUE)

setwd(datapath)
source("100.common-variables.r")
source("101.common-functions.r")
source("102.gamlss-recode.r")
source("ZZZ_function.R")
source("300.variables.r")
source("301.functions.r")

# Dummy variables required by zzz_cent default arguments
m0 <- NULL; m2 <- NULL
res <- list(mu=NULL, sigma=NULL, nu=NULL)
data4 <- data.frame(Age=NULL)

# ============================================================
# Sheet / feature config (same as calibration script)
# ============================================================
var <- c('Global.table',
         'aseg.vol.table',
         'lh.aparc.volume.table','rh.aparc.volume.table',
         'lh.aparc.thickness.table','rh.aparc.thickness.table',
         'lh.aparc.area.table','rh.aparc.area.table')

# ============================================================
# Helper: compute Z scores given a GAMLSS model and data
# ============================================================
compute_Z <- function(model1, all_data1) {
  # Replicate exactly what Calibration-normative-model-using-new-dataset.R does
  # for computing mu/sigma on the link scale, then add site random effects.

  data2 <- all_data1

  # --- mu on link scale (population, no ranef) ---
  tem_rnd <- matrix(0, dim(data2)[1], 1)
  cat('    [compute_Z] model.frame mu...\n')
  Model.Frame <- model.frame(formula = model1$mu.formula, data = data2)
  cat('    [compute_Z] model.matrix mu...\n')
  Model.Matrix <- model.matrix(model1$mu.formula, Model.Frame)
  cat('    [compute_Z] mu colnames:', paste(head(colnames(Model.Matrix),5), collapse=', '), '\n')
  cat('    [compute_Z] mu.coefficients names:', paste(head(names(model1$mu.coefficients),5), collapse=', '), '\n')
  Fit.fix <- matrix(model1$mu.coefficients[colnames(Model.Matrix)], ncol=1,
                    dimnames=list(colnames(Model.Matrix), 'Beta'))
  if (!is.null(model1$mu.coefSmo[[1]])) {
    Fit.fix[length(Fit.fix)] <- 0
    for (iz in 1:dim(data2)[1]) {
      s <- as.character(data2$Site_ZZZ[iz])
      if (s %in% names(model1$mu.coefSmo[[1]]$coef)) {
        tem_rnd[iz] <- model1$mu.coefSmo[[1]]$coef[s]
      } else {
        tem_rnd[iz] <- mean(model1$mu.coefSmo[[1]]$coef)
      }
    }
  }
  mu <- exp(as.vector(Model.Matrix %*% Fit.fix) + as.vector(tem_rnd))

  # --- sigma on link scale (population, no ranef) ---
  tem_rnd <- matrix(0, dim(data2)[1], 1)
  cat('    [compute_Z] model.frame sigma...\n')
  Model.Frame <- model.frame(formula = model1$sigma.formula, data = data2)
  cat('    [compute_Z] model.matrix sigma...\n')
  Model.Matrix <- model.matrix(model1$sigma.formula, Model.Frame)
  cat('    [compute_Z] sigma colnames:', paste(head(colnames(Model.Matrix),5), collapse=', '), '\n')
  cat('    [compute_Z] sigma.coefficients names:', paste(head(names(model1$sigma.coefficients),5), collapse=', '), '\n')
  Fit.fix <- matrix(model1$sigma.coefficients[colnames(Model.Matrix)], ncol=1,
                    dimnames=list(colnames(Model.Matrix), 'Beta'))
  if (!is.null(model1$sigma.coefSmo[[1]])) {
    Fit.fix[length(Fit.fix)] <- 0
    for (iz in 1:dim(data2)[1]) {
      s <- as.character(data2$Site_ZZZ[iz])
      if (s %in% names(model1$sigma.coefSmo[[1]]$coef)) {
        tem_rnd[iz] <- model1$sigma.coefSmo[[1]]$coef[s]
      } else {
        tem_rnd[iz] <- mean(model1$sigma.coefSmo[[1]]$coef)
      }
    }
  }
  sigma <- exp(as.vector(Model.Matrix %*% Fit.fix) + as.vector(tem_rnd))

  # --- nu ---
  cat('    [compute_Z] model.frame nu...\n')
  Model.Frame <- model.frame(formula = model1$nu.formula, data = data2)
  Model.Matrix <- model.matrix(model1$nu.formula, Model.Frame)
  Fit.fix <- matrix(model1$nu.coefficients[colnames(Model.Matrix)], ncol=1,
                    dimnames=list(colnames(Model.Matrix), 'Beta'))
  if (!is.null(model1$nu.coefSmo[[1]])) Fit.fix[length(Fit.fix)] <- 0
  nu <- as.vector(Model.Matrix[, colnames(Model.Matrix)[1]] * Fit.fix[1])

  cat('    [compute_Z] calling zzz_cent...\n')
  Z <- zzz_cent(obj=model1, type='z-scores', mu=mu, sigma=sigma, nu=nu,
                xname='Age', xvalues=all_data1$Age, yval=all_data1$feature,
                calibration=FALSE, lpar=3)
  return(as.numeric(Z))
}

# ============================================================
# Main loop: collect Z scores before & after calibration
# ============================================================
all_results <- list()

for (sheet in var) {

  setwd(datapath)
  MRI <- read_excel(MR_path, sheet=sheet)
  MRI <- as.data.frame(MRI)

  meta_cols <- c('ID','Freesufer_Path1','Freesufer_Path2','Freesufer_Path3',
                 'Freesurfer_Path1','Freesurfer_Path2','Freesurfer_Path3',
                 'euler_number_l','euler_number_r','lhSurfaceHoles','rhSurfaceHoles')

  if ('ID' %in% colnames(MRI)) {
    rownames(MRI) <- as.character(MRI$ID)
  } else if (all(c('Freesufer_Path2','Freesufer_Path3') %in% colnames(MRI))) {
    rownames(MRI) <- paste0(MRI$Freesufer_Path2, MRI$Freesufer_Path3)
  } else if (all(c('Freesurfer_Path2','Freesurfer_Path3') %in% colnames(MRI))) {
    rownames(MRI) <- paste0(MRI$Freesurfer_Path2, MRI$Freesurfer_Path3)
  }

  tem_feature <- setdiff(colnames(MRI), meta_cols)
  str_label <- sheet

  if (sheet == 'aseg.vol.table') {
    if (all(c('Left.Cerebellum.White.Matter','Right.Cerebellum.White.Matter',
              'Left.Cerebellum.Cortex','Right.Cerebellum.Cortex') %in% colnames(MRI))) {
      MRI[,'cerebellum_WM']    <- MRI$Left.Cerebellum.White.Matter + MRI$Right.Cerebellum.White.Matter
      MRI[,'cerebellum_GM']    <- MRI$Left.Cerebellum.Cortex + MRI$Right.Cerebellum.Cortex
      MRI[,'cerebellum_total'] <- MRI[,'cerebellum_WM'] + MRI[,'cerebellum_GM']
    }
    if (all(c('CC_Anterior','CC_Central','CC_Mid_Anterior','CC_Mid_Posterior','CC_Posterior') %in% colnames(MRI))) {
      MRI[,'CC'] <- MRI$CC_Anterior + MRI$CC_Central + MRI$CC_Mid_Anterior +
        MRI$CC_Mid_Posterior + MRI$CC_Posterior
    }
    tem_feature <- setdiff(colnames(MRI), meta_cols)
  }

  if (sheet == 'Global.table') {
    for (col in setdiff(colnames(MRI), meta_cols)) MRI[,col] <- as.numeric(MRI[,col])
    for (ii in 1:nrow(MRI)) {
      if (all(c('lhMeanThickness','rhMeanThickness','lhVertex','rhVertex') %in% colnames(MRI))) {
        MRI[ii,'mean_thickness'] <- (MRI[ii,'lhMeanThickness']*MRI[ii,'lhVertex'] +
                                       MRI[ii,'rhMeanThickness']*MRI[ii,'rhVertex']) /
          (MRI[ii,'lhVertex'] + MRI[ii,'rhVertex'])
      }
      if (all(c('lh_totaISA2','rh_totaISA2') %in% colnames(MRI))) {
        MRI[ii,'total_surface_arrea'] <- MRI[ii,'lh_totaISA2'] + MRI[ii,'rh_totaISA2']
      }
    }
    tem_feature <- intersect(c('GMV','sGMV','WMV','Ventricles','TCV','mean_thickness','total_surface_arrea'), colnames(MRI))
    str_label   <- 'Global_feature'
  }

  for (i in tem_feature) {
    cat('Processing:', sheet, '-', i, '\n')

    # --- original (uncalibrated) model rds ---
    raw_rdsfile <- paste0(i, '_normative_model.rds')
    if (i == 'Brain.Stem')        raw_rdsfile <- 'BrainStem_normative_model.rds'
    if (i == 'cerebellum_total')  raw_rdsfile <- 'Cerebellum_normative_model.rds'
    if (i == 'mean_thickness')    raw_rdsfile <- 'Cortical_thickness_normative_model.rds'
    if (i == 'total_surface_arrea') raw_rdsfile <- 'Total_surface_area_normative_model.rds'
    raw_full <- file.path(modelpath_raw, str_label, raw_rdsfile)

    # --- calibrated model rds ---
    cal_rdsfile <- paste0(str_label, '_', i, '_loop_our_model_new_site_calibrated.rds')
    cal_full    <- file.path(modelpath_cal, str_label, cal_rdsfile)

    if (!file.exists(raw_full)) { cat('  raw model not found, skip\n'); next }
    if (!file.exists(cal_full)) { cat('  calibrated model not found, skip\n'); next }

    raw_results <- readRDS(raw_full)
    m2_raw <- raw_results$m2
    m2_cal <- readRDS(cal_full)

    # --- clinical data ---
    setwd(datapath)
    data1 <- read.csv(clinical_path, header=TRUE)
    data1 <- data1[toupper(trimws(as.character(data1$Diagnosis))) == 'TD', ]
    # Only use Center to define sites (ignore Province and Manufacturer)
    data1$Site_ZZZ <- data1$Center

    if ('ID' %in% colnames(data1)) {
      rownames(data1) <- as.character(data1$ID)
    } else if (all(c('Freesufer_Path2','Freesufer_Path3') %in% colnames(data1))) {
      rownames(data1) <- paste0(data1$Freesufer_Path2, data1$Freesufer_Path3)
    } else if (all(c('Freesurfer_Path2','Freesurfer_Path3') %in% colnames(data1))) {
      rownames(data1) <- paste0(data1$Freesurfer_Path2, data1$Freesurfer_Path3)
    }

    inter_row <- intersect(rownames(data1), rownames(MRI))
    if (length(inter_row) == 0) { cat('  no matched subjects, skip\n'); next }

    data1 <- cbind(data1[inter_row, ], MRI[inter_row, i])
    colnames(data1)[ncol(data1)] <- 'tem_feature'
    data1$Site_ZZZ <- as.factor(data1$Site_ZZZ)
    data1$Sex      <- factor(data1$Sex, levels=c('Female','Male'))
    data1 <- data1[order(data1$Age), ]
    data1[,'feature'] <- data1$tem_feature
    data1 <- data1[!is.na(data1$tem_feature), ]
    data1 <- data1[data1$feature > (mean(data1$feature) - 3*sd(data1$feature)) &
                     data1$feature < (mean(data1$feature) + 3*sd(data1$feature)), ]
    data1 <- data1[, c('Age','Sex','Site_ZZZ','feature','tem_feature')]
    if (nrow(data1) == 0) next

    # --- compute Z scores ---
    Z_raw <- tryCatch(compute_Z(m2_raw, data1), error=function(e){ cat('  Z_raw failed:', e$message,'\n'); NULL })
    Z_cal <- tryCatch(compute_Z(m2_cal, data1), error=function(e){ cat('  Z_cal failed:', e$message,'\n'); NULL })
    if (is.null(Z_raw) || is.null(Z_cal)) next

    all_results[[length(all_results)+1]] <- data.frame(
      sheet    = sheet,
      str      = str_label,
      feature  = i,
      subj_id  = rownames(data1),
      Age      = data1$Age,
      Sex      = as.character(data1$Sex),
      Site     = as.character(data1$Site_ZZZ),
      Z_raw    = Z_raw,
      Z_cal    = Z_cal
    )
  }
}

# ============================================================
# Merge all results
# ============================================================
df_all <- do.call(rbind, all_results)
write.csv(df_all, file.path(savepath, 'validation_Z_scores.csv'), row.names=FALSE)
cat('Total features validated:', length(unique(df_all$feature)), '\n')
cat('Total subject-feature records:', nrow(df_all), '\n')

# ============================================================
# Summary stats per feature
# ============================================================
summary_df <- df_all %>%
  group_by(str, feature) %>%
  summarise(
    n          = n(),
    mean_Z_raw = mean(Z_raw, na.rm=TRUE),
    sd_Z_raw   = sd(Z_raw,   na.rm=TRUE),
    mean_Z_cal = mean(Z_cal, na.rm=TRUE),
    sd_Z_cal   = sd(Z_cal,   na.rm=TRUE),
    .groups='drop'
  ) %>%
  mutate(
    mean_improvement = abs(mean_Z_raw) - abs(mean_Z_cal),   # positive = calibration pulled mean closer to 0
    sd_improvement   = abs(sd_Z_raw - 1) - abs(sd_Z_cal - 1) # positive = sd closer to 1 after calibration
  )
write.csv(summary_df, file.path(savepath, 'validation_summary.csv'), row.names=FALSE)

# ============================================================
# Plot 1: Mean Z score before vs after calibration per feature
# (ideal = 0)
# ============================================================
plot_mean <- summary_df %>%
  select(str, feature, mean_Z_raw, mean_Z_cal) %>%
  pivot_longer(cols=c(mean_Z_raw, mean_Z_cal),
               names_to='Model', values_to='mean_Z') %>%
  mutate(Model = recode(Model,
                        'mean_Z_raw' = 'Before calibration',
                        'mean_Z_cal' = 'After calibration'))

p1 <- ggplot(plot_mean, aes(x=reorder(feature, mean_Z), y=mean_Z, color=Model, group=Model)) +
  geom_hline(yintercept=0, linetype='dashed', color='grey50', linewidth=0.6) +
  geom_point(size=1.8, alpha=0.85) +
  facet_wrap(~str, scales='free_y', ncol=2) +
  scale_color_manual(values=c('Before calibration'='#B0B0B0', 'After calibration'='#E84935')) +
  coord_flip() +
  labs(title='Mean Z score per feature (TD validation set)',
       subtitle='Ideal value = 0',
       x='', y='Mean Z score', color='') +
  theme_bw(base_size=9) +
  theme(legend.position='bottom',
        strip.text=element_text(size=7),
        axis.text.y=element_text(size=5))

ggsave(file.path(savepath, 'plot1_mean_Z_per_feature.pdf'), p1,
       width=12, height=14)

# ============================================================
# Plot 2: SD of Z score before vs after calibration per feature
# (ideal = 1)
# ============================================================
plot_sd <- summary_df %>%
  select(str, feature, sd_Z_raw, sd_Z_cal) %>%
  pivot_longer(cols=c(sd_Z_raw, sd_Z_cal),
               names_to='Model', values_to='sd_Z') %>%
  mutate(Model = recode(Model,
                        'sd_Z_raw' = 'Before calibration',
                        'sd_Z_cal' = 'After calibration'))

p2 <- ggplot(plot_sd, aes(x=reorder(feature, sd_Z), y=sd_Z, color=Model)) +
  geom_hline(yintercept=1, linetype='dashed', color='grey50', linewidth=0.6) +
  geom_point(size=1.8, alpha=0.85) +
  facet_wrap(~str, scales='free_y', ncol=2) +
  scale_color_manual(values=c('Before calibration'='#B0B0B0', 'After calibration'='#E84935')) +
  coord_flip() +
  labs(title='SD of Z score per feature (TD validation set)',
       subtitle='Ideal value = 1',
       x='', y='SD of Z score', color='') +
  theme_bw(base_size=9) +
  theme(legend.position='bottom',
        strip.text=element_text(size=7),
        axis.text.y=element_text(size=5))

ggsave(file.path(savepath, 'plot2_sd_Z_per_feature.pdf'), p2,
       width=12, height=14)

# ============================================================
# Plot 3: Overall Z score distribution before vs after
# violin + box, collapsed across all features
# ============================================================
df_dist <- df_all %>%
  select(feature, Z_raw, Z_cal) %>%
  pivot_longer(cols=c(Z_raw, Z_cal),
               names_to='Model', values_to='Z') %>%
  mutate(Model = recode(Model,
                        'Z_raw' = 'Before calibration',
                        'Z_cal' = 'After calibration'))

p3 <- ggplot(df_dist, aes(x=Model, y=Z, fill=Model)) +
  geom_violin(alpha=0.5, trim=TRUE, linewidth=0.3) +
  geom_boxplot(width=0.15, outlier.size=0.4, outlier.alpha=0.3, linewidth=0.4) +
  geom_hline(yintercept=0, linetype='dashed', color='grey40', linewidth=0.5) +
  scale_fill_manual(values=c('Before calibration'='#B0B0B0', 'After calibration'='#E84935')) +
  coord_cartesian(ylim=c(-5, 5)) +
  labs(title='Overall Z score distribution (all features, TD validation set)',
       subtitle='Ideal: median = 0, IQR symmetric around 0',
       x='', y='Z score') +
  theme_bw(base_size=11) +
  theme(legend.position='none')

ggsave(file.path(savepath, 'plot3_overall_Z_distribution.pdf'), p3,
       width=6, height=5)

# ============================================================
# Plot 4: QQ plot - Z scores vs standard normal
# ============================================================
set.seed(42)
df_qq <- df_all %>%
  select(Z_raw, Z_cal) %>%
  pivot_longer(everything(), names_to='Model', values_to='Z') %>%
  mutate(Model = recode(Model,
                        'Z_raw' = 'Before calibration',
                        'Z_cal' = 'After calibration')) %>%
  filter(!is.na(Z), is.finite(Z)) %>%
  group_by(Model) %>%
  arrange(Z, .by_group=TRUE) %>%
  mutate(theoretical = qnorm(ppoints(n()))) %>%
  ungroup()

p4 <- ggplot(df_qq, aes(x=theoretical, y=Z, color=Model)) +
  geom_point(size=0.6, alpha=0.4) +
  geom_abline(slope=1, intercept=0, linetype='dashed', color='grey30', linewidth=0.6) +
  scale_color_manual(values=c('Before calibration'='#B0B0B0', 'After calibration'='#E84935')) +
  coord_cartesian(xlim=c(-4,4), ylim=c(-6,6)) +
  labs(title='QQ plot: Z scores vs standard normal',
       subtitle='Points on diagonal = perfect normality',
       x='Theoretical quantiles', y='Sample Z scores', color='') +
  theme_bw(base_size=11) +
  theme(legend.position='bottom')

ggsave(file.path(savepath, 'plot4_QQ_plot.pdf'), p4,
       width=6, height=5)

# ============================================================
# Plot 5: Mean Z score by site - check site bias
# ============================================================
site_df <- df_all %>%
  group_by(Site) %>%
  summarise(
    mean_Z_raw = mean(Z_raw, na.rm=TRUE),
    mean_Z_cal = mean(Z_cal, na.rm=TRUE),
    n = n(),
    .groups='drop'
  ) %>%
  pivot_longer(cols=c(mean_Z_raw, mean_Z_cal),
               names_to='Model', values_to='mean_Z') %>%
  mutate(Model = recode(Model,
                        'mean_Z_raw' = 'Before calibration',
                        'mean_Z_cal' = 'After calibration'))

p5 <- ggplot(site_df, aes(x=reorder(Site, mean_Z), y=mean_Z,
                           color=Model, group=Model)) +
  geom_hline(yintercept=0, linetype='dashed', color='grey50', linewidth=0.6) +
  geom_line(linewidth=0.5, alpha=0.6) +
  geom_point(size=3) +
  scale_color_manual(values=c('Before calibration'='#B0B0B0', 'After calibration'='#E84935')) +
  coord_flip() +
  labs(title='Mean Z score by site (TD validation set)',
       subtitle='Ideal = 0 for all sites; calibration should reduce site bias',
       x='Site', y='Mean Z score', color='') +
  theme_bw(base_size=11) +
  theme(legend.position='bottom')

ggsave(file.path(savepath, 'plot5_site_bias.pdf'), p5,
       width=8, height=max(4, length(unique(site_df$Site))*0.5 + 2))

cat('\nValidation complete. Results saved to:', savepath, '\n')
cat('Files generated:\n')
cat('  validation_Z_scores.csv   - per-subject Z scores for all features\n')
cat('  validation_summary.csv    - mean/SD per feature before & after calibration\n')
cat('  plot1_mean_Z_per_feature.pdf\n')
cat('  plot2_sd_Z_per_feature.pdf\n')
cat('  plot3_overall_Z_distribution.pdf\n')
cat('  plot4_QQ_plot.pdf\n')
cat('  plot5_site_bias.pdf\n') 