rm(list=ls())
library(readxl)

datapath='/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/Source-codes'
clinical_datapath='/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Datasets/Datasets-new/Clinical_vars.csv'
MR_datapath='/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Datasets/Datasets-new/MR_measures.xlsx'
modelpath<-'/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Models/GAMLSS/DK'
savepath='/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Results/Newsite_area_volume'

if (!(dir.exists(savepath)))
{dir.create(savepath, recursive = TRUE)}

setwd(datapath)
source("100.common-variables.r")
source("101.common-functions.r")
source("ZZZ_function.R")
source("300.variables.r")
source("301.functions.r")

set_rownames_from_cols <- function(df) {
  if ('ID' %in% colnames(df)) {
    rownames(df) <- as.character(df$ID)
    return(df)
  }
  if ('case_dir' %in% colnames(df)) {
    rownames(df) <- as.character(df$case_dir)
    return(df)
  }
  if (ncol(df) >= 1) {
    first_col <- df[[1]]
    if ((is.character(first_col) || is.factor(first_col)) && sum(is.na(first_col)) == 0) {
      if (length(unique(first_col)) == nrow(df)) {
        rownames(df) <- as.character(first_col)
        return(df)
      }
    }
  }
  if (all(c('Freesufer_Path2','Freesufer_Path3') %in% colnames(df))) {
    rownames(df) <- paste0(df$Freesufer_Path2, df$Freesufer_Path3)
    return(df)
  }
  if (all(c('Freesurfer_Path2','Freesurfer_Path3') %in% colnames(df))) {
    rownames(df) <- paste0(df$Freesurfer_Path2, df$Freesurfer_Path3)
    return(df)
  }
  df
}

get_global_covariates <- function(MR_datapath) {
  g <- read_excel(MR_datapath, sheet = "Global.table")
  g <- as.data.frame(g)
  g <- set_rownames_from_cols(g)
  if (all(c('lh_totaISA2', 'rh_totaISA2') %in% colnames(g))) {
    g$total_surface_arrea <- g$lh_totaISA2 + g$rh_totaISA2
  }
  area_candidates <- c("total_surface_arrea", "total_surface_area", "Total_surface_area", "TotalSurfaceArea")
  vol_candidates <- c("TCV", "eTIV", "ICV", "EstimatedTotalIntraCranialVol", "EstimatedTotalIntraCranialVolume")
  area_col <- intersect(area_candidates, colnames(g))
  vol_col <- intersect(vol_candidates, colnames(g))
  list(
    data = g,
    area_col = if (length(area_col) > 0) area_col[1] else NULL,
    vol_col = if (length(vol_col) > 0) vol_col[1] else NULL
  )
}

global_cov <- get_global_covariates(MR_datapath)

var<-c('Global.table',
       'aseg.vol.table',
       'lh.aparc.volume.table','rh.aparc.volume.table',
       'lh.aparc.thickness.table','rh.aparc.thickness.table',
       'lh.aparc.area.table','rh.aparc.area.table')

for(sheet in var)
{
  setwd(datapath)
  MRI <- read_excel(MR_datapath,sheet=sheet)
  MRI <- as.data.frame(MRI)
  MRI <- set_rownames_from_cols(MRI)
  meta_cols <- c('ID','case_dir','Freesufer_Path1','Freesufer_Path2','Freesufer_Path3','Freesurfer_Path1','Freesurfer_Path2','Freesurfer_Path3','euler_number_l','euler_number_r','lhSurfaceHoles','rhSurfaceHoles')
  tem_feature <- setdiff(colnames(MRI), meta_cols)
  str <- sheet

  if(sheet=="aseg.vol.table")
  {
    if(all(c('Left.Cerebellum.White.Matter','Right.Cerebellum.White.Matter','Left.Cerebellum.Cortex','Right.Cerebellum.Cortex') %in% colnames(MRI)))
    {
      MRI[,'cerebellum_WM']<-MRI$Left.Cerebellum.White.Matter+MRI$Right.Cerebellum.White.Matter
      MRI[,'cerebellum_GM']<-MRI$Left.Cerebellum.Cortex+MRI$Right.Cerebellum.Cortex
      MRI[,'cerebellum_total']<-MRI[,'cerebellum_WM']+MRI[,'cerebellum_GM']
    }
    if(all(c('CC_Anterior','CC_Central','CC_Mid_Anterior','CC_Mid_Posterior','CC_Posterior') %in% colnames(MRI)))
    {
      MRI[,'CC']<-MRI$CC_Anterior+MRI$CC_Central+MRI$CC_Mid_Anterior+MRI$CC_Mid_Posterior+MRI$CC_Posterior
    }
    tem_feature<-setdiff(colnames(MRI),meta_cols)
    str=sheet
  }

  if(sheet=='Global.table')
  {
    for(col in setdiff(colnames(MRI),meta_cols))
    {
      MRI[,col]<-as.numeric(MRI[,col])
    }

    for(i in 1:dim(MRI)[1])
    {
      if(all(c('lhMeanThickness','rhMeanThickness','lhVertex','rhVertex') %in% colnames(MRI)))
      {
        MRI[i,'mean_thickness']<-
          (MRI[i,'lhMeanThickness']*MRI[i,'lhVertex']+
             MRI[i,'rhMeanThickness']*MRI[i,'rhVertex'])/(MRI[i,'lhVertex']+MRI[i,'rhVertex'])
      } else if(all(c('lhMeanThickness','rhMeanThickness') %in% colnames(MRI)))
      {
        MRI[i,'mean_thickness']<-(MRI[i,'lhMeanThickness']+MRI[i,'rhMeanThickness'])/2
      }
      if(all(c('lh_totaISA2','rh_totaISA2') %in% colnames(MRI)))
      {
        MRI[i,'total_surface_arrea']<-MRI[i,'lh_totaISA2']+MRI[i,'rh_totaISA2']
      }
    }

    tem_feature<-intersect(c('GMV','sGMV','WMV','Ventricles','TCV','mean_thickness','total_surface_arrea'),colnames(MRI))
    str='Global_feature'
  }

  if (!(dir.exists(paste0(savepath,'/',str))))
  {dir.create(paste0(savepath,'/',str), recursive = TRUE)}

  setwd(paste0(savepath,'/',str))

  Z_data<-list()
  Quant_data<-list()

  for(i in tem_feature[1:length(tem_feature)])
  {
    print(i)
    setwd(paste0(modelpath,'/',str))

    rdsfile<-paste0(i,'_normative_model.rds')
    if(i=='Brain.Stem'){rdsfile='BrainStem_normative_model.rds'}
    if(i=='cerebellum_total'){rdsfile='Cerebellum_normative_model.rds'}
    if(i=='mean_thickness'){rdsfile='Cortical_thickness_normative_model.rds'}
    if(i=='total_surface_arrea'){rdsfile='Total_surface_area_normative_model.rds'}

    if(file.exists(rdsfile)){
      print('file exist')
      results<-readRDS(rdsfile)
    } else {next}

    cal_outfile <- paste0(savepath,'/',str,'/',str,'_',i,'_loop_our_model_new_site_calibrated.rds')
    if(file.exists(cal_outfile)){print(paste0('calibrated file already exists, skipping: ',i)); next}

    m2<-results$m2
    m0<-results$m0

    setwd(datapath)
    data1<-read.csv(clinical_datapath,header=TRUE)
    data1 <- data1[toupper(trimws(as.character(data1$Diagnosis))) == 'TD',]
    data1$Site_ZZZ<-data1$Center

    library(dplyr)
    data1 <- set_rownames_from_cols(data1)

    setwd(paste0(savepath,'/',str))
    inter_row<-intersect(rownames(data1),rownames(MRI))
    data1=cbind(data1[inter_row,],MRI[inter_row,i])
    if(dim(data1)[1]==0){next}
    colnames(data1)[dim(data1)[2]]=c('tem_feature')

    if (!is.null(global_cov$area_col)) {
      data1$Area_cov <- as.numeric(global_cov$data[inter_row, global_cov$area_col])
    }
    if (!is.null(global_cov$vol_col)) {
      data1$Volume_cov <- as.numeric(global_cov$data[inter_row, global_cov$vol_col])
    }

    data1$Site_ZZZ<-as.factor(data1$Site_ZZZ)
    data1$Sex<-as.factor(data1$Sex)
    data1$Sex<-factor(data1$Sex,levels=c('Female','Male'))

    data1<-data1[order(data1$Age),]
    data1[,'feature']<-data1$tem_feature

    data1<-data1[!is.na(data1$tem_feature),]
    data1<-data1[data1$feature>(mean(data1$feature)-3*sd(data1$feature))&
                   data1$feature<(mean(data1$feature)+3*sd(data1$feature)),]

    keep_cols <- c('Age','Sex','Site_ZZZ','tem_feature','feature')
    extra_cov <- intersect(c('Area_cov','Volume_cov'), colnames(data1))
    data1 <- data1[, c(keep_cols, extra_cov)]

    data2<-data1

    tem_rnd<-matrix(0,dim(data2)[1],1)
    Model.Frame<-model.frame(formula = m2$mu.formula,data=data2)
    Model.Matrix<-model.matrix(m2$mu.formula,Model.Frame)
    Fit.fix<-matrix(m2$mu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))

    if(!is.null(m2$mu.coefSmo[[1]]))
    {Fit.fix[length(Fit.fix)]=0;
    for(iz in 1:dim(data2)[1]){
      tem_rnd[iz]<-0
    }
    }else{print('no random effects for this term')}

    mufix<-as.vector(Model.Matrix %*% Fit.fix)

    tem_rnd<-matrix(0,dim(data2)[1],1)
    Model.Frame<-model.frame(formula = m2$sigma.formula,data=data2)
    Model.Matrix<-model.matrix(m2$sigma.formula,Model.Frame)
    Fit.fix<-matrix(m2$sigma.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
    if(!is.null(m2$sigma.coefSmo[[1]]))
    {Fit.fix[length(Fit.fix)]=0;
    for(iz in 1:dim(data2)[1]){
      tem_rnd[iz]<-0
    }
    }else{print('no random effects for this term')}

    sigmafix<-as.vector(Model.Matrix %*% Fit.fix)

    data2[,'mufix']<-mufix
    data2[,'sigmafix']<-sigmafix

    data2 <- na.omit(data2)

    con=gamlss.control()
    m2_cal <- tryCatch({
      cal_cov <- intersect(c('Area_cov','Volume_cov'), colnames(data2))
      if (i %in% c('TCV')) {
        cal_cov <- setdiff(cal_cov, 'Volume_cov')
      }
      if (i %in% c('total_surface_arrea', 'total_surface_area', 'Total_surface_area', 'TotalSurfaceArea')) {
        cal_cov <- setdiff(cal_cov, 'Area_cov')
      }
      mu_rhs <- c('offset(mufix)', cal_cov, 'random(Site_ZZZ)')
      sigma_rhs <- c('offset(sigmafix)', cal_cov, 'random(Site_ZZZ)')
      mu_cal_formula <- as.formula(paste('feature ~', paste(mu_rhs, collapse = ' + ')))
      sigma_cal_formula <- as.formula(paste('feature ~', paste(sigma_rhs, collapse = ' + ')))
      gamlss(formula=mu_cal_formula,
             sigma.formula = sigma_cal_formula,
             control=con,
             family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
             data=data2)
    }, error = function(e) {
      message(paste0('WARNING: gamlss calibration failed for ', i, ', keeping original model. Error: ', e$message))
      NULL
    })

    if(is.null(m2_cal)){
      setwd(paste0(savepath,'/',str))
      saveRDS(m2, paste0(str,'_',i,'_loop_our_model_new_site_calibrated.rds'))
      next
    }

    model_cal<-m2
    for(new_site in names(m2_cal$mu.coefSmo[[1]][1]$coef))
    {
      if(!is.null(m2$mu.coefSmo[[1]][1]$coef))
      {
        model_cal$mu.coefSmo[[1]][1]$coef[new_site]<-m2_cal$mu.coefSmo[[1]][1]$coef[new_site]
      }
      if(!is.null(m2$sigma.coefSmo[[1]][1]$coef))
      {
        model_cal$sigma.coefSmo[[1]][1]$coef[new_site]<-m2_cal$sigma.coefSmo[[1]][1]$coef[new_site]
      }
    }

    setwd(paste0(savepath,'/',str))
    saveRDS(model_cal, paste0(str,'_',i,'_loop_our_model_new_site_calibrated.rds'))
  }
}
