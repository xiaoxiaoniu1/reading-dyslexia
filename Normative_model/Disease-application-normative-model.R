rm(list=ls())
#set your own directory
datapath='/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/Source-codes'  # Change the directory where you save the Source-codes  
clinical_datapath='/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Datasets/Datasets-diseases/Clinical_vars_DD.csv'  # Change the directory where you save the clinical variables  
MR_datapath='/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Datasets/Datasets-diseases/MR_measures_DD.xlsx'  # Change the directory where you save the MR measures  
savepath1='/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Results/Newsite' #Chnage the path you have saved your normative models in "Normative-model-fit.R"; here we used our models for example
savepath='/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Results/Disease';# Create and determine the directory where you would save the results  
if (!(dir.exists(savepath)))
{dir.create(savepath, recursive = TRUE)}
setwd(datapath)
source("100.common-variables.r")
source("101.common-functions.r")
source("102.gamlss-recode.r")
source("ZZZ_function.R")
source("300.variables.r")
source("301.functions.r")

library(readxl)
library(dplyr)
library(stringr)


var<-c('Global.table',
       'aseg.vol.table',
       'lh.aparc.volume.table','rh.aparc.volume.table',
       'lh.aparc.thickness.table','rh.aparc.thickness.table',
       'lh.aparc.area.table','rh.aparc.area.table')

for(sheet in var) 
{ 
  
  setwd(datapath)
  MRI <- read_excel(MR_datapath,sheet=sheet)

  MRI<-as.data.frame(MRI)
  
  if('Freesurfer_Path2' %in% colnames(MRI)) {
    rownames(MRI)<-paste0(MRI$Freesurfer_Path2,MRI$Freesurfer_Path3)
  } else if('Freesufer_Path2' %in% colnames(MRI)) {
    rownames(MRI)<-paste0(MRI$Freesufer_Path2,MRI$Freesufer_Path3)
  } else if('case_dir' %in% colnames(MRI)) {
    rownames(MRI)<-MRI$case_dir
  } else if('Measure.volume' %in% colnames(MRI)) {
    rownames(MRI)<-MRI$Measure.volume
  } else {
    rownames(MRI)<-as.character(MRI[,1])
  }
  if(str_detect(sheet,'aparc'))
  {
    tem_feature<-colnames(MRI)[c(2:35)];
  }
  
  str=sheet;
  
  
  if(sheet=="aseg.vol.table")
  {
    MRI[,'cerebellum_WM']<-MRI$Left.Cerebellum.White.Matter+MRI$Right.Cerebellum.White.Matter
    MRI[,'cerebellum_GM']<-MRI$Left.Cerebellum.Cortex+MRI$Right.Cerebellum.Cortex
    MRI[,'cerebellum_total']<-MRI[,'cerebellum_WM']+MRI[,'cerebellum_GM'];
    MRI[,'CC']<-MRI$CC_Anterior+MRI$CC_Central+MRI$CC_Mid_Anterior+
      MRI$CC_Mid_Posterior+MRI$CC_Posterior
    
    tem_feature<-colnames(MRI)[c(6:9,12:14,16:17,24:31,69:71)];
    
    
    str=sheet;
    
  }
  
  if(sheet=='Global.table')
  {
    
    
    for(col in colnames(MRI)[2:19])
    {
      MRI[,col]<-as.numeric(MRI[,col])
    }
    
    for(i in 1:dim(MRI)[1])
    {
      MRI[i,'mean_thickness']<-
        (MRI[i,'lhMeanThickness']*MRI[i,'lhVertex']+
           MRI[i,'rhMeanThickness']*MRI[i,'rhVertex'])/(MRI[i,'lhVertex']+MRI[i,'rhVertex'])
      
      MRI[i,'total_surface_arrea']<-MRI[i,'lh_totaISA2']+MRI[i,'rh_totaISA2'];
    }
    
    tem_feature<-colnames(MRI)[c(2,3,4,5,13,23,24)];
    
    str='Global_feature';
    
  }
  
  
  if (!(dir.exists(paste0(savepath,'/',str))))
  {dir.create(paste0(savepath,'/',str), recursive = TRUE)}
  
  setwd(paste0(savepath,'/',str))
  
  Z_data<-list();
  Quant_data<-list()
  
  
  for(i in tem_feature[1:length(tem_feature)])
  {
    
    print(i)
    setwd(paste0(savepath1,'/',str))
    
    
    rdsfile<-paste0(str,'_',i,'_loop_our_model_new_site_calibrated.rds')
    
    if(i=='Brain.Stem'){rdsfile=paste0(str,'_BrainStem_loop_our_model_new_site_calibrated.rds')}
    if(i=='cerebellum_total'){rdsfile=paste0(str,'_Cerebellum_loop_our_model_new_site_calibrated.rds')}
    if(i=='mean_thickness'){rdsfile=paste0(str,'_mean_thickness_loop_our_model_new_site_calibrated.rds')}
    if(i=='total_surface_arrea'){rdsfile=paste0(str,'_total_surface_arrea_loop_our_model_new_site_calibrated.rds')}

    if(file.exists(rdsfile)){
      
    results<-readRDS(rdsfile)
    
    #for each feature, we should load the oirginal clinical information
    setwd(datapath)  
    data1<-read.csv(clinical_datapath,header=TRUE);
    
    # 如果存在 group_d_or_c，则按其编码过滤目标人群；
    # 当前 Clinical_vars_control.csv 不含该列，因此默认对表中的所有被试（TD 和 DD）应用校准后的模型
    if('group_d_or_c' %in% colnames(data1)) {
      data1 <- data1[data1$group_d_or_c == 0, ]
      cat('Filtered to', nrow(data1), 'disease subjects\n')
    }
    
    data1$Site_ZZZ<-data1$Center
    if('Freesurfer_Path2' %in% colnames(data1)) {
      rownames(data1)<-paste0(data1$Freesurfer_Path2,data1$Freesurfer_Path3)
    } else if('Freesufer_Path2' %in% colnames(data1)) {
      rownames(data1)<-paste0(data1$Freesufer_Path2,data1$Freesufer_Path3)
    } else {
      rownames(data1)<-as.character(data1$ID)
    }
    data1$Site_ZZZ<-as.factor(data1$Site_ZZZ)
    
    setwd(paste0(savepath,'/',str))
   
    inter_row<-intersect(rownames(data1),rownames(MRI))
    cat('inter_row length:', length(inter_row), '\n')
    cat('data1 NA in Age:', sum(is.na(data1$Age)), '\n')
    cat('data1 NA in Sex:', sum(is.na(data1$Sex)), '\n')
    cat('data1 Site_ZZZ levels:', length(unique(data1$Site_ZZZ)), '\n')
    if(length(inter_row)==0){cat('WARNING: no intersection, skipping\n'); next}
    data1=cbind(data1[inter_row,],MRI[inter_row,i])
  
    colnames(data1)[dim(data1)[2]]=c('tem_feature')
    
    
    data1$Site_ZZZ<-as.factor(data1$Site_ZZZ)
    
    data1$Sex<-as.factor(data1$Sex)
    
    data1$Sex<-factor(data1$Sex,levels=c('Female','Male'))
    
    
    data1<-data1[order(data1$Age),]
    data1[,'feature']<-data1$tem_feature
    
    

    data1<-data1[,c('Age','Sex','Site_ZZZ','feature','tem_feature','Diagnosis')]
    
    all_data1<-data1
    all_data1$Sex<-as.factor(all_data1$Sex)
    all_data1$Site_ZZZ<-as.factor(all_data1$Site_ZZZ)
    
  
 
    m2<-results
    # 读取原始模型文件以获取 m0（包含 bfpNA powers）
    orig_model_base <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Models/GAMLSS/DK"
    orig_rdsfile <- file.path(orig_model_base, str, paste0(i, "_normative_model.rds"))
    if(i=="Brain.Stem") orig_rdsfile <- file.path(orig_model_base, str, "BrainStem_normative_model.rds")
    if(i=="cerebellum_total") orig_rdsfile <- file.path(orig_model_base, str, "Cerebellum_normative_model.rds")
    if(i=="mean_thickness") orig_rdsfile <- file.path(orig_model_base, str, "Cortical_thickness_normative_model.rds")
    if(i=="total_surface_arrea") orig_rdsfile <- file.path(orig_model_base, str, "Total_surface_area_normative_model.rds")
    if(file.exists(orig_rdsfile)) {
      orig_model <- readRDS(orig_rdsfile)
      m0 <- orig_model$m0
    } else {
      m0 <- results
    }
    model1<-m2;
    
    # 从模型的 mu.x 推断 bfpNA 的列数，覆盖全局 bfpNA 函数
    mu_n_powers <- max(2, length(grep("^bfpNA", colnames(model1$mu.x))))
    bfpNA <- local({
      n <- mu_n_powers
      function(x, powers=NULL, shift=0, scale=1) {
        if(is.null(powers) || length(powers)==0) powers <- seq_len(n)
        nobs <- length(x)
        X <- matrix(0, nrow=nobs, ncol=length(powers))
        for(j in seq_along(powers)) X[,j] <- x^powers[j]
        X
      }
    })
    
    fix_formula <- function(f, m0) {
      f_str <- paste(deparse(f), collapse=" ")
      if(!is.null(m0$mu.coefSmo[[1]]$power)) {
        old_str <- paste0("c(m0$mu.coefSmo[[1]]$power)")
        new_str <- paste0("c(", paste(m0$mu.coefSmo[[1]]$power, collapse=","), ")")
        f_str <- gsub(old_str, new_str, f_str, fixed=TRUE)
      }
      if(!is.null(m0$sigma.coefSmo[[1]]$power)) {
        old_str <- paste0("c(m0$sigma.coefSmo[[1]]$power)")
        new_str <- paste0("c(", paste(m0$sigma.coefSmo[[1]]$power, collapse=","), ")")
        f_str <- gsub(old_str, new_str, f_str, fixed=TRUE)
      }
      as.formula(f_str)
    }
    # 将 m0 注入公式环境，让 bfpNA 能正确解析 powers
    environment(model1$mu.formula) <- environment()
    environment(model1$sigma.formula) <- environment()
    environment(model1$nu.formula) <- environment()
    
    Z_score_sum<-NULL;
    Quant_score_sum<-NULL
    
    cat('all_data1 nrow:', nrow(all_data1), '\n')
    cat('Age range:', min(all_data1$Age,na.rm=TRUE), '-', max(all_data1$Age,na.rm=TRUE), '\n')
    cat('tem_feature NA:', sum(is.na(all_data1$tem_feature)), '\n')
    cat('mu.formula:', deparse(model1$mu.formula), '\n')
    cat('mu.coefSmo[[1]]$power:', model1$mu.coefSmo[[1]]$power, '\n')
    
    tem_rnd<-matrix(0,dim(all_data1)[1],1);
    tryCatch({
      Model.Frame<-model.frame(formula = model1$mu.formula,data=all_data1);
    }, error = function(e) {
      cat('Error in model.frame:', e$message, '\n')
      stop(e)
    })
    Model.Matrix<-model.matrix(model1$mu.formula,Model.Frame)
    Fit.fix<-matrix(model1$mu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
    # 修复：bfpNA 返回多列，列名带数字后缀，需要用部分匹配填充 NA
    if(any(is.na(Fit.fix))) {
      for(cn in rownames(Fit.fix)[is.na(Fit.fix)]) {
        base_cn <- sub('[0-9]+$', '', cn)
        matched <- names(model1$mu.coefficients)[startsWith(names(model1$mu.coefficients), base_cn)]
        if(length(matched)>0) Fit.fix[cn,] <- model1$mu.coefficients[matched[1]]
      }
    }
    if(!is.null(model1$mu.coefSmo[[1]]))
    {Fit.fix[length(Fit.fix)]=0;
    for(iz in 1:dim(all_data1)[1]){
      if(all_data1$Site_ZZZ[iz] %in% names(model1$mu.coefSmo[[1]]$coef)){
        tem_rnd[iz]<-model1$mu.coefSmo[[1]]$coef[as.character(all_data1$Site_ZZZ[iz])]
      }else{tem_rnd[iz]<-mean(model1$mu.coefSmo[[1]]$coef)}
    }
    }else{print('no random effects for this term')};
    
    
    mu<-exp(as.vector(Model.Matrix %*% Fit.fix)+as.vector(tem_rnd))
    
    

    tem_rnd<-matrix(0,dim(all_data1)[1],1);
    Model.Frame<-model.frame(formula = model1$sigma.formula,data=all_data1);
    Model.Matrix<-model.matrix(model1$sigma.formula,Model.Frame)
    Fit.fix<-matrix(model1$sigma.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
    if(any(is.na(Fit.fix))) {
      for(cn in rownames(Fit.fix)[is.na(Fit.fix)]) {
        base_cn <- sub('[0-9]+$', '', cn)
        matched <- names(model1$sigma.coefficients)[startsWith(names(model1$sigma.coefficients), base_cn)]
        if(length(matched)>0) Fit.fix[cn,] <- model1$sigma.coefficients[matched[1]]
      }
    }
    if(!is.null(model1$sigma.coefSmo[[1]]))
    {Fit.fix[length(Fit.fix)]=0;
    for(iz in 1:dim(all_data1)[1]){
      if(all_data1$Site_ZZZ[iz] %in% names(model1$sigma.coefSmo[[1]]$coef)){
        tem_rnd[iz]<-model1$sigma.coefSmo[[1]]$coef[as.character(all_data1$Site_ZZZ[iz])]
      }else{tem_rnd[iz]<-mean(model1$sigma.coefSmo[[1]]$coef)}
    }
    }else{print('no random effects for this term')};
    
    sigma<-exp(as.vector(Model.Matrix %*% Fit.fix)+as.vector(tem_rnd))
    cat('sigma NA count:', sum(is.na(sigma)), 'min:', min(sigma, na.rm=TRUE), '\n')
    cat('sigma Fit.fix NA:', sum(is.na(Fit.fix)), '\n')
    cat('sigma colnames(Model.Matrix):', colnames(Model.Matrix), '\n')
    cat('sigma.coefficients names:', names(model1$sigma.coefficients), '\n')

    

    Model.Frame<-model.frame(formula = model1$nu.formula,data=all_data1);
    Model.Matrix<-model.matrix(model1$nu.formula,Model.Frame)
    Fit.fix<-matrix(model1$nu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
    if(!is.null(model1$nu.coefSmo[[1]]))
    {Fit.fix[length(Fit.fix)]=0;
    }else{print('no random effects for this term')};
    
    nu<-as.vector(Model.Matrix[,colnames(Model.Matrix)[1]]*Fit.fix[1])
    

    
    if(length(mu)!=dim(all_data1)[1])
    {
      print("Error, Please Check Data!!!")
    }
    
    
    Z_score_sum<-zzz_cent(obj=model1,type=c("z-scores"),mu=mu,sigma=sigma,nu=nu,
                          xname = 'Age',xvalues=all_data1$Age,yval=all_data1$tem_feature,
                          calibration=FALSE,lpar=3)
    
    Quant_score_sum<-zzz_cent(obj=model1,type=c("z-scores"),mu=mu,sigma=sigma,nu=nu,
                              xname = 'Age',xvalues=all_data1$Age,yval=all_data1$tem_feature,
                              calibration=FALSE,lpar=3,cdf=TRUE)
    
    
    Z_score_sum<-data.frame(Z_score_sum);
    colnames(Z_score_sum)<-c('Z_score');
    rownames(Z_score_sum)<-rownames(all_data1)
    
    Quant_score_sum<-data.frame(Quant_score_sum);
    colnames(Quant_score_sum)<-c('Quant_score');
    rownames(Quant_score_sum)<-rownames(all_data1)
    
    Z_data[[i]]<-Z_score_sum
    Quant_data[[i]]<-Quant_score_sum
    
    
    results_new<-list();
    results_new$m2<-m2
    results_new$Zscore<-Z_data
    results_new$Quant_data<-Quant_data
    results_new$data1<-data1
    results_new$all_data<-all_data1
    results_new$str<-str
    results_new$i<-i

    
    setwd(paste0(savepath,'/',str))
    
    saveRDS(results_new,paste0(str,'_',i,'_model_new.rds'))
    
  }
  
}

}
