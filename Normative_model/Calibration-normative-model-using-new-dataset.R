rm(list=ls())

#set your own directory
datapath='/data/home/tqi/data1/share/after_freesurfer/CODE/Normative_model/Source-codes'  # Change the directory where you save the Source-codes  
clinical_datapath='/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Datasets/Datasets-new/Clinical_vars_HC.csv'  # Change the directory where you save the clinical variables  
MR_datapath='/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Datasets/Datasets-new/MR_measures_HC.xlsx'  # Change the directory where you save the MR measures  
modelpath<-'/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Models/GAMLSS/DK' #Chnage the path you have saved your normative models in "Normative-model-fit.R"
savepath='/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/Normative_model/Results/Newsite';# Create and determine the directory where you would save the results  



setwd(datapath)
source("100.common-variables.r")
source("101.common-functions.r")
source("ZZZ_function.R")
source("300.variables.r")
source("301.functions.r")

split_ratio <- 0.7
split_seed <- 20260415

normalize_subject_id <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- basename(x)
  x <- gsub("[^a-z0-9]", "", x)
  x
}

get_re_coef <- function(x) {
  if (is.null(x) || length(x) < 1) return(NULL)
  if (!is.null(x[[1]]$coef)) return(x[[1]]$coef)
  if (is.list(x[[1]]) && length(x[[1]]) >= 1 && !is.null(x[[1]][[1]]$coef)) return(x[[1]][[1]]$coef)
  NULL
}

has_re_coef <- function(x) {
  !is.null(get_re_coef(x))
}

debug_sep <- function(title="") {
  cat("\n", paste(rep("=", 80), collapse=""), "\n", sep="")
  if (nzchar(title)) cat(title, "\n")
}

debug_kv <- function(name, value) {
  cat(sprintf("[DEBUG] %s: %s\n", name, paste(capture.output(print(value)), collapse=" ")))
}

debug_num_summary <- function(x, label) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    cat(sprintf("[DEBUG] %s: no finite values\n", label))
    return(invisible(NULL))
  }
  cat(sprintf(
    "[DEBUG] %s: n=%d min=%.6f q25=%.6f median=%.6f mean=%.6f q75=%.6f max=%.6f sd=%.6f\n",
    label, length(x), min(x), stats::quantile(x, 0.25), stats::median(x), mean(x),
    stats::quantile(x, 0.75), max(x), stats::sd(x)
  ))
}

debug_factor_table <- function(x, label, max_n=20) {
  tab <- sort(table(as.character(x)), decreasing=TRUE)
  if (length(tab) == 0) {
    cat(sprintf("[DEBUG] %s: empty\n", label))
    return(invisible(NULL))
  }
  show_tab <- utils::head(tab, max_n)
  cat(sprintf("[DEBUG] %s (top %d):\n", label, length(show_tab)))
  print(show_tab)
  if (length(tab) > max_n) {
    cat(sprintf("[DEBUG] %s: ... %d more levels omitted\n", label, length(tab) - max_n))
  }
}

append_error_log <- function(error_csv_path, row) {
  row <- as.data.frame(row, stringsAsFactors=FALSE)
  if (file.exists(error_csv_path)) {
    old <- utils::read.csv(error_csv_path, stringsAsFactors=FALSE)
    out <- dplyr::bind_rows(old, row)
  } else {
    out <- row
  }
  utils::write.csv(out, error_csv_path, row.names=FALSE)
}

calc_mu_sigma_nu <- function(model_obj, df) {
  df$Site_ZZZ <- as.character(df$Site_ZZZ)
  df$Site_ZZZ[is.na(df$Site_ZZZ) | df$Site_ZZZ == ""] <- "UnknownSite"
  df$Site_ZZZ <- as.factor(df$Site_ZZZ)
  df$Sex <- as.character(df$Sex)
  df$Sex <- factor(df$Sex, levels=c("Female","Male"))

  get_re <- function(x) {
    if (is.null(x) || length(x) < 1) return(NULL)
    if (!is.null(x[[1]]$coef)) return(x[[1]]$coef)
    if (is.list(x[[1]]) && length(x[[1]]) >= 1 && !is.null(x[[1]][[1]]$coef)) return(x[[1]][[1]]$coef)
    NULL
  }

  mu_re <- get_re(model_obj$mu.coefSmo)
  sigma_re <- get_re(model_obj$sigma.coefSmo)

  Model.Frame <- model.frame(formula = model_obj$mu.formula, data=df)
  Model.Matrix <- model.matrix(model_obj$mu.formula, Model.Frame)
  Fit.fix <- matrix(model_obj$mu.coefficients[colnames(Model.Matrix)], ncol=1, dimnames=list(colnames(Model.Matrix),"Beta"))
  if (any(is.na(Fit.fix))) {
    for (cn in rownames(Fit.fix)[is.na(Fit.fix)]) {
      base_cn <- sub("[0-9]+$", "", cn)
      matched <- names(model_obj$mu.coefficients)[startsWith(names(model_obj$mu.coefficients), base_cn)]
      if (length(matched) > 0) Fit.fix[cn,] <- model_obj$mu.coefficients[matched[1]]
    }
  }
  tem_rnd <- rep(0, nrow(df))
  if (!is.null(mu_re)) {
    Fit.fix[length(Fit.fix)] <- 0
    for (iz in seq_len(nrow(df))) {
      site <- as.character(df$Site_ZZZ[iz])
      if (site %in% names(mu_re)) tem_rnd[iz] <- mu_re[site] else tem_rnd[iz] <- mean(mu_re)
    }
  }
  mu <- exp(as.vector(Model.Matrix %*% Fit.fix) + as.vector(tem_rnd))

  Model.Frame <- model.frame(formula = model_obj$sigma.formula, data=df)
  Model.Matrix <- model.matrix(model_obj$sigma.formula, Model.Frame)
  Fit.fix <- matrix(model_obj$sigma.coefficients[colnames(Model.Matrix)], ncol=1, dimnames=list(colnames(Model.Matrix),"Beta"))
  if (any(is.na(Fit.fix))) {
    for (cn in rownames(Fit.fix)[is.na(Fit.fix)]) {
      base_cn <- sub("[0-9]+$", "", cn)
      matched <- names(model_obj$sigma.coefficients)[startsWith(names(model_obj$sigma.coefficients), base_cn)]
      if (length(matched) > 0) Fit.fix[cn,] <- model_obj$sigma.coefficients[matched[1]]
    }
  }
  tem_rnd <- rep(0, nrow(df))
  if (!is.null(sigma_re)) {
    Fit.fix[length(Fit.fix)] <- 0
    for (iz in seq_len(nrow(df))) {
      site <- as.character(df$Site_ZZZ[iz])
      if (site %in% names(sigma_re)) tem_rnd[iz] <- sigma_re[site] else tem_rnd[iz] <- mean(sigma_re)
    }
  }
  sigma <- exp(as.vector(Model.Matrix %*% Fit.fix) + as.vector(tem_rnd))

  Model.Frame <- model.frame(formula = model_obj$nu.formula, data=df)
  Model.Matrix <- model.matrix(model_obj$nu.formula, Model.Frame)
  Fit.fix <- matrix(model_obj$nu.coefficients[colnames(Model.Matrix)], ncol=1, dimnames=list(colnames(Model.Matrix),"Beta"))
  if (any(is.na(Fit.fix))) {
    for (cn in rownames(Fit.fix)[is.na(Fit.fix)]) {
      base_cn <- sub("[0-9]+$", "", cn)
      matched <- names(model_obj$nu.coefficients)[startsWith(names(model_obj$nu.coefficients), base_cn)]
      if (length(matched) > 0) Fit.fix[cn,] <- model_obj$nu.coefficients[matched[1]]
    }
  }
  nu <- as.vector(Model.Matrix[, colnames(Model.Matrix)[1]] * Fit.fix[1])

  list(mu=mu, sigma=sigma, nu=nu)
}

var<-c('Global.table',
       'aseg.vol.table',
       'lh.aparc.volume.table','rh.aparc.volume.table',
       'lh.aparc.thickness.table','rh.aparc.thickness.table',
       'lh.aparc.area.table','rh.aparc.area.table')


for(sheet in var) #for local feature
{ 
  

    setwd(datapath)
    MRI <- readxl::read_excel(MR_datapath,sheet=sheet)
    
    MRI<-MRI[MRI$Freesurfer_Path2!='Epilepsy_dicom_info_nii'&
               MRI$Freesurfer_Path2!='ET'&
               MRI$Freesurfer_Path2!='Guojibu_HC_nii',]
    MRI<-MRI[!is.na(MRI$Freesurfer_Path3),]
    MRI<-as.data.frame(MRI)
    
    rownames(MRI)<-paste0(MRI$Freesurfer_Path2,MRI$Freesurfer_Path3)
    
  
    tem_feature<-colnames(MRI)[c(2:35)];
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
  #tem_feature<-colnames(MRI)[c(5,13,23,24)];

  str='Global_feature';#str_lab='Global_feature';
  
}

  
if (!(dir.exists(paste0(savepath,'/',str))))
{dir.create(paste0(savepath,'/',str))}

setwd(paste0(savepath,'/',str))
error_log_csv <- file.path(paste0(savepath,'/',str), paste0(str, "_feature_errors.csv"))
if (file.exists(error_log_csv)) file.remove(error_log_csv)
cat(sprintf("[DEBUG] error log csv: %s\n", error_log_csv))



Z_data<-list();
Quant_data<-list()

setwd(datapath)
clinical0<-read.csv(clinical_datapath,header=TRUE, stringsAsFactors=FALSE)
clinical0$Site_ZZZ<-paste0(clinical0$Province,clinical0$Center,clinical0$Manufacturer)

library(dplyr)

if (!("ID" %in% colnames(clinical0))) {
  stop("Clinical file missing column: ID")
}
if (!("Freesurfer_Path1" %in% colnames(MRI))) {
  stop(paste("MRI sheet missing column: Freesurfer_Path1 for", sheet))
}
clinical0$SubjectID <- normalize_subject_id(clinical0$ID)
MRI$SubjectID <- normalize_subject_id(MRI$Freesurfer_Path1)
rownames(clinical0) <- clinical0$SubjectID
rownames(MRI) <- MRI$SubjectID

inter_row_sheet<-intersect(rownames(clinical0),rownames(MRI))
cat("\n--- Alignment diagnostics ---\n")
cat("sheet:", sheet, "\n")
cat("clinical rows:", nrow(clinical0), " unique subject ids:", length(unique(rownames(clinical0))), "\n")
cat("MRI rows:", nrow(MRI), " unique subject ids:", length(unique(rownames(MRI))), "\n")
cat("intersection:", length(inter_row_sheet), "\n")
cat("clinical head ids:", paste(utils::head(rownames(clinical0), 3), collapse=" | "), "\n")
cat("MRI head ids:", paste(utils::head(rownames(MRI), 3), collapse=" | "), "\n")
if (length(inter_row_sheet) > 0) {
  cat("intersection head ids:", paste(utils::head(inter_row_sheet, 3), collapse=" | "), "\n")
}
if (length(inter_row_sheet) < 2) {
  only_clin <- setdiff(rownames(clinical0), rownames(MRI))
  only_mri <- setdiff(rownames(MRI), rownames(clinical0))
  cat("clinical-only head ids:", paste(utils::head(only_clin, 5), collapse=" | "), "\n")
  cat("MRI-only head ids:", paste(utils::head(only_mri, 5), collapse=" | "), "\n")
  stop("Too few subjects after aligning clinical and MRI data")
}

set.seed(split_seed)
split_df <- clinical0[inter_row_sheet, , drop=FALSE]
split_df$Site_ZZZ <- paste0(split_df$Province, split_df$Center, split_df$Manufacturer)
cal_ids <- c()
for (site_i in unique(split_df$Site_ZZZ)) {
  site_ids <- rownames(split_df)[split_df$Site_ZZZ == site_i]
  site_n <- max(1, floor(length(site_ids) * split_ratio))
  cal_ids <- c(cal_ids, sample(site_ids, size=site_n, replace=FALSE))
}
cal_ids <- unique(cal_ids)
split_map <- setNames(rep("holdout", length(inter_row_sheet)), inter_row_sheet)
split_map[cal_ids] <- "calibration"
write.csv(data.frame(SubjectKey=inter_row_sheet, Split=split_map[inter_row_sheet]),
          file=paste0(str,"_split_",split_ratio,"_",split_seed,".csv"),
          row.names=FALSE)


for(i in tem_feature[1:length(tem_feature)])
{

feature_step <- "feature initialization"
debug_sep(sprintf("[FEATURE START] sheet=%s | feature=%s", sheet, i))

tryCatch({
  cat(sprintf("[DEBUG] Enter feature loop for %s / %s\n", sheet, i))
  feature_step <- "load normative model"
  setwd(paste0(modelpath,'/',str))

  rdsfile<-paste0(i,'_normative_model.rds')

  if(i=='Brain.Stem'){rdsfile='BrainStem_normative_model.rds'}
  if(i=='cerebellum_total'){rdsfile='Cerebellum_normative_model.rds'}
  if(i=='mean_thickness'){rdsfile='Cortical_thickness_normative_model.rds'}
  if(i=='total_surface_arrea'){rdsfile='Total_surface_area_normative_model.rds'}

  debug_kv("model_rds", file.path(paste0(modelpath,'/',str), rdsfile))

  if(file.exists(rdsfile)){
    print('file exist')
    results<-readRDS(rdsfile)
    cat("[DEBUG] RDS loaded successfully\n")
  } else {
    stop(paste("Normative model RDS not found:", rdsfile))
  }
  
  
m2<-results$m2;
m0<-results$m0
feature_step <- "prepare feature data"
setwd(paste0(savepath,'/',str))

feature_values <- MRI[inter_row_sheet,i,drop=TRUE]
raw_feature_values <- as.numeric(feature_values)
debug_num_summary(raw_feature_values, "raw_feature_values")
data1<-data.frame(clinical0[inter_row_sheet,,drop=FALSE], tem_feature=raw_feature_values, row.names=inter_row_sheet, stringsAsFactors=FALSE)



data1$Site_ZZZ<-as.factor(data1$Site_ZZZ)

data1$Sex<-as.factor(data1$Sex)

data1$Sex<-factor(data1$Sex,levels=c('Female','Male'))


data1<-data1[order(data1$Age),]
data1[,'SubjectKey']<-rownames(data1)
data1[,'Split']<-split_map[data1$SubjectKey]
data1[,'feature']<-data1$tem_feature

cat(sprintf("[DEBUG] before filtering: n=%d\n", nrow(data1)))
debug_factor_table(data1$Split, "split counts before filtering")
debug_factor_table(data1$Site_ZZZ, "site counts before filtering")

#remove the extreme values
n_before_na <- nrow(data1)
data1<-data1[!is.na(data1$tem_feature),]
cat(sprintf("[DEBUG] removed missing feature rows: %d\n", n_before_na - nrow(data1)))
feature_mean <- mean(data1$feature)
feature_sd <- sd(data1$feature)
debug_kv("feature_mean_before_outlier_filter", feature_mean)
debug_kv("feature_sd_before_outlier_filter", feature_sd)
outlier_lower <- feature_mean - 3 * feature_sd
outlier_upper <- feature_mean + 3 * feature_sd
if (is.finite(outlier_lower) && is.finite(outlier_upper)) {
  data1<-data1[data1$feature > outlier_lower & data1$feature < outlier_upper,]
} else {
  cat("[DEBUG] outlier bounds are not finite; skipping 3SD filter\n")
}
cat(sprintf("[DEBUG] after filtering: n=%d\n", nrow(data1)))
debug_num_summary(data1$feature, "feature_after_filter")
debug_factor_table(data1$Split, "split counts after filtering")
debug_factor_table(data1$Site_ZZZ, "site counts after filtering")

#select the columns
data1<-data1[,c('SubjectKey','Split','Age','Sex','Site_ZZZ','tem_feature','feature')]

data2<-data1;

feature_step <- "build mu/sigma offsets"
#for mu
tem_rnd<-matrix(0,dim(data2)[1],1);
Model.Frame<-model.frame(formula = m2$mu.formula,data=data2);
Model.Matrix<-model.matrix(m2$mu.formula,Model.Frame)
Fit.fix<-matrix(m2$mu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))

if(!is.null(m2$mu.coefSmo[[1]]))
{Fit.fix[length(Fit.fix)]=0;
for(iz in 1:dim(data2)[1]){
  tem_rnd[iz]<-0
}
}else{print('no random effects for this term')};

mufix<-as.vector(Model.Matrix %*% Fit.fix)
debug_num_summary(mufix, "mufix")


#for sigma
tem_rnd<-matrix(0,dim(data2)[1],1);
Model.Frame<-model.frame(formula = m2$sigma.formula,data=data2);
Model.Matrix<-model.matrix(m2$sigma.formula,Model.Frame)
Fit.fix<-matrix(m2$sigma.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
if(!is.null(m2$sigma.coefSmo[[1]]))
{Fit.fix[length(Fit.fix)]=0;
for(iz in 1:dim(data2)[1]){
  tem_rnd[iz]<-0
}
}else{print('no random effects for this term')};

sigmafix<-as.vector(Model.Matrix %*% Fit.fix)
debug_num_summary(sigmafix, "sigmafix")

data2[,'mufix']<-mufix
data2[,'sigmafix']<-sigmafix

data2 <- na.omit(data2)
cat(sprintf("[DEBUG] rows after na.omit on design data: %d\n", nrow(data2)))
debug_num_summary(data2$feature, "data2$feature")
debug_num_summary(data2$mufix, "data2$mufix")
debug_num_summary(data2$sigmafix, "data2$sigmafix")

con=gamlss.control()
data2$Site_ZZZ <- as.character(data2$Site_ZZZ)
data2$Site_ZZZ[is.na(data2$Site_ZZZ) | data2$Site_ZZZ == ""] <- "UnknownSite"
data2$Site_ZZZ <- as.factor(data2$Site_ZZZ)
data2$Sex <- as.character(data2$Sex)
data2$Sex <- factor(data2$Sex, levels=c('Female','Male'))

data2$Split <- as.character(data2$Split)
data2$Split[is.na(data2$Split) | data2$Split == ""] <- "holdout"
data2_cal <- droplevels(data2[data2$Split == "calibration", ])
data2_holdout <- droplevels(data2[data2$Split == "holdout", ])

cat(sprintf("[DEBUG] calibration rows: %d | holdout rows: %d\n", nrow(data2_cal), nrow(data2_holdout)))
debug_factor_table(data2_cal$Site_ZZZ, "calibration site counts")
debug_factor_table(data2_holdout$Site_ZZZ, "holdout site counts")
debug_num_summary(data2_cal$feature, "calibration feature")
debug_num_summary(data2_holdout$feature, "holdout feature")

if (nrow(data2_cal) < 2) stop("Too few calibration subjects after filtering")

feature_step <- "fit calibration gamlss"
n_site_cal <- length(unique(as.character(data2_cal$Site_ZZZ)))
cat(sprintf("[DEBUG] unique calibration sites: %d\n", n_site_cal))
if (n_site_cal < 2) {
  cat("[DEBUG] calibration model mode: offset-only\n")
  m2_cal<-gamlss(formula=feature~offset(mufix),
                 sigma.formula = feature~offset(sigmafix),
                 control=con,
                 family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                 data=data2_cal)
} else {
  cat("[DEBUG] calibration model mode: random(site) on mu and sigma\n")
  m2_cal<-gamlss(formula=feature~offset(mufix)+random(Site_ZZZ),
                 sigma.formula = feature~offset(sigmafix)+random(Site_ZZZ),
                 control=con,
                 family = GG(mu.link='log',sigma.link = 'log',nu.link = 'identity'),
                 data=data2_cal)
}
cat("[DEBUG] calibration gamlss fit completed\n")
debug_kv("m2_cal$converged", m2_cal$converged)
debug_kv("m2_cal$family", m2_cal$family[1])
debug_num_summary(as.numeric(fitted(m2_cal, what="mu")), "fitted mu (calibration)")
debug_num_summary(as.numeric(fitted(m2_cal, what="sigma")), "fitted sigma (calibration)")


model_cal<-m2;
mu_cal_re <- get_re_coef(m2_cal$mu.coefSmo)
sigma_cal_re <- get_re_coef(m2_cal$sigma.coefSmo)
if (!is.null(mu_cal_re) && has_re_coef(m2$mu.coefSmo))
{
  for(new_site in names(mu_cal_re))
  {
    model_cal$mu.coefSmo[[1]][1]$coef[new_site]<-mu_cal_re[new_site]
  }
}
if (!is.null(sigma_cal_re) && has_re_coef(m2$sigma.coefSmo))
{
  for(new_site in names(sigma_cal_re))
  {
    model_cal$sigma.coefSmo[[1]][1]$coef[new_site]<-sigma_cal_re[new_site]
  }
}

setwd(paste0(savepath,'/',str))
feature_step <- "save calibrated model RDS"
saveRDS(model_cal, paste0(str,'_',i,'_loop_our_model_new_site_calibrated.rds'))
cat("[DEBUG] calibrated model RDS saved\n")

if (nrow(data2_cal) >= 2) {
  feature_step <- "export calibration z-scores"
  par_uncal_cal <- calc_mu_sigma_nu(m2, data2_cal)
  z_uncal_cal <- zzz_cent(obj=m2,type=c("z-scores"),mu=par_uncal_cal$mu,sigma=par_uncal_cal$sigma,nu=par_uncal_cal$nu,
                          xname = "Age",xvalues=data2_cal$Age,yval=data2_cal$tem_feature,
                          calibration=FALSE,lpar=3)
  par_cal_cal <- calc_mu_sigma_nu(model_cal, data2_cal)
  z_cal_cal <- zzz_cent(obj=model_cal,type=c("z-scores"),mu=par_cal_cal$mu,sigma=par_cal_cal$sigma,nu=par_cal_cal$nu,
                        xname = "Age",xvalues=data2_cal$Age,yval=data2_cal$tem_feature,
                        calibration=FALSE,lpar=3)

  summary_cal <- rbind(
    data.frame(Model="uncalibrated", N=length(z_uncal_cal), Z_mean=mean(z_uncal_cal, na.rm=TRUE), Z_sd=sd(z_uncal_cal, na.rm=TRUE), Z_var=var(z_uncal_cal, na.rm=TRUE)),
    data.frame(Model="calibrated", N=length(z_cal_cal), Z_mean=mean(z_cal_cal, na.rm=TRUE), Z_sd=sd(z_cal_cal, na.rm=TRUE), Z_var=var(z_cal_cal, na.rm=TRUE))
  )
  debug_num_summary(z_uncal_cal, "z_uncal_cal")
  debug_num_summary(z_cal_cal, "z_cal_cal")
  write.csv(summary_cal, file=paste0(str,"_",i,"_calibration_zscore_summary.csv"), row.names=FALSE)

  z_by_site_cal <- rbind(
    data.frame(Model="uncalibrated", Site_ZZZ=data2_cal$Site_ZZZ, Z=as.numeric(z_uncal_cal)),
    data.frame(Model="calibrated", Site_ZZZ=data2_cal$Site_ZZZ, Z=as.numeric(z_cal_cal))
  )
  z_site_sum_cal <- z_by_site_cal %>%
    group_by(Model, Site_ZZZ) %>%
    summarize(N=dplyr::n(), Z_mean=mean(Z, na.rm=TRUE), Z_sd=sd(Z, na.rm=TRUE), Z_var=var(Z, na.rm=TRUE), .groups="drop")
  write.csv(as.data.frame(z_site_sum_cal), file=paste0(str,"_",i,"_calibration_zscore_by_site.csv"), row.names=FALSE)
}

if (nrow(data2_holdout) >= 2) {
  feature_step <- "export holdout z-scores"
  par_uncal <- calc_mu_sigma_nu(m2, data2_holdout)
  z_uncal <- zzz_cent(obj=m2,type=c("z-scores"),mu=par_uncal$mu,sigma=par_uncal$sigma,nu=par_uncal$nu,
                      xname = "Age",xvalues=data2_holdout$Age,yval=data2_holdout$tem_feature,
                      calibration=FALSE,lpar=3)
  par_cal <- calc_mu_sigma_nu(model_cal, data2_holdout)
  z_cal <- zzz_cent(obj=model_cal,type=c("z-scores"),mu=par_cal$mu,sigma=par_cal$sigma,nu=par_cal$nu,
                    xname = "Age",xvalues=data2_holdout$Age,yval=data2_holdout$tem_feature,
                    calibration=FALSE,lpar=3)

  summary_all <- rbind(
    data.frame(Model="uncalibrated", N=length(z_uncal), Z_mean=mean(z_uncal, na.rm=TRUE), Z_sd=sd(z_uncal, na.rm=TRUE), Z_var=var(z_uncal, na.rm=TRUE)),
    data.frame(Model="calibrated", N=length(z_cal), Z_mean=mean(z_cal, na.rm=TRUE), Z_sd=sd(z_cal, na.rm=TRUE), Z_var=var(z_cal, na.rm=TRUE))
  )
  debug_num_summary(z_uncal, "z_uncal_holdout")
  debug_num_summary(z_cal, "z_cal_holdout")
  write.csv(summary_all, file=paste0(str,"_",i,"_holdout_zscore_summary.csv"), row.names=FALSE)

  z_by_site <- rbind(
    data.frame(Model="uncalibrated", Site_ZZZ=data2_holdout$Site_ZZZ, Z=as.numeric(z_uncal)),
    data.frame(Model="calibrated", Site_ZZZ=data2_holdout$Site_ZZZ, Z=as.numeric(z_cal))
  )
  z_site_sum <- z_by_site %>%
    group_by(Model, Site_ZZZ) %>%
    summarize(N=dplyr::n(), Z_mean=mean(Z, na.rm=TRUE), Z_sd=sd(Z, na.rm=TRUE), Z_var=var(Z, na.rm=TRUE), .groups="drop")
  write.csv(as.data.frame(z_site_sum), file=paste0(str,"_",i,"_holdout_zscore_by_site.csv"), row.names=FALSE)
}
#plot calibrate models
feature_step <- "build original reference centiles"
#both male and female
#original model
model1<-m2;
num_length=5000
if(has_re_coef(model1$mu.coefSmo))
{
data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
} else
{
data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}
data4 <- do.call( what=expand.grid, args=data3 )

data4[,'feature']<-0
tem_rnd<-matrix(0,dim(data4)[1],1);
Model.Frame<-model.frame(formula = model1$mu.formula,data=data4);
Model.Matrix<-model.matrix(model1$mu.formula,Model.Frame)
Fit.fix<-matrix(model1$mu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
mu_re <- get_re_coef(model1$mu.coefSmo)
if(!is.null(mu_re))
{Fit.fix[length(Fit.fix)]=0;
for(iz in 1:dim(data4)[1]){
  if(data4$Site_ZZZ[iz] %in% names(mu_re)){
    tem_rnd[iz]<-mu_re[as.character(data4$Site_ZZZ[iz])]
  }else{tem_rnd[iz]<-mean(mu_re)}
}
}else{print('no random effects for this term')};


mu0<-exp(as.vector(Model.Matrix %*% Fit.fix)+as.vector(tem_rnd))



if(has_re_coef(model1$sigma.coefSmo))
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
} else
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}
data4 <- do.call( what=expand.grid, args=data3 )

data4[,'feature']<-0
tem_rnd<-matrix(0,dim(data4)[1],1);
Model.Frame<-model.frame(formula = model1$sigma.formula,data=data4);
Model.Matrix<-model.matrix(model1$sigma.formula,Model.Frame)
Fit.fix<-matrix(model1$sigma.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
sigma_re <- get_re_coef(model1$sigma.coefSmo)
if(!is.null(sigma_re))
{Fit.fix[length(Fit.fix)]=0;
for(iz in 1:dim(data4)[1]){
  if(data4$Site_ZZZ[iz] %in% names(sigma_re)){
    tem_rnd[iz]<-sigma_re[as.character(data4$Site_ZZZ[iz])]
  }else{tem_rnd[iz]<-mean(sigma_re)}
}
}else{print('no random effects for this term')};

sigma0<-exp(as.vector(Model.Matrix %*% Fit.fix)+as.vector(tem_rnd))



if(has_re_coef(model1$nu.coefSmo))
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
} else
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}

data4 <- do.call( what=expand.grid, args=data3 )

Model.Frame<-model.frame(formula = model1$nu.formula,data=data4);
Model.Matrix<-model.matrix(model1$nu.formula,Model.Frame)
Fit.fix<-matrix(model1$nu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
nu_re <- get_re_coef(model1$nu.coefSmo)
if(!is.null(nu_re))
{Fit.fix[length(Fit.fix)]=0;
}else{print('no random effects for this term')};

nu0<-as.vector(Model.Matrix[,colnames(Model.Matrix)[1]]*Fit.fix[1])



tem_par<-mu0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;

if(Seg>1)
{
for(Seg1 in c(2:Seg))
{
  par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
}
par=par/Seg
}
mu=par


tem_par<-sigma0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;
if(Seg>1)
{
for(Seg1 in c(2:Seg))
{
  par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
}
par=par/Seg
}
sigma=par


tem_par<-nu0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;
if(Seg>1)
{
for(Seg1 in c(2:Seg))
{
  par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
}
par=par/Seg
}
nu=par


p2<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
             cent = c(0.5, 2.5, 50, 97.5,99.5),xname = 'Age',xvalues=data4$Age[1:num_length],
             calibration=FALSE,lpar=3)
p2[,'sigma']<-sigma


#calibrate model
feature_step <- "build calibrated centiles"
model1<-model_cal;
num_length=5000
if(has_re_coef(model1$mu.coefSmo))
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$mu.coefSmo[[1]]$coef))
} else
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}
data4 <- do.call( what=expand.grid, args=data3 )

data4[,'feature']<-0
tem_rnd<-matrix(0,dim(data4)[1],1);
Model.Frame<-model.frame(formula = model1$mu.formula,data=data4);
Model.Matrix<-model.matrix(model1$mu.formula,Model.Frame)
Fit.fix<-matrix(model1$mu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
mu_re <- get_re_coef(model1$mu.coefSmo)
if(!is.null(mu_re))
{Fit.fix[length(Fit.fix)]=0;
for(iz in 1:dim(data4)[1]){
  if(data4$Site_ZZZ[iz] %in% names(mu_re)){
    tem_rnd[iz]<-mu_re[as.character(data4$Site_ZZZ[iz])]
  }else{tem_rnd[iz]<-mean(mu_re)}
}
}else{print('no random effects for this term')};


mu0<-exp(as.vector(Model.Matrix %*% Fit.fix)+as.vector(tem_rnd))

#mu0 <- predict(model1, newdata = data4, type = "response", what = "mu")


if(has_re_coef(model1$sigma.coefSmo))
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$sigma.coefSmo[[1]]$coef))
} else
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}
data4 <- do.call( what=expand.grid, args=data3 )

data4[,'feature']<-0
tem_rnd<-matrix(0,dim(data4)[1],1);
Model.Frame<-model.frame(formula = model1$sigma.formula,data=data4);
Model.Matrix<-model.matrix(model1$sigma.formula,Model.Frame)
Fit.fix<-matrix(model1$sigma.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
sigma_re <- get_re_coef(model1$sigma.coefSmo)
if(!is.null(sigma_re))
{Fit.fix[length(Fit.fix)]=0;
for(iz in 1:dim(data4)[1]){
  if(data4$Site_ZZZ[iz] %in% names(sigma_re)){
    tem_rnd[iz]<-sigma_re[as.character(data4$Site_ZZZ[iz])]
  }else{tem_rnd[iz]<-mean(sigma_re)}
}
}else{print('no random effects for this term')};

sigma0<-exp(as.vector(Model.Matrix %*% Fit.fix)+as.vector(tem_rnd))
#sigma0 <- predict(model1, newdata = data4, type = "response", what = "sigma")


if(has_re_coef(model1$nu.coefSmo))
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=names(model1$nu.coefSmo[[1]]$coef))
} else
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}

data4 <- do.call( what=expand.grid, args=data3 )

Model.Frame<-model.frame(formula = model1$nu.formula,data=data4);
Model.Matrix<-model.matrix(model1$nu.formula,Model.Frame)
Fit.fix<-matrix(model1$nu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
nu_re <- get_re_coef(model1$nu.coefSmo)
if(!is.null(nu_re))
{Fit.fix[length(Fit.fix)]=0;
}else{print('no random effects for this term')};

nu0<-as.vector(Model.Matrix[,colnames(Model.Matrix)[1]]*Fit.fix[1])

#nu0 <- predict(model1, newdata = data4, type = "response", what = "nu")

tem_par<-mu0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;

if(Seg>1)
{
  for(Seg1 in c(2:Seg))
  {
    par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
  }
  par=par/Seg
}
mu=par


tem_par<-sigma0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;
if(Seg>1)
{
  for(Seg1 in c(2:Seg))
  {
    par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
  }
  par=par/Seg
}
sigma=par


tem_par<-nu0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;
if(Seg>1)
{
  for(Seg1 in c(2:Seg))
  {
    par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
  }
  par=par/Seg
}
nu=par


p2_cal<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
             cent = c(0.5, 2.5, 50, 97.5,99.5),xname = 'Age',xvalues=data4$Age[1:num_length],
             calibration=FALSE,lpar=3)
p2_cal[,'sigma']<-sigma





#site_specific model
feature_step <- "build site-specific centiles"
model1<-model_cal;
num_length=5000
if(has_re_coef(model1$mu.coefSmo))
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=unique(data1$Site_ZZZ))
} else
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}
data4 <- do.call( what=expand.grid, args=data3 )

data4[,'feature']<-0
tem_rnd<-matrix(0,dim(data4)[1],1);
Model.Frame<-model.frame(formula = model1$mu.formula,data=data4);
Model.Matrix<-model.matrix(model1$mu.formula,Model.Frame)
Fit.fix<-matrix(model1$mu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
mu_re <- get_re_coef(model1$mu.coefSmo)
if(!is.null(mu_re))
{Fit.fix[length(Fit.fix)]=0;
for(iz in 1:dim(data4)[1]){
  if(data4$Site_ZZZ[iz] %in% names(mu_re)){
    tem_rnd[iz]<-mu_re[as.character(data4$Site_ZZZ[iz])]
  }else{tem_rnd[iz]<-mean(mu_re)}
}
}else{print('no random effects for this term')};


mu0<-exp(as.vector(Model.Matrix %*% Fit.fix)+as.vector(tem_rnd))




if(has_re_coef(model1$sigma.coefSmo))
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=unique(data1$Site_ZZZ))
} else
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}
data4 <- do.call( what=expand.grid, args=data3 )

data4[,'feature']<-0
tem_rnd<-matrix(0,dim(data4)[1],1);
Model.Frame<-model.frame(formula = model1$sigma.formula,data=data4);
Model.Matrix<-model.matrix(model1$sigma.formula,Model.Frame)
Fit.fix<-matrix(model1$sigma.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
sigma_re <- get_re_coef(model1$sigma.coefSmo)
if(!is.null(sigma_re))
{Fit.fix[length(Fit.fix)]=0;
for(iz in 1:dim(data4)[1]){
  if(data4$Site_ZZZ[iz] %in% names(sigma_re)){
    tem_rnd[iz]<-sigma_re[as.character(data4$Site_ZZZ[iz])]
  }else{tem_rnd[iz]<-mean(sigma_re)}
}
}else{print('no random effects for this term')};

sigma0<-exp(as.vector(Model.Matrix %*% Fit.fix)+as.vector(tem_rnd))



if(has_re_coef(model1$nu.coefSmo))
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male"),Site_ZZZ=unique(data1$Site_ZZZ))
} else
{
  data3 <- list(Age=seq(min(data1$Age),max(data1$Age),length.out=num_length),Sex=c("Female","Male")) 
}

data4 <- do.call( what=expand.grid, args=data3 )

Model.Frame<-model.frame(formula = model1$nu.formula,data=data4);
Model.Matrix<-model.matrix(model1$nu.formula,Model.Frame)
Fit.fix<-matrix(model1$nu.coefficients[colnames(Model.Matrix)],ncol=1,dimnames = list(colnames(Model.Matrix),'Beta'))
nu_re <- get_re_coef(model1$nu.coefSmo)
if(!is.null(nu_re))
{Fit.fix[length(Fit.fix)]=0;
}else{print('no random effects for this term')};

nu0<-as.vector(Model.Matrix[,colnames(Model.Matrix)[1]]*Fit.fix[1])



tem_par<-mu0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;

if(Seg>1)
{
  for(Seg1 in c(2:Seg))
  {
    par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
  }
  par=par/Seg
}
mu=par


tem_par<-sigma0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;
if(Seg>1)
{
  for(Seg1 in c(2:Seg))
  {
    par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
  }
  par=par/Seg
}
sigma=par


tem_par<-nu0
par<-tem_par[1:num_length]
Seg=length(tem_par)/num_length;
if(Seg>1)
{
  for(Seg1 in c(2:Seg))
  {
    par=par+tem_par[((Seg1-1)*num_length+1):(Seg1*num_length)]
  }
  par=par/Seg
}
nu=par


p2_site<-zzz_cent(obj=model1,type=c("centiles"),mu=mu,sigma=sigma,nu=nu,
                 cent = c(0.5, 2.5, 50, 97.5,99.5),xname = 'Age',xvalues=data4$Age[1:num_length],
                 calibration=FALSE,lpar=3)
p2_site[,'sigma']<-sigma


library(reshape2);
colnames(p2)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma');

colnames(p2_cal)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma');

colnames(p2_site)<-c('Age','lower99CI','lower95CI','median','upper95CI','upper99CI','sigma');


if(!(sheet %in% c('lh.aparc.thickness.table','rh.aparc.thickness.table')))
   {
     scale1=10000;
     ylab1='×10^4 mm3';
   }

if(sheet %in% c('lh.aparc.thickness.table','rh.aparc.thickness.table'))
   {
     scale1=1;
     ylab1='mm';
   }



# png(filename = paste0(str,'_',i,'_all_without_sex_stratified.png'), 
#     width = 1480,           
#     height = 740,          
#     units = "px",          
#     bg = "white",          
#     res = 300)     
# 
# p3<-ggplot()+
#   #geom_line(data=mydata,aes(x=Age,y=Value,group=Percentile,color=Percentile))+
#   geom_point(data=data1[data1$Sex=='Female',],aes(x=Age,y=tem_feature/scale1),
#              colour=c('grey'),shape=16,size=3,alpha = 0.1)+
#   geom_point(data=data1[data1$Sex=='Male',],aes(x=Age,y=tem_feature/scale1),
#              colour=c('grey'),shape=17,size=3,alpha = 0.1)+
#   geom_line(data=p2,aes(x=Age,y=median/scale1),color=c('#262626'),linewidth=2,linetype=c('solid'))+
#   geom_line(data=p2_cal,aes(x=Age,y=median/scale1),color=c('#4FBBD8'),linewidth=2,linetype=c('dashed'))+
#   geom_line(data=p2_site,aes(x=Age,y=median/scale1),color=c('#E84935'),linewidth=2,linetype=c('dotted'))+
#   # geom_line(data=p2,aes(x=Age,y=lower99CI/scale1),color=c('#262626'),linewidth=1,linetype=c('dashed'))+
#   # geom_line(data=p2,aes(x=Age,y=lower95CI/scale1),color=c('#262626'),linewidth=1,linetype=c('dotted'))+
#   # geom_line(data=p2,aes(x=Age,y=upper95CI/scale1),color=c('#262626'),linewidth=1,linetype=c('dotted'))+
#   # geom_line(data=p2,aes(x=Age,y=upper99CI/scale1),color=c('#262626'),linewidth=1,linetype=c('dashed'))+
#   labs(title=paste0(i,' ',ylab1),x='',y='')+
#   theme_bw()+
#   theme(
#     axis.title = element_text(family = "serif",size=12,color = "black"),
#     axis.text.x = element_text(
#       size = 12,              
#       color = "black",          
#       family = "serif"     
#     ),
#     axis.text.y = element_text(
#       size = 10,              
#       color = "black",         
#       #face = "bold" ,          
#       family = "serif"
#     )
#   )+
#   # scale_x_log10(breaks = c(6,18,35,80),  
#   #               labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
#   scale_x_continuous(breaks = c(6,18,35,80),  
#                      labels = c("6 yr", "18 yr", "35 yr", "80 yr")) 
# 
# print(p3)  
# dev.off()

cat(sprintf("[FEATURE DONE] sheet=%s | feature=%s\n", sheet, i))

}, error = function(e) {
  debug_sep(sprintf("[FEATURE ERROR] sheet=%s | feature=%s", sheet, i))
  err_msg <- conditionMessage(e)
  cal_n <- if (exists("data2_cal")) nrow(data2_cal) else NA_integer_
  holdout_n <- if (exists("data2_holdout")) nrow(data2_holdout) else NA_integer_
  total_n <- if (exists("data1")) nrow(data1) else NA_integer_
  unique_sites_cal <- if (exists("data2_cal") && "Site_ZZZ" %in% colnames(data2_cal)) length(unique(as.character(data2_cal$Site_ZZZ))) else NA_integer_
  feature_mean_fail <- if (exists("data1") && "feature" %in% colnames(data1)) mean(data1$feature, na.rm=TRUE) else NA_real_
  feature_sd_fail <- if (exists("data1") && "feature" %in% colnames(data1)) sd(data1$feature, na.rm=TRUE) else NA_real_

  cat(sprintf("[ERROR] failed step: %s\n", feature_step))
  cat(sprintf("[ERROR] message: %s\n", err_msg))
  cat(sprintf("[ERROR] calibration rows at failure: %s\n", cal_n))
  cat(sprintf("[ERROR] holdout rows at failure: %s\n", holdout_n))

  if (exists("data1") && "Site_ZZZ" %in% colnames(data1)) {
    debug_factor_table(data1$Site_ZZZ, "site counts at failure")
  }
  if (exists("data1") && "feature" %in% colnames(data1)) {
    debug_num_summary(data1$feature, "feature values at failure")
  }

  append_error_log(
    error_log_csv,
    data.frame(
      sheet=sheet,
      feature=i,
      failed_step=feature_step,
      error_message=err_msg,
      n_total=total_n,
      n_calibration=cal_n,
      n_holdout=holdout_n,
      n_site_calibration=unique_sites_cal,
      feature_mean=feature_mean_fail,
      feature_sd=feature_sd_fail,
      timestamp=as.character(Sys.time()),
      stringsAsFactors=FALSE
    )
  )
  cat(sprintf("[ERROR] appended to error csv: %s\n", error_log_csv))
  cat(sprintf("[FEATURE SKIP] sheet=%s | feature=%s\n", sheet, i))
  return(NULL)
})


}

}



