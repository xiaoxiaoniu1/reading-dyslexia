if(!requireNamespace("pacman",quietly=TRUE)) install.packages("pacman",repos="https://cloud.r-project.org/")
pacman::p_load(readxl,dplyr,mgcv,ggplot2)

`%||%` <- function(x,y) if(is.null(x)||length(x)==0||all(is.na(x))) y else x

args <- commandArgs(trailingOnly=TRUE)
get_arg <- function(key,default=NULL){
  i <- match(key,args)
  if(is.na(i) || i==length(args)) default else args[i+1]
}

result_dir <- get_arg("--result-dir","/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_DGLM_stratified")
demo_file  <- get_arg("--demo-file","/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/all_data_cqt_mean_1.5.xlsx")
mind_dir   <- get_arg("--mind-combat-dir","/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_combat")
out_dir    <- get_arg("--out-dir",file.path(result_dir,"VIOLIN_Followup_stratified"))
label <- get_arg("--analysis-label",basename(normalizePath(result_dir,winslash="/",mustWork=FALSE)))
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)

safe <- function(x) gsub("[^A-Za-z0-9_]+","_",x)
cat("result_dir:",result_dir,"\nout_dir:",out_dir,"\n")

df <- read_excel(demo_file,sheet="Sheet1") |>
  mutate(
    original_project = as.character(`original-project`),
    id_old = as.character(id_old),
    Age_Raw = as.numeric(age_month),
    file_base = paste0(original_project,"_",id_old,"_MIND_DK318_combat"),
    degree_file = file.path(mind_dir,paste0(file_base,"_degree.csv")),
    Diagnosis_Fac = factor(ifelse(group_d_or_c==0,"TD","DD"),levels=c("TD","DD")),
    AgeGroup_Fac = factor(ifelse(group_age==1,"Adult","Child"),levels=c("Child","Adult")),
    Sex = factor(ifelse(sex==1,"Male","Female"),levels=c("Male","Female")),
    Site = factor(site),
    has_file = file.exists(degree_file)
  ) |>
  filter(!is.na(original_project),!is.na(id_old),original_project!="",id_old!="",!is.na(Age_Raw),has_file)
if(!nrow(df)) stop("No valid degree files found.")

df$Age_Centered <- df$Age_Raw - mean(df$Age_Raw,na.rm=TRUE)

read_deg <- function(fp){
  d <- read.csv(fp,check.names=FALSE)
  v <- as.numeric(d$degree)
  names(v) <- as.character(d$ROI)
  v
}

tmp <- read_deg(df$degree_file[1])
rois <- names(tmp)
deg <- matrix(NA_real_,nrow=nrow(df),ncol=length(rois),dimnames=list(df$file_base,rois))
for(i in seq_len(nrow(df))){
  v <- read_deg(df$degree_file[i])
  if(!identical(names(v),rois)) stop("ROI mismatch: ",df$degree_file[i])
  deg[i,] <- v
}

fit_adjust_model <- function(d,response){
  rhs <- c("s(Age_Centered,k=5)")
  if(nlevels(droplevels(d$Diagnosis_Fac))>1) rhs <- c(rhs,"Diagnosis_Fac")
  if(nlevels(droplevels(d$AgeGroup_Fac))>1) rhs <- c(rhs,"AgeGroup_Fac")
  if(nlevels(droplevels(d$Sex))>1) rhs <- c(rhs,"Sex")
  if(nlevels(droplevels(d$Site))>1) rhs <- c(rhs,"Site")
  tryCatch(
    gam(as.formula(paste(response,"~",paste(rhs,collapse="+"))),data=d,method="REML"),
    error=function(e) NULL
  )
}

disp_y <- function(d){
  f <- fit_adjust_model(d,"y")
  if(is.null(f)) abs(d$y - mean(d$y,na.rm=TRUE)) else abs(residuals(f,type="response"))
}

mk <- function(csv,pcol,kind,subset_var,subset_level,x_var,effect){
  data.frame(
    csv=csv,
    pcol=pcol,
    kind=kind,
    subset_var=subset_var,
    subset_level=subset_level,
    x_var=x_var,
    facet_var="Sex",
    family="interaction_followup",
    threshold=if(grepl("_FDR$",pcol)) "FDR" else "uncorrected",
    effect=effect,
    tag=file.path(subset_level,"interaction_followup",if(grepl("_FDR$",pcol)) "FDR" else "uncorrected",effect),
    stringsAsFactors=FALSE
  )
}

tasks <- bind_rows(
  mk("DGLM_DK318_Child_degree_results.csv","p_mean_Diagnosis_Sex_Interaction","mean","AgeGroup_Fac","Child","Diagnosis_Fac","mean_diagnosis_by_sex"),
  mk("DGLM_DK318_Child_degree_results.csv","p_mean_Diagnosis_Sex_Interaction_FDR","mean","AgeGroup_Fac","Child","Diagnosis_Fac","mean_diagnosis_by_sex"),
  mk("DGLM_DK318_Child_degree_results.csv","p_disp_Diagnosis_Sex_Interaction","disp","AgeGroup_Fac","Child","Diagnosis_Fac","disp_diagnosis_by_sex"),
  mk("DGLM_DK318_Child_degree_results.csv","p_disp_Diagnosis_Sex_Interaction_FDR","disp","AgeGroup_Fac","Child","Diagnosis_Fac","disp_diagnosis_by_sex"),
  mk("DGLM_DK318_Adult_degree_results.csv","p_mean_Diagnosis_Sex_Interaction","mean","AgeGroup_Fac","Adult","Diagnosis_Fac","mean_diagnosis_by_sex"),
  mk("DGLM_DK318_Adult_degree_results.csv","p_mean_Diagnosis_Sex_Interaction_FDR","mean","AgeGroup_Fac","Adult","Diagnosis_Fac","mean_diagnosis_by_sex"),
  mk("DGLM_DK318_Adult_degree_results.csv","p_disp_Diagnosis_Sex_Interaction","disp","AgeGroup_Fac","Adult","Diagnosis_Fac","disp_diagnosis_by_sex"),
  mk("DGLM_DK318_Adult_degree_results.csv","p_disp_Diagnosis_Sex_Interaction_FDR","disp","AgeGroup_Fac","Adult","Diagnosis_Fac","disp_diagnosis_by_sex"),
  mk("DGLM_DK318_DD_Adult_vs_Child_degree_results.csv","p_mean_AgeGroup_Sex_Interaction","mean","Diagnosis_Fac","DD","AgeGroup_Fac","mean_agegroup_by_sex"),
  mk("DGLM_DK318_DD_Adult_vs_Child_degree_results.csv","p_mean_AgeGroup_Sex_Interaction_FDR","mean","Diagnosis_Fac","DD","AgeGroup_Fac","mean_agegroup_by_sex"),
  mk("DGLM_DK318_DD_Adult_vs_Child_degree_results.csv","p_disp_AgeGroup_Sex_Interaction","disp","Diagnosis_Fac","DD","AgeGroup_Fac","disp_agegroup_by_sex"),
  mk("DGLM_DK318_DD_Adult_vs_Child_degree_results.csv","p_disp_AgeGroup_Sex_Interaction_FDR","disp","Diagnosis_Fac","DD","AgeGroup_Fac","disp_agegroup_by_sex"),
  mk("DGLM_DK318_TD_Adult_vs_Child_degree_results.csv","p_mean_AgeGroup_Sex_Interaction","mean","Diagnosis_Fac","TD","AgeGroup_Fac","mean_agegroup_by_sex"),
  mk("DGLM_DK318_TD_Adult_vs_Child_degree_results.csv","p_mean_AgeGroup_Sex_Interaction_FDR","mean","Diagnosis_Fac","TD","AgeGroup_Fac","mean_agegroup_by_sex"),
  mk("DGLM_DK318_TD_Adult_vs_Child_degree_results.csv","p_disp_AgeGroup_Sex_Interaction","disp","Diagnosis_Fac","TD","AgeGroup_Fac","disp_agegroup_by_sex"),
  mk("DGLM_DK318_TD_Adult_vs_Child_degree_results.csv","p_disp_AgeGroup_Sex_Interaction_FDR","disp","Diagnosis_Fac","TD","AgeGroup_Fac","disp_agegroup_by_sex")
)

write.csv(tasks,file.path(out_dir,"violin_followup_task_catalog.csv"),row.names=FALSE)

validate_plot_data <- function(d,x_var,facet_var,min_n=3){
  d[[x_var]] <- droplevels(factor(d[[x_var]]))
  if(nlevels(d[[x_var]]) < 2) return(FALSE)
  tab <- table(as.character(d[[facet_var]]), d[[x_var]])
  if(!all(dim(tab) >= 2)) return(FALSE)
  all(tab >= min_n)
}

build_violin_plot <- function(d,task,feature){
  x_lab <- if(task$x_var=="Diagnosis_Fac") "Diagnosis" else "Age Group"
  y_lab <- if(task$kind=="disp") "|Residual| of Degree" else "Degree"
  fill_vals <- if(task$x_var=="Diagnosis_Fac") c(TD="#1f77b4",DD="#d62728") else c(Child="#2ca02c",Adult="#9467bd")
  fill_vals <- fill_vals[names(fill_vals) %in% levels(droplevels(d[[task$x_var]]))]

  ggplot(d,aes(x=.data[[task$x_var]],y=PlotY,fill=.data[[task$x_var]])) +
    geom_violin(trim=FALSE,alpha=.45,color=NA,width=.95) +
    geom_boxplot(width=.18,outlier.shape=NA,alpha=.85,color="#333333",fill="white") +
    geom_jitter(width=.10,height=0,size=1.5,alpha=.55,color="#333333") +
    scale_fill_manual(values=fill_vals,drop=FALSE) +
    facet_wrap(stats::as.formula(paste("~",task$facet_var)),nrow=1,scales="free_y") +
    labs(
      title=paste(feature,"|",task$subset_level,"|",task$effect),
      subtitle=paste(task$pcol,"|",label),
      x=x_lab,
      y=y_lab,
      fill=x_lab
    ) +
    theme_bw(base_size=14) +
    theme(
      legend.position="top",
      panel.grid.minor=element_blank(),
      strip.background=element_rect(fill="#f2f2f2",color="#d9d9d9")
    )
}

run_one <- function(d_all,feature,task){
  d <- d_all
  d$y <- deg[,feature]
  d <- d |>
    filter(is.finite(y), .data[[task$subset_var]] == task$subset_level)
  if(nrow(d) < 24 || var(d$y,na.rm=TRUE) < 1e-12) return(list(ok=FALSE,reason="data"))

  d$PlotY <- if(task$kind=="disp") disp_y(d) else d$y
  d <- d[is.finite(d$PlotY),,drop=FALSE]
  if(!validate_plot_data(d,task$x_var,task$facet_var,min_n=3)) return(list(ok=FALSE,reason="group_count"))

  od <- file.path(out_dir,task$tag)
  dir.create(od,recursive=TRUE,showWarnings=FALSE)
  out_file <- file.path(od,paste0("VIOLIN_",safe(task$pcol),"_",safe(feature),".png"))
  g <- build_violin_plot(d,task,feature)
  ggsave(out_file,g,width=10,height=4.8,dpi=300)
  list(ok=TRUE,main=out_file,diff=NA_character_)
}

manifest <- list()
for(i in seq_len(nrow(tasks))){
  task <- as.list(tasks[i,])
  csvp <- file.path(result_dir,task$csv)
  if(!file.exists(csvp)) {
    cat("skip missing",csvp,"\n")
    next
  }
  tab <- read.csv(csvp,check.names=FALSE)
  if(!(task$pcol %in% names(tab)) || !("feature" %in% names(tab))) {
    cat("skip missing col",task$pcol,"in",csvp,"\n")
    next
  }
  sig <- tab |> filter(!is.na(.data[[task$pcol]]) & .data[[task$pcol]] < 0.05)
  cat("task",task$pcol,"| subset",task$subset_level,":",nrow(sig),"significant ROIs\n")
  if(!nrow(sig)) next

  rows <- vector("list",nrow(sig))
  for(j in seq_len(nrow(sig))){
    ft <- as.character(sig$feature[j])
    if(!(ft %in% colnames(deg))){
      rows[[j]] <- data.frame(family=task$family,threshold=task$threshold,effect=task$effect,task=task$pcol,subset=task$subset_level,kind=task$kind,feature=ft,p=sig[[task$pcol]][j],ok=FALSE,reason="missing_feature",main=NA,diff=NA)
      next
    }
    rr <- run_one(df,ft,task)
    rows[[j]] <- data.frame(family=task$family,threshold=task$threshold,effect=task$effect,task=task$pcol,subset=task$subset_level,kind=task$kind,feature=ft,p=sig[[task$pcol]][j],ok=rr$ok,reason=rr$reason %||% NA,main=rr$main %||% NA,diff=rr$diff %||% NA)
  }

  task_out <- bind_rows(rows)
  dir.create(file.path(out_dir,task$tag),recursive=TRUE,showWarnings=FALSE)
  write.csv(task_out,file.path(out_dir,task$tag,paste0("summary_",safe(task$pcol),"_",safe(task$subset_level),".csv")),row.names=FALSE)
  manifest[[length(manifest)+1]] <- task_out
}

manifest_df <- bind_rows(manifest)
write.csv(manifest_df,file.path(out_dir,"violin_followup_manifest.csv"),row.names=FALSE)
if(nrow(manifest_df)) {
  write.csv(manifest_df |> group_by(family,threshold,subset,kind,effect) |> summarise(n_sig=n(),n_ok=sum(ok,na.rm=TRUE),.groups="drop"),file.path(out_dir,"violin_followup_overview_by_effect.csv"),row.names=FALSE)
  write.csv(manifest_df |> group_by(task,subset) |> summarise(n_sig=n(),n_ok=sum(ok,na.rm=TRUE),.groups="drop"),file.path(out_dir,"violin_followup_task_overview.csv"),row.names=FALSE)
}
cat("done. success:",sum(manifest_df$ok,na.rm=TRUE),"\n")
