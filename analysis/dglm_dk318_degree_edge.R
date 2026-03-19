# ============================================================
# DGLM (Dispersion Generalized Linear Model) Analysis
# for DK-318 MIND degree and edges
# 
# Interaction check using DGLM with:
#  - Mean model: y ~ Diagnosis * AgeGroup + Sex
#  - Variance model: ~ Diagnosis * AgeGroup + Sex
# ============================================================

# 自动安装缺失的包
packages <- c("readxl", "dplyr", "purrr", "dglm", "stringr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    cat("Installing", pkg, "...\n")
    install.packages(pkg, repos = "https://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

library(readxl)
library(dplyr)
library(purrr)
library(dglm)
library(stringr)

# ----------------------------
# 1) Paths
# ----------------------------
demo_file    <- "/data1/tqi/share/after_freesurfer/FILE/all_data_cqt.xlsx"
mind_csv_dir <- "/data1/tqi/share/after_freesurfer/fs_subjects_all/MIND_out_combat_degree/"

out_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/MIND_DGLM"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

out_degree_csv <- file.path(out_dir, "DGLM_DK318_degree_results.csv")
out_edge_csv   <- file.path(out_dir, "DGLM_DK318_edge_results.csv")

# ----------------------------
# 2) Load demo + recode factors
# ----------------------------
df <- read_excel(demo_file, sheet = "Sheet1")

df <- df %>%
  mutate(
    original_project = as.character(`original-project`),
    id_old = as.character(id_old),
    subj_prefix = paste0(original_project, "_", id_old, "_MIND_DK318"),
    file_base = paste0(subj_prefix, "_combat_labeled"),
    mind_file = file.path(mind_csv_dir, paste0(file_base, ".csv")),
    
    Diagnosis = ifelse(group_d_or_c == 0, "TD", "DD"),
    AgeGroup  = ifelse(group_age == 1, "Adult", "Child"),
    Sex       = ifelse(sex == 1, "Male", "Female"),
    
    Diagnosis = factor(Diagnosis, levels = c("TD", "DD")),
    AgeGroup  = factor(AgeGroup,  levels = c("Child", "Adult")),
    Sex       = factor(Sex),
    
    has_file  = file.exists(mind_file)
  ) %>%
  filter(!is.na(original_project), !is.na(id_old), original_project != "", id_old != "")

cat("Subjects in demo:", nrow(df), "\n")
cat("Subjects with MIND csv:", sum(df$has_file), "\n")

df2 <- df %>% filter(has_file)

# ----------------------------
# 3) Read matrix helpers
# ----------------------------
read_mind_csv <- function(fp) {
  mat <- as.matrix(read.csv(fp, row.names = 1, check.names = FALSE))
  storage.mode(mat) <- "numeric"
  mat
}

calc_degree <- function(mat) {
  diag(mat) <- NA
  rowMeans(mat, na.rm = TRUE)
}

tmp <- read_mind_csv(df2$mind_file[1])
n_roi <- nrow(tmp)
cat("ROI number:", n_roi, "\n")

upper_idx <- which(upper.tri(matrix(1, n_roi, n_roi)), arr.ind = TRUE)
n_edge <- nrow(upper_idx)
cat("Number of edges:", n_edge, "\n")

# ----------------------------
# 4) Build matrices
# ----------------------------
n_subj <- nrow(df2)
roi_names <- rownames(tmp)
if (is.null(roi_names) || any(is.na(roi_names))) {
  roi_names <- paste0("ROI_", seq_len(n_roi))
}

degree_mat <- matrix(NA_real_, nrow = n_subj, ncol = n_roi,
                     dimnames = list(df2$file_base, roi_names))
edge_mat <- matrix(NA_real_, nrow = n_subj, ncol = n_edge,
                   dimnames = list(df2$file_base, paste0("E_", seq_len(n_edge))))

for (s in seq_len(n_subj)) {
  m <- read_mind_csv(df2$mind_file[s])
  degree_mat[s, ] <- calc_degree(m)
  edge_mat[s, ] <- m[upper.tri(m)]
}

cat("Degree matrix:", paste(dim(degree_mat), collapse = " x "), "\n")
cat("Edge matrix:", paste(dim(edge_mat), collapse = " x "), "\n")

# ============================================================
# 5) DGLM function
# ============================================================
run_dglm_per_feature <- function(y_mat, df_sub, feature_names,
                                 min_n = 10, var_eps = 1e-12) {
  
  res <- lapply(seq_len(ncol(y_mat)), function(k) {
    y <- y_mat[, k]
    ok <- is.finite(y)
    y2 <- y[ok]
    d2 <- df_sub[ok, , drop = FALSE]
    d2$Diagnosis <- droplevels(d2$Diagnosis)
    d2$AgeGroup  <- droplevels(d2$AgeGroup)
    d2$Sex       <- droplevels(d2$Sex)

    na_row <- data.frame(
      feature = feature_names[k], n = length(y2),
      p_mean_Diagnosis = NA_real_, p_mean_AgeGroup = NA_real_, p_mean_Interaction = NA_real_,
      p_disp_Diagnosis = NA_real_, p_disp_AgeGroup = NA_real_, p_disp_Interaction = NA_real_
    )

    if (length(y2) < min_n) return(na_row)

    v <- var(y2)
    if (!is.finite(v) || v < var_eps) return(na_row)

    fit_err <- NULL
    fit <- tryCatch(
      dglm(y2 ~ Diagnosis * AgeGroup + Sex,
           dformula = ~ Diagnosis * AgeGroup + Sex,
           family = gaussian(link = "identity"),
           data = d2),
      error = function(e) { fit_err <<- conditionMessage(e); NULL }
    )

    if (is.null(fit)) {
      if (!is.null(fit_err) && k == 1) cat("[dglm error on feature 1]:", fit_err, "\n")
      return(na_row)
    }
    
    # Helper: extract p-values from a coefficient table
    # Auto-detects p-value column: "Pr(>|t|)" for mean model, "Pr(>|Chi|)" for dispersion model
    extract_ps <- function(ct) {
      p_diag <- NA_real_; p_age <- NA_real_; p_int <- NA_real_
      if (is.null(ct)) return(list(p_diag = p_diag, p_age = p_age, p_int = p_int))
      candidates <- c("Pr(>|t|)", "Pr(>|Chi|)", "Pr(>Chi)", "p.value")
      p_col <- candidates[candidates %in% colnames(ct)][1]
      if (is.na(p_col)) return(list(p_diag = p_diag, p_age = p_age, p_int = p_int))
      rn <- rownames(ct)
      diag_rows <- grep("^Diagnosis", rn)
      for (r in diag_rows) { if (!grepl(":", rn[r])) { p_diag <- ct[r, p_col]; break } }
      age_rows <- grep("^AgeGroup", rn)
      for (r in age_rows)  { if (!grepl(":", rn[r])) { p_age  <- ct[r, p_col]; break } }
      int_rows <- grep("Diagnosis.*:.*AgeGroup|AgeGroup.*:.*Diagnosis", rn)
      if (length(int_rows) > 0) p_int <- ct[int_rows[1], p_col]
      list(p_diag = p_diag, p_age = p_age, p_int = p_int)
    }

    # Mean model coefficients (glm -> "Pr(>|t|)")
    # summary(fit) gives mean submodel; fit$dispersion.fit is the dispersion glm object
    mean_ct <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
    mean_ps <- extract_ps(mean_ct)
    disp_ct <- tryCatch(summary(fit$dispersion.fit)$coefficients, error = function(e) NULL)
    disp_ps <- extract_ps(disp_ct)

    data.frame(
      feature               = feature_names[k],
      n                     = length(y2),
      p_mean_Diagnosis      = mean_ps$p_diag,
      p_mean_AgeGroup       = mean_ps$p_age,
      p_mean_Interaction    = mean_ps$p_int,
      p_disp_Diagnosis      = disp_ps$p_diag,
      p_disp_AgeGroup       = disp_ps$p_age,
      p_disp_Interaction    = disp_ps$p_int
    )
  })

  out <- bind_rows(res)
  out$p_mean_Diagnosis_FDR   <- p.adjust(out$p_mean_Diagnosis,   method = "fdr")
  out$p_mean_AgeGroup_FDR    <- p.adjust(out$p_mean_AgeGroup,    method = "fdr")
  out$p_mean_Interaction_FDR <- p.adjust(out$p_mean_Interaction, method = "fdr")
  out$p_disp_Diagnosis_FDR   <- p.adjust(out$p_disp_Diagnosis,   method = "fdr")
  out$p_disp_AgeGroup_FDR    <- p.adjust(out$p_disp_AgeGroup,    method = "fdr")
  out$p_disp_Interaction_FDR <- p.adjust(out$p_disp_Interaction, method = "fdr")

  out
}

# ============================================================
# 5b) Quick diagnostic: run dglm on first ROI
# ============================================================
cat("\n=== Diagnostic ===\n")
cat("Diagnosis counts:", paste(names(table(df2$Diagnosis)), table(df2$Diagnosis), sep="=", collapse=", "), "\n")
cat("AgeGroup  counts:", paste(names(table(df2$AgeGroup)),  table(df2$AgeGroup),  sep="=", collapse=", "), "\n")
cat("Sex       counts:", paste(names(table(df2$Sex)),        table(df2$Sex),        sep="=", collapse=", "), "\n")
diag_y  <- degree_mat[, 1]
diag_ok <- is.finite(diag_y)
diag_d  <- df2[diag_ok, ]
diag_d$Diagnosis <- droplevels(diag_d$Diagnosis)
diag_d$AgeGroup  <- droplevels(diag_d$AgeGroup)
diag_d$Sex       <- droplevels(diag_d$Sex)
cat("First ROI: n valid =", sum(diag_ok), ", var =", round(var(diag_y[diag_ok]), 8), "\n")
diag_fit <- tryCatch(
  dglm(diag_y[diag_ok] ~ Diagnosis * AgeGroup + Sex,
       dformula = ~ Diagnosis * AgeGroup + Sex,
       family = gaussian(link = "identity"),
       data = diag_d),
  error = function(e) { cat("dglm ERROR:", conditionMessage(e), "\n"); NULL }
)
if (!is.null(diag_fit)) {
  mean_ct2 <- summary(diag_fit)$coefficients
  disp_ct2 <- tryCatch(summary(diag_fit$dispersion.fit)$coefficients, error = function(e) NULL)
  cat("Mean model coef colnames: ", paste(colnames(mean_ct2), collapse=", "), "\n")
  cat("Disp model coef colnames: ", paste(colnames(disp_ct2), collapse=", "), "\n")
  p_col2 <- c("Pr(>|t|)","Pr(>|Chi|)","Pr(>Chi)","p.value")
  p_col2 <- p_col2[p_col2 %in% colnames(mean_ct2)][1]
  cat("DiagnosisDD mean p-value: ", mean_ct2["DiagnosisDD", p_col2], "\n")
  cat("dglm fit names: ", paste(names(diag_fit), collapse=", "), "\n")
} else {
  cat(">>> dglm FAILED on first ROI. All results will be NA. Fix the error above first.\n")
}

# ============================================================
# 6) Run DGLM for degree
# ============================================================
cat("\n=== Running DGLM (degree) ===\n")
degree_res <- run_dglm_per_feature(degree_mat, df2, colnames(degree_mat))

cat("Significant mean Interaction (degree, FDR<0.05):",
    sum(degree_res$p_mean_Interaction_FDR < 0.05, na.rm = TRUE), "\n")
cat("Significant dispersion Interaction (degree, FDR<0.05):",
    sum(degree_res$p_disp_Interaction_FDR < 0.05, na.rm = TRUE), "\n")

write.csv(degree_res, out_degree_csv, row.names = FALSE)
cat("Saved:", out_degree_csv, "\n")

# ============================================================
# 7) Run DGLM for edges (chunked)
# ============================================================
run_dglm_edgewise_chunked <- function(edge_mat, df_sub, upper_idx, roi_names,
                                      chunk_size = 2000, min_n = 10, var_eps = 1e-12) {
  n_edge <- ncol(edge_mat)
  chunks <- split(seq_len(n_edge), ceiling(seq_len(n_edge) / chunk_size))
  
  out_list <- vector("list", length(chunks))
  
  for (ci in seq_along(chunks)) {
    idxs <- chunks[[ci]]
    cat("Edge chunk", ci, "/", length(chunks), "\n")
    
    chunk_res <- lapply(idxs, function(k) {
      y <- edge_mat[, k]
      ok <- is.finite(y)
      y2 <- y[ok]
      d2 <- df_sub[ok, , drop = FALSE]
      # Drop unused factor levels after subsetting
      d2$Diagnosis <- droplevels(d2$Diagnosis)
      d2$AgeGroup  <- droplevels(d2$AgeGroup)
      d2$Sex       <- droplevels(d2$Sex)

      v <- ifelse(length(y2) > 1, var(y2), NA_real_)

      na_edge_row <- data.frame(
        edge_id = k, i = upper_idx[k, 1], j = upper_idx[k, 2],
        region_i = roi_names[upper_idx[k, 1]], region_j = roi_names[upper_idx[k, 2]],
        n = length(y2),
        p_mean_Diagnosis = NA_real_, p_mean_AgeGroup = NA_real_, p_mean_Interaction = NA_real_,
        p_disp_Diagnosis = NA_real_, p_disp_AgeGroup = NA_real_, p_disp_Interaction = NA_real_
      )

      if (length(y2) < min_n || !is.finite(v) || v < var_eps) return(na_edge_row)

      fit <- tryCatch(
        dglm(y2 ~ Diagnosis * AgeGroup + Sex,
             dformula = ~ Diagnosis * AgeGroup + Sex,
             family = gaussian(link = "identity"),
             data = d2),
        error = function(e) NULL
      )

      if (is.null(fit)) return(na_edge_row)
      
      # Helper: extract p-values; auto-detects column name
      extract_ps_edge <- function(ct) {
        p_diag <- NA_real_; p_age <- NA_real_; p_int <- NA_real_
        if (is.null(ct)) return(list(p_diag = p_diag, p_age = p_age, p_int = p_int))
        candidates <- c("Pr(>|t|)", "Pr(>|Chi|)", "Pr(>Chi)", "p.value")
        p_col <- candidates[candidates %in% colnames(ct)][1]
        if (is.na(p_col)) return(list(p_diag = p_diag, p_age = p_age, p_int = p_int))
        rn <- rownames(ct)
        diag_rows <- grep("^Diagnosis", rn)
        for (r in diag_rows) { if (!grepl(":", rn[r])) { p_diag <- ct[r, p_col]; break } }
        age_rows <- grep("^AgeGroup", rn)
        for (r in age_rows)  { if (!grepl(":", rn[r])) { p_age  <- ct[r, p_col]; break } }
        int_rows <- grep("Diagnosis.*:.*AgeGroup|AgeGroup.*:.*Diagnosis", rn)
        if (length(int_rows) > 0) p_int <- ct[int_rows[1], p_col]
        list(p_diag = p_diag, p_age = p_age, p_int = p_int)
      }

      mean_ct <- tryCatch(summary(fit)$coefficients, error = function(e) NULL)
      mean_ps <- extract_ps_edge(mean_ct)
      disp_ct <- tryCatch(summary(fit$dispersion.fit)$coefficients, error = function(e) NULL)
      disp_ps <- extract_ps_edge(disp_ct)

      data.frame(
        edge_id               = k,
        i                     = upper_idx[k, 1],
        j                     = upper_idx[k, 2],
        region_i              = roi_names[upper_idx[k, 1]],
        region_j              = roi_names[upper_idx[k, 2]],
        n                     = length(y2),
        p_mean_Diagnosis      = mean_ps$p_diag,
        p_mean_AgeGroup       = mean_ps$p_age,
        p_mean_Interaction    = mean_ps$p_int,
        p_disp_Diagnosis      = disp_ps$p_diag,
        p_disp_AgeGroup       = disp_ps$p_age,
        p_disp_Interaction    = disp_ps$p_int
      )
    })

    out_list[[ci]] <- bind_rows(chunk_res)
  }

  out <- bind_rows(out_list)
  out$p_mean_Diagnosis_FDR   <- p.adjust(out$p_mean_Diagnosis,   method = "fdr")
  out$p_mean_AgeGroup_FDR    <- p.adjust(out$p_mean_AgeGroup,    method = "fdr")
  out$p_mean_Interaction_FDR <- p.adjust(out$p_mean_Interaction, method = "fdr")
  out$p_disp_Diagnosis_FDR   <- p.adjust(out$p_disp_Diagnosis,   method = "fdr")
  out$p_disp_AgeGroup_FDR    <- p.adjust(out$p_disp_AgeGroup,    method = "fdr")
  out$p_disp_Interaction_FDR <- p.adjust(out$p_disp_Interaction, method = "fdr")

  out
}

cat("\n=== Running DGLM (edges) ===\n")
edge_res <- run_dglm_edgewise_chunked(edge_mat, df2, upper_idx, roi_names, chunk_size = 2000)

cat("Significant mean Interaction (edges, FDR<0.05):",
    sum(edge_res$p_mean_Interaction_FDR < 0.05, na.rm = TRUE), "\n")
cat("Significant dispersion Interaction (edges, FDR<0.05):",
    sum(edge_res$p_disp_Interaction_FDR < 0.05, na.rm = TRUE), "\n")

write.csv(edge_res, out_edge_csv, row.names = FALSE)
cat("Saved:", out_edge_csv, "\n")

# ============================================================
# 8) Significant results summary tables
# ============================================================
threshold <- 0.05

# --- Degree significant results ---
deg_sig_cols <- c("p_mean_Diagnosis_FDR", "p_mean_AgeGroup_FDR", "p_mean_Interaction_FDR",
                   "p_disp_Diagnosis_FDR", "p_disp_AgeGroup_FDR", "p_disp_Interaction_FDR")
deg_sig_labels <- c("mean_Diagnosis", "mean_AgeGroup", "mean_Interaction",
                    "disp_Diagnosis", "disp_AgeGroup", "disp_Interaction")

cat("\n==============================\n")
cat("Significant Degree Results (FDR <", threshold, ")\n")
cat("==============================\n")
for (ci in seq_along(deg_sig_cols)) {
  col <- deg_sig_cols[ci]
  sig <- degree_res[!is.na(degree_res[[col]]) & degree_res[[col]] < threshold, ]
  cat(sprintf("\n[%s] n_sig = %d\n", deg_sig_labels[ci], nrow(sig)))
  if (nrow(sig) > 0) {
    p_raw_col <- sub("_FDR", "", col)
    print(data.frame(
      feature   = sig$feature,
      n         = sig$n,
      p_raw     = round(sig[[p_raw_col]], 4),
      p_FDR     = round(sig[[col]], 4)
    ), row.names = FALSE)
  }
}

# Save significant degree table (any column significant)
deg_any_sig <- degree_res[
  rowSums(sapply(deg_sig_cols, function(col) {
    !is.na(degree_res[[col]]) & degree_res[[col]] < threshold
  })) > 0, ]
out_deg_sig_csv <- file.path(out_dir, "DGLM_DK318_degree_significant.csv")
write.csv(deg_any_sig, out_deg_sig_csv, row.names = FALSE)
cat(sprintf("\nDegree: %d ROIs with at least one significant effect. Saved: %s\n",
            nrow(deg_any_sig), out_deg_sig_csv))

# --- Edge significant results ---
edge_sig_cols <- c("p_mean_Diagnosis_FDR", "p_mean_AgeGroup_FDR", "p_mean_Interaction_FDR",
                   "p_disp_Diagnosis_FDR", "p_disp_AgeGroup_FDR", "p_disp_Interaction_FDR")
edge_sig_labels <- c("mean_Diagnosis", "mean_AgeGroup", "mean_Interaction",
                     "disp_Diagnosis", "disp_AgeGroup", "disp_Interaction")

cat("\n==============================\n")
cat("Significant Edge Results (FDR <", threshold, ")\n")
cat("==============================\n")
for (ci in seq_along(edge_sig_cols)) {
  col <- edge_sig_cols[ci]
  sig <- edge_res[!is.na(edge_res[[col]]) & edge_res[[col]] < threshold, ]
  cat(sprintf("\n[%s] n_sig = %d\n", edge_sig_labels[ci], nrow(sig)))
  if (nrow(sig) > 0 && nrow(sig) <= 20) {
    p_raw_col <- sub("_FDR", "", col)
    print(data.frame(
      region_i  = sig$region_i,
      region_j  = sig$region_j,
      n         = sig$n,
      p_raw     = round(sig[[p_raw_col]], 4),
      p_FDR     = round(sig[[col]], 4)
    ), row.names = FALSE)
  } else if (nrow(sig) > 20) {
    p_raw_col <- sub("_FDR", "", col)
    top <- sig[order(sig[[col]]), ][1:20, ]
    cat("  (showing top 20 by FDR)\n")
    print(data.frame(
      region_i  = top$region_i,
      region_j  = top$region_j,
      n         = top$n,
      p_raw     = round(top[[p_raw_col]], 4),
      p_FDR     = round(top[[col]], 4)
    ), row.names = FALSE)
  }
}

# Save significant edge table (any column significant)
edge_any_sig <- edge_res[
  rowSums(sapply(edge_sig_cols, function(col) {
    !is.na(edge_res[[col]]) & edge_res[[col]] < threshold
  })) > 0, ]
out_edge_sig_csv <- file.path(out_dir, "DGLM_DK318_edge_significant.csv")
write.csv(edge_any_sig, out_edge_sig_csv, row.names = FALSE)
cat(sprintf("\nEdge: %d edges with at least one significant effect. Saved: %s\n",
            nrow(edge_any_sig), out_edge_sig_csv))

cat("\nDONE.\n")
