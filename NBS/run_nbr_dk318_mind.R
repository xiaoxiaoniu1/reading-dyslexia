# ============================================================
# NBR / NBS analysis for DK318 MIND matrices
#
# Updated goals:
#   1) Full factorial NBS model:
#        Diagnosis * AgeGroup * Sex
#   2) Four planned Diagnosis simple effects:
#        Child Female: DD vs TD
#        Child Male:   DD vs TD
#        Adult Female: DD vs TD
#        Adult Male:   DD vs TD
#   3) Post-NBS summaries and figures:
#        component summary, significant edge list, node degree,
#        component strength violin plot, chord plot, MNI projection plot,
#        lobe-level heatmap, BrainNet Viewer .node/.edge files
#
# Input and output paths from the original script are kept unchanged.
# ROI coordinate/name paths are added as requested.
# ============================================================

packages <- c(
  "readxl", "dplyr", "stringr", "parallel", "NBR",
  "ggplot2", "tidyr", "tibble", "purrr", "circlize",
  "igraph", "emmeans"
)

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(parallel)
  library(NBR)
  library(ggplot2)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(circlize)
  library(igraph)
  library(emmeans)
})

# ----------------------------
# 1) Paths and parameters
# ----------------------------
demo_file <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/all_data_cqt_mean_1.5.xlsx"
mind_dir  <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_combat"
out_dir   <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NBS"

roi_coord_file <- "/data/home/tqi/data1/share/after_freesurfer/FILE/DK-318/DK318_MNI_Coordinates.csv"
roi_name_file  <- "/data/home/tqi/data1/share/after_freesurfer/FILE/DK-318/DK318_roi_names.csv"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
brainnet_dir <- file.path(out_dir, "brainnet_files")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(brainnet_dir, showWarnings = FALSE, recursive = TRUE)

# Test setting requested by the user. For formal analysis, increase to 5000 or 10000.
nperm_full <- 10L
nperm_posthoc <- 10L
thrP_nbs <- 0.01
component_p_alpha <- 0.05
component_p_metric <- "ncompFWE"  # choices: "ncompFWE" or "strnFWE"
full_model <- "~ Diagnosis * AgeGroup * Sex"

# Full model uses NBR::nbr_lm.
run_full_nbr <- TRUE

# Planned contrasts use a model-matrix linear contrast equivalent to emmeans.
# This avoids calling emmeans separately for 50,403 edges.
run_posthoc_contrast_nbs <- TRUE

# Optional validation: compare the first edge with emmeans output.
validate_emmeans_one_edge <- FALSE

# Plot limits. Chord and MNI plots are capped for readability.
plot_max_edges <- 400L

set.seed(20250101)

available_cores <- parallel::detectCores(logical = FALSE)
cores_to_use <- if (is.na(available_cores)) NULL else max(1L, available_cores - 1L)

# ----------------------------
# 2) Helper functions
# ----------------------------
read_mind_csv <- function(fp) {
  mat <- as.matrix(read.csv(fp, row.names = 1, check.names = FALSE))
  storage.mode(mat) <- "numeric"
  mat
}

vectorize_upper_triangle <- function(mat) {
  mat[upper.tri(mat, diag = FALSE)]
}

safe_colnames <- function(mat) {
  rn <- rownames(mat)
  cn <- colnames(mat)
  if (is.null(rn) || is.null(cn)) {
    stop("Matrix must contain both row names and column names.")
  }
  if (!identical(rn, cn)) {
    stop("Matrix row names and column names are not identical.")
  }
  rn
}

sanitize_name <- function(x) {
  x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

save_text_summary <- function(object, file) {
  txt <- capture.output({
    cat("===== class =====\n")
    print(class(object))
    cat("\n===== names =====\n")
    print(names(object))
    cat("\n===== structure =====\n")
    str(object, max.level = 2)
    cat("\n===== print(object) =====\n")
    print(object)
  })
  writeLines(txt, con = file)
}

infer_dk_lobe <- function(roi) {
  base <- as.character(roi)
  base <- gsub("^(lh|rh)[_\\.]", "", base)
  base <- gsub("_part[0-9]+$", "", base)

  frontal <- c(
    "caudalmiddlefrontal", "rostralmiddlefrontal", "superiorfrontal",
    "lateralorbitofrontal", "medialorbitofrontal", "frontalpole",
    "parsopercularis", "parsorbitalis", "parstriangularis", "precentral"
  )
  temporal <- c(
    "bankssts", "entorhinal", "fusiform", "inferiortemporal",
    "middletemporal", "parahippocampal", "superiortemporal",
    "temporalpole", "transversetemporal"
  )
  parietal <- c(
    "inferiorparietal", "postcentral", "precuneus",
    "superiorparietal", "supramarginal"
  )
  occipital <- c("cuneus", "lateraloccipital", "lingual", "pericalcarine")
  cingulate <- c(
    "caudalanteriorcingulate", "rostralanteriorcingulate",
    "posteriorcingulate", "isthmuscingulate"
  )

  out <- rep("Other", length(base))
  out[base %in% frontal] <- "Frontal"
  out[base %in% temporal] <- "Temporal"
  out[base %in% parietal] <- "Parietal"
  out[base %in% occipital] <- "Occipital"
  out[base %in% cingulate] <- "Cingulate"
  out[base %in% "insula"] <- "Insula"
  out
}

infer_canonical_network <- function(roi) {
  base <- as.character(roi)
  base <- gsub("^(lh|rh)[_\\.]", "", base)
  base <- gsub("_part[0-9]+$", "", base)

  vis <- c("cuneus", "lateraloccipital", "lingual", "pericalcarine")
  somatomotor <- c("precentral", "postcentral", "paracentral", "transversetemporal")
  dorsal_attention <- c("superiorparietal", "inferiorparietal", "precuneus", "frontaleyefields")
  ventral_attention <- c("bankssts", "supramarginal", "insula", "parsopercularis", "caudalmiddlefrontal")
  limbic <- c("entorhinal", "temporalpole", "parahippocampal", "fusiform", "medialorbitofrontal")
  default_mode <- c(
    "posteriorcingulate", "isthmuscingulate", "rostralanteriorcingulate",
    "caudalanteriorcingulate", "medialorbitofrontal", "inferiorparietal",
    "precuneus", "middletemporal", "superiorfrontal", "rostralmiddlefrontal"
  )
  frontoparietal <- c(
    "rostralmiddlefrontal", "caudalmiddlefrontal", "superiorfrontal",
    "inferiorparietal", "parsopercularis", "parsorbitalis", "parstriangularis",
    "lateralorbitofrontal", "frontalpole"
  )

  out <- rep("Other", length(base))
  out[base %in% vis] <- "VIS"
  out[base %in% somatomotor] <- "SMN"
  out[base %in% ventral_attention] <- "VAN"
  out[base %in% limbic] <- "Limbic"
  out[base %in% default_mode] <- "DMN"
  out[base %in% dorsal_attention] <- "DAN"
  out[base %in% frontoparietal] <- "FPN"
  out
}

get_network_palette <- function() {
  c(
    VIS = 1,
    SMN = 2,
    DAN = 3,
    VAN = 4,
    Limbic = 5,
    DMN = 6,
    FPN = 7,
    Other = 8
  )
}

scale_to_range <- function(x, to = c(1, 6), default = mean(to)) {
  x <- as.numeric(x)
  if (length(x) == 0L) return(numeric(0))
  if (all(is.na(x))) return(rep(default, length(x)))
  rng <- range(x, na.rm = TRUE)
  if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0) {
    return(rep(default, length(x)))
  }
  scaled <- (x - rng[1]) / diff(rng)
  to[1] + scaled * diff(to)
}


prepare_roi_metadata <- function(roi_names, roi_coord_file, roi_name_file, nnodes) {
  roi_meta <- data.frame(
    ROI_index = seq_along(roi_names),
    ROI = as.character(roi_names),
    stringsAsFactors = FALSE
  )

  if (file.exists(roi_name_file)) {
    roi_name_df <- read.csv(roi_name_file, check.names = FALSE, stringsAsFactors = FALSE)
    roi_name_col <- names(roi_name_df)[1]
    roi_name_df <- data.frame(
      ROI_index_from_file = seq_len(nrow(roi_name_df)),
      ROI_from_file = as.character(roi_name_df[[roi_name_col]]),
      stringsAsFactors = FALSE
    )
    write.csv(roi_name_df, file.path(tab_dir, "DK318_roi_names_from_server.csv"), row.names = FALSE)

    if (nrow(roi_name_df) == nnodes && !identical(roi_name_df$ROI_from_file, roi_meta$ROI)) {
      warning("ROI names in DK318_roi_names.csv are not identical to matrix row names. Matrix order is used as the primary order.")
    }
  } else {
    warning("ROI name file does not exist: ", roi_name_file)
  }

  if (file.exists(roi_coord_file)) {
    coord <- read.csv(roi_coord_file, check.names = FALSE, stringsAsFactors = FALSE)
    needed <- c("label", "x", "y", "z", "hemi")
    if (!all(needed %in% names(coord))) {
      stop("DK318_MNI_Coordinates.csv must contain columns: label, x, y, z, hemi")
    }
    coord <- coord %>%
      mutate(
        hemi_prefix = case_when(
          tolower(hemi) %in% c("left", "lh", "l") ~ "lh",
          tolower(hemi) %in% c("right", "rh", "r") ~ "rh",
          TRUE ~ as.character(hemi)
        ),
        ROI = paste0(hemi_prefix, "_", label),
        dk_region = gsub("_part[0-9]+$", "", label),
        lobe = infer_dk_lobe(ROI),
        network = infer_canonical_network(ROI)
      ) %>%
      select(ROI, x, y, z, hemi, label, dk_region, lobe, network, n_vertices)

    roi_meta <- roi_meta %>% left_join(coord, by = "ROI")

    if (any(is.na(roi_meta$x)) && nrow(coord) == nnodes) {
      warning("Some ROI coordinates were not matched by name. Coordinates will be filled by file order for unmatched rows.")
      miss <- which(is.na(roi_meta$x))
      roi_meta$x[miss] <- coord$x[miss]
      roi_meta$y[miss] <- coord$y[miss]
      roi_meta$z[miss] <- coord$z[miss]
      roi_meta$hemi[miss] <- coord$hemi[miss]
      roi_meta$label[miss] <- coord$label[miss]
      roi_meta$dk_region[miss] <- coord$dk_region[miss]
      roi_meta$lobe[miss] <- coord$lobe[miss]
      roi_meta$network[miss] <- coord$network[miss]
      roi_meta$n_vertices[miss] <- coord$n_vertices[miss]
    }
  } else {
    warning("ROI coordinate file does not exist: ", roi_coord_file)
    roi_meta$x <- NA_real_
    roi_meta$y <- NA_real_
    roi_meta$z <- NA_real_
    roi_meta$hemi <- NA_character_
    roi_meta$label <- roi_meta$ROI
    roi_meta$dk_region <- gsub("_part[0-9]+$", "", gsub("^(lh|rh)_", "", roi_meta$ROI))
    roi_meta$lobe <- infer_dk_lobe(roi_meta$ROI)
    roi_meta$network <- infer_canonical_network(roi_meta$ROI)
    roi_meta$n_vertices <- NA_real_
  }

  roi_meta$lobe[is.na(roi_meta$lobe)] <- "Other"
  roi_meta$network[is.na(roi_meta$network)] <- infer_canonical_network(roi_meta$ROI[is.na(roi_meta$network)])
  roi_meta$network[is.na(roi_meta$network)] <- "Other"
  roi_meta
}

make_edge_annotation <- function(roi_meta, nnodes) {
  tri_pos <- which(upper.tri(matrix(0, nnodes, nnodes), diag = FALSE), arr.ind = TRUE)
  edge_annot <- data.frame(
    edge_index = seq_len(nrow(tri_pos)),
    i = tri_pos[, 1],
    j = tri_pos[, 2],
    ROI1 = roi_meta$ROI[tri_pos[, 1]],
    ROI2 = roi_meta$ROI[tri_pos[, 2]],
    edge_label = paste(roi_meta$ROI[tri_pos[, 1]], roi_meta$ROI[tri_pos[, 2]], sep = "__"),
    stringsAsFactors = FALSE
  )

  roi1_meta <- roi_meta %>%
    rename_with(~ paste0(.x, "1"), -ROI_index) %>%
    rename(i = ROI_index, ROI1_meta = ROI1)
  roi2_meta <- roi_meta %>%
    rename_with(~ paste0(.x, "2"), -ROI_index) %>%
    rename(j = ROI_index, ROI2_meta = ROI2)

  edge_annot <- edge_annot %>%
    left_join(roi1_meta, by = "i") %>%
    left_join(roi2_meta, by = "j") %>%
    mutate(
      lobe_pair = ifelse(lobe1 <= lobe2, paste(lobe1, lobe2, sep = "--"), paste(lobe2, lobe1, sep = "--")),
      dk_pair = ifelse(dk_region1 <= dk_region2, paste(dk_region1, dk_region2, sep = "--"), paste(dk_region2, dk_region1, sep = "--"))
    )

  edge_annot
}

run_nbr_analysis <- function(net, nnodes, idata, mod, result_prefix,
                             nperm, thrP = 0.01, cores = NULL) {
  cat("\n==============================\n")
  cat("Running:", result_prefix, "\n")
  cat("Model:", mod, "\n")
  cat("Permutations:", nperm, "\n")
  cat("thrP:", thrP, "\n")
  cat("==============================\n")

  nbr_res <- NBR::nbr_lm(
    net = net,
    nnodes = nnodes,
    idata = idata,
    mod = mod,
    diag = FALSE,
    nperm = nperm,
    thrP = thrP,
    thrT = NULL,
    cores = cores,
    nudist = TRUE,
    verbose = TRUE,
    na.action = na.exclude
  )

  saveRDS(nbr_res, file = file.path(out_dir, paste0(result_prefix, "_nbr_result.rds")))
  save_text_summary(nbr_res, file = file.path(out_dir, paste0(result_prefix, "_nbr_summary.txt")))

  edge_res <- NBR::edge_lm(
    net = net,
    nnodes = nnodes,
    idata = idata,
    mod = mod,
    diag = FALSE,
    padj = "fdr",
    cores = cores,
    verbose = TRUE,
    na.action = na.exclude
  )

  write.csv(edge_res, file.path(out_dir, paste0(result_prefix, "_edge_lm_fdr.csv")), row.names = FALSE)

  invisible(list(nbr = nbr_res, edge = edge_res))
}

as_component_edges <- function(comp_obj, fwe_obj, effect_name, edge_annot,
                               p_metric = "ncompFWE", fwe_alpha = 0.05) {
  if (is.null(comp_obj) || length(comp_obj) == 0L) {
    return(list(summary = data.frame(), edges = data.frame()))
  }

  comp_df <- as.data.frame(comp_obj)
  if (nrow(comp_df) == 0L) {
    return(list(summary = data.frame(), edges = data.frame()))
  }

  names(comp_df) <- c("edge_index", "i", "j", "Component", "strn")
  comp_df <- comp_df %>% mutate(across(c(edge_index, i, j, Component), as.integer))

  fwe_df <- as.data.frame(fwe_obj)
  if (nrow(fwe_df) == 0L) {
    fwe_df <- comp_df %>%
      group_by(Component) %>%
      summarise(ncomp = n(), strn = sum(abs(strn), na.rm = TRUE), .groups = "drop") %>%
      mutate(ncompFWE = NA_real_, strnFWE = NA_real_)
  }

  fwe_df <- fwe_df %>%
    mutate(
      Component = as.integer(Component),
      effect = effect_name,
      component_sig = .data[[p_metric]] < fwe_alpha
    )

  edge_df <- comp_df %>%
    left_join(fwe_df, by = "Component") %>%
    left_join(edge_annot, by = "edge_index", suffix = c("", "_annot")) %>%
    mutate(effect = effect_name)

  list(summary = fwe_df, edges = edge_df)
}

export_full_nbr_components <- function(nbr_res, result_prefix, edge_annot,
                                       p_metric = "ncompFWE", fwe_alpha = 0.05) {
  if (is.null(nbr_res$components) || is.null(nbr_res$fwe)) {
    warning("NBR result does not contain expected components/fwe fields.")
    return(invisible(NULL))
  }

  terms <- names(nbr_res$components)
  all_summary <- list()
  all_edges <- list()

  for (term in terms) {
    out <- as_component_edges(
      comp_obj = nbr_res$components[[term]],
      fwe_obj = nbr_res$fwe[[term]],
      effect_name = term,
      edge_annot = edge_annot,
      p_metric = p_metric,
      fwe_alpha = fwe_alpha
    )

    term_safe <- sanitize_name(term)
    write.csv(out$summary, file.path(tab_dir, paste0(result_prefix, "_", term_safe, "_component_summary.csv")), row.names = FALSE)
    write.csv(out$edges, file.path(tab_dir, paste0(result_prefix, "_", term_safe, "_suprathreshold_edges.csv")), row.names = FALSE)

    all_summary[[term]] <- out$summary
    all_edges[[term]] <- out$edges
  }

  all_summary_df <- bind_rows(all_summary)
  all_edges_df <- bind_rows(all_edges)
  write.csv(all_summary_df, file.path(tab_dir, paste0(result_prefix, "_ALL_component_summary.csv")), row.names = FALSE)
  write.csv(all_edges_df, file.path(tab_dir, paste0(result_prefix, "_ALL_suprathreshold_edges.csv")), row.names = FALSE)

  invisible(list(summary = all_summary_df, edges = all_edges_df))
}

# ----------------------------
# 3) emmeans-equivalent linear contrasts and NBS
# ----------------------------
make_design_and_contrast <- function(idata, model_formula, age_group, sex) {
  idata2 <- idata %>%
    mutate(
      Diagnosis = factor(Diagnosis, levels = c("TD", "DD")),
      AgeGroup = factor(AgeGroup, levels = c("Child", "Adult")),
      Sex = factor(Sex, levels = c("Female", "Male"))
    )

  f <- as.formula(model_formula)
  X <- model.matrix(f, data = idata2)
  if (qr(X)$rank < ncol(X)) {
    stop("The design matrix is rank deficient. Check whether all Diagnosis x AgeGroup x Sex cells exist.")
  }

  nd <- data.frame(
    Diagnosis = factor(c("TD", "DD"), levels = levels(idata2$Diagnosis)),
    AgeGroup = factor(c(age_group, age_group), levels = levels(idata2$AgeGroup)),
    Sex = factor(c(sex, sex), levels = levels(idata2$Sex))
  )
  Xnew <- model.matrix(f, data = nd)

  missing_cols <- setdiff(colnames(X), colnames(Xnew))
  if (length(missing_cols) > 0L) {
    for (cc in missing_cols) Xnew <- cbind(Xnew, 0)
  }
  Xnew <- Xnew[, colnames(X), drop = FALSE]

  cvec <- Xnew[2, ] - Xnew[1, ]
  names(cvec) <- colnames(X)

  list(idata = idata2, X = X, cvec = cvec)
}

compute_contrast_stats_from_design <- function(net, X, cvec) {
  XtX <- crossprod(X)
  XtX_inv <- solve(XtX)
  B <- XtX_inv %*% crossprod(X, net)
  resid <- net - X %*% B
  df_resid <- nrow(X) - qr(X)$rank
  sigma2 <- colSums(resid * resid, na.rm = TRUE) / df_resid
  cvar <- as.numeric(t(cvec) %*% XtX_inv %*% cvec)
  estimate <- as.numeric(t(cvec) %*% B)
  se <- sqrt(sigma2 * cvar)
  t_value <- estimate / se
  p_value <- 2 * pt(abs(t_value), df = df_resid, lower.tail = FALSE)
  list(estimate = estimate, se = se, t = t_value, p = p_value, df = df_resid)
}

find_components_from_edges <- function(edge_indices, t_values, tri_pos, nnodes, tcrit) {
  if (length(edge_indices) == 0L) {
    return(list(edges = data.frame(), summary = data.frame()))
  }

  pairs <- tri_pos[edge_indices, , drop = FALSE]
  g <- igraph::graph_from_data_frame(
    data.frame(from = as.character(pairs[, 1]), to = as.character(pairs[, 2])),
    directed = FALSE,
    vertices = data.frame(name = as.character(seq_len(nnodes)))
  )
  memb <- igraph::components(g)$membership
  raw_comp <- memb[as.character(pairs[, 1])]
  comp <- as.integer(factor(raw_comp))
  strn_signed <- (abs(t_values[edge_indices]) - tcrit) * sign(t_values[edge_indices])

  edge_df <- data.frame(
    edge_index = edge_indices,
    i = pairs[, 1],
    j = pairs[, 2],
    Component = comp,
    strn = strn_signed,
    stringsAsFactors = FALSE
  )

  summ <- edge_df %>%
    group_by(Component) %>%
    summarise(
      ncomp = n(),
      strn = sum(abs(strn), na.rm = TRUE),
      mean_signed_strn = mean(strn, na.rm = TRUE),
      .groups = "drop"
    )

  list(edges = edge_df, summary = summ)
}

run_contrast_nbs <- function(net, nnodes, idata, model_formula,
                             contrast_name, age_group, sex,
                             nperm, thrP, edge_annot,
                             p_metric = "ncompFWE", fwe_alpha = 0.05) {
  cat("\n==============================\n")
  cat("Running planned contrast NBS:", contrast_name, "\n")
  cat("Contrast: DD - TD at AgeGroup =", age_group, ", Sex =", sex, "\n")
  cat("Model:", model_formula, "\n")
  cat("Permutations:", nperm, "\n")
  cat("thrP:", thrP, "\n")
  cat("==============================\n")

  design <- make_design_and_contrast(idata, model_formula, age_group, sex)
  obs <- compute_contrast_stats_from_design(net, design$X, design$cvec)
  tri_pos <- which(upper.tri(matrix(0, nnodes, nnodes), diag = FALSE), arr.ind = TRUE)
  tcrit <- qt(1 - thrP / 2, df = obs$df)

  edge_stats <- data.frame(
    edge_index = seq_len(ncol(net)),
    contrast = contrast_name,
    AgeGroup = age_group,
    Sex = sex,
    estimate = obs$estimate,
    SE = obs$se,
    df = obs$df,
    t = obs$t,
    p = obs$p,
    p_fdr = p.adjust(obs$p, method = "fdr"),
    stringsAsFactors = FALSE
  ) %>% left_join(edge_annot, by = "edge_index")

  supra_edges <- which(obs$p < thrP)
  obs_comp <- find_components_from_edges(supra_edges, obs$t, tri_pos, nnodes, tcrit)

  null_dist <- data.frame(
    perm = seq_len(nperm),
    max_ncomp = 0,
    max_strn = 0
  )

  for (pp in seq_len(nperm)) {
    if (pp %% 10 == 0L) cat("Permutation", pp, "of", nperm, "\n")
    perm_idata <- design$idata[sample(seq_len(nrow(design$idata))), , drop = FALSE]
    Xp <- model.matrix(as.formula(model_formula), data = perm_idata)
    Xp <- Xp[, colnames(design$X), drop = FALSE]
    perm <- compute_contrast_stats_from_design(net, Xp, design$cvec)
    perm_edges <- which(perm$p < thrP)
    perm_comp <- find_components_from_edges(perm_edges, perm$t, tri_pos, nnodes, tcrit)
    if (nrow(perm_comp$summary) > 0L) {
      null_dist$max_ncomp[pp] <- max(perm_comp$summary$ncomp, na.rm = TRUE)
      null_dist$max_strn[pp] <- max(perm_comp$summary$strn, na.rm = TRUE)
    }
  }

  comp_summary <- obs_comp$summary
  comp_edges <- obs_comp$edges

  if (nrow(comp_summary) > 0L) {
    comp_summary <- comp_summary %>%
      mutate(
        ncompFWE = sapply(ncomp, function(x) mean(null_dist$max_ncomp > x)),
        strnFWE = sapply(strn, function(x) mean(null_dist$max_strn > x)),
        contrast = contrast_name,
        AgeGroup = age_group,
        Sex = sex,
        component_sig = .data[[p_metric]] < fwe_alpha
      )

    comp_edges <- comp_edges %>%
      left_join(comp_summary, by = "Component", suffix = c("", "_component")) %>%
      left_join(edge_stats %>% select(edge_index, estimate, SE, df, t, p, p_fdr, starts_with("ROI"), starts_with("x"), starts_with("y"), starts_with("z"), starts_with("hemi"), starts_with("label"), starts_with("dk_region"), starts_with("lobe"), edge_label, lobe_pair, dk_pair), by = "edge_index") %>%
      mutate(
        contrast = contrast_name,
        AgeGroup = age_group,
        Sex = sex,
        direction = case_when(
          t > 0 ~ "DD > TD",
          t < 0 ~ "DD < TD",
          TRUE ~ "zero"
        )
      )
  } else {
    comp_summary <- data.frame()
    comp_edges <- data.frame()
  }

  result <- list(
    contrast = contrast_name,
    AgeGroup = age_group,
    Sex = sex,
    edge_stats = edge_stats,
    components = comp_edges,
    fwe = comp_summary,
    nudist = null_dist,
    model_formula = model_formula,
    contrast_vector = design$cvec,
    thrP = thrP,
    nperm = nperm
  )

  safe <- sanitize_name(contrast_name)
  saveRDS(result, file.path(out_dir, paste0("DK318_MIND_posthoc_", safe, "_contrast_nbs_result.rds")))
  write.csv(edge_stats, file.path(tab_dir, paste0("DK318_MIND_posthoc_", safe, "_edgewise_emmeans_equivalent.csv")), row.names = FALSE)
  write.csv(comp_summary, file.path(tab_dir, paste0("DK318_MIND_posthoc_", safe, "_component_summary.csv")), row.names = FALSE)
  write.csv(comp_edges, file.path(tab_dir, paste0("DK318_MIND_posthoc_", safe, "_suprathreshold_edges.csv")), row.names = FALSE)
  write.csv(null_dist, file.path(tab_dir, paste0("DK318_MIND_posthoc_", safe, "_null_distribution.csv")), row.names = FALSE)

  invisible(result)
}

validate_first_edge_with_emmeans <- function(net, idata, model_formula) {
  tmp <- data.frame(y = net[, 1], idata)
  fit <- lm(as.formula(paste0("y ", model_formula)), data = tmp)
  emm <- emmeans::emmeans(fit, ~ Diagnosis | AgeGroup * Sex)
  con <- emmeans::contrast(emm, method = list(DD_vs_TD = c(-1, 1)))
  print(summary(con))
}

# ----------------------------
# 4) Plotting and post-NBS outputs
# ----------------------------
get_significant_edges <- function(comp_edges, p_metric = "ncompFWE", fwe_alpha = 0.05) {
  if (is.null(comp_edges) || nrow(comp_edges) == 0L) return(data.frame())
  comp_edges %>% filter(.data[[p_metric]] < fwe_alpha)
}

make_node_degree_table <- function(edges, roi_meta) {
  if (nrow(edges) == 0L) return(data.frame())

  d1 <- edges %>% transmute(ROI_index = i, ROI = ROI1)
  d2 <- edges %>% transmute(ROI_index = j, ROI = ROI2)

  bind_rows(d1, d2) %>%
    count(ROI_index, ROI, name = "degree") %>%
    left_join(roi_meta, by = c("ROI_index", "ROI")) %>%
    arrange(desc(degree), ROI_index)
}

plot_component_strength <- function(edges_sig, net, pheno, effect_name) {
  if (nrow(edges_sig) == 0L) return(invisible(NULL))

  comps <- sort(unique(edges_sig$Component))
  strength_list <- lapply(comps, function(cc) {
    idx <- unique(edges_sig$edge_index[edges_sig$Component == cc])
    if (length(idx) == 0L) return(NULL)
    vals <- rowMeans(net[, idx, drop = FALSE], na.rm = TRUE)
    data.frame(
      subj_prefix = pheno$subj_prefix,
      Diagnosis = pheno$Diagnosis,
      AgeGroup = pheno$AgeGroup,
      Sex = pheno$Sex,
      Group8 = interaction(pheno$AgeGroup, pheno$Sex, pheno$Diagnosis, sep = "_", drop = TRUE),
      Component = paste0("Component_", cc),
      component_strength = vals,
      stringsAsFactors = FALSE
    )
  })

  strength_df <- bind_rows(strength_list)
  if (nrow(strength_df) == 0L) return(invisible(NULL))

  safe <- sanitize_name(effect_name)
  write.csv(strength_df, file.path(tab_dir, paste0(safe, "_subject_component_strength.csv")), row.names = FALSE)

  p <- ggplot(strength_df, aes(x = Group8, y = component_strength)) +
    geom_violin(trim = FALSE, alpha = 0.35) +
    geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.6) +
    geom_jitter(width = 0.12, size = 0.8, alpha = 0.55) +
    facet_wrap(~ Component, scales = "free_y") +
    labs(
      title = paste0(effect_name, " component strength"),
      x = "AgeGroup_Sex_Diagnosis",
      y = "Mean MIND value within NBS component"
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(fig_dir, paste0(safe, "_component_strength_violin.png")), p, width = 11, height = 6, dpi = 300)
  invisible(strength_df)
}

plot_lobe_heatmap <- function(edges_sig, edge_annot, effect_name) {
  if (nrow(edges_sig) == 0L) return(invisible(NULL))

  possible <- edge_annot %>%
    count(lobe_pair, name = "possible_edges")

  obs <- edges_sig %>%
    count(lobe_pair, name = "observed_edges") %>%
    left_join(possible, by = "lobe_pair") %>%
    mutate(
      normalized_proportion = observed_edges / possible_edges,
      lobe_a = sapply(strsplit(lobe_pair, "--", fixed = TRUE), `[`, 1),
      lobe_b = sapply(strsplit(lobe_pair, "--", fixed = TRUE), `[`, 2)
    )

  safe <- sanitize_name(effect_name)
  write.csv(obs, file.path(tab_dir, paste0(safe, "_lobe_level_edge_counts.csv")), row.names = FALSE)

  p1 <- ggplot(obs, aes(x = lobe_a, y = lobe_b, fill = observed_edges)) +
    geom_tile() +
    geom_text(aes(label = observed_edges), size = 3) +
    coord_equal() +
    labs(title = paste0(effect_name, " lobe-level raw edge count"), x = NULL, y = NULL, fill = "Edges") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  p2 <- ggplot(obs, aes(x = lobe_a, y = lobe_b, fill = normalized_proportion)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.3f", normalized_proportion)), size = 3) +
    coord_equal() +
    labs(title = paste0(effect_name, " lobe-level normalized proportion"), x = NULL, y = NULL, fill = "Proportion") +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(file.path(fig_dir, paste0(safe, "_lobe_heatmap_raw_count.png")), p1, width = 7, height = 6, dpi = 300)
  ggsave(file.path(fig_dir, paste0(safe, "_lobe_heatmap_normalized.png")), p2, width = 7, height = 6, dpi = 300)

  invisible(obs)
}

plot_mni_projection <- function(edges_sig, roi_meta, effect_name, max_edges = 400L) {
  if (nrow(edges_sig) == 0L) return(invisible(NULL))
  if (any(is.na(roi_meta$x))) {
    warning("MNI coordinates are missing. MNI projection plot skipped for ", effect_name)
    return(invisible(NULL))
  }

  edges_plot <- edges_sig %>%
    mutate(abs_t = abs(t)) %>%
    arrange(desc(abs_t)) %>%
    slice_head(n = max_edges)

  deg <- make_node_degree_table(edges_plot, roi_meta)
  if (nrow(deg) == 0L) return(invisible(NULL))

  nodes <- deg %>%
    select(ROI_index, ROI, x, y, z, degree)

  line_base <- edges_plot %>%
    transmute(
      edge_index,
      Component,
      t,
      abs_t = abs(t),
      x1 = x1, y1 = y1, z1 = z1,
      x2 = x2, y2 = y2, z2 = z2
    )

  line_df <- bind_rows(
    line_base %>% transmute(projection = "Axial X-Y", x = x1, y = y1, xend = x2, yend = y2, abs_t = abs_t),
    line_base %>% transmute(projection = "Coronal X-Z", x = x1, y = z1, xend = x2, yend = z2, abs_t = abs_t),
    line_base %>% transmute(projection = "Sagittal Y-Z", x = y1, y = z1, xend = y2, yend = z2, abs_t = abs_t)
  )

  node_df <- bind_rows(
    nodes %>% transmute(projection = "Axial X-Y", x = x, y = y, degree = degree),
    nodes %>% transmute(projection = "Coronal X-Z", x = x, y = z, degree = degree),
    nodes %>% transmute(projection = "Sagittal Y-Z", x = y, y = z, degree = degree)
  )

  safe <- sanitize_name(effect_name)
  p <- ggplot() +
    geom_segment(
      data = line_df,
      aes(x = x, y = y, xend = xend, yend = yend, linewidth = abs_t),
      alpha = 0.25
    ) +
    geom_point(data = node_df, aes(x = x, y = y, size = degree), alpha = 0.85) +
    facet_wrap(~ projection) +
    coord_equal() +
    labs(
      title = paste0(effect_name, " MNI projection"),
      x = "MNI coordinate",
      y = "MNI coordinate",
      linewidth = "|t|",
      size = "Degree"
    ) +
    theme_bw(base_size = 10)

  ggsave(file.path(fig_dir, paste0(safe, "_MNI_projection.png")), p, width = 10, height = 5, dpi = 300)
  invisible(p)
}

plot_chord <- function(edges_sig, effect_name, max_edges = 400L) {
  if (nrow(edges_sig) == 0L) return(invisible(NULL))

  edges_plot <- edges_sig %>%
    mutate(abs_t = abs(t)) %>%
    arrange(desc(abs_t)) %>%
    slice_head(n = max_edges)

  chord_df <- edges_plot %>%
    transmute(from = ROI1, to = ROI2, value = abs_t)

  safe <- sanitize_name(effect_name)
  png(file.path(fig_dir, paste0(safe, "_chord_plot.png")), width = 2700, height = 2700, res = 300)
  tryCatch({
    circlize::circos.clear()
    circlize::chordDiagram(chord_df, directional = 0, annotationTrack = "grid", transparency = 0.65)
    title(effect_name)
    circlize::circos.clear()
  }, error = function(e) {
    message("Chord plot failed for ", effect_name, ": ", conditionMessage(e))
  })
  dev.off()

  invisible(NULL)
}

write_brainnet_files <- function(edges_sig, roi_meta, effect_name) {
  if (nrow(edges_sig) == 0L) return(invisible(NULL))
  if (any(is.na(roi_meta$x))) {
    warning("MNI coordinates are missing. BrainNet files skipped for ", effect_name)
    return(invisible(NULL))
  }

  safe <- sanitize_name(effect_name)
  comps <- sort(unique(edges_sig$Component))
  network_palette <- get_network_palette()

  component_summary_export <- edges_sig %>%
    distinct(Component, ncomp, strn, ncompFWE, strnFWE, contrast, AgeGroup, Sex, component_sig)
  write.csv(
    component_summary_export,
    file = file.path(brainnet_dir, paste0(safe, "_BrainNet_component_summary.csv")),
    row.names = FALSE
  )

  for (cc in comps) {
    e <- edges_sig %>% filter(Component == cc)
    if (nrow(e) == 0L) next

    node_ids <- sort(unique(c(e$i, e$j)))
    nodes <- roi_meta %>%
      filter(ROI_index %in% node_ids) %>%
      arrange(ROI_index)

    deg <- make_node_degree_table(e, roi_meta) %>%
      select(ROI_index, degree)
    nodes <- nodes %>%
      left_join(deg, by = "ROI_index") %>%
      mutate(
        degree = ifelse(is.na(degree), 0, degree),
        network = ifelse(is.na(network) | network == "", "Other", network),
        brainnet_color = unname(network_palette[network]),
        brainnet_color = ifelse(is.na(brainnet_color), network_palette[["Other"]], brainnet_color),
        brainnet_size = scale_to_range(degree, to = c(2, 8), default = 4)
      )

    idx_map <- setNames(seq_along(nodes$ROI_index), nodes$ROI_index)
    edge_mat <- matrix(0, nrow = nrow(nodes), ncol = nrow(nodes))
    edge_mat_td_gt_dd <- matrix(0, nrow = nrow(nodes), ncol = nrow(nodes))
    edge_mat_td_lt_dd <- matrix(0, nrow = nrow(nodes), ncol = nrow(nodes))

    for (kk in seq_len(nrow(e))) {
      a <- idx_map[as.character(e$i[kk])]
      b <- idx_map[as.character(e$j[kk])]
      edge_weight <- abs(e$t[kk])
      edge_mat[a, b] <- edge_weight
      edge_mat[b, a] <- edge_weight

      if (e$t[kk] < 0) {
        edge_mat_td_gt_dd[a, b] <- edge_weight
        edge_mat_td_gt_dd[b, a] <- edge_weight
      }
      if (e$t[kk] > 0) {
        edge_mat_td_lt_dd[a, b] <- edge_weight
        edge_mat_td_lt_dd[b, a] <- edge_weight
      }
    }

    node_table <- data.frame(
      x = nodes$x,
      y = nodes$y,
      z = nodes$z,
      color = nodes$brainnet_color,
      size = nodes$brainnet_size,
      label = nodes$ROI,
      stringsAsFactors = FALSE
    )

    node_annotation <- nodes %>%
      transmute(
        ROI_index,
        ROI,
        hemi,
        label,
        dk_region,
        lobe,
        network,
        x,
        y,
        z,
        degree,
        brainnet_color,
        brainnet_size
      )

    edge_annotation <- e %>%
      transmute(
        edge_index,
        Component,
        ROI1,
        ROI2,
        i,
        j,
        t,
        abs_t = abs(t),
        p,
        p_fdr,
        direction,
        td_direction = case_when(
          t < 0 ~ "TD > DD",
          t > 0 ~ "TD < DD",
          TRUE ~ "zero"
        ),
        network1,
        network2,
        lobe1,
        lobe2,
        dk_region1,
        dk_region2
      )

    write.table(
      node_table,
      file = file.path(brainnet_dir, paste0(safe, "_component_", cc, ".node")),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    write.table(
      edge_mat,
      file = file.path(brainnet_dir, paste0(safe, "_component_", cc, ".edge")),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    write.table(
      edge_mat_td_gt_dd,
      file = file.path(brainnet_dir, paste0(safe, "_component_", cc, "_TD_gt_DD.edge")),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    write.table(
      edge_mat_td_lt_dd,
      file = file.path(brainnet_dir, paste0(safe, "_component_", cc, "_TD_lt_DD.edge")),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    write.csv(
      node_annotation,
      file = file.path(brainnet_dir, paste0(safe, "_component_", cc, "_node_annotation.csv")),
      row.names = FALSE
    )
    write.csv(
      edge_annotation,
      file = file.path(brainnet_dir, paste0(safe, "_component_", cc, "_edge_annotation.csv")),
      row.names = FALSE
    )
  }

  write.csv(
    data.frame(
      network = names(network_palette),
      brainnet_color_code = unname(network_palette),
      stringsAsFactors = FALSE
    ),
    file = file.path(brainnet_dir, paste0(safe, "_BrainNet_network_color_key.csv")),
    row.names = FALSE
  )

  invisible(NULL)
}

make_post_nbs_outputs <- function(result, roi_meta, edge_annot, net, pheno,
                                  p_metric = "ncompFWE", fwe_alpha = 0.05,
                                  max_edges = 400L) {
  effect_name <- result$contrast
  comp_edges <- result$components
  if (is.null(comp_edges) || nrow(comp_edges) == 0L) {
    cat("No suprathreshold NBS component for ", effect_name, "\n")
    return(invisible(NULL))
  }

  sig_edges <- get_significant_edges(comp_edges, p_metric = p_metric, fwe_alpha = fwe_alpha)
  safe <- sanitize_name(effect_name)
  write.csv(sig_edges, file.path(tab_dir, paste0(safe, "_FWE_significant_edges.csv")), row.names = FALSE)

  if (nrow(sig_edges) == 0L) {
    cat("No FWE-significant component for ", effect_name, "\n")
    return(invisible(NULL))
  }

  degree <- make_node_degree_table(sig_edges, roi_meta)
  write.csv(degree, file.path(tab_dir, paste0(safe, "_node_degree.csv")), row.names = FALSE)

  plot_component_strength(sig_edges, net, pheno, effect_name)
  plot_lobe_heatmap(sig_edges, edge_annot, effect_name)
  plot_mni_projection(sig_edges, roi_meta, effect_name, max_edges = max_edges)
  plot_chord(sig_edges, effect_name, max_edges = max_edges)
  write_brainnet_files(sig_edges, roi_meta, effect_name)

  invisible(NULL)
}

# ----------------------------
# 5) Read phenotype table
# ----------------------------
df <- read_excel(demo_file, sheet = "Sheet1")

df <- df %>%
  mutate(
    original_project = as.character(`original-project`),
    id_old = as.character(id_old),
    Age = as.numeric(age_month),
    subj_prefix = paste0(original_project, "_", id_old),
    file_base = paste0(subj_prefix, "_MIND_DK318_combat"),
    mind_file = file.path(mind_dir, paste0(file_base, ".csv")),
    Diagnosis = ifelse(group_d_or_c == 0, "TD", "DD"),
    AgeGroup = ifelse(group_age == 1, "Adult", "Child"),
    Sex = ifelse(sex == 1, "Male", "Female"),
    Diagnosis = factor(Diagnosis, levels = c("TD", "DD")),
    AgeGroup = factor(AgeGroup, levels = c("Child", "Adult")),
    Sex = factor(Sex, levels = c("Female", "Male")),
    has_file = file.exists(mind_file)
  ) %>%
  filter(!is.na(original_project), !is.na(id_old), original_project != "", id_old != "") %>%
  filter(has_file) %>%
  select(original_project, id_old, subj_prefix, file_base, mind_file,
         Diagnosis, AgeGroup, Sex, Age)

if (nrow(df) == 0L) {
  stop("No valid subjects with DK318 MIND matrix files were found.")
}

cat("Subjects included:", nrow(df), "\n")
cat("Diagnosis counts:\n")
print(table(df$Diagnosis, useNA = "ifany"))
cat("AgeGroup counts:\n")
print(table(df$AgeGroup, useNA = "ifany"))
cat("Sex counts:\n")
print(table(df$Sex, useNA = "ifany"))
cat("Diagnosis x AgeGroup x Sex counts:\n")
print(ftable(df$Diagnosis, df$AgeGroup, df$Sex))

cell_counts <- as.data.frame(ftable(df$Diagnosis, df$AgeGroup, df$Sex))
colnames(cell_counts) <- c("Diagnosis", "AgeGroup", "Sex", "n")
write.csv(cell_counts, file.path(tab_dir, "Diagnosis_AgeGroup_Sex_cell_counts.csv"), row.names = FALSE)
if (any(cell_counts$n == 0L)) {
  warning("At least one Diagnosis x AgeGroup x Sex cell has n = 0. Full factorial model may be rank deficient.")
}

write.csv(df, file.path(out_dir, "included_subjects_and_covariates.csv"), row.names = FALSE)

# ----------------------------
# 6) Build edge matrix
# ----------------------------
first_mat <- read_mind_csv(df$mind_file[1])
roi_names <- safe_colnames(first_mat)
nnodes <- nrow(first_mat)

if (nnodes != 318L) {
  warning("Detected nnodes = ", nnodes, ", not 318. Please double-check the input matrices.")
}

nedge <- nnodes * (nnodes - 1L) / 2L
tri_pos <- which(upper.tri(matrix(0, nnodes, nnodes), diag = FALSE), arr.ind = TRUE)
edge_labels <- apply(tri_pos, 1, function(idx) {
  paste(roi_names[idx[1]], roi_names[idx[2]], sep = "__")
})

net2d <- matrix(NA_real_, nrow = nrow(df), ncol = nedge,
                dimnames = list(df$file_base, edge_labels))

for (i in seq_len(nrow(df))) {
  mat <- read_mind_csv(df$mind_file[i])

  if (!all(dim(mat) == c(nnodes, nnodes))) {
    stop("Matrix dimension mismatch in: ", df$mind_file[i],
         " | expected ", nnodes, "x", nnodes,
         " got ", paste(dim(mat), collapse = "x"))
  }

  these_rois <- safe_colnames(mat)
  if (!identical(these_rois, roi_names)) {
    stop("ROI order mismatch in: ", df$mind_file[i])
  }

  if (!isTRUE(all.equal(mat, t(mat), tolerance = 1e-8))) {
    warning("Matrix is not perfectly symmetric: ", df$mind_file[i])
  }

  diag(mat) <- 0
  net2d[i, ] <- vectorize_upper_triangle(mat)
}

cat("Network matrix dimension:", paste(dim(net2d), collapse = " x "), "\n")

roi_meta <- prepare_roi_metadata(roi_names, roi_coord_file, roi_name_file, nnodes)
edge_annot <- make_edge_annotation(roi_meta, nnodes)

write.csv(roi_meta, file.path(out_dir, "roi_metadata_with_mni_coordinates.csv"), row.names = FALSE)
write.csv(edge_annot, file.path(out_dir, "edge_annotation_roi_mni_lobe.csv"), row.names = FALSE)
saveRDS(net2d, file.path(out_dir, "dk318_mind_upper_triangle_matrix.rds"))

# ----------------------------
# 7) Full factorial NBR model
# ----------------------------
idata_full <- df %>% select(Diagnosis, AgeGroup, Sex)

if (validate_emmeans_one_edge) {
  validate_first_edge_with_emmeans(net2d, idata_full, full_model)
}

if (run_full_nbr) {
  res_full <- run_nbr_analysis(
    net = net2d,
    nnodes = nnodes,
    idata = idata_full,
    mod = full_model,
    result_prefix = "DK318_MIND_full_Diagnosis_by_AgeGroup_by_Sex",
    nperm = nperm_full,
    thrP = thrP_nbs,
    cores = cores_to_use
  )

  export_full_nbr_components(
    nbr_res = res_full$nbr,
    result_prefix = "DK318_MIND_full_Diagnosis_by_AgeGroup_by_Sex",
    edge_annot = edge_annot,
    p_metric = component_p_metric,
    fwe_alpha = component_p_alpha
  )
}

# ----------------------------
# 8) Planned post-hoc simple-effect NBS contrasts
# ----------------------------
posthoc_plan <- data.frame(
  contrast = c(
    "Child_Female_DD_vs_TD",
    "Child_Male_DD_vs_TD",
    "Adult_Female_DD_vs_TD",
    "Adult_Male_DD_vs_TD"
  ),
  AgeGroup = c("Child", "Child", "Adult", "Adult"),
  Sex = c("Female", "Male", "Female", "Male"),
  stringsAsFactors = FALSE
)
write.csv(posthoc_plan, file.path(tab_dir, "planned_posthoc_contrasts.csv"), row.names = FALSE)

posthoc_results <- list()

if (run_posthoc_contrast_nbs) {
  for (kk in seq_len(nrow(posthoc_plan))) {
    res <- run_contrast_nbs(
      net = net2d,
      nnodes = nnodes,
      idata = idata_full,
      model_formula = full_model,
      contrast_name = posthoc_plan$contrast[kk],
      age_group = posthoc_plan$AgeGroup[kk],
      sex = posthoc_plan$Sex[kk],
      nperm = nperm_posthoc,
      thrP = thrP_nbs,
      edge_annot = edge_annot,
      p_metric = component_p_metric,
      fwe_alpha = component_p_alpha
    )
    posthoc_results[[posthoc_plan$contrast[kk]]] <- res

    make_post_nbs_outputs(
      result = res,
      roi_meta = roi_meta,
      edge_annot = edge_annot,
      net = net2d,
      pheno = df,
      p_metric = component_p_metric,
      fwe_alpha = component_p_alpha,
      max_edges = plot_max_edges
    )
  }
}

posthoc_summary <- bind_rows(lapply(posthoc_results, function(x) x$fwe))
write.csv(posthoc_summary, file.path(tab_dir, "DK318_MIND_posthoc_ALL_component_summary.csv"), row.names = FALSE)

# ----------------------------
# 9) Session info
# ----------------------------
analysis_settings <- data.frame(
  setting = c(
    "demo_file", "mind_dir", "out_dir", "roi_coord_file", "roi_name_file",
    "model", "thrP_nbs", "nperm_full", "nperm_posthoc", "component_p_metric", "component_p_alpha"
  ),
  value = c(
    demo_file, mind_dir, out_dir, roi_coord_file, roi_name_file,
    full_model, as.character(thrP_nbs), as.character(nperm_full), as.character(nperm_posthoc),
    component_p_metric, as.character(component_p_alpha)
  )
)
write.csv(analysis_settings, file.path(out_dir, "analysis_settings.csv"), row.names = FALSE)

session_info_file <- file.path(out_dir, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), con = session_info_file)

cat("\nAll analyses finished. Results saved to:\n", out_dir, "\n")
cat("Tables saved to:\n", tab_dir, "\n")
cat("Figures saved to:\n", fig_dir, "\n")
cat("BrainNet files saved to:\n", brainnet_dir, "\n")
