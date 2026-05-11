# ============================================================
# NBR / NBS analysis for DK318 MIND matrices
#
# Statistics + tables + BrainNet-export workflow:
#   1) Overall Diagnosis effect
#   2) Diagnosis-related interactions
#   3) Planned simple effects within AgeGroup x Sex strata
#
# Input and output paths from the original script are kept unchanged.
# Plotting / figure generation code has been removed.
# BrainNet export is retained for downstream network visualization.
# ============================================================

packages <- c(
  "readxl", "dplyr", "stringr", "parallel", "NBR",
  "tidyr", "tibble", "purrr"
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
  library(tidyr)
  library(tibble)
  library(purrr)
})

get_env_flag <- function(var_name, default = TRUE) {
  raw <- Sys.getenv(var_name, unset = if (default) "true" else "false")
  raw <- trimws(tolower(raw))

  if (raw %in% c("true", "t", "1", "yes", "y", "on")) return(TRUE)
  if (raw %in% c("false", "f", "0", "no", "n", "off")) return(FALSE)

  stop(
    "Invalid boolean value for ", var_name, ": '", raw,
    "'. Use true/false, 1/0, yes/no, or on/off."
  )
}

# ----------------------------
# 1) Paths and parameters
# ----------------------------
demo_file <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/all_data_cqt_mean_1.5.xlsx"
mind_dir  <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_combat"
out_dir   <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NBS_1000"

roi_coord_file <- "/data/home/tqi/data1/share/after_freesurfer/FILE/DK-318/DK318_MNI_Coordinates.csv"
roi_name_file  <- "/data/home/tqi/data1/share/after_freesurfer/FILE/DK-318/DK318_roi_names.csv"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
tab_dir <- file.path(out_dir, "tables")
brainnet_dir <- file.path(out_dir, "brainnet_files")
dir.create(tab_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(brainnet_dir, showWarnings = FALSE, recursive = TRUE)

# ----------------------------
# NBR analysis parameters
# ----------------------------
analysis_mode <- Sys.getenv("NBR_ANALYSIS_MODE", unset = "debug")
analysis_mode <- trimws(tolower(analysis_mode))

if (!(analysis_mode %in% c("debug", "formal"))) {
  stop("Invalid NBR_ANALYSIS_MODE: ", analysis_mode, ". Use 'debug' or 'formal'.")
}

nperm_debug <- 20L
nperm_formal <- 1000L

nperm_nbr <- if (analysis_mode == "debug") nperm_debug else nperm_formal
#这里是单向的p值，但是我们要算整体，所以应该除2
thrP_nbs <- 0.005
component_p_alpha <- 0.005
component_p_metric <- "ncompFWE"

overall_model <- "~ D + A + S + Age_c"
interaction_model <- "~ D + A + S + DA + DS + AS + DAS + Age_c"

run_overall_nbr <- get_env_flag("NBR_RUN_OVERALL", default = TRUE)
run_interaction_nbr <- get_env_flag("NBR_RUN_INTERACTION", default = TRUE)
run_age_sex_simple_effect_nbr <- get_env_flag("NBR_RUN_AGE_SEX_SIMPLE", default = TRUE)
run_age_only_simple_effect_nbr <- get_env_flag("NBR_RUN_AGE_ONLY_SIMPLE", default = TRUE)


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

flatten_list_columns <- function(df) {
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  df[] <- lapply(df, function(col) {
    if (!is.list(col)) {
      if (is.factor(col)) return(as.character(col))
      return(col)
    }

    vapply(col, function(x) {
      if (length(x) == 0L || all(is.na(x))) return(NA_character_)
      paste(as.character(x), collapse = ";")
    }, character(1))
  })
  df
}

save_csv <- function(df, file, row.names = FALSE) {
  utils::write.csv(flatten_list_columns(df), file = file, row.names = row.names)
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
  # 基于Yeo7网络映射的DK318脑区分类
  # Based on Yeo7 network mapping from CSV files
  
  # 读取Yeo7映射文件
  lh_mapping <- read.csv("/data/home/tqi/data1/share/after_freesurfer/FILE/DK-318/lh_DK318_to_Yeo7_mapping.csv", 
                         stringsAsFactors = FALSE)
  rh_mapping <- read.csv("/data/home/tqi/data1/share/after_freesurfer/FILE/DK-318/rh_DK318_to_Yeo7_mapping.csv", 
                         stringsAsFactors = FALSE)
  
  # 添加半球前缀
  lh_mapping$full_name <- paste0("lh_", lh_mapping$DK318_region)
  rh_mapping$full_name <- paste0("rh_", rh_mapping$DK318_region)
  
  # 合并左右半球
  all_mapping <- rbind(lh_mapping, rh_mapping)
  
  # 创建网络名称映射（Yeo7标准名称到缩写）
  network_abbrev <- c(
    "Visual" = "VIS",
    "Somatomotor" = "SMN",
    "Dorsal Attention" = "DAN",
    "Ventral Attention" = "VAN",
    "Salience/Ventral Attention" = "VAN",
    "Limbic" = "Limbic",
    "Frontoparietal" = "FPN",
    "Control" = "FPN",
    "Default" = "DMN",
    "Unknown" = "Other"
  )
  
  # 为每个ROI查找对应的网络
  out <- sapply(roi, function(r) {
    # 查找匹配的脑区
    idx <- which(all_mapping$full_name == r)
    if (length(idx) > 0) {
      network <- all_mapping$primary_Yeo7_network[idx[1]]
      abbrev <- network_abbrev[network]
      if (is.na(abbrev)) return("Other")
      return(as.character(abbrev))
    } else {
      return("Other")
    }
  })
  
  return(out)
}

get_network_palette <- function() {
  # 按照Yeo7标准顺序
  # 1=Visual, 2=Somatomotor, 3=Dorsal Attention, 4=Ventral Attention,
  # 5=Limbic, 6=Frontoparietal, 7=Default
  c(
    VIS = 1,
    SMN = 2,
    DAN = 3,
    VAN = 4,
    Limbic = 5,
    FPN = 6,      # Frontoparietal = 6 (按Yeo7标准)
    DMN = 7,      # Default Mode = 7 (按Yeo7标准)
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
    save_csv(roi_name_df, file.path(tab_dir, "DK318_roi_names_from_server.csv"), row.names = FALSE)

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
      )

    if (!("n_vertices" %in% names(coord))) {
      coord$n_vertices <- NA_real_
    }

    coord <- coord %>%
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
  roi_meta <- roi_meta %>%
    mutate(
      ROI = as.character(ROI),
      lobe = as.character(lobe),
      dk_region = as.character(dk_region),
      network = as.character(network)
    )

  make_pair <- function(x, y) {
    x <- as.character(x)
    y <- as.character(y)
    ifelse(x <= y, paste(x, y, sep = "--"), paste(y, x, sep = "--"))
  }

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
    rename(i = ROI_index)
  roi2_meta <- roi_meta %>%
    rename_with(~ paste0(.x, "2"), -ROI_index) %>%
    rename(j = ROI_index)

  edge_annot <- edge_annot %>%
    left_join(roi1_meta, by = "i") %>%
    left_join(roi2_meta, by = "j") %>%
    mutate(
      lobe_pair = make_pair(lobe1, lobe2),
      dk_pair = make_pair(dk_region1, dk_region2),
      network_pair = make_pair(network1, network2)
    )

  edge_annot
}

resolve_edge_roi_columns <- function(edges) {
  out <- edges

  if (!("ROI1" %in% names(out))) {
    for (candidate in c("ROI1.x", "ROI1.y")) {
      if (candidate %in% names(out)) {
        out$ROI1 <- out[[candidate]]
        break
      }
    }
  }

  if (!("ROI2" %in% names(out))) {
    for (candidate in c("ROI2.x", "ROI2.y")) {
      if (candidate %in% names(out)) {
        out$ROI2 <- out[[candidate]]
        break
      }
    }
  }

  out
}

check_nbr_supports_alternative <- function() {
  fml <- names(formals(NBR::nbr_lm))
  if (!("alternative" %in% fml)) {
    stop(
      "Current NBR::nbr_lm does not expose an 'alternative' argument. ",
      "Please update NBR from CRAN. Directional NBS requires alternative = 'greater' or 'lower'."
    )
  }
}

compute_edge_stats_lm <- function(net, idata, mod, term_name) {
  X <- model.matrix(as.formula(mod), data = idata)

  if (!(term_name %in% colnames(X))) {
    stop("Term ", term_name, " is not found in model matrix columns.")
  }

  if (qr(X)$rank < ncol(X)) {
    stop("Design matrix is rank deficient for model: ", mod)
  }

  XtX_inv <- solve(crossprod(X))
  B <- XtX_inv %*% crossprod(X, net)

  resid <- net - X %*% B
  df_resid <- nrow(X) - qr(X)$rank
  sigma2 <- colSums(resid * resid, na.rm = TRUE) / df_resid

  term_pos <- which(colnames(X) == term_name)
  se <- sqrt(sigma2 * XtX_inv[term_pos, term_pos])

  estimate <- as.numeric(B[term_pos, ])
  t_value <- estimate / se
  p_value <- 2 * pt(abs(t_value), df = df_resid, lower.tail = FALSE)

  data.frame(
    edge_index = seq_len(ncol(net)),
    term = term_name,
    estimate = estimate,
    SE = se,
    df = df_resid,
    t = t_value,
    p = p_value,
    p_fdr = p.adjust(p_value, method = "fdr"),
    stringsAsFactors = FALSE
  )
}

run_nbr_model <- function(net, nnodes, idata, mod, result_prefix,
                          alternative, nperm, thrP, cores = NULL) {
  cat("\n==============================\n")
  cat("Running NBR:", result_prefix, "\n")
  cat("Model:", mod, "\n")
  cat("Alternative:", alternative, "\n")
  cat("Permutations:", nperm, "\n")
  cat("thrP:", thrP, "\n")
  cat("==============================\n")

  res <- NBR::nbr_lm(
    net = net,
    nnodes = nnodes,
    idata = idata,
    mod = mod,
    alternative = alternative,
    diag = FALSE,
    nperm = nperm,
    thrP = thrP,
    thrT = NULL,
    cores = cores,
    nudist = TRUE,
    verbose = TRUE,
    na.action = na.exclude
  )

  saveRDS(res, file = file.path(out_dir, paste0(result_prefix, "_nbr_result.rds")))
  save_text_summary(res, file = file.path(out_dir, paste0(result_prefix, "_nbr_summary.txt")))

  invisible(res)
}

extract_nbr_term <- function(nbr_res, term_name, edge_annot, edge_stats,
                             analysis_name,
                             p_metric = "ncompFWE",
                             alpha = 0.05) {
  if (is.null(nbr_res$components[[term_name]]) ||
      is.null(nbr_res$fwe[[term_name]])) {
    return(list(
      summary = data.frame(),
      edges = data.frame(),
      sig_edges = data.frame(),
      result_for_plot = list(contrast = analysis_name, components = data.frame())
    ))
  }

  comp <- as.data.frame(nbr_res$components[[term_name]])
  if (nrow(comp) == 0L) {
    return(list(
      summary = data.frame(),
      edges = data.frame(),
      sig_edges = data.frame(),
      result_for_plot = list(contrast = analysis_name, components = data.frame())
    ))
  }

  if (ncol(comp) < 5L) {
    stop("Unexpected NBR component table format for term: ", term_name)
  }

  names(comp)[1:5] <- c("edge_index", "i", "j", "Component", "strn")

  comp <- comp %>%
    mutate(
      edge_index = as.integer(edge_index),
      i = as.integer(i),
      j = as.integer(j),
      Component = as.integer(Component)
    )

  fwe <- as.data.frame(nbr_res$fwe[[term_name]])
  if (nrow(fwe) > 0L) {
    fwe <- fwe %>% mutate(Component = as.integer(Component))
  }

  edges <- comp %>%
    left_join(fwe, by = "Component") %>%
    left_join(edge_stats, by = "edge_index") %>%
    left_join(edge_annot, by = c("edge_index", "i", "j")) %>%
    mutate(
      analysis = analysis_name,
      term = term_name,
      component_sig = ifelse(!is.na(.data[[p_metric]]), .data[[p_metric]] < alpha, FALSE)
    )

  sig_edges <- edges %>%
    filter(!is.na(.data[[p_metric]]) & .data[[p_metric]] < alpha)

  summary <- fwe %>%
    mutate(
      analysis = analysis_name,
      term = term_name,
      component_sig = ifelse(!is.na(.data[[p_metric]]), .data[[p_metric]] < alpha, FALSE)
    )

  result_for_plot <- list(
    contrast = analysis_name,
    components = edges,
    fwe = summary
  )

  list(
    summary = summary,
    edges = edges,
    sig_edges = sig_edges,
    result_for_plot = result_for_plot
  )
}

add_direction_labels <- function(edges, analysis_type) {
  if (nrow(edges) == 0L) return(edges)

  edges %>%
    mutate(
      direction = case_when(
        analysis_type == "diagnosis" & t > 0 ~ "TD > DD",
        analysis_type == "diagnosis" & t < 0 ~ "DD > TD",
        analysis_type == "interaction" & t > 0 ~ "positive interaction contrast",
        analysis_type == "interaction" & t < 0 ~ "negative interaction contrast",
        TRUE ~ "zero"
      )
    )
}

export_nbr_null_distribution <- function(nbr_res, term_name, analysis_name) {
  if (is.null(nbr_res$nudist)) return(invisible(NULL))

  nd <- nbr_res$nudist
  nd_df <- tryCatch(as.data.frame(nd), error = function(e) data.frame())
  if (nrow(nd_df) == 0L) return(invisible(NULL))

  safe <- sanitize_name(paste0(analysis_name, "_", term_name))
  save_csv(nd_df, file.path(tab_dir, paste0(safe, "_null_distribution_raw.csv")), row.names = FALSE)

  invisible(nd_df)
}

# ----------------------------
# 3) Table exports
# ----------------------------
get_significant_edges <- function(comp_edges, p_metric = "ncompFWE", fwe_alpha = 0.05) {
  if (is.null(comp_edges) || nrow(comp_edges) == 0L) return(data.frame())
  comp_edges %>% filter(.data[[p_metric]] < fwe_alpha)
}

make_node_degree_table <- function(edges, roi_meta) {
  if (nrow(edges) == 0L) return(data.frame())

  edges <- resolve_edge_roi_columns(edges)

  required_cols <- c("i", "j", "ROI1", "ROI2")
  missing_cols <- setdiff(required_cols, names(edges))
  if (length(missing_cols) > 0L) {
    stop(
      "Cannot create node degree table because required columns are missing: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  d1 <- edges %>% transmute(ROI_index = i, ROI = ROI1)
  d2 <- edges %>% transmute(ROI_index = j, ROI = ROI2)

  bind_rows(d1, d2) %>%
    count(ROI_index, ROI, name = "degree") %>%
    left_join(roi_meta, by = c("ROI_index", "ROI")) %>%
    arrange(desc(degree), ROI_index)
}

write_analysis_tables <- function(result,
                                  roi_meta,
                                  p_metric = "ncompFWE",
                                  fwe_alpha = 0.05) {
  effect_name <- result$contrast
  comp_edges <- result$components
  if (is.null(comp_edges) || nrow(comp_edges) == 0L) {
    cat("No suprathreshold NBS component for ", effect_name, "\n")
    return(invisible(NULL))
  }

  safe <- sanitize_name(effect_name)

  save_csv(
    comp_edges,
    file.path(tab_dir, paste0(safe, "_suprathreshold_component_edges.csv")),
    row.names = FALSE
  )

  if (!is.null(result$fwe) && nrow(result$fwe) > 0L) {
    save_csv(
      result$fwe,
      file.path(tab_dir, paste0(safe, "_component_summary.csv")),
      row.names = FALSE
    )
  }

  sig_edges <- get_significant_edges(comp_edges, p_metric = p_metric, fwe_alpha = fwe_alpha)
  save_csv(
    sig_edges,
    file.path(tab_dir, paste0(safe, "_edges_in_FWE_significant_components.csv")),
    row.names = FALSE
  )

  if (nrow(sig_edges) == 0L) {
    cat("No FWE-significant component for ", effect_name, "\n")
    return(invisible(NULL))
  }

  degree <- make_node_degree_table(sig_edges, roi_meta)
  save_csv(
    degree,
    file.path(tab_dir, paste0(safe, "_node_degree.csv")),
    row.names = FALSE
  )

  write_brainnet_files(
    edges_sig = sig_edges,
    roi_meta = roi_meta,
    effect_name = effect_name
  )

  invisible(list(sig_edges = sig_edges, degree = degree))
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
    mutate(
      contrast = if ("contrast" %in% names(.)) as.character(contrast) else effect_name,
      AgeGroup = if ("AgeGroup" %in% names(.)) as.character(AgeGroup) else NA_character_,
      Sex = if ("Sex" %in% names(.)) as.character(Sex) else NA_character_
    ) %>%
    distinct(Component, ncomp, strn, ncompFWE, strnFWE, contrast, AgeGroup, Sex, component_sig)
  save_csv(
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
    edge_mat_dd_gt_td <- matrix(0, nrow = nrow(nodes), ncol = nrow(nodes))

    for (kk in seq_len(nrow(e))) {
      a <- idx_map[as.character(e$i[kk])]
      b <- idx_map[as.character(e$j[kk])]
      edge_weight <- abs(e$t[kk])
      edge_mat[a, b] <- edge_weight
      edge_mat[b, a] <- edge_weight

      if (e$t[kk] > 0) {
        edge_mat_td_gt_dd[a, b] <- edge_weight
        edge_mat_td_gt_dd[b, a] <- edge_weight
      }
      if (e$t[kk] < 0) {
        edge_mat_dd_gt_td[a, b] <- edge_weight
        edge_mat_dd_gt_td[b, a] <- edge_weight
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
          t > 0 ~ "TD > DD",
          t < 0 ~ "DD > TD",
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
      edge_mat_dd_gt_td,
      file = file.path(brainnet_dir, paste0(safe, "_component_", cc, "_DD_gt_TD.edge")),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    save_csv(
      node_annotation,
      file = file.path(brainnet_dir, paste0(safe, "_component_", cc, "_node_annotation.csv")),
      row.names = FALSE
    )
    save_csv(
      edge_annotation,
      file = file.path(brainnet_dir, paste0(safe, "_component_", cc, "_edge_annotation.csv")),
      row.names = FALSE
    )
  }

  save_csv(
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

# ----------------------------
# Effect coding and continuous age covariate
# ----------------------------
df <- df %>%
  mutate(
    Diagnosis = factor(Diagnosis, levels = c("TD", "DD")),
    AgeGroup = factor(AgeGroup, levels = c("Child", "Adult")),
    Sex = factor(Sex, levels = c("Female", "Male")),
    Age_c = Age - mean(Age, na.rm = TRUE),
    D = ifelse(Diagnosis == "TD", 0.5, -0.5),
    A = ifelse(AgeGroup == "Adult", 0.5, -0.5),
    S = ifelse(Sex == "Male", 0.5, -0.5),
    DA = D * A,
    DS = D * S,
    AS = A * S,
    DAS = D * A * S
  )

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
save_csv(cell_counts, file.path(tab_dir, "Diagnosis_AgeGroup_Sex_cell_counts.csv"), row.names = FALSE)
if (any(cell_counts$n == 0L)) {
  warning("At least one Diagnosis x AgeGroup x Sex cell has n = 0. Interaction model may be rank deficient.")
}

save_csv(df, file.path(out_dir, "included_subjects_and_covariates.csv"), row.names = FALSE)

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

save_csv(roi_meta, file.path(out_dir, "roi_metadata_with_mni_coordinates.csv"), row.names = FALSE)
save_csv(edge_annot, file.path(out_dir, "edge_annotation_roi_mni_lobe.csv"), row.names = FALSE)
saveRDS(net2d, file.path(out_dir, "dk318_mind_upper_triangle_matrix.rds"))

# ----------------------------
# 7) NBR models: overall Diagnosis effect and Diagnosis-related interactions
# ----------------------------
check_nbr_supports_alternative()

idata_nbr <- df %>%
  select(
    subj_prefix,
    Diagnosis, AgeGroup, Sex, Age, Age_c,
    D, A, S, DA, DS, AS, DAS
  )

saveRDS(idata_nbr, file.path(out_dir, "idata_nbr_effect_coded.rds"))
save_csv(idata_nbr, file.path(out_dir, "idata_nbr_effect_coded.csv"), row.names = FALSE)

overall_results <- list()
interaction_results <- list()
simple_results <- list()
nbr_result_objects <- list()

# ----------------------------
# 7.1 Overall Diagnosis effect
# ----------------------------
if (run_overall_nbr) {
  res_overall_DD_gt_TD <- run_nbr_model(
    net = net2d,
    nnodes = nnodes,
    idata = idata_nbr,
    mod = overall_model,
    result_prefix = "NBR_overall_DD_gt_TD",
    alternative = "lower",
    nperm = nperm_nbr,
    thrP = thrP_nbs,
    cores = cores_to_use
  )

  res_overall_TD_gt_DD <- run_nbr_model(
    net = net2d,
    nnodes = nnodes,
    idata = idata_nbr,
    mod = overall_model,
    result_prefix = "NBR_overall_TD_gt_DD",
    alternative = "greater",
    nperm = nperm_nbr,
    thrP = thrP_nbs,
    cores = cores_to_use
  )

  nbr_result_objects[["overall_DD_gt_TD_D"]] <- list(res = res_overall_DD_gt_TD, term = "D", analysis = "overall_DD_gt_TD")
  nbr_result_objects[["overall_TD_gt_DD_D"]] <- list(res = res_overall_TD_gt_DD, term = "D", analysis = "overall_TD_gt_DD")

  edge_stats_overall <- compute_edge_stats_lm(
    net = net2d,
    idata = idata_nbr,
    mod = overall_model,
    term_name = "D"
  )

  overall_results[["overall_DD_gt_TD"]] <- extract_nbr_term(
    nbr_res = res_overall_DD_gt_TD,
    term_name = "D",
    edge_annot = edge_annot,
    edge_stats = edge_stats_overall,
    analysis_name = "overall_DD_gt_TD",
    p_metric = component_p_metric,
    alpha = component_p_alpha
  )

  overall_results[["overall_TD_gt_DD"]] <- extract_nbr_term(
    nbr_res = res_overall_TD_gt_DD,
    term_name = "D",
    edge_annot = edge_annot,
    edge_stats = edge_stats_overall,
    analysis_name = "overall_TD_gt_DD",
    p_metric = component_p_metric,
    alpha = component_p_alpha
  )
}

# ----------------------------
# 7.2 Diagnosis-related interaction effects
# ----------------------------
if (run_interaction_nbr) {
  res_interaction_pos <- run_nbr_model(
    net = net2d,
    nnodes = nnodes,
    idata = idata_nbr,
    mod = interaction_model,
    result_prefix = "NBR_interactions_positive_direction",
    alternative = "greater",
    nperm = nperm_nbr,
    thrP = thrP_nbs,
    cores = cores_to_use
  )

  res_interaction_neg <- run_nbr_model(
    net = net2d,
    nnodes = nnodes,
    idata = idata_nbr,
    mod = interaction_model,
    result_prefix = "NBR_interactions_negative_direction",
    alternative = "lower",
    nperm = nperm_nbr,
    thrP = thrP_nbs,
    cores = cores_to_use
  )

  interaction_terms <- c("DA", "DS", "DAS")

  for (tt in interaction_terms) {
    nbr_result_objects[[paste0("interaction_", tt, "_positive")]] <- list(
      res = res_interaction_pos,
      term = tt,
      analysis = paste0("interaction_", tt, "_positive")
    )
    nbr_result_objects[[paste0("interaction_", tt, "_negative")]] <- list(
      res = res_interaction_neg,
      term = tt,
      analysis = paste0("interaction_", tt, "_negative")
    )

    edge_stats_tt <- compute_edge_stats_lm(
      net = net2d,
      idata = idata_nbr,
      mod = interaction_model,
      term_name = tt
    )

    analysis_pos <- paste0("interaction_", tt, "_positive")
    analysis_neg <- paste0("interaction_", tt, "_negative")

    interaction_results[[analysis_pos]] <- extract_nbr_term(
      nbr_res = res_interaction_pos,
      term_name = tt,
      edge_annot = edge_annot,
      edge_stats = edge_stats_tt,
      analysis_name = analysis_pos,
      p_metric = component_p_metric,
      alpha = component_p_alpha
    )

    interaction_results[[analysis_neg]] <- extract_nbr_term(
      nbr_res = res_interaction_neg,
      term_name = tt,
      edge_annot = edge_annot,
      edge_stats = edge_stats_tt,
      analysis_name = analysis_neg,
      p_metric = component_p_metric,
      alpha = component_p_alpha
    )
  }
}

# ----------------------------
# 8) Planned simple-effect NBR contrasts
# ----------------------------
run_simple_effect_nbr <- function(age_group, direction,
                                  net2d, df_nbr, nnodes,
                                  nperm, thrP, cores = NULL,
                                  sex_value = NULL) {
  keep <- df_nbr$AgeGroup == age_group
  label_parts <- c(age_group)

  if (!is.null(sex_value)) {
    keep <- keep & df_nbr$Sex == sex_value
    label_parts <- c(label_parts, sex_value)
  }

  dat <- df_nbr[keep, , drop = FALSE]
  net_sub <- net2d[keep, , drop = FALSE]

  subgroup_label <- paste(label_parts, collapse = " ")

  if (nrow(dat) == 0L) {
    stop("No subjects found for ", subgroup_label)
  }

  if (length(unique(dat$Diagnosis)) < 2L) {
    stop("Only one Diagnosis level exists in ", subgroup_label)
  }

  dat <- dat %>%
    mutate(
      Diagnosis = factor(Diagnosis, levels = c("TD", "DD")),
      Sex = factor(Sex, levels = c("Female", "Male")),
      D = ifelse(Diagnosis == "TD", 0.5, -0.5),
      S = ifelse(Sex == "Male", 0.5, -0.5),
      DS = D * S
    )

  age_sd <- sd(dat$Age_c, na.rm = TRUE)
  age_has_var <- is.finite(age_sd) && age_sd > 0
  sex_has_two_levels <- length(unique(stats::na.omit(as.character(dat$Sex)))) >= 2L

  if (is.null(sex_value)) {
    if (sex_has_two_levels && age_has_var) {
      mod <- "~ D + S + DS + Age_c"
    } else if (sex_has_two_levels && !age_has_var) {
      mod <- "~ D + S + DS"
    } else if (!sex_has_two_levels && age_has_var) {
      mod <- "~ D + Age_c"
    } else {
      mod <- "~ D"
    }
  } else {
    if (age_has_var) {
      mod <- "~ D + Age_c"
    } else {
      mod <- "~ D"
    }
  }

  alt <- if (direction == "TD_gt_DD") {
    "greater"
  } else if (direction == "DD_gt_TD") {
    "lower"
  } else {
    stop("direction must be TD_gt_DD or DD_gt_TD")
  }

  prefix_core <- paste(label_parts, collapse = "_")
  prefix <- paste0("NBR_", prefix_core, "_", direction)

  res <- run_nbr_model(
    net = net_sub,
    nnodes = nnodes,
    idata = dat,
    mod = mod,
    result_prefix = prefix,
    alternative = alt,
    nperm = nperm,
    thrP = thrP,
    cores = cores
  )

  edge_stats <- compute_edge_stats_lm(
    net = net_sub,
    idata = dat,
    mod = mod,
    term_name = "D"
  )

  analysis_name <- paste0(prefix_core, "_", direction)

  out <- extract_nbr_term(
    nbr_res = res,
    term_name = "D",
    edge_annot = edge_annot,
    edge_stats = edge_stats,
    analysis_name = analysis_name,
    p_metric = component_p_metric,
    alpha = component_p_alpha
  )

  list(
    result = res,
    extracted = out,
    data = dat,
    net = net_sub,
    model = mod,
    age_group = age_group,
    sex = sex_value,
    direction = direction,
    analysis_name = analysis_name,
    result_prefix = prefix
  )
}

age_only_posthoc_plan <- expand.grid(
  AgeGroup = c("Child", "Adult"),
  direction = c("DD_gt_TD", "TD_gt_DD"),
  stringsAsFactors = FALSE
) %>%
  arrange(AgeGroup, direction)

age_sex_posthoc_plan <- expand.grid(
  AgeGroup = c("Child", "Adult"),
  Sex = c("Female", "Male"),
  direction = c("DD_gt_TD", "TD_gt_DD"),
  stringsAsFactors = FALSE
) %>%
  arrange(AgeGroup, Sex, direction)

save_csv(age_only_posthoc_plan, file.path(tab_dir, "planned_simple_effects_age_only.csv"), row.names = FALSE)
save_csv(age_sex_posthoc_plan, file.path(tab_dir, "planned_simple_effects_age_sex.csv"), row.names = FALSE)

if (run_age_only_simple_effect_nbr) {
  for (kk in seq_len(nrow(age_only_posthoc_plan))) {
    nm <- paste(
      age_only_posthoc_plan$AgeGroup[kk],
      age_only_posthoc_plan$direction[kk],
      sep = "_"
    )

    simple_results[[nm]] <- run_simple_effect_nbr(
      age_group = age_only_posthoc_plan$AgeGroup[kk],
      direction = age_only_posthoc_plan$direction[kk],
      net2d = net2d,
      df_nbr = idata_nbr,
      nnodes = nnodes,
      nperm = nperm_nbr,
      thrP = thrP_nbs,
      cores = cores_to_use,
      sex_value = NULL
    )

    nbr_result_objects[[paste0("simple_", nm, "_D")]] <- list(
      res = simple_results[[nm]]$result,
      term = "D",
      analysis = nm
    )

    saveRDS(simple_results[[nm]], file.path(out_dir, paste0(simple_results[[nm]]$result_prefix, "_full_output.rds")))
  }
}

if (run_age_sex_simple_effect_nbr) {
  for (kk in seq_len(nrow(age_sex_posthoc_plan))) {
    nm <- paste(
      age_sex_posthoc_plan$AgeGroup[kk],
      age_sex_posthoc_plan$Sex[kk],
      age_sex_posthoc_plan$direction[kk],
      sep = "_"
    )

    simple_results[[nm]] <- run_simple_effect_nbr(
      age_group = age_sex_posthoc_plan$AgeGroup[kk],
      sex_value = age_sex_posthoc_plan$Sex[kk],
      direction = age_sex_posthoc_plan$direction[kk],
      net2d = net2d,
      df_nbr = idata_nbr,
      nnodes = nnodes,
      nperm = nperm_nbr,
      thrP = thrP_nbs,
      cores = cores_to_use
    )

    nbr_result_objects[[paste0("simple_", nm, "_D")]] <- list(
      res = simple_results[[nm]]$result,
      term = "D",
      analysis = nm
    )

    saveRDS(simple_results[[nm]], file.path(out_dir, paste0(simple_results[[nm]]$result_prefix, "_full_output.rds")))
  }
}

# ----------------------------
# 9) Export NBR results to tables
# ----------------------------
all_extracted <- list()

if (length(overall_results) > 0L) {
  all_extracted <- c(all_extracted, overall_results)
}

if (length(interaction_results) > 0L) {
  all_extracted <- c(all_extracted, interaction_results)
}

if (length(simple_results) > 0L) {
  simple_extracted <- lapply(simple_results, function(x) x$extracted)
  all_extracted <- c(all_extracted, simple_extracted)
}

if (length(all_extracted) > 0L) {
  all_component_summary <- bind_rows(lapply(all_extracted, function(x) x$summary))
  all_component_edges <- bind_rows(lapply(all_extracted, function(x) x$edges))
  all_sig_edges <- bind_rows(lapply(all_extracted, function(x) x$sig_edges))
} else {
  all_component_summary <- data.frame()
  all_component_edges <- data.frame()
  all_sig_edges <- data.frame()
}

save_csv(
  all_component_summary,
  file.path(tab_dir, "NBR_ALL_component_summary.csv"),
  row.names = FALSE
)

save_csv(
  all_component_edges,
  file.path(tab_dir, "NBR_ALL_suprathreshold_component_edges.csv"),
  row.names = FALSE
)

save_csv(
  all_sig_edges,
  file.path(tab_dir, "NBR_ALL_edges_in_FWE_significant_components.csv"),
  row.names = FALSE
)

for (obj in nbr_result_objects) {
  export_nbr_null_distribution(
    nbr_res = obj$res,
    term_name = obj$term,
    analysis_name = obj$analysis
  )
}

for (nm in names(all_extracted)) {
  x <- all_extracted[[nm]]

  if (is.null(x$result_for_plot) ||
      is.null(x$result_for_plot$components) ||
      nrow(x$result_for_plot$components) == 0L) {
    next
  }

  analysis_label <- x$result_for_plot$contrast
  analysis_type <- if (grepl("^interaction_", analysis_label)) "interaction" else "diagnosis"
  x$result_for_plot$components <- add_direction_labels(
    x$result_for_plot$components,
    analysis_type = analysis_type
  )

  tryCatch(
    {
      write_analysis_tables(
        result = x$result_for_plot,
        roi_meta = roi_meta,
        p_metric = component_p_metric,
        fwe_alpha = component_p_alpha
      )
    },
    error = function(e) {
      warning(
        "Failed to export per-analysis tables for ", analysis_label,
        ": ", conditionMessage(e)
      )
    }
  )
}

# ----------------------------
# 10) Session info
# ----------------------------
analysis_settings <- data.frame(
  setting = c(
    "demo_file", "mind_dir", "out_dir", "roi_coord_file", "roi_name_file",
    "analysis_mode", "overall_model", "interaction_model",
    "thrP_nbs", "nperm_nbr", "component_p_metric", "component_p_alpha",
    "run_overall_nbr", "run_interaction_nbr", "run_age_only_simple_effect_nbr", "run_age_sex_simple_effect_nbr",
    "TIV_covariate", "site_covariate"
  ),
  value = c(
    demo_file, mind_dir, out_dir, roi_coord_file, roi_name_file,
    analysis_mode, overall_model, interaction_model,
    as.character(thrP_nbs), as.character(nperm_nbr),
    component_p_metric, as.character(component_p_alpha),
    as.character(run_overall_nbr), as.character(run_interaction_nbr),
    as.character(run_age_only_simple_effect_nbr), as.character(run_age_sex_simple_effect_nbr),
    "not included; unavailable",
    "not included; handled by upstream ComBat"
  )
)
save_csv(analysis_settings, file.path(out_dir, "analysis_settings.csv"), row.names = FALSE)

session_info_file <- file.path(out_dir, "sessionInfo.txt")
writeLines(capture.output(sessionInfo()), con = session_info_file)

cat("\nAll analyses finished. Results saved to:\n", out_dir, "\n")
cat("Tables saved to:\n", tab_dir, "\n")
cat("BrainNet files saved to:\n", brainnet_dir, "\n")
