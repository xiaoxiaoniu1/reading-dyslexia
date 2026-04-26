if(!requireNamespace("pacman", quietly = TRUE)) {
  install.packages("pacman", repos = "https://cloud.r-project.org/")
}
pacman::p_load(readxl, dplyr, ggplot2)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x
}

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(key, default = NULL) {
  i <- match(key, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}

get_flag <- function(key, default = TRUE) {
  if (!(key %in% args)) return(default)
  i <- match(key, args)
  if (i == length(args)) return(TRUE)
  val <- tolower(args[i + 1])
  if (val %in% c("false", "0", "no", "n")) FALSE else TRUE
}

model_arg <- tolower(get_arg("--model", "default"))
sex_arg <- tolower(get_arg("--sex", "all"))
result_dir_arg <- get_arg("--result-dir", NULL)

include_network <- get_flag("--include-network", TRUE)
p_threshold <- as.numeric(get_arg("--p-threshold", "0.05"))

result_base_dir <- "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5"
default_roots <- list(
  default = file.path(result_base_dir, "MIND_DK318_DGLM"),
  main = file.path(result_base_dir, "MIND_DK318_DGLM_main"),
  sensitivity = file.path(result_base_dir, "MIND_DK318_DGLM_sensitivity"),
  main_sex = file.path(result_base_dir, "MIND_DK318_DGLM_main_sex"),
  sensitivity_sex = file.path(result_base_dir, "MIND_DK318_DGLM_sensitivity_sex")
)

demo_file <- get_arg(
  "--demo-file",
  file.path(result_base_dir, "all_data_cqt_mean_1.5.xlsx")
)
mind_dir <- get_arg(
  "--mind-combat-dir",
  file.path(result_base_dir, "MIND_DK318_combat")
)

degree_csv_name <- "DGLM_DK318_degree_results.csv"
network_csv_name <- "DGLM_DK318_degree_network.csv"

interaction_effects <- c(
  "Diagnosis_AgeGroup",
  "Diagnosis_Sex",
  "AgeGroup_Sex",
  "Diagnosis_AgeGroup_Sex"
)
domains <- c("mean", "disp")

safe <- function(x) {
  x <- gsub("[^A-Za-z0-9_]+", "_", as.character(x))
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

resolve_degree_dir <- function(path) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  if (file.exists(file.path(path, degree_csv_name))) return(path)
  nested <- file.path(path, "degree")
  if (file.exists(file.path(nested, degree_csv_name))) return(nested)
  NULL
}

find_network_dirs <- function(degree_dir) {
  net_root <- file.path(degree_dir, "network_FDR")
  if (!dir.exists(net_root)) return(character(0))
  out <- list.dirs(net_root, recursive = FALSE, full.names = TRUE)
  out[file.exists(file.path(out, network_csv_name))]
}

sex_dirs_for_root <- function(root, sex = "all") {
  labels <- if (sex == "all") c("male", "female") else sex
  out <- character(0)

  for (s in labels) {
    d <- file.path(root, paste0(s, "_sex"), "degree")
    if (file.exists(file.path(d, degree_csv_name))) out <- c(out, d)

    d2 <- file.path(root, paste0(s, "_sex"))
    if (file.exists(file.path(d2, degree_csv_name))) out <- c(out, d2)
  }

  unique(out)
}

discover_degree_dirs <- function(model = "default", sex = "all", result_dir = NULL) {
  if (!is.null(result_dir) && !is.na(result_dir)) {
    d <- resolve_degree_dir(result_dir)
    if (!is.null(d)) return(d)

    out <- sex_dirs_for_root(result_dir, sex)
    if (length(out)) return(out)

    stop("Cannot find ", degree_csv_name, " under: ", result_dir)
  }

  models <- if (model == "all") names(default_roots) else model
  out <- character(0)

  for (m in models) {
    if (!(m %in% names(default_roots))) next
    root <- default_roots[[m]]

    if (grepl("_sex$", m)) {
      out <- c(out, sex_dirs_for_root(root, sex))
    } else {
      d <- resolve_degree_dir(root)
      if (!is.null(d)) out <- c(out, d)
    }
  }

  unique(out)
}

parse_sex_filter <- function(active_result_dir) {
  low <- tolower(normalizePath(active_result_dir, winslash = "/", mustWork = FALSE))
  if (grepl("male_sex", low)) return("Male")
  if (grepl("female_sex", low)) return("Female")
  NA_character_
}

include_age_for_dir <- function(active_result_dir) {
  low <- tolower(normalizePath(active_result_dir, winslash = "/", mustWork = FALSE))
  if (grepl("dglm_main", low) && !grepl("sensitivity", low)) return(FALSE)
  TRUE
}

parse_p_col <- function(col) {
  if (!startsWith(col, "p_")) return(NULL)
  rest <- sub("^p_", "", col)

  for (domain in domains) {
    prefix <- paste0(domain, "_")
    if (!startsWith(rest, prefix)) next
    tail <- substring(rest, nchar(prefix) + 1)

    for (eff in interaction_effects[order(nchar(interaction_effects), decreasing = TRUE)]) {
      if (tail == eff) {
        return(list(domain = domain, effect = eff, correction = "uncorrected"))
      }
      if (tail == paste0(eff, "_FDR")) {
        return(list(domain = domain, effect = eff, correction = "FDR"))
      }
      net_prefix <- paste0(eff, "_networkFDR_")
      if (startsWith(tail, net_prefix)) {
        return(list(domain = domain, effect = eff, correction = "networkFDR"))
      }
    }
  }

  NULL
}

build_tasks_from_table <- function(tab) {
  rows <- list()

  for (col in names(tab)) {
    parsed <- parse_p_col(col)
    if (is.null(parsed)) next

    effect <- parsed$effect
    domain <- parsed$domain
    correction <- parsed$correction

    if (effect == "Diagnosis_AgeGroup") {
      x_var <- "Diagnosis_Fac"
      facet_row <- NA_character_
      facet_col <- "AgeGroup_Fac"
    } else if (effect == "Diagnosis_Sex") {
      x_var <- "Diagnosis_Fac"
      facet_row <- NA_character_
      facet_col <- "Sex"
    } else if (effect == "AgeGroup_Sex") {
      x_var <- "AgeGroup_Fac"
      facet_row <- NA_character_
      facet_col <- "Sex"
    } else if (effect == "Diagnosis_AgeGroup_Sex") {
      x_var <- "Diagnosis_Fac"
      facet_row <- "Sex"
      facet_col <- "AgeGroup_Fac"
    } else {
      next
    }

    rows[[length(rows) + 1]] <- data.frame(
      pcol = col,
      domain = domain,
      effect = effect,
      correction = correction,
      x_var = x_var,
      facet_row = facet_row,
      facet_col = facet_col,
      stringsAsFactors = FALSE
    )
  }

  if (length(rows) == 0) return(data.frame())
  bind_rows(rows)
}

cat("model:", model_arg, "\n")
cat("sex:", sex_arg, "\n")
cat("result_dir override:", result_dir_arg %||% "none", "\n")

degree_dirs <- discover_degree_dirs(model_arg, sex_arg, result_dir_arg)
if (!length(degree_dirs)) stop("No degree result directories found.")

result_dirs <- degree_dirs
if (include_network) {
  result_dirs <- unique(c(result_dirs, unlist(lapply(degree_dirs, find_network_dirs))))
}

cat("Selected result directories:\n")
print(result_dirs)

# ----------------------------
# Load demo and degree matrix
# ----------------------------
df <- read_excel(demo_file, sheet = "Sheet1") %>%
  mutate(
    original_project = as.character(`original-project`),
    id_old = as.character(id_old),
    Age_Raw = as.numeric(age_month),
    file_base = paste0(original_project, "_", id_old, "_MIND_DK318_combat"),
    degree_file = file.path(mind_dir, paste0(file_base, "_degree.csv")),
    Diagnosis_Fac = factor(ifelse(group_d_or_c == 0, "TD", "DD"), levels = c("TD", "DD")),
    AgeGroup_Fac = factor(ifelse(group_age == 1, "Adult", "Child"), levels = c("Child", "Adult")),
    Sex = factor(ifelse(sex == 1, "Male", "Female"), levels = c("Female", "Male")),
    has_file = file.exists(degree_file)
  ) %>%
  filter(
    !is.na(original_project),
    !is.na(id_old),
    original_project != "",
    id_old != "",
    !is.na(Age_Raw),
    has_file
  ) %>%
  group_by(AgeGroup_Fac) %>%
  mutate(age_within_group = Age_Raw - mean(Age_Raw, na.rm = TRUE)) %>%
  ungroup()

read_deg <- function(fp) {
  d <- read.csv(fp, check.names = FALSE)
  stopifnot(all(c("ROI", "degree") %in% colnames(d)))
  v <- as.numeric(d$degree)
  names(v) <- as.character(d$ROI)
  v
}

tmp <- read_deg(df$degree_file[1])
rois <- names(tmp)

deg <- matrix(
  NA_real_,
  nrow = nrow(df),
  ncol = length(rois),
  dimnames = list(df$file_base, rois)
)

for (i in seq_len(nrow(df))) {
  v <- read_deg(df$degree_file[i])
  if (!identical(names(v), rois)) stop("ROI mismatch: ", df$degree_file[i])
  deg[i, ] <- v
}

# ----------------------------
# Plot helpers
# ----------------------------
fit_dispersion_proxy <- function(d, include_age = TRUE, sex_filtered = FALSE) {
  rhs <- if (sex_filtered) {
    "Diagnosis_Fac * AgeGroup_Fac"
  } else {
    "Diagnosis_Fac * AgeGroup_Fac * Sex"
  }
  if (include_age) rhs <- paste(rhs, "+ age_within_group")
  form <- as.formula(paste("y ~", rhs))

  fit <- tryCatch(lm(form, data = d), error = function(e) NULL)
  if (is.null(fit)) return(abs(d$y - mean(d$y, na.rm = TRUE)))
  abs(residuals(fit))
}

validate_plot_data <- function(d, x_var, facet_row = NA, facet_col = NA, min_n = 2) {
  if (!nrow(d) || !(x_var %in% names(d))) return(FALSE)

  d[[x_var]] <- droplevels(factor(d[[x_var]]))
  if (nlevels(d[[x_var]]) < 2) return(FALSE)

  panel_id <- rep("Overall", nrow(d))
  if (!is.na(facet_row) && facet_row %in% names(d)) {
    panel_id <- paste(panel_id, as.character(d[[facet_row]]), sep = "__")
  }
  if (!is.na(facet_col) && facet_col %in% names(d)) {
    panel_id <- paste(panel_id, as.character(d[[facet_col]]), sep = "__")
  }

  split_ok <- split(d[[x_var]], panel_id)
  all(vapply(split_ok, function(x) {
    tab <- table(droplevels(factor(x)))
    length(tab) >= 2 && all(tab >= min_n)
  }, logical(1)))
}

fill_values_for_x <- function(x_var, levels_x) {
  palette <- c(
    TD = "#4E79A7",
    DD = "#E15759",
    Child = "#59A14F",
    Adult = "#F28E2B",
    Female = "#B07AA1",
    Male = "#76B7B2"
  )

  out <- palette[levels_x]
  missing <- is.na(out)
  if (any(missing)) {
    fallback <- c("#4E79A7", "#E15759", "#59A14F", "#F28E2B", "#B07AA1", "#76B7B2")
    out[missing] <- fallback[seq_len(sum(missing))]
  }
  out
}

x_label <- function(x_var) {
  if (x_var == "Diagnosis_Fac") return("Diagnosis")
  if (x_var == "AgeGroup_Fac") return("Age Group")
  if (x_var == "Sex") return("Sex")
  x_var
}

effect_label <- function(effect) {
  gsub("_", " × ", effect)
}

build_violin_plot <- function(plot_df, task, feature, p_value, estimate_value, active_label) {
  x_var <- task$x_var
  x_lab <- x_label(x_var)
  y_lab <- if (task$domain == "disp") "|Residual| of Degree" else "Degree"
  fill_vals <- fill_values_for_x(x_var, levels(droplevels(plot_df[[x_var]])))

  subtitle <- paste0(
    task$domain, " ", effect_label(task$effect),
    " | ", task$correction,
    " | p=", signif(p_value, 3)
  )
  if (is.finite(estimate_value)) {
    subtitle <- paste0(subtitle, " | estimate=", signif(estimate_value, 3))
  }
  subtitle <- paste0(subtitle, " | ", active_label)

  g <- ggplot(plot_df, aes(x = .data[[x_var]], y = PlotY, fill = .data[[x_var]])) +
    geom_violin(trim = FALSE, alpha = 0.45, color = NA, width = 0.95) +
    geom_boxplot(width = 0.18, outlier.shape = NA, alpha = 0.85, color = "#333333", fill = "white") +
    geom_jitter(width = 0.10, height = 0, size = 1.6, alpha = 0.55, color = "#333333") +
    scale_fill_manual(values = fill_vals, drop = FALSE) +
    labs(
      title = feature,
      subtitle = subtitle,
      x = x_lab,
      y = y_lab,
      fill = x_lab
    ) +
    theme_bw(base_size = 14) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "#f2f2f2", color = "#d9d9d9"),
      axis.text.x = element_text(size = 11)
    )

  row_var <- if (!is.na(task$facet_row) && task$facet_row %in% names(plot_df)) task$facet_row else NULL
  col_var <- if (!is.na(task$facet_col) && task$facet_col %in% names(plot_df)) task$facet_col else NULL

  if (!is.null(row_var) || !is.null(col_var)) {
    row_part <- row_var %||% "."
    col_part <- col_var %||% "."
    g <- g + facet_grid(stats::as.formula(paste(row_part, "~", col_part)), scales = "free_y")
  }

  g
}

run_one_plot <- function(feature, task, out_dir, p_value, estimate_value,
                         active_label, sex_filter = NA_character_, include_age = TRUE) {
  if (!(feature %in% colnames(deg))) {
    return(list(ok = FALSE, reason = "missing_feature", main = NA_character_))
  }

  d <- df
  if (!is.na(sex_filter)) d <- d %>% filter(Sex == sex_filter)

  d$y <- deg[match(d$file_base, rownames(deg)), feature]
  d <- d %>% filter(is.finite(y))

  if (nrow(d) < 12 || var(d$y, na.rm = TRUE) < 1e-12) {
    return(list(ok = FALSE, reason = "data", main = NA_character_))
  }

  sex_filtered <- !is.na(sex_filter)

  d$PlotY <- if (task$domain == "disp") {
    fit_dispersion_proxy(d, include_age = include_age, sex_filtered = sex_filtered)
  } else {
    d$y
  }
  d <- d[is.finite(d$PlotY), , drop = FALSE]

  if (!validate_plot_data(d, task$x_var, task$facet_row, task$facet_col, min_n = 2)) {
    return(list(ok = FALSE, reason = "group_count", main = NA_character_))
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  out_file <- file.path(
    out_dir,
    paste0("VIOLIN_", safe(task$pcol), "_", safe(feature), ".png")
  )

  g <- build_violin_plot(
    plot_df = d,
    task = task,
    feature = feature,
    p_value = p_value,
    estimate_value = estimate_value,
    active_label = active_label
  )

  has_row <- !is.na(task$facet_row) && task$facet_row %in% names(d)
  has_col <- !is.na(task$facet_col) && task$facet_col %in% names(d)

  width <- if (has_row && has_col) 11 else if (has_col) 8.8 else 7.0
  height <- if (has_row && has_col) 6.8 else if (has_col) 4.9 else 4.6

  ggsave(out_file, g, width = width, height = height, dpi = 300)

  list(ok = TRUE, reason = NA_character_, main = out_file)
}

process_result_dir <- function(active_result_dir) {
  active_result_dir <- normalizePath(active_result_dir, winslash = "/", mustWork = FALSE)
  active_label <- basename(active_result_dir)
  is_network <- basename(dirname(active_result_dir)) == "network_FDR"

  csv_name <- if (is_network) network_csv_name else degree_csv_name
  csv_path <- file.path(active_result_dir, csv_name)

  if (!file.exists(csv_path)) {
    cat("skip missing", csv_path, "\n")
    return(data.frame())
  }

  tab <- read.csv(csv_path, check.names = FALSE)
  if (!("feature" %in% names(tab))) {
    cat("skip missing feature column:", csv_path, "\n")
    return(data.frame())
  }

  tasks <- build_tasks_from_table(tab)
  if (!nrow(tasks)) {
    cat("no interaction p columns found:", csv_path, "\n")
    return(data.frame())
  }

  if (is_network) {
    tasks <- tasks %>%
      filter(correction == "networkFDR")
  }

  if (!nrow(tasks)) {
    cat("no plottable violin tasks found:", csv_path, "\n")
    return(data.frame())
  }

  sex_filter <- parse_sex_filter(active_result_dir)
  include_age <- include_age_for_dir(active_result_dir)

  root_out <- file.path(active_result_dir, "violin_followup")
  dir.create(root_out, recursive = TRUE, showWarnings = FALSE)

  manifest <- list()

  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i, ]
    pcol <- task$pcol
    if (!(pcol %in% names(tab))) next

    pvals <- suppressWarnings(as.numeric(tab[[pcol]]))
    sig_idx <- which(is.finite(pvals) & pvals < p_threshold)

    task_out_dir <- file.path(root_out, task$domain, task$effect, task$correction, safe(pcol))
    cat("task", pcol, ":", length(sig_idx), "sig ROI in", task_out_dir, "\n")

    if (length(sig_idx) == 0) next

    estimate_col <- paste0("estimate_", task$domain, "_", task$effect)
    estimate_vals <- if (estimate_col %in% names(tab)) {
      suppressWarnings(as.numeric(tab[[estimate_col]]))
    } else {
      rep(NA_real_, nrow(tab))
    }

    rows <- vector("list", length(sig_idx))

    for (ii in seq_along(sig_idx)) {
      j <- sig_idx[ii]
      ft <- as.character(tab$feature[j])
      rr <- run_one_plot(
        feature = ft,
        task = task,
        out_dir = task_out_dir,
        p_value = pvals[j],
        estimate_value = estimate_vals[j],
        active_label = active_label,
        sex_filter = sex_filter,
        include_age = include_age
      )

      rows[[ii]] <- data.frame(
        result_dir = active_result_dir,
        task = pcol,
        domain = task$domain,
        effect = task$effect,
        correction = task$correction,
        feature = ft,
        p = pvals[j],
        estimate = estimate_vals[j],
        sex_filter = sex_filter %||% NA_character_,
        include_age = include_age,
        ok = rr$ok,
        reason = rr$reason %||% NA_character_,
        main = rr$main %||% NA_character_,
        stringsAsFactors = FALSE
      )
    }

    task_result <- bind_rows(rows)
    write.csv(
      task_result,
      file.path(task_out_dir, paste0("summary_", safe(pcol), ".csv")),
      row.names = FALSE
    )

    manifest[[length(manifest) + 1]] <- task_result
  }

  manifest_df <- if (length(manifest)) bind_rows(manifest) else data.frame()

  if (nrow(manifest_df)) {
    write.csv(
      manifest_df,
      file.path(active_result_dir, "violin_followup_manifest.csv"),
      row.names = FALSE
    )
    write.csv(
      manifest_df %>%
        group_by(task, domain, effect, correction) %>%
        summarise(n_sig = n(), n_ok = sum(ok, na.rm = TRUE), .groups = "drop"),
      file.path(active_result_dir, "violin_followup_task_overview.csv"),
      row.names = FALSE
    )
  }

  manifest_df
}

all_manifest <- bind_rows(lapply(result_dirs, process_result_dir))

cat("done. success:", sum(all_manifest$ok, na.rm = TRUE), "\n")
