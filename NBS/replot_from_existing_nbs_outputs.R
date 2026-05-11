#!/usr/bin/env Rscript

pkgs <- c("dplyr","ggplot2","stringr","circlize")
for (p in pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = "https://cran.r-project.org")
suppressPackageStartupMessages({library(dplyr); library(ggplot2); library(stringr); library(circlize)})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) normalizePath(args[1], winslash = "/", mustWork = FALSE) else "/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/MIND_DK318_NBS_100"
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
brainnet_dir <- file.path(out_dir, "brainnet_files")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(brainnet_dir, showWarnings = FALSE, recursive = TRUE)

san <- function(x) gsub("^_|_$", "", gsub("_+", "_", gsub("[^A-Za-z0-9_\\-]+", "_", x)))
palette_net <- c(VIS=1,SMN=2,DAN=3,VAN=4,Limbic=5,DMN=6,FPN=7,Other=8)
scale_rng <- function(x, to=c(2,8), def=4) {
  if (!length(x)) return(numeric(0)); x <- as.numeric(x)
  if (all(is.na(x))) return(rep(def, length(x)))
  r <- range(x, na.rm = TRUE); if (!is.finite(r[1]) || !is.finite(r[2]) || diff(r)==0) return(rep(def, length(x)))
  to[1] + (x-r[1])/diff(r)*diff(to)
}
first_col <- function(df, nm) { hit <- nm[nm %in% names(df)][1]; if (is.na(hit)) stop("Missing columns: ", paste(nm, collapse=", ")); hit }
fix_edges <- function(x) {
  if (!nrow(x)) return(x)
  if (!"ROI1" %in% names(x)) x$ROI1 <- x[[first_col(x, c("ROI1","ROI1.x","ROI1.y"))]]
  if (!"ROI2" %in% names(x)) x$ROI2 <- x[[first_col(x, c("ROI2","ROI2.x","ROI2.y"))]]
  if (!"direction" %in% names(x)) x$direction <- NA_character_
  x
}
node_degree <- function(edges, roi_meta) {
  edges <- fix_edges(edges)
  bind_rows(transmute(edges, ROI_index=i, ROI=ROI1), transmute(edges, ROI_index=j, ROI=ROI2)) %>%
    count(ROI_index, ROI, name="degree") %>%
    left_join(roi_meta, by=c("ROI_index","ROI")) %>%
    arrange(desc(degree), ROI_index)
}
write_component_top_roi_table <- function(edges, roi_meta, nm, top_n=10L) {
  if (!nrow(edges)) return()
  edges <- fix_edges(edges)
  out <- bind_rows(lapply(sort(unique(edges$Component)), function(cc) {
    e_comp <- edges %>% filter(Component == cc)
    nd <- node_degree(e_comp, roi_meta)
    if (!nrow(nd)) return(NULL)
    nd %>%
      mutate(
        Component = cc,
        rank_within_component = row_number()
      ) %>%
      slice_head(n = top_n) %>%
      select(Component, rank_within_component, ROI_index, ROI, degree, everything())
  }))
  if (!nrow(out)) return()
  write.csv(out, file.path(tab_dir, paste0(san(nm), "_component_top10_rois_by_degree.csv")), row.names=FALSE)
}
prepare_grouped_df <- function(dd, nm) {
  m_overall <- str_match(nm, "^overall_(DD_gt_TD|TD_gt_DD)$")
  m_age_sex <- str_match(nm, "^(Child|Adult)_(Female|Male)_(DD_gt_TD|TD_gt_DD)$")
  m_age_only <- str_match(nm, "^(Child|Adult)_(DD_gt_TD|TD_gt_DD)$")
  xv <- "Group8"; xl <- "AgeGroup_Sex_Diagnosis"
  if (!all(is.na(m_overall))) {
    dd <- dd %>%
      mutate(Group2=factor(Diagnosis, levels=c("TD","DD"), labels=c("overall_TD", "overall_DD")))
    xv <- "Group2"; xl <- "Diagnosis in overall sample"
  } else if (!all(is.na(m_age_sex))) {
    dd <- dd %>%
      filter(AgeGroup==m_age_sex[,2], Sex==m_age_sex[,3]) %>%
      mutate(Group2=factor(Diagnosis, levels=c("TD","DD"), labels=c(paste0(m_age_sex[,2], "_", m_age_sex[,3], "_TD"), paste0(m_age_sex[,2], "_", m_age_sex[,3], "_DD"))))
    xv <- "Group2"; xl <- "Diagnosis within matched AgeGroup and Sex"
  } else if (!all(is.na(m_age_only))) {
    dd <- dd %>%
      filter(AgeGroup==m_age_only[,2]) %>%
      mutate(Group2=factor(Diagnosis, levels=c("TD","DD"), labels=c(paste0(m_age_only[,2], "_TD"), paste0(m_age_only[,2], "_DD"))))
    xv <- "Group2"; xl <- "Diagnosis within matched AgeGroup"
  }
  list(data=dd, xvar=xv, xlabel=xl)
}
group_palette <- function(groups) {
  u <- unique(as.character(groups))
  cols <- ifelse(grepl("(^|_)TD$", u), "#4C9BE8", ifelse(grepl("(^|_)DD$", u), "#E6864A", "#7A7A7A"))
  setNames(cols, u)
}
group_labels <- function(groups) {
  u <- unique(as.character(groups))
  labs <- ifelse(grepl("(^|_)TD$", u), "TD", ifelse(grepl("(^|_)DD$", u), "DD", u))
  setNames(labs, u)
}
plot_strength <- function(edges, net, pheno, nm) {
  if (!nrow(edges) || is.null(net) || !nrow(net)) return()
  dd <- bind_rows(lapply(sort(unique(edges$Component)), function(cc) {
    idx <- unique(edges$edge_index[edges$Component==cc]); idx <- idx[idx>=1 & idx<=ncol(net)]; if (!length(idx)) return(NULL)
    data.frame(subj_prefix=pheno$subj_prefix, Diagnosis=pheno$Diagnosis, AgeGroup=pheno$AgeGroup, Sex=pheno$Sex,
      Group8=interaction(pheno$AgeGroup, pheno$Sex, pheno$Diagnosis, sep="_", drop=TRUE),
      Component=paste0("Component_", cc), component_strength=rowMeans(net[, idx, drop=FALSE], na.rm=TRUE))
  }))
  if (!nrow(dd)) return()
  gp <- prepare_grouped_df(dd, nm)
  dd <- gp$data; xv <- gp$xvar; xl <- gp$xlabel
  pal <- group_palette(dd[[xv]])
  glab <- group_labels(dd[[xv]])
  write.csv(dd, file.path(tab_dir, paste0(san(nm), "_subject_component_strength.csv")), row.names=FALSE)
  p <- ggplot(dd, aes(x=.data[[xv]], y=component_strength, fill=.data[[xv]], color=.data[[xv]])) +
    geom_violin(trim=FALSE, alpha=.30, show.legend=TRUE) +
    geom_boxplot(width=.15, outlier.shape=NA, alpha=.55, show.legend=FALSE) +
    geom_jitter(width=.12, size=.8, alpha=.65, show.legend=FALSE) +
    facet_wrap(~Component, scales="free_y") +
    scale_fill_manual(values=pal, breaks=names(glab), labels=unname(glab), name=NULL) +
    scale_color_manual(values=pal, breaks=names(glab), labels=unname(glab), name=NULL) +
    guides(color="none", fill=guide_legend(override.aes=list(alpha=.9))) +
    labs(title=paste0(nm, " component strength"), x=xl, y="Mean MIND value within NBS component") +
    theme_bw(10) +
    theme(
      axis.text.x=element_text(angle=45, hjust=1),
      legend.position=c(.98, .98),
      legend.justification=c(1, 1),
      legend.direction="horizontal",
      legend.background=element_rect(fill=scales::alpha("white", .75), colour=NA),
      legend.key=element_blank(),
      legend.margin=margin(1, 1, 1, 1),
      legend.box.margin=margin(0, 0, 0, 0)
    )
  ggsave(file.path(fig_dir, paste0(san(nm), "_component_strength_violin.png")), p, width=11, height=6, dpi=300)
}
plot_top_roi_strength <- function(edges, net, pheno, nm, roi_meta, top_n=3L) {
  if (!nrow(edges) || is.null(net) || !nrow(net)) return()
  nd <- node_degree(edges, roi_meta)
  top_roi <- nd %>% slice_head(n=top_n)
  if (!nrow(top_roi)) return()
  dd <- bind_rows(lapply(seq_len(nrow(top_roi)), function(k) {
    roi_idx <- top_roi$ROI_index[k]
    roi_nm <- top_roi$ROI[k]
    roi_lab <- roi_nm
    idx <- unique(edges$edge_index[edges$i==roi_idx | edges$j==roi_idx])
    idx <- idx[idx>=1 & idx<=ncol(net)]
    if (!length(idx)) return(NULL)
    data.frame(subj_prefix=pheno$subj_prefix, Diagnosis=pheno$Diagnosis, AgeGroup=pheno$AgeGroup, Sex=pheno$Sex,
      Group8=interaction(pheno$AgeGroup, pheno$Sex, pheno$Diagnosis, sep="_", drop=TRUE),
      ROI=factor(roi_lab, levels=top_roi$ROI),
      ROI_name=roi_nm, ROI_index=roi_idx, degree=top_roi$degree[k], roi_strength=rowMeans(net[, idx, drop=FALSE], na.rm=TRUE))
  }))
  if (!nrow(dd)) return()
  gp <- prepare_grouped_df(dd, nm)
  dd <- gp$data; xv <- gp$xvar; xl <- gp$xlabel
  pal <- group_palette(dd[[xv]])
  glab <- group_labels(dd[[xv]])
  write.csv(dd, file.path(tab_dir, paste0(san(nm), "_top3_roi_subject_strength.csv")), row.names=FALSE)
  p <- ggplot(dd, aes(x=ROI, y=roi_strength, fill=.data[[xv]], color=.data[[xv]])) +
    geom_violin(position=position_dodge(width=.85), trim=FALSE, alpha=.28, show.legend=TRUE) +
    geom_boxplot(position=position_dodge(width=.85), width=.15, outlier.shape=NA, alpha=.55, show.legend=FALSE) +
    geom_jitter(position=position_jitterdodge(jitter.width=.12, dodge.width=.85), size=.8, alpha=.6, show.legend=FALSE) +
    scale_fill_manual(values=pal, breaks=names(glab), labels=unname(glab), name=NULL) +
    scale_color_manual(values=pal, breaks=names(glab), labels=unname(glab), name=NULL) +
    guides(color="none", fill=guide_legend(override.aes=list(alpha=.9))) +
    labs(title=paste0(nm, " top ", top_n, " ROI strength"), subtitle="Top ROIs ranked by number of significant component edges", x="Top ROIs in NBS component", y="Mean MIND value across component edges incident to ROI", fill=xl, color=xl) +
    theme_bw(10) +
    theme(
      axis.text.x=element_text(angle=25, hjust=1),
      legend.position=c(.98, .98),
      legend.justification=c(1, 1),
      legend.direction="horizontal",
      legend.background=element_rect(fill=scales::alpha("white", .75), colour=NA),
      legend.key=element_blank(),
      legend.margin=margin(1, 1, 1, 1),
      legend.box.margin=margin(0, 0, 0, 0)
    )
  ggsave(file.path(fig_dir, paste0(san(nm), "_top3_roi_strength_violin.png")), p, width=10, height=6, dpi=300)
}
plot_lobe <- function(edges, edge_annot, nm) {
  if (!nrow(edges)) return()
  
  # Create network_pair if not exists
  if (!"network_pair" %in% names(edges)) {
    edges <- edges %>% mutate(network_pair = ifelse(network1 <= network2, paste(network1, network2, sep="--"), paste(network2, network1, sep="--")))
  }
  if (!"network_pair" %in% names(edge_annot)) {
    edge_annot <- edge_annot %>% mutate(network_pair = ifelse(network1 <= network2, paste(network1, network2, sep="--"), paste(network2, network1, sep="--")))
  }
  
  # Network-level statistics
  obs_net <- edges %>% count(network_pair, name="observed_edges") %>% 
    left_join(count(edge_annot, network_pair, name="possible_edges"), by="network_pair") %>%
    mutate(normalized_proportion=observed_edges/possible_edges, 
           network_a=sapply(strsplit(network_pair,"--",fixed=TRUE), `[`, 1), 
           network_b=sapply(strsplit(network_pair,"--",fixed=TRUE), `[`, 2))
  
  # Define network order (reverse for y-axis to match the example)
  network_order <- c("VIS", "DMN", "SMN", "Limbic", "VAN", "DAN", "FPN")
  
  # Create complete grid with all possible combinations
  all_combinations <- expand.grid(
    network_a = network_order,
    network_b = network_order,
    stringsAsFactors = FALSE
  )
  
  # Only keep upper triangle (including diagonal)
  all_combinations <- all_combinations %>%
    mutate(
      idx_a = match(network_a, network_order),
      idx_b = match(network_b, network_order)
    ) %>%
    filter(idx_a >= idx_b) %>%
    select(-idx_a, -idx_b)
  
  # Merge with observed data
  plot_data <- all_combinations %>%
    left_join(obs_net, by = c("network_a", "network_b")) %>%
    mutate(
      observed_edges = ifelse(is.na(observed_edges), 0, observed_edges),
      normalized_proportion = ifelse(is.na(normalized_proportion), 0, normalized_proportion),
      network_a = factor(network_a, levels = network_order),
      network_b = factor(network_b, levels = rev(network_order))
    )
  
  write.csv(obs_net, file.path(tab_dir, paste0(san(nm), "_network_level_edge_counts.csv")), row.names=FALSE)
  
  # Plot 1: Raw count with circles
  p1 <- ggplot(plot_data, aes(x=network_a, y=network_b)) + 
    geom_point(aes(size=observed_edges, fill=observed_edges), shape=21, color="white", stroke=0.5) + 
    scale_size_continuous(range = c(0, 20), limits = c(0, NA), guide = "none") +
    scale_fill_gradient(low="white", high="#D62728", limits = c(0, NA), name="Edges") +
    coord_equal() + 
    labs(title=paste0(nm, " network-level raw edge count"), x=NULL, y=NULL) + 
    theme_minimal(base_size=12) + 
    theme(
      axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=11),
      axis.text.y=element_text(size=11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill="white", color=NA),
      plot.background = element_rect(fill="white", color=NA),
      legend.position = "right",
      plot.title = element_text(hjust=0.5, size=13)
    )
  
  # Plot 2: Normalized proportion with circles and color scale
  p2 <- ggplot(plot_data, aes(x=network_a, y=network_b)) + 
    geom_point(aes(size=normalized_proportion, fill=normalized_proportion), shape=21, color="white", stroke=0.5) + 
    scale_size_continuous(range = c(0, 20), limits = c(0, NA), guide = "none") +
    scale_fill_gradient2(
      low="#3182BD", mid="#FFFFBF", high="#D62728", 
      midpoint=0.05, limits = c(0, NA), 
      name="Proportion"
    ) +
    coord_equal() + 
    labs(title=paste0(nm, " network-level normalized proportion"), x=NULL, y=NULL) + 
    theme_minimal(base_size=12) + 
    theme(
      axis.text.x=element_text(angle=45, hjust=1, vjust=1, size=11),
      axis.text.y=element_text(size=11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_rect(fill="white", color=NA),
      plot.background = element_rect(fill="white", color=NA),
      legend.position = "right",
      plot.title = element_text(hjust=0.5, size=13)
    )
  
  ggsave(file.path(fig_dir, paste0(san(nm), "_lobe_heatmap_raw_count.png")), p1, width=8, height=7, dpi=300, bg="white")
  ggsave(file.path(fig_dir, paste0(san(nm), "_lobe_heatmap_normalized.png")), p2, width=8, height=7, dpi=300, bg="white")
}
plot_adj <- function(edges, roi_meta, nm) {
  n <- nrow(roi_meta); mt <- mb <- matrix(0, n, n)
  for (k in seq_len(nrow(edges))) { i <- edges$i[k]; j <- edges$j[k]; mt[i,j] <- edges$t[k]; mt[j,i] <- edges$t[k]; mb[i,j] <- 1; mb[j,i] <- 1 }
  rownames(mt) <- colnames(mt) <- roi_meta$ROI; rownames(mb) <- colnames(mb) <- roi_meta$ROI
  write.csv(mt, file.path(tab_dir, paste0(san(nm), "_NBS_t_matrix.csv"))); write.csv(mb, file.path(tab_dir, paste0(san(nm), "_NBS_binary_mask.csv")))
  df <- as.data.frame(as.table(mt)); names(df) <- c("ROI1","ROI2","t")
  p <- ggplot(df, aes(ROI1, ROI2, fill=t)) + geom_tile() + labs(title=paste0(nm, " NBS t-value matrix"), x=NULL, y=NULL, fill="t") + theme_bw(8) + theme(axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks=element_blank())
  ggsave(file.path(fig_dir, paste0(san(nm), "_NBS_t_matrix_heatmap.png")), p, width=7, height=7, dpi=300)
}
plot_chord2 <- function(edges, nm, max_edges=400L) {
  if (!nrow(edges)) return(); edges <- fix_edges(edges)
  cd <- edges %>% mutate(abs_t=abs(t)) %>% arrange(desc(abs_t)) %>% slice_head(n=max_edges) %>% transmute(from=ROI1, to=ROI2, value=abs_t)
  png(file.path(fig_dir, paste0(san(nm), "_chord_plot.png")), width=2700, height=2700, res=300)
  tryCatch({circlize::circos.clear(); circlize::chordDiagram(cd, directional=0, annotationTrack="grid", transparency=.65); title(nm); circlize::circos.clear()}, error=function(e) message("Chord plot failed for ", nm, ": ", conditionMessage(e)))
  dev.off()
}
write_brainnet <- function(edges, roi_meta, nm) {
  if (!nrow(edges) || any(is.na(roi_meta$x))) return(); edges <- fix_edges(edges)
  sm <- unique(edges[, intersect(c("Component","ncomp","ncompFWE","strn","strnFWE","component_sig"), names(edges)), drop=FALSE])
  if (ncol(sm)) write.csv(sm, file.path(brainnet_dir, paste0(san(nm), "_BrainNet_component_summary.csv")), row.names=FALSE)
  for (cc in sort(unique(edges$Component))) {
    e <- filter(edges, Component==cc); ids <- sort(unique(c(e$i, e$j))); nodes <- roi_meta %>% filter(ROI_index %in% ids) %>% arrange(ROI_index)
    dg <- node_degree(e, roi_meta) %>% select(ROI_index, degree)
    nodes <- left_join(nodes, dg, by="ROI_index") %>% mutate(degree=ifelse(is.na(degree),0,degree), network=ifelse(is.na(network)|network=="","Other",network), brainnet_color=unname(palette_net[network]), brainnet_color=ifelse(is.na(brainnet_color), palette_net[["Other"]], brainnet_color), brainnet_size=scale_rng(degree))
    mp <- setNames(seq_along(nodes$ROI_index), nodes$ROI_index); em <- em1 <- em2 <- matrix(0, nrow(nodes), nrow(nodes))
    for (k in seq_len(nrow(e))) { a <- mp[as.character(e$i[k])]; b <- mp[as.character(e$j[k])]; w <- abs(e$t[k]); em[a,b] <- em[b,a] <- w; if (e$t[k] > 0) em1[a,b] <- em1[b,a] <- w; if (e$t[k] < 0) em2[a,b] <- em2[b,a] <- w }
    write.table(data.frame(x=nodes$x,y=nodes$y,z=nodes$z,color=nodes$brainnet_color,size=nodes$brainnet_size,label=nodes$ROI), file.path(brainnet_dir, paste0(san(nm), "_component_", cc, ".node")), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(em,  file.path(brainnet_dir, paste0(san(nm), "_component_", cc, ".edge")), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(em1, file.path(brainnet_dir, paste0(san(nm), "_component_", cc, "_TD_gt_DD.edge")), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.table(em2, file.path(brainnet_dir, paste0(san(nm), "_component_", cc, "_DD_gt_TD.edge")), sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    write.csv(nodes %>% select(ROI_index, ROI, hemi, label, dk_region, lobe, network, x, y, z, degree, brainnet_color, brainnet_size), file.path(brainnet_dir, paste0(san(nm), "_component_", cc, "_node_annotation.csv")), row.names=FALSE)
    ea <- e[, intersect(c("edge_index","Component","ROI1","ROI2","i","j","t","p","p_fdr","direction","network1","network2","lobe1","lobe2","dk_region1","dk_region2"), names(e)), drop=FALSE]; ea$abs_t <- abs(ea$t)
    ea$td_direction <- ifelse(ea$t>0, "TD > DD", ifelse(ea$t<0, "DD > TD", "zero"))
    write.csv(ea, file.path(brainnet_dir, paste0(san(nm), "_component_", cc, "_edge_annotation.csv")), row.names=FALSE)
  }
  write.csv(data.frame(network=names(palette_net), brainnet_color_code=unname(palette_net)), file.path(brainnet_dir, paste0(san(nm), "_BrainNet_network_color_key.csv")), row.names=FALSE)
}

roi_meta <- read.csv(file.path(out_dir, "roi_metadata_with_mni_coordinates.csv"), check.names=FALSE)
edge_annot <- read.csv(file.path(out_dir, "edge_annotation_roi_mni_lobe.csv"), check.names=FALSE)

# Add network_pair if not exists
if ("network1" %in% names(edge_annot) && "network2" %in% names(edge_annot) && !"network_pair" %in% names(edge_annot)) {
  edge_annot <- edge_annot %>% mutate(network_pair = ifelse(network1 <= network2, paste(network1, network2, sep="--"), paste(network2, network1, sep="--")))
}

idata <- read.csv(file.path(out_dir, "idata_nbr_effect_coded.csv"), check.names=FALSE)
net2d <- readRDS(file.path(out_dir, "dk318_mind_upper_triangle_matrix.rds"))
idata <- idata %>% mutate(Diagnosis=factor(Diagnosis, levels=c("TD","DD")), AgeGroup=factor(AgeGroup, levels=c("Child","Adult")), Sex=factor(Sex, levels=c("Female","Male")))
files <- list.files(tab_dir, pattern="_edges_in_FWE_significant_components\\.csv$", full.names=TRUE)
files <- files[!grepl("NBR_ALL_edges_in_FWE_significant_components\\.csv$", files)]
all_edges_file <- file.path(tab_dir, "NBR_ALL_edges_in_FWE_significant_components.csv")
all_edges <- if (file.exists(all_edges_file)) read.csv(all_edges_file, check.names=FALSE) else data.frame()
analysis_tables <- list()
if (length(files)) {
  for (fp in files) {
    nm <- sub("_edges_in_FWE_significant_components\\.csv$", "", basename(fp))
    analysis_tables[[nm]] <- read.csv(fp, check.names=FALSE)
  }
}
if (nrow(all_edges) && "analysis" %in% names(all_edges)) {
  split_edges <- split(all_edges, all_edges$analysis)
  for (nm in names(split_edges)) {
    if (!nm %in% names(analysis_tables)) analysis_tables[[nm]] <- split_edges[[nm]]
  }
}
if (!length(analysis_tables)) stop("No per-analysis edge tables found in ", tab_dir, " and no usable NBR_ALL_edges_in_FWE_significant_components.csv found.")
message("Found ", length(analysis_tables), " analyses.")
for (nm in names(analysis_tables)) {
  message("Replotting: ", nm)
  e <- analysis_tables[[nm]]; if (!nrow(e)) next; e <- fix_edges(e)
  if (!"direction" %in% names(e) || all(is.na(e$direction))) e$direction <- ifelse(grepl("^interaction_", nm), ifelse(e$t>0, "positive interaction contrast", ifelse(e$t<0, "negative interaction contrast", "zero")), ifelse(e$t>0, "TD > DD", ifelse(e$t<0, "DD > TD", "zero")))
  m_age_sex <- str_match(nm, "^(Child|Adult)_(Female|Male)_(DD_gt_TD|TD_gt_DD)$")
  m_age_only <- str_match(nm, "^(Child|Adult)_(DD_gt_TD|TD_gt_DD)$")

  if (!all(is.na(m_age_sex))) {
    keep <- idata$AgeGroup == m_age_sex[,2] & idata$Sex == m_age_sex[,3]
    net <- net2d[keep,,drop=FALSE]
    pheno <- idata[keep,,drop=FALSE]
  } else if (!all(is.na(m_age_only))) {
    keep <- idata$AgeGroup == m_age_only[,2]
    net <- net2d[keep,,drop=FALSE]
    pheno <- idata[keep,,drop=FALSE]
  } else {
    net <- net2d
    pheno <- idata
  }

  write.csv(node_degree(e, roi_meta), file.path(tab_dir, paste0(san(nm), "_node_degree.csv")), row.names=FALSE)
  write_component_top_roi_table(e, roi_meta, nm)
  plot_strength(e, net, pheno, nm); plot_top_roi_strength(e, net, pheno, nm, roi_meta); plot_adj(e, roi_meta, nm); plot_lobe(e, edge_annot, nm); plot_chord2(e, nm); write_brainnet(e, roi_meta, nm)
}
message("Done: ", fig_dir, " and ", brainnet_dir)
