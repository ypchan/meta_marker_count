# -*- coding: UTF-8 -*-

# ============================================================
# Fig_community_assembly_network.R
#
# Community assembly process and co-occurrence network structure
# for Bacteria, Archaea and Fungi along sediment depth.
#
# ASCII-only script for RStudio and server use.
# ============================================================

# ============================================================
# Input file: microbial abundance table
# ============================================================
#
# Default file:
#   all.marker_rpm.genus.long.tsv
#
# Required columns:
#   sample_id
#   domain
#   rank
#   lineage
#   marker_rpm
#   depth
#
# Optional columns:
#   year
#   month
#   site_id
#   site_type
#   marker
#   taxon_marker_reads
#
# Example:
# sample_id        year month depth marker domain   rank  lineage                               marker_rpm
# 201704_MF1_0002 2017 4     0-2   16S    Archaea  genus d__Archaea;p__Asgardarchaeota;g__X     0.354
# 201704_MF1_0002 2017 4     0-2   ITS    Fungi    genus k__Fungi;p__Ascomycota;g__Y             0.014
#
# Notes:
#   1. lineage is used as the unique taxon identifier.
#   2. marker_rpm is used as abundance.
#   3. depth is used to group community assembly and network analyses.
#   4. For UNITE fungi, SH-aware lineage labels are recommended.
#
# Main outputs:
#   community_assembly_network/assembly_pairwise_RCbray.tsv
#   community_assembly_network/assembly_process_summary.tsv
#   community_assembly_network/network_edges.tsv
#   community_assembly_network/network_nodes.tsv
#   community_assembly_network/network_stats.tsv
#   community_assembly_network/community_assembly_process.pdf
#   community_assembly_network/network_structure_by_depth.pdf
#
# ============================================================

options(stringsAsFactors = FALSE)

# ----------------------------
# 0. Packages
# ----------------------------

packages <- c(
  "tidyverse",
  "data.table",
  "vegan",
  "igraph",
  "ggraph",
  "patchwork",
  "ggsci",
  "scales"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(tidyverse)
library(data.table)
library(vegan)
library(igraph)
library(ggraph)
library(patchwork)
library(ggsci)
library(scales)

# ----------------------------
# 1. Parameters
# ----------------------------

abund_file <- "all.marker_rpm.genus.long.tsv"

output_dir <- "community_assembly_network"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

domain_keep <- c("Bacteria", "Archaea", "Fungi")
target_rank <- "genus"

depth_col <- "depth"

# Set to NULL to use the order appearing in the table.
depth_order <- NULL

# Example:
# depth_order <- c("0-2", "6-8", "10-12", "12-14", "20-22", "28-30")
# depth_order <- c("00-10", "10-20", "20-30", "30-40", "40-50", "50-60")

# Feature filtering
min_detected_samples <- 3
min_total_rpm <- 0

# Community assembly null model
run_assembly_analysis <- TRUE
assembly_n_perm <- 499
assembly_seed <- 123
assembly_min_samples <- 6

# RCbray thresholds
rc_high <- 0.95
rc_low <- -0.95

# Network analysis
run_network_analysis <- TRUE
network_min_samples <- 6
network_top_n_taxa_per_domain_depth <- 80
network_min_prevalence <- 0.20
network_cor_method <- "spearman"
network_r_cutoff <- 0.60
network_q_cutoff <- 0.05
network_p_adjust <- "BH"

# Plot
node_size <- 1.8
edge_width <- 0.35

# ----------------------------
# 2. Helper functions
# ----------------------------

theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = "black"),
      axis.text = element_text(color = "black", size = base_size),
      axis.title = element_text(color = "black", size = base_size + 1),
      axis.line = element_line(linewidth = 0.35, color = "black"),
      axis.ticks = element_line(linewidth = 0.35, color = "black"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size + 1),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(0.35, "cm"),
      plot.title = element_text(face = "bold", size = base_size + 2),
      plot.margin = margin(4, 4, 4, 4)
    )
}

rank_prefix <- c(
  domain = "d__",
  phylum = "p__",
  class = "c__",
  order = "o__",
  family = "f__",
  genus = "g__",
  species = "s__"
)

extract_rank_label <- function(lineage, rank) {
  prefix <- rank_prefix[[rank]]
  
  if (is.null(prefix)) {
    return(lineage)
  }
  
  parts <- stringr::str_split(lineage, ";", simplify = FALSE)[[1]]
  hit <- parts[stringr::str_starts(parts, fixed(prefix))]
  
  if (length(hit) == 0) {
    out <- tail(parts, 1)
  } else {
    out <- hit[1]
  }
  
  out <- stringr::str_remove(out, paste0("^", prefix))
  out <- stringr::str_replace_all(out, "_", " ")
  out
}

make_safe_id <- function(x) {
  make.unique(gsub("[^A-Za-z0-9_.]+", "_", x))
}

make_community_matrix <- function(df, samples, features) {
  df %>%
    group_by(sample_id, feature_id) %>%
    summarise(abundance = sum(marker_rpm, na.rm = TRUE), .groups = "drop") %>%
    complete(
      sample_id = samples,
      feature_id = features,
      fill = list(abundance = 0)
    ) %>%
    pivot_wider(
      names_from = feature_id,
      values_from = abundance,
      values_fill = 0
    ) %>%
    arrange(match(sample_id, samples))
}

filter_features_for_group <- function(df) {
  df %>%
    group_by(feature_id, feature_label, domain) %>%
    summarise(
      total_rpm = sum(marker_rpm, na.rm = TRUE),
      detected_samples = n_distinct(sample_id[marker_rpm > 0]),
      prevalence = detected_samples / n_distinct(df$sample_id),
      .groups = "drop"
    ) %>%
    filter(
      total_rpm > min_total_rpm,
      detected_samples >= min_detected_samples
    )
}

classify_rc_process <- function(rc) {
  case_when(
    is.na(rc) ~ "NA",
    rc > rc_high ~ "Dispersal limitation",
    rc < rc_low ~ "Homogenizing dispersal",
    TRUE ~ "Drift or weak selection"
  )
}

calc_rcbray_null <- function(mat, n_perm = 499, seed = 123) {
  set.seed(seed)
  
  sample_ids <- rownames(mat)
  
  if (nrow(mat) < 3 || ncol(mat) < 2) {
    return(tibble())
  }
  
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  
  keep_rows <- rowSums(mat) > 0
  keep_cols <- colSums(mat) > 0
  
  mat <- mat[keep_rows, keep_cols, drop = FALSE]
  sample_ids <- rownames(mat)
  
  if (nrow(mat) < 3 || ncol(mat) < 2) {
    return(tibble())
  }
  
  obs <- as.vector(vegan::vegdist(mat, method = "bray"))
  
  pair_index <- combn(sample_ids, 2)
  pair_df <- tibble(
    sample1 = pair_index[1, ],
    sample2 = pair_index[2, ],
    obs_bray = obs
  )
  
  null_less <- rep(0, length(obs))
  null_equal <- rep(0, length(obs))
  
  for (i in seq_len(n_perm)) {
    null_mat <- mat
    
    for (j in seq_len(ncol(null_mat))) {
      null_mat[, j] <- sample(null_mat[, j], size = nrow(null_mat), replace = FALSE)
    }
    
    null_dist <- as.vector(vegan::vegdist(null_mat, method = "bray"))
    
    null_less <- null_less + as.integer(null_dist < obs)
    null_equal <- null_equal + as.integer(null_dist == obs)
  }
  
  rc <- ((null_less + 0.5 * null_equal) / n_perm) * 2 - 1
  
  pair_df %>%
    mutate(
      rc_bray = rc,
      assembly_process = classify_rc_process(rc_bray)
    )
}

pairwise_cor_network <- function(mat, method = "spearman") {
  mat <- as.matrix(mat)
  mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  
  feature_ids <- colnames(mat)
  
  if (length(feature_ids) < 2) {
    return(tibble())
  }
  
  results <- list()
  idx <- 1
  
  for (i in seq_len(length(feature_ids) - 1)) {
    for (j in (i + 1):length(feature_ids)) {
      x <- mat[, i]
      y <- mat[, j]
      
      ok <- is.finite(x) & is.finite(y)
      
      if (sum(ok) < network_min_samples) {
        next
      }
      
      if (length(unique(x[ok])) < 2 || length(unique(y[ok])) < 2) {
        next
      }
      
      ct <- suppressWarnings(
        cor.test(x[ok], y[ok], method = method, exact = FALSE)
      )
      
      results[[idx]] <- tibble(
        from = feature_ids[i],
        to = feature_ids[j],
        rho = unname(ct$estimate),
        p_value = ct$p.value,
        n = sum(ok)
      )
      
      idx <- idx + 1
    }
  }
  
  bind_rows(results)
}

network_stats_from_edges <- function(edges, nodes, depth_value, domain_value) {
  if (nrow(edges) == 0) {
    return(tibble(
      depth = depth_value,
      domain = domain_value,
      nodes = nrow(nodes),
      edges = 0,
      positive_edges = 0,
      negative_edges = 0,
      density = NA_real_,
      mean_degree = NA_real_,
      modularity = NA_real_,
      largest_component = NA_integer_
    ))
  }
  
  g <- igraph::graph_from_data_frame(
    d = edges %>% select(from, to, rho, q_value, sign),
    vertices = nodes,
    directed = FALSE
  )
  
  comp <- igraph::components(g)
  
  if (igraph::ecount(g) > 0 && igraph::vcount(g) > 2) {
    clu <- suppressWarnings(igraph::cluster_louvain(g))
    mod <- igraph::modularity(clu)
  } else {
    mod <- NA_real_
  }
  
  tibble(
    depth = depth_value,
    domain = domain_value,
    nodes = igraph::vcount(g),
    edges = igraph::ecount(g),
    positive_edges = sum(E(g)$sign == "Positive"),
    negative_edges = sum(E(g)$sign == "Negative"),
    density = igraph::edge_density(g, loops = FALSE),
    mean_degree = mean(igraph::degree(g)),
    modularity = mod,
    largest_component = max(comp$csize)
  )
}

plot_one_network <- function(edges, nodes, title_text) {
  if (nrow(edges) == 0 || nrow(nodes) == 0) {
    return(
      ggplot() +
        annotate("text", x = 0, y = 0, label = paste0(title_text, "\nNo edges")) +
        theme_void()
    )
  }
  
  g <- igraph::graph_from_data_frame(
    d = edges,
    vertices = nodes,
    directed = FALSE
  )
  
  ggraph(g, layout = "fr") +
    geom_edge_link(
      aes(color = sign),
      alpha = 0.45,
      linewidth = edge_width
    ) +
    geom_node_point(
      aes(color = domain),
      size = node_size,
      alpha = 0.9
    ) +
    scale_edge_color_manual(
      values = c(
        Positive = "#F4A6A6",
        Negative = "#74CBE8"
      ),
      name = "Interaction"
    ) +
    scale_color_npg(name = "Domain") +
    labs(title = title_text) +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      legend.position = "bottom",
      plot.margin = margin(3, 3, 3, 3)
    )
}

# ----------------------------
# 3. Read abundance table
# ----------------------------

if (!file.exists(abund_file)) {
  stop("Abundance file not found: ", abund_file)
}

abund_raw <- data.table::fread(
  abund_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

names(abund_raw) <- trimws(names(abund_raw))

required_cols <- c("sample_id", "domain", "rank", "lineage", "marker_rpm", depth_col)
missing_cols <- setdiff(required_cols, names(abund_raw))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

abund_dat <- abund_raw %>%
  mutate(
    sample_id = as.character(sample_id),
    domain = as.character(domain),
    rank = as.character(rank),
    lineage = as.character(lineage),
    marker_rpm = as.numeric(marker_rpm),
    depth_group = as.character(.data[[depth_col]])
  ) %>%
  filter(
    domain %in% domain_keep,
    rank == target_rank,
    !is.na(marker_rpm),
    marker_rpm > 0
  ) %>%
  mutate(
    feature_id = lineage,
    feature_safe_id = make_safe_id(lineage),
    feature_label = map_chr(lineage, extract_rank_label, rank = target_rank)
  )

if (is.null(depth_order)) {
  depth_order <- unique(abund_dat$depth_group)
}

abund_dat <- abund_dat %>%
  mutate(
    depth_group = factor(depth_group, levels = depth_order),
    domain = factor(domain, levels = domain_keep)
  )

# ----------------------------
# 4. Community assembly analysis
# ----------------------------

assembly_pairwise_all <- tibble()
assembly_summary <- tibble()

if (run_assembly_analysis) {
  message("[INFO] Running RCbray null-model community assembly analysis")
  
  assembly_results <- list()
  idx <- 1
  
  for (dm in domain_keep) {
    for (dp in depth_order) {
      sub <- abund_dat %>%
        filter(domain == dm, depth_group == dp)
      
      samples <- sort(unique(sub$sample_id))
      
      if (length(samples) < assembly_min_samples) {
        next
      }
      
      feat_info <- filter_features_for_group(sub)
      features <- feat_info$feature_id
      
      if (length(features) < 2) {
        next
      }
      
      mat_df <- make_community_matrix(sub, samples, features)
      mat <- mat_df %>%
        select(-sample_id) %>%
        as.matrix()
      
      rownames(mat) <- mat_df$sample_id
      
      one <- calc_rcbray_null(
        mat = mat,
        n_perm = assembly_n_perm,
        seed = assembly_seed
      ) %>%
        mutate(
          domain = dm,
          depth = as.character(dp)
        )
      
      assembly_results[[idx]] <- one
      idx <- idx + 1
    }
  }
  
  assembly_pairwise_all <- bind_rows(assembly_results)
  
  if (nrow(assembly_pairwise_all) > 0) {
    assembly_summary <- assembly_pairwise_all %>%
      group_by(domain, depth, assembly_process) %>%
      summarise(n_pairs = n(), .groups = "drop") %>%
      group_by(domain, depth) %>%
      mutate(
        total_pairs = sum(n_pairs),
        percentage = n_pairs / total_pairs * 100
      ) %>%
      ungroup()
    
    write_tsv(
      assembly_pairwise_all,
      file.path(output_dir, "assembly_pairwise_RCbray.tsv")
    )
    
    write_tsv(
      assembly_summary,
      file.path(output_dir, "assembly_process_summary.tsv")
    )
    
    p_assembly <- assembly_summary %>%
      mutate(
        depth = factor(depth, levels = depth_order),
        domain = factor(domain, levels = domain_keep),
        assembly_process = factor(
          assembly_process,
          levels = c(
            "Dispersal limitation",
            "Homogenizing dispersal",
            "Drift or weak selection",
            "NA"
          )
        )
      ) %>%
      ggplot(aes(x = depth, y = percentage, fill = assembly_process)) +
      geom_col(width = 0.75, color = "black", linewidth = 0.2) +
      facet_wrap(~ domain, nrow = 1) +
      scale_fill_manual(
        values = c(
          "Dispersal limitation" = "#D95F02",
          "Homogenizing dispersal" = "#1B9E77",
          "Drift or weak selection" = "grey70",
          "NA" = "grey90"
        ),
        name = "Assembly process"
      ) +
      scale_y_continuous(
        limits = c(0, 100),
        expand = expansion(mult = c(0, 0.02))
      ) +
      labs(
        title = "Community assembly process along sediment depth",
        x = "Sediment depth",
        y = "Pairwise proportion (%)"
      ) +
      theme_nature(base_size = 8) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      )
    
    ggsave(
      file.path(output_dir, "community_assembly_process.pdf"),
      p_assembly,
      width = 180,
      height = 90,
      units = "mm",
      device = "pdf"
    )
    
    ggsave(
      file.path(output_dir, "community_assembly_process.png"),
      p_assembly,
      width = 180,
      height = 90,
      units = "mm",
      dpi = 600
    )
  }
}

# ----------------------------
# 5. Network analysis
# ----------------------------

network_edges_all <- tibble()
network_nodes_all <- tibble()
network_stats_all <- tibble()
network_plots <- list()

if (run_network_analysis) {
  message("[INFO] Running co-occurrence network analysis")
  
  edge_results <- list()
  node_results <- list()
  stat_results <- list()
  plot_results <- list()
  
  idx <- 1
  pidx <- 1
  
  for (dm in domain_keep) {
    for (dp in depth_order) {
      sub <- abund_dat %>%
        filter(domain == dm, depth_group == dp)
      
      samples <- sort(unique(sub$sample_id))
      
      if (length(samples) < network_min_samples) {
        next
      }
      
      feat_info <- sub %>%
        group_by(feature_id, feature_label, domain) %>%
        summarise(
          total_rpm = sum(marker_rpm, na.rm = TRUE),
          detected_samples = n_distinct(sample_id[marker_rpm > 0]),
          prevalence = detected_samples / length(samples),
          .groups = "drop"
        ) %>%
        filter(prevalence >= network_min_prevalence) %>%
        arrange(desc(total_rpm)) %>%
        slice_head(n = network_top_n_taxa_per_domain_depth)
      
      features <- feat_info$feature_id
      
      if (length(features) < 3) {
        next
      }
      
      mat_df <- make_community_matrix(sub, samples, features)
      
      mat <- mat_df %>%
        select(-sample_id) %>%
        as.matrix()
      
      colnames(mat) <- features
      
      edges <- pairwise_cor_network(
        mat = mat,
        method = network_cor_method
      )
      
      if (nrow(edges) > 0) {
        edges <- edges %>%
          mutate(
            q_value = p.adjust(p_value, method = network_p_adjust),
            sign = ifelse(rho > 0, "Positive", "Negative")
          ) %>%
          filter(
            abs(rho) >= network_r_cutoff,
            q_value <= network_q_cutoff
          )
      }
      
      nodes <- feat_info %>%
        transmute(
          name = feature_id,
          label = feature_label,
          domain = as.character(domain),
          total_rpm,
          detected_samples,
          prevalence
        )
      
      edges <- edges %>%
        mutate(
          depth = as.character(dp),
          domain = dm
        )
      
      nodes <- nodes %>%
        mutate(
          depth = as.character(dp),
          domain_group = dm
        )
      
      stats <- network_stats_from_edges(
        edges = edges,
        nodes = nodes,
        depth_value = as.character(dp),
        domain_value = dm
      )
      
      edge_results[[idx]] <- edges
      node_results[[idx]] <- nodes
      stat_results[[idx]] <- stats
      
      p_net <- plot_one_network(
        edges = edges,
        nodes = nodes,
        title_text = paste(dm, dp, sep = " | ")
      )
      
      plot_results[[pidx]] <- p_net
      
      idx <- idx + 1
      pidx <- pidx + 1
    }
  }
  
  network_edges_all <- bind_rows(edge_results)
  network_nodes_all <- bind_rows(node_results)
  network_stats_all <- bind_rows(stat_results)
  
  write_tsv(
    network_edges_all,
    file.path(output_dir, "network_edges.tsv")
  )
  
  write_tsv(
    network_nodes_all,
    file.path(output_dir, "network_nodes.tsv")
  )
  
  write_tsv(
    network_stats_all,
    file.path(output_dir, "network_stats.tsv")
  )
  
  if (length(plot_results) > 0) {
    p_network_all <- wrap_plots(plot_results, ncol = length(domain_keep)) +
      plot_annotation(
        title = "Co-occurrence network structure along sediment depth",
        theme = theme(
          plot.title = element_text(face = "bold", size = 12)
        )
      )
    
    ggsave(
      file.path(output_dir, "network_structure_by_depth.pdf"),
      p_network_all,
      width = 220,
      height = max(120, 45 * ceiling(length(plot_results) / length(domain_keep))),
      units = "mm",
      device = "pdf"
    )
    
    ggsave(
      file.path(output_dir, "network_structure_by_depth.png"),
      p_network_all,
      width = 220,
      height = max(120, 45 * ceiling(length(plot_results) / length(domain_keep))),
      units = "mm",
      dpi = 600
    )
  }
  
  if (nrow(network_stats_all) > 0) {
    p_network_stats <- network_stats_all %>%
      mutate(
        depth = factor(depth, levels = depth_order),
        domain = factor(domain, levels = domain_keep)
      ) %>%
      pivot_longer(
        cols = c(edges, positive_edges, negative_edges, density, mean_degree, modularity),
        names_to = "metric",
        values_to = "value"
      ) %>%
      ggplot(aes(x = depth, y = value, group = domain, color = domain)) +
      geom_line(linewidth = 0.45) +
      geom_point(size = 1.6) +
      facet_wrap(~ metric, scales = "free_y", ncol = 3) +
      scale_color_npg() +
      labs(
        title = "Network properties along sediment depth",
        x = "Sediment depth",
        y = NULL,
        color = "Domain"
      ) +
      theme_nature(base_size = 8) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "bottom"
      )
    
    ggsave(
      file.path(output_dir, "network_statistics_depth_trend.pdf"),
      p_network_stats,
      width = 180,
      height = 110,
      units = "mm",
      device = "pdf"
    )
    
    ggsave(
      file.path(output_dir, "network_statistics_depth_trend.png"),
      p_network_stats,
      width = 180,
      height = 110,
      units = "mm",
      dpi = 600
    )
  }
}

message("[DONE] Results written to: ", output_dir)