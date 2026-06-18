# ============================================================
# Fig_env_microbiome_structure.R
#
# Environmental drivers of bacterial, archaeal and fungal
# alpha diversity and community composition
# ============================================================
#
# Required input file 1: environmental metadata table
# ------------------------------------------------------------
# Default file name:
#   env_metadata.tsv
#
# Required columns:
#   sample_id
#
# Other numeric columns are treated as environmental variables.
#
# Example:
#   sample_id        TN     SOC    AP     pH     CEC    MMT    MMP
#   201704_MF1_0002 1.25   12.3   8.5    6.82   10.1   0.24   0.13
#   201704_MF1_0608 1.10   10.8   7.2    6.75   9.8    0.20   0.10
#
#
# Required input file 2: microbial abundance table
# ------------------------------------------------------------
# Default file name:
#   all.marker_rpm.genus.long.tsv
#
# Required columns:
#   sample_id
#   domain
#   lineage
#   marker_rpm
#
# Recommended columns:
#   taxon_marker_reads
#   rank
#   marker
#
# Notes:
#   1. lineage is used as the unique taxon identifier.
#   2. taxon_marker_reads is preferred for Chao1.
#   3. marker_rpm is used for Bray-Curtis community composition.
#   4. The script supports Bacteria, Archaea and Fungi.
#
# Main outputs:
#   env_microbiome_structure/env_alpha_spearman.tsv
#   env_microbiome_structure/env_community_mantel.tsv
#   env_microbiome_structure/env_microbiome_structure.pdf
#   env_microbiome_structure/env_microbiome_structure.png
#
# ============================================================

# ----------------------------
# 0. Packages
# ----------------------------

packages <- c(
  "tidyverse",
  "data.table",
  "vegan",
  "patchwork",
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
library(patchwork)
library(scales)

# ----------------------------
# 1. Parameters
# ----------------------------

env_file <- "env_metadata.tsv"
abund_file <- "all.marker_rpm.genus.long.tsv"

output_dir <- "env_microbiome_structure"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

domain_keep <- c("Bacteria", "Archaea", "Fungi")

# Use NULL to automatically use all numeric environmental columns.
env_vars <- NULL

# Example:
# env_vars <- c("TN", "SOC", "AP", "pH", "CEC", "MMT", "MMP")

target_rank <- "genus"

min_total_abundance <- 0
min_detected_samples <- 3

alpha_indices <- c("Chao1", "Shannon", "Pielou")

cor_method <- "spearman"
p_adjust_method <- "BH"

mantel_method <- "spearman"
mantel_permutations <- 999

# Plot settings
show_tile_values <- TRUE
show_stars <- TRUE

# ----------------------------
# 2. Theme
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

p_to_star <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ ""
  )
}

safe_cor_test <- function(x, y, method = "spearman") {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]
  
  n <- length(x)
  
  if (n < 4 || length(unique(x)) < 2 || length(unique(y)) < 2) {
    return(tibble(r = NA_real_, p = NA_real_, n = n))
  }
  
  ct <- suppressWarnings(cor.test(x, y, method = method, exact = FALSE))
  
  tibble(
    r = unname(ct$estimate),
    p = ct$p.value,
    n = n
  )
}

# ----------------------------
# 3. Read input files
# ----------------------------

if (!file.exists(env_file)) {
  stop("Environmental metadata table not found: ", env_file)
}

if (!file.exists(abund_file)) {
  stop("Microbial abundance table not found: ", abund_file)
}

env_raw <- data.table::fread(
  env_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

names(env_raw) <- trimws(names(env_raw))

if (!"sample_id" %in% names(env_raw)) {
  stop("Environmental table must contain sample_id column.")
}

env_raw <- env_raw %>%
  mutate(sample_id = as.character(sample_id))

if (is.null(env_vars)) {
  env_vars <- env_raw %>%
    select(-sample_id) %>%
    select(where(is.numeric)) %>%
    names()
}

if (length(env_vars) == 0) {
  stop("No numeric environmental variables found.")
}

missing_env_vars <- setdiff(env_vars, names(env_raw))
if (length(missing_env_vars) > 0) {
  stop("Missing environmental variables: ", paste(missing_env_vars, collapse = ", "))
}

env_dat <- env_raw %>%
  select(sample_id, all_of(env_vars)) %>%
  mutate(across(all_of(env_vars), as.numeric))

abund_raw <- data.table::fread(
  abund_file,
  sep = "\t",
  header = TRUE,
  data.table = FALSE,
  check.names = FALSE
)

names(abund_raw) <- trimws(names(abund_raw))

required_abund_cols <- c("sample_id", "domain", "lineage", "marker_rpm")
missing_abund_cols <- setdiff(required_abund_cols, names(abund_raw))

if (length(missing_abund_cols) > 0) {
  stop("Missing required columns in abundance table: ", paste(missing_abund_cols, collapse = ", "))
}

abund_dat <- abund_raw %>%
  mutate(
    sample_id = as.character(sample_id),
    domain = as.character(domain),
    lineage = as.character(lineage),
    marker_rpm = as.numeric(marker_rpm)
  ) %>%
  filter(domain %in% domain_keep)

if ("rank" %in% names(abund_dat)) {
  abund_dat <- abund_dat %>%
    filter(rank == target_rank)
}

if (!"taxon_marker_reads" %in% names(abund_dat)) {
  message("[WARN] taxon_marker_reads column not found. Chao1 will be calculated from marker_rpm, which is less ideal.")
  abund_dat <- abund_dat %>%
    mutate(taxon_marker_reads = marker_rpm)
} else {
  abund_dat <- abund_dat %>%
    mutate(taxon_marker_reads = as.numeric(taxon_marker_reads))
}

common_samples <- intersect(env_dat$sample_id, abund_dat$sample_id)

if (length(common_samples) < 4) {
  stop("Too few shared samples between environmental table and abundance table.")
}

env_dat <- env_dat %>%
  filter(sample_id %in% common_samples)

abund_dat <- abund_dat %>%
  filter(sample_id %in% common_samples)

# ----------------------------
# 4. Build domain-specific abundance matrices
# ----------------------------

make_domain_matrix <- function(df, domain_name, value_col = "marker_rpm") {
  sub <- df %>%
    filter(domain == domain_name) %>%
    group_by(sample_id, lineage) %>%
    summarise(value = sum(.data[[value_col]], na.rm = TRUE), .groups = "drop")
  
  taxa_keep <- sub %>%
    group_by(lineage) %>%
    summarise(
      total_abundance = sum(value, na.rm = TRUE),
      detected_samples = n_distinct(sample_id[value > 0]),
      .groups = "drop"
    ) %>%
    filter(
      total_abundance > min_total_abundance,
      detected_samples >= min_detected_samples
    ) %>%
    pull(lineage)
  
  sub %>%
    filter(lineage %in% taxa_keep) %>%
    complete(
      sample_id = common_samples,
      lineage = taxa_keep,
      fill = list(value = 0)
    ) %>%
    pivot_wider(
      names_from = lineage,
      values_from = value,
      values_fill = 0
    ) %>%
    arrange(match(sample_id, common_samples))
}

rpm_mats <- setNames(
  lapply(domain_keep, function(d) make_domain_matrix(abund_dat, d, "marker_rpm")),
  domain_keep
)

count_mats <- setNames(
  lapply(domain_keep, function(d) make_domain_matrix(abund_dat, d, "taxon_marker_reads")),
  domain_keep
)

# ----------------------------
# 5. Calculate alpha diversity
# ----------------------------

calc_alpha_one_domain <- function(mat_df, domain_name) {
  sample_ids <- mat_df$sample_id
  mat <- mat_df %>%
    select(-sample_id) %>%
    as.matrix()
  
  mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  
  if (ncol(mat) == 0) {
    return(tibble())
  }
  
  richness <- specnumber(mat)
  shannon <- diversity(mat, index = "shannon")
  pielou <- ifelse(richness > 1, shannon / log(richness), NA_real_)
  
  chao1 <- rep(NA_real_, nrow(mat))
  
  for (i in seq_len(nrow(mat))) {
    xi <- mat[i, ]
    xi <- xi[xi > 0]
    if (length(xi) > 0) {
      est <- suppressWarnings(vegan::estimateR(round(xi)))
      chao1[i] <- unname(est["S.chao1"])
    }
  }
  
  tibble(
    sample_id = sample_ids,
    domain = domain_name,
    Chao1 = chao1,
    Shannon = shannon,
    Pielou = pielou
  )
}

alpha_df <- map2_dfr(
  count_mats,
  names(count_mats),
  calc_alpha_one_domain
)

write_tsv(
  alpha_df,
  file.path(output_dir, "alpha_diversity_by_domain.tsv")
)

# ----------------------------
# 6. Spearman correlation between environment and alpha diversity
# ----------------------------

alpha_env <- alpha_df %>%
  left_join(env_dat, by = "sample_id")

alpha_cor <- expand_grid(
  domain = domain_keep,
  env_var = env_vars,
  alpha_index = alpha_indices
) %>%
  mutate(result = pmap(
    list(domain, env_var, alpha_index),
    function(domain, env_var, alpha_index) {
      sub <- alpha_env %>% filter(domain == !!domain)
      safe_cor_test(
        x = sub[[env_var]],
        y = sub[[alpha_index]],
        method = cor_method
      )
    }
  )) %>%
  unnest(result) %>%
  group_by(domain) %>%
  mutate(q = p.adjust(p, method = p_adjust_method)) %>%
  ungroup() %>%
  mutate(
    star = p_to_star(p),
    label = ifelse(
      is.na(r),
      "",
      paste0(sprintf("%.2f", r), ifelse(star == "", "", paste0("\n", star)))
    )
  )

write_tsv(
  alpha_cor,
  file.path(output_dir, "env_alpha_spearman.tsv")
)

# ----------------------------
# 7. Mantel test between environment and community composition
# ----------------------------

mantel_one <- function(domain_name, env_var) {
  mat_df <- rpm_mats[[domain_name]]
  
  sample_ids <- mat_df$sample_id
  
  mat <- mat_df %>%
    select(-sample_id) %>%
    as.matrix()
  
  mode(mat) <- "numeric"
  mat[is.na(mat)] <- 0
  
  env_sub <- env_dat %>%
    filter(sample_id %in% sample_ids) %>%
    arrange(match(sample_id, sample_ids))
  
  keep <- is.finite(env_sub[[env_var]])
  
  mat <- mat[keep, , drop = FALSE]
  x <- env_sub[[env_var]][keep]
  
  if (length(x) < 4 || ncol(mat) == 0 || length(unique(x)) < 2) {
    return(tibble(
      domain = domain_name,
      env_var = env_var,
      mantel_r = NA_real_,
      mantel_p = NA_real_,
      n = length(x)
    ))
  }
  
  community_dist <- vegan::vegdist(mat, method = "bray")
  env_dist <- dist(scale(x), method = "euclidean")
  
  mt <- vegan::mantel(
    community_dist,
    env_dist,
    method = mantel_method,
    permutations = mantel_permutations
  )
  
  tibble(
    domain = domain_name,
    env_var = env_var,
    mantel_r = unname(mt$statistic),
    mantel_p = mt$signif,
    n = length(x)
  )
}

mantel_df <- expand_grid(
  domain = domain_keep,
  env_var = env_vars
) %>%
  mutate(result = map2(domain, env_var, mantel_one)) %>%
  unnest(result) %>%
  mutate(
    p_group = case_when(
      is.na(mantel_p) ~ "NA",
      mantel_p < 0.01 ~ "< 0.01",
      mantel_p < 0.05 ~ "0.01 - 0.05",
      TRUE ~ ">= 0.05"
    ),
    r_group = case_when(
      is.na(mantel_r) ~ "NA",
      mantel_r < 0.2 ~ "< 0.2",
      mantel_r < 0.4 ~ "0.2 - 0.4",
      TRUE ~ "> 0.4"
    )
  )

write_tsv(
  mantel_df,
  file.path(output_dir, "env_community_mantel.tsv")
)

# ----------------------------
# 8. Heatmap panel
# ----------------------------

env_vars_factor <- rev(env_vars)

alpha_cor <- alpha_cor %>%
  mutate(
    env_var = factor(env_var, levels = env_vars_factor),
    alpha_index = factor(alpha_index, levels = alpha_indices),
    domain = factor(domain, levels = domain_keep)
  )

p_heatmap <- ggplot(
  alpha_cor,
  aes(x = alpha_index, y = env_var, fill = r)
) +
  geom_tile(color = "grey85", linewidth = 0.25) +
  {
    if (show_tile_values) {
      geom_text(aes(label = label), size = 2.7, lineheight = 0.75)
    }
  } +
  facet_grid(. ~ domain) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman's R"
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Environmental variables and alpha diversity"
  ) +
  theme_nature(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position = "left",
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.25)
  )

# ----------------------------
# 9. Mantel network panel
# ----------------------------

env_nodes <- tibble(
  node = env_vars,
  type = "environment",
  x = 0,
  y = seq_along(env_vars_factor)[match(env_vars, env_vars_factor)]
)

domain_nodes <- tibble(
  domain = domain_keep,
  node = paste0(domain_keep, "\ncommunity\ncomposition"),
  type = "community",
  x = 1,
  y = seq(
    from = max(env_nodes$y) - 0.6,
    to = min(env_nodes$y) + 0.6,
    length.out = length(domain_keep)
  )
)

edge_df <- mantel_df %>%
  left_join(env_nodes %>% select(env_var = node, x_start = x, y_start = y), by = "env_var") %>%
  left_join(domain_nodes %>% select(domain, x_end = x, y_end = y), by = "domain") %>%
  mutate(
    p_group = factor(p_group, levels = c("< 0.01", "0.01 - 0.05", ">= 0.05", "NA")),
    r_group = factor(r_group, levels = c("< 0.2", "0.2 - 0.4", "> 0.4", "NA")),
    edge_alpha = ifelse(is.na(mantel_p), 0.15, ifelse(mantel_p < 0.05, 0.95, 0.35))
  )

p_mantel <- ggplot() +
  geom_curve(
    data = edge_df,
    aes(
      x = x_start,
      y = y_start,
      xend = x_end,
      yend = y_end,
      color = p_group,
      linewidth = r_group,
      alpha = edge_alpha
    ),
    curvature = 0.18,
    lineend = "round"
  ) +
  geom_text(
    data = env_nodes,
    aes(x = x - 0.03, y = y, label = node),
    hjust = 1,
    size = 3.2
  ) +
  geom_text(
    data = domain_nodes,
    aes(x = x + 0.03, y = y, label = node),
    hjust = 0,
    size = 3.2
  ) +
  scale_color_manual(
    values = c(
      "< 0.01" = "#D95F02",
      "0.01 - 0.05" = "#1B9E77",
      ">= 0.05" = "grey75",
      "NA" = "grey90"
    ),
    name = "Mantel's P"
  ) +
  scale_linewidth_manual(
    values = c(
      "< 0.2" = 0.35,
      "0.2 - 0.4" = 0.9,
      "> 0.4" = 2.0,
      "NA" = 0.25
    ),
    name = "Mantel's R"
  ) +
  scale_alpha_identity() +
  coord_cartesian(
    xlim = c(-0.25, 1.35),
    ylim = range(env_nodes$y) + c(-0.6, 0.6),
    clip = "off"
  ) +
  labs(title = "Environmental variables and community composition") +
  theme_void(base_size = 8) +
  theme(
    plot.title = element_text(face = "bold", size = 10, hjust = 0),
    legend.position = "right",
    plot.margin = margin(8, 20, 8, 8)
  )

# ----------------------------
# 10. Combined figure
# ----------------------------

final_fig <- p_heatmap + p_mantel +
  plot_layout(widths = c(1.45, 1.05), guides = "collect") +
  plot_annotation(
    tag_levels = "a",
    theme = theme(
      plot.tag = element_text(face = "bold", size = 12),
      legend.position = "right"
    )
  )

final_fig

ggsave(
  file.path(output_dir, "env_microbiome_structure.pdf"),
  final_fig,
  width = 210,
  height = 95,
  units = "mm",
  device = cairo_pdf
)

ggsave(
  file.path(output_dir, "env_microbiome_structure.png"),
  final_fig,
  width = 210,
  height = 95,
  units = "mm",
  dpi = 600
)

message("[DONE] Results written to: ", output_dir)