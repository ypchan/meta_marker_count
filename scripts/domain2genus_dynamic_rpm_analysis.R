#!/usr/bin/env Rscript

# ============================================================
# domain_genus_dynamic_rpm_analysis.R
#
# Purpose:
#   Analyze domain-to-genus temporal dynamics using RPM-normalized
#   abundances across fixed stations, habitat groups and depth layers.
#
# Main questions:
#   1. Which taxa show strong temporal instability?
#   2. Which time intervals are major disturbance points?
#   3. Which domain, site group and depth layer are most affected?
#   4. Which genus-level taxa drive the observed changes?
#
# Input:
#   Long-format TSV table with columns:
#     sample_id, date, site_id, site_group, depth_bin, domain, genus, rpm
#
# Output:
#   TSV tables and PDF figures.
#
# Author:
#   For mangrove sediment domain/genus RPM dynamics
# ============================================================


# ------------------------------------------------------------
# 0. Packages
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(ggsci)
  library(patchwork)
  library(scales)
})


# ------------------------------------------------------------
# 1. User parameters
# ------------------------------------------------------------

# Input file
input_file <- "domain_genus_rpm.tsv"

# Output directory
outdir <- "domain_genus_dynamic_results"

# Site and depth order
site_levels <- c("MF1", "MF2", "MG1", "MG2", "TD1", "TD2")
site_group_levels <- c("MF", "MG", "TD")

depth_levels <- c(
  "00-10", "10-20", "20-30",
  "30-40", "40-50", "50-60"
)

domain_levels <- c("Bacteria", "Archaea", "Fungi")

# Filtering parameters
# A genus must appear in at least this proportion of samples.
min_prevalence <- 0.05

# A genus must have at least this mean RPM across all samples.
min_mean_rpm <- 0.10

# log2FC threshold used to define obvious shock.
# abs(log2FC) >= 2 means at least 4-fold change.
shock_log2fc <- 2

# Strong shock threshold.
# abs(log2FC) >= 3 means at least 8-fold change.
strong_shock_log2fc <- 3

# Pseudo-count strategy for RPM log2FC.
# If TRUE, pseudo = 1% quantile of non-zero RPM.
use_data_driven_pseudo <- TRUE

# If the data-driven pseudo fails, use this value.
fallback_pseudo <- 0.1

# Number of top taxa to show in plots
top_n_taxa <- 20

# Figure size
fig_width <- 12
fig_height <- 8


# ------------------------------------------------------------
# 2. Helper functions
# ------------------------------------------------------------

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

message2 <- function(...) {
  message(format(Sys.time(), "[%Y-%m-%d %H:%M:%S] "), ...)
}

safe_filename <- function(x) {
  x %>%
    stringr::str_replace_all("[^A-Za-z0-9_\\-]+", "_") %>%
    stringr::str_replace_all("_+", "_") %>%
    stringr::str_replace_all("^_|_$", "")
}

parse_sample_date <- function(x) {
  # This function accepts:
  #   2023-03-01
  #   2023-03
  #   202303
  #   2023/03/01
  x <- as.character(x)
  
  out <- suppressWarnings(as.Date(x))
  
  need_parse <- is.na(out)
  
  if (any(need_parse)) {
    x2 <- x[need_parse]
    
    # Format: YYYY-MM
    idx_ym_dash <- grepl("^\\d{4}-\\d{2}$", x2)
    if (any(idx_ym_dash)) {
      out[need_parse][idx_ym_dash] <- as.Date(paste0(x2[idx_ym_dash], "-01"))
    }
    
    # Format: YYYYMM
    idx_ym <- grepl("^\\d{6}$", x2)
    if (any(idx_ym)) {
      out[need_parse][idx_ym] <- as.Date(
        paste0(substr(x2[idx_ym], 1, 4), "-", substr(x2[idx_ym], 5, 6), "-01")
      )
    }
    
    # Format: YYYY/MM/DD
    idx_slash <- grepl("^\\d{4}/\\d{2}/\\d{2}$", x2)
    if (any(idx_slash)) {
      out[need_parse][idx_slash] <- as.Date(gsub("/", "-", x2[idx_slash]))
    }
  }
  
  out
}

check_required_columns <- function(df, required_cols) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(
      "Missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
}


# ------------------------------------------------------------
# 3. Read and check input
# ------------------------------------------------------------

message2("Reading input: ", input_file)

taxa_raw <- readr::read_tsv(input_file, show_col_types = FALSE)

required_cols <- c(
  "sample_id", "date", "site_id", "site_group",
  "depth_bin", "domain", "genus", "rpm"
)

check_required_columns(taxa_raw, required_cols)

# Basic standardization
taxa <- taxa_raw %>%
  mutate(
    sample_id = as.character(sample_id),
    date = parse_sample_date(date),
    year_month = format(date, "%Y-%m"),
    
    site_id = as.character(site_id),
    site_group = as.character(site_group),
    depth_bin = as.character(depth_bin),
    domain = as.character(domain),
    genus = as.character(genus),
    
    genus = if_else(
      is.na(genus) | genus == "" | genus == "NA",
      "Unclassified_genus",
      genus
    ),
    
    rpm = as.numeric(rpm),
    rpm = replace_na(rpm, 0),
    rpm = if_else(rpm < 0, 0, rpm)
  )

if (any(is.na(taxa$date))) {
  stop(
    "Some dates cannot be parsed. Please use YYYY-MM-DD, YYYY-MM, or YYYYMM.",
    call. = FALSE
  )
}

# Apply factor order
taxa <- taxa %>%
  mutate(
    site_group = factor(site_group, levels = site_group_levels),
    site_id = factor(site_id, levels = site_levels),
    depth_bin = factor(depth_bin, levels = depth_levels),
    domain = factor(domain, levels = domain_levels),
    taxon_id = paste(domain, genus, sep = "|")
  )

# Warn if unknown levels exist
if (any(is.na(taxa$site_group))) {
  warning("Some site_group values are not in site_group_levels.")
}
if (any(is.na(taxa$site_id))) {
  warning("Some site_id values are not in site_levels.")
}
if (any(is.na(taxa$depth_bin))) {
  warning("Some depth_bin values are not in depth_levels.")
}
if (any(is.na(taxa$domain))) {
  warning("Some domain values are not in domain_levels.")
}

message2("Input rows: ", nrow(taxa))
message2("Samples: ", n_distinct(taxa$sample_id))
message2("Taxa: ", n_distinct(taxa$taxon_id))
message2("Dates: ", paste(sort(unique(taxa$year_month)), collapse = ", "))


# ------------------------------------------------------------
# 4. Aggregate duplicated rows
# ------------------------------------------------------------
# If the same sample_id + domain + genus appears multiple times,
# sum RPM to avoid duplicate inflation.

taxa_sample <- taxa %>%
  group_by(
    sample_id, date, year_month,
    site_group, site_id, depth_bin,
    domain, genus, taxon_id
  ) %>%
  summarise(
    rpm = sum(rpm, na.rm = TRUE),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 5. Build sample metadata and zero-fill selected taxa
# ------------------------------------------------------------
# Important:
#   Many abundance tables only include detected taxa.
#   For temporal log2FC, undetected taxa should be treated as RPM = 0,
#   but only within actually sampled samples.
#
# Strategy:
#   1. Use observed sample_id rows as the true sampling design.
#   2. Calculate prevalence and mean RPM using sample number as denominator.
#   3. Keep non-rare taxa.
#   4. Complete sample × taxon table with RPM = 0.

sample_meta <- taxa_sample %>%
  distinct(sample_id, date, year_month, site_group, site_id, depth_bin)

n_samples <- nrow(sample_meta)

taxa_pre_filter <- taxa_sample %>%
  group_by(domain, genus, taxon_id) %>%
  summarise(
    positive_samples = n_distinct(sample_id[rpm > 0]),
    prevalence = positive_samples / n_samples,
    total_rpm = sum(rpm, na.rm = TRUE),
    mean_rpm = total_rpm / n_samples,
    max_rpm = max(rpm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(total_rpm))

readr::write_tsv(
  taxa_pre_filter,
  file.path(outdir, "taxa_pre_filter_summary.tsv")
)

taxa_keep <- taxa_pre_filter %>%
  filter(
    prevalence >= min_prevalence,
    mean_rpm >= min_mean_rpm
  )

message2("Taxa before filtering: ", nrow(taxa_pre_filter))
message2("Taxa retained for dynamics: ", nrow(taxa_keep))

readr::write_tsv(
  taxa_keep,
  file.path(outdir, "taxa_retained_for_dynamics.tsv")
)

# Complete selected taxa across observed samples
taxa_complete <- sample_meta %>%
  tidyr::crossing(
    taxa_keep %>% select(domain, genus, taxon_id)
  ) %>%
  left_join(
    taxa_sample %>% select(sample_id, taxon_id, rpm),
    by = c("sample_id", "taxon_id")
  ) %>%
  mutate(
    rpm = replace_na(rpm, 0),
    log1p_rpm = log1p(rpm)
  )


# ------------------------------------------------------------
# 6. Choose pseudo value for log2FC
# ------------------------------------------------------------

nonzero_rpm <- taxa_complete$rpm[taxa_complete$rpm > 0]

if (use_data_driven_pseudo && length(nonzero_rpm) > 10) {
  pseudo <- as.numeric(quantile(nonzero_rpm, probs = 0.01, na.rm = TRUE))
  if (is.na(pseudo) || pseudo <= 0) {
    pseudo <- fallback_pseudo
  }
} else {
  pseudo <- fallback_pseudo
}

message2("Pseudo value for log2FC: ", signif(pseudo, 4))

writeLines(
  paste0("pseudo\t", pseudo),
  con = file.path(outdir, "pseudo_value.tsv")
)


# ------------------------------------------------------------
# 7. Domain-level total RPM dynamics
# ------------------------------------------------------------

domain_total <- taxa_sample %>%
  group_by(
    sample_id, date, year_month,
    site_group, site_id, depth_bin,
    domain
  ) %>%
  summarise(
    total_rpm = sum(rpm, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(
  domain_total,
  file.path(outdir, "domain_total_rpm_by_sample.tsv")
)

# Aggregated by site_group and depth
domain_total_group_depth <- domain_total %>%
  group_by(date, year_month, site_group, depth_bin, domain) %>%
  summarise(
    mean_total_rpm = mean(total_rpm, na.rm = TRUE),
    sd_total_rpm = sd(total_rpm, na.rm = TRUE),
    n_sites = n_distinct(site_id),
    .groups = "drop"
  )

readr::write_tsv(
  domain_total_group_depth,
  file.path(outdir, "domain_total_rpm_by_group_depth.tsv")
)


# ------------------------------------------------------------
# 8. Calculate adjacent-time log2FC using RPM
# ------------------------------------------------------------
# For each taxon in each fixed unit:
#   fixed unit = site_id × depth_bin
#
# log2FC compares adjacent sampled time points:
#   log2FC = log2((RPM_t + pseudo) / (RPM_previous + pseudo))
#
# Positive log2FC: increase
# Negative log2FC: decrease
# abs(log2FC): change intensity

taxa_fc <- taxa_complete %>%
  arrange(taxon_id, site_id, depth_bin, date) %>%
  group_by(taxon_id, domain, genus, site_group, site_id, depth_bin) %>%
  mutate(
    prev_date = lag(date),
    prev_year_month = lag(year_month),
    prev_rpm = lag(rpm),
    
    log2fc = log2((rpm + pseudo) / (prev_rpm + pseudo)),
    abs_log2fc = abs(log2fc),
    
    direction = case_when(
      is.na(log2fc) ~ NA_character_,
      log2fc > 0 ~ "increase",
      log2fc < 0 ~ "decrease",
      TRUE ~ "stable"
    ),
    
    interval = if_else(
      is.na(prev_year_month),
      NA_character_,
      paste0(prev_year_month, "_to_", year_month)
    ),
    
    is_shock = abs_log2fc >= shock_log2fc,
    is_strong_shock = abs_log2fc >= strong_shock_log2fc
  ) %>%
  ungroup() %>%
  filter(!is.na(log2fc))

readr::write_tsv(
  taxa_fc,
  file.path(outdir, "taxa_adjacent_time_log2fc.tsv")
)


# ------------------------------------------------------------
# 9. Taxon-level volatility ranking
# ------------------------------------------------------------
# Main metrics:
#
# prevalence:
#   fraction of samples where the taxon was detected.
#
# mean_rpm:
#   average RPM across all samples.
#
# volatility:
#   median absolute log2FC.
#   This measures typical temporal instability.
#
# max_shock:
#   maximum absolute log2FC.
#   This identifies extreme one-time jumps.
#
# affected_units:
#   number of site-depth-time observations with abs(log2FC) >= shock threshold.
#
# affected_ratio:
#   affected_units / total comparisons.
#
# strong_shock_units:
#   number of observations with abs(log2FC) >= strong threshold.

taxa_volatility <- taxa_fc %>%
  group_by(domain, genus, taxon_id) %>%
  summarise(
    mean_rpm = mean(rpm, na.rm = TRUE),
    median_rpm = median(rpm, na.rm = TRUE),
    max_rpm = max(rpm, na.rm = TRUE),
    
    prevalence = mean(rpm > 0, na.rm = TRUE),
    
    volatility = median(abs_log2fc, na.rm = TRUE),
    mean_abs_log2fc = mean(abs_log2fc, na.rm = TRUE),
    max_shock = max(abs_log2fc, na.rm = TRUE),
    
    affected_units = sum(is_shock, na.rm = TRUE),
    strong_shock_units = sum(is_strong_shock, na.rm = TRUE),
    n_comparisons = n(),
    affected_ratio = affected_units / n_comparisons,
    
    n_increase = sum(log2fc >= shock_log2fc, na.rm = TRUE),
    n_decrease = sum(log2fc <= -shock_log2fc, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(desc(volatility), desc(max_shock), desc(affected_ratio))

readr::write_tsv(
  taxa_volatility,
  file.path(outdir, "taxa_volatility_ranking.tsv")
)


# ------------------------------------------------------------
# 10. Directional synchrony
# ------------------------------------------------------------
# This asks:
#   when a taxon changes strongly, does it increase/decrease
#   synchronously across many site-depth units?
#
# sync_ratio close to 1:
#   changes are mostly in the same direction.
#
# sync_ratio close to 0.5:
#   mixed increase/decrease, less coherent.

taxon_interval_sync <- taxa_fc %>%
  filter(is_shock) %>%
  group_by(domain, genus, taxon_id, interval) %>%
  summarise(
    shocked_units = n(),
    up_units = sum(log2fc > 0, na.rm = TRUE),
    down_units = sum(log2fc < 0, na.rm = TRUE),
    sync_ratio = pmax(up_units, down_units) / shocked_units,
    dominant_direction = if_else(up_units >= down_units, "increase", "decrease"),
    median_log2fc = median(log2fc, na.rm = TRUE),
    median_abs_log2fc = median(abs_log2fc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(shocked_units), desc(sync_ratio), desc(median_abs_log2fc))

readr::write_tsv(
  taxon_interval_sync,
  file.path(outdir, "taxon_interval_directional_synchrony.tsv")
)

taxa_volatility2 <- taxa_volatility %>%
  left_join(
    taxon_interval_sync %>%
      group_by(taxon_id) %>%
      slice_max(
        order_by = shocked_units * sync_ratio * median_abs_log2fc,
        n = 1,
        with_ties = FALSE
      ) %>%
      ungroup() %>%
      select(
        taxon_id,
        most_synchronous_interval = interval,
        sync_shocked_units = shocked_units,
        sync_ratio,
        dominant_direction
      ),
    by = "taxon_id"
  )

readr::write_tsv(
  taxa_volatility2,
  file.path(outdir, "taxa_volatility_ranking_with_synchrony.tsv")
)


# ------------------------------------------------------------
# 11. Event score by time interval
# ------------------------------------------------------------
# Event score identifies the strongest disturbance time intervals.
#
# event_score:
#   median abs(log2FC) across taxa and site-depth units.
#
# shocked_taxa:
#   number of taxa with obvious changes.
#
# affected_units:
#   number of taxon × site × depth comparisons with obvious changes.

global_event_score <- taxa_fc %>%
  group_by(domain, interval) %>%
  summarise(
    event_score = median(abs_log2fc, na.rm = TRUE),
    mean_abs_log2fc = mean(abs_log2fc, na.rm = TRUE),
    q90_abs_log2fc = quantile(abs_log2fc, 0.90, na.rm = TRUE),
    
    shocked_taxa = n_distinct(taxon_id[is_shock]),
    strong_shocked_taxa = n_distinct(taxon_id[is_strong_shock]),
    
    affected_units = sum(is_shock, na.rm = TRUE),
    strong_affected_units = sum(is_strong_shock, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(desc(event_score))

readr::write_tsv(
  global_event_score,
  file.path(outdir, "global_event_score_by_domain_interval.tsv")
)


event_score_group_depth <- taxa_fc %>%
  group_by(domain, interval, site_group, depth_bin) %>%
  summarise(
    event_score = median(abs_log2fc, na.rm = TRUE),
    mean_abs_log2fc = mean(abs_log2fc, na.rm = TRUE),
    q90_abs_log2fc = quantile(abs_log2fc, 0.90, na.rm = TRUE),
    
    shocked_taxa = n_distinct(taxon_id[is_shock]),
    affected_units = sum(is_shock, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  arrange(desc(event_score))

readr::write_tsv(
  event_score_group_depth,
  file.path(outdir, "event_score_by_domain_interval_group_depth.tsv")
)


# ------------------------------------------------------------
# 12. Taxa contributing to each event interval
# ------------------------------------------------------------
# This table tells you:
#   during each interval, which taxa changed most strongly?

taxon_event_contribution <- taxa_fc %>%
  group_by(domain, genus, taxon_id, interval) %>%
  summarise(
    median_log2fc = median(log2fc, na.rm = TRUE),
    median_abs_log2fc = median(abs_log2fc, na.rm = TRUE),
    max_abs_log2fc = max(abs_log2fc, na.rm = TRUE),
    
    shocked_units = sum(is_shock, na.rm = TRUE),
    strong_shocked_units = sum(is_strong_shock, na.rm = TRUE),
    n_units = n(),
    shocked_ratio = shocked_units / n_units,
    
    .groups = "drop"
  ) %>%
  filter(shocked_units > 0) %>%
  arrange(interval, domain, desc(shocked_units), desc(median_abs_log2fc))

readr::write_tsv(
  taxon_event_contribution,
  file.path(outdir, "taxon_contribution_by_event_interval.tsv")
)


# ------------------------------------------------------------
# 13. Candidate dynamic taxa
# ------------------------------------------------------------
# Candidate taxa are taxa that are:
#   not too rare
#   have clear shocks
#   affect multiple site-depth-time units

candidate_taxa <- taxa_volatility2 %>%
  filter(
    prevalence >= min_prevalence,
    mean_rpm >= min_mean_rpm,
    max_shock >= shock_log2fc,
    affected_ratio >= 0.05
  ) %>%
  arrange(desc(max_shock), desc(volatility), desc(affected_ratio))

readr::write_tsv(
  candidate_taxa,
  file.path(outdir, "candidate_dynamic_taxa.tsv")
)

candidate_taxa_strict <- taxa_volatility2 %>%
  filter(
    prevalence >= 0.10,
    mean_rpm >= 1,
    max_shock >= strong_shock_log2fc,
    affected_ratio >= 0.10
  ) %>%
  arrange(desc(max_shock), desc(volatility), desc(affected_ratio))

readr::write_tsv(
  candidate_taxa_strict,
  file.path(outdir, "candidate_dynamic_taxa_strict.tsv")
)


# ------------------------------------------------------------
# 14. Plot style
# ------------------------------------------------------------

theme_nature_like <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90", linewidth = 0.25),
      panel.border = element_rect(color = "grey75", linewidth = 0.4),
      strip.background = element_rect(fill = "grey95", color = "grey75", linewidth = 0.4),
      strip.text = element_text(color = "black", face = "bold"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(color = "black"),
      legend.key = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0),
      plot.subtitle = element_text(color = "grey30")
    )
}


# ------------------------------------------------------------
# 15. Figure 1: domain total dynamics
# ------------------------------------------------------------

p_domain_total <- ggplot(
  domain_total_group_depth,
  aes(x = date, y = mean_total_rpm, color = domain, group = domain)
) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.8) +
  facet_grid(depth_bin ~ site_group, scales = "free_y") +
  scale_color_npg() +
  scale_y_continuous(trans = "log1p") +
  theme_nature_like(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7)
  ) +
  labs(
    title = "Domain-level temporal dynamics",
    subtitle = "Mean total RPM by habitat group and depth layer",
    x = NULL,
    y = "Mean total RPM, log1p scale",
    color = "Domain"
  )

ggsave(
  file.path(outdir, "Fig1_domain_total_dynamics.pdf"),
  p_domain_total,
  width = fig_width,
  height = fig_height
)


# ------------------------------------------------------------
# 16. Figure 2: top volatile taxa
# ------------------------------------------------------------

plot_top_volatile <- taxa_volatility2 %>%
  group_by(domain) %>%
  slice_max(order_by = volatility, n = top_n_taxa, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    genus_label = paste(domain, genus, sep = " | "),
    genus_label = fct_reorder(genus_label, volatility)
  )

p_top_volatile <- ggplot(
  plot_top_volatile,
  aes(x = volatility, y = genus_label)
) +
  geom_point(
    aes(size = affected_ratio, color = domain),
    alpha = 0.85
  ) +
  scale_color_npg() +
  scale_size_continuous(range = c(1.8, 7)) +
  theme_nature_like(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 7),
    panel.grid.major.y = element_blank()
  ) +
  labs(
    title = "Most volatile genus-level taxa",
    subtitle = "Volatility is median absolute adjacent-time log2FC",
    x = "Volatility score: median |log2FC|",
    y = NULL,
    color = "Domain",
    size = "Affected ratio"
  )

ggsave(
  file.path(outdir, "Fig2_top_volatile_taxa.pdf"),
  p_top_volatile,
  width = 9,
  height = 11
)


# ------------------------------------------------------------
# 17. Figure 3: global event score
# ------------------------------------------------------------

# Keep interval order by time
interval_order <- taxa_fc %>%
  distinct(prev_date, date, interval) %>%
  arrange(prev_date, date) %>%
  pull(interval)

global_event_score <- global_event_score %>%
  mutate(interval = factor(interval, levels = interval_order))

p_global_event <- ggplot(
  global_event_score,
  aes(x = interval, y = event_score, color = domain, group = domain)
) +
  geom_line(linewidth = 0.8) +
  geom_point(aes(size = shocked_taxa), alpha = 0.9) +
  scale_color_npg() +
  theme_nature_like(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(
    title = "Disturbance intensity by time interval",
    subtitle = "Event score is median absolute adjacent-time log2FC",
    x = NULL,
    y = "Event score",
    color = "Domain",
    size = "Shocked taxa"
  )

ggsave(
  file.path(outdir, "Fig3_global_event_score.pdf"),
  p_global_event,
  width = 9,
  height = 5
)


# ------------------------------------------------------------
# 18. Figure 4: event score by site group and depth
# ------------------------------------------------------------

event_score_group_depth <- event_score_group_depth %>%
  mutate(interval = factor(interval, levels = interval_order))

p_event_group_depth <- ggplot(
  event_score_group_depth,
  aes(x = interval, y = event_score, color = domain, group = domain)
) +
  geom_line(linewidth = 0.6) +
  geom_point(aes(size = shocked_taxa), alpha = 0.85) +
  facet_grid(depth_bin ~ site_group, scales = "free_y") +
  scale_color_npg() +
  theme_nature_like(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
  ) +
  labs(
    title = "Disturbance intensity across habitat groups and depth layers",
    subtitle = "Event score separated by domain, habitat group, and depth",
    x = NULL,
    y = "Event score",
    color = "Domain",
    size = "Shocked taxa"
  )

ggsave(
  file.path(outdir, "Fig4_event_score_group_depth.pdf"),
  p_event_group_depth,
  width = 13,
  height = 9
)


# ------------------------------------------------------------
# 19. Figure 5: heatmap-like tile plot of volatile taxa
# ------------------------------------------------------------
# This ggplot-based heatmap avoids additional heatmap package dependencies.
# Each domain is plotted separately.

plot_volatile_heatmap <- function(target_domain, n_top = 50) {
  
  message2("Plotting heatmap for ", target_domain)
  
  top_taxa <- taxa_volatility2 %>%
    filter(domain == target_domain) %>%
    slice_max(order_by = volatility, n = n_top, with_ties = FALSE) %>%
    pull(taxon_id)
  
  if (length(top_taxa) == 0) {
    warning("No taxa found for domain: ", target_domain)
    return(NULL)
  }
  
  df <- taxa_complete %>%
    filter(domain == target_domain, taxon_id %in% top_taxa) %>%
    group_by(taxon_id, genus, site_group, depth_bin, date, year_month) %>%
    summarise(
      rpm = mean(rpm, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    group_by(taxon_id) %>%
    mutate(
      z = as.numeric(scale(log1p(rpm))),
      z = replace_na(z, 0),
      z = pmax(pmin(z, 3), -3),
      col_id = paste(site_group, depth_bin, year_month, sep = " | ")
    ) %>%
    ungroup()
  
  row_order <- taxa_volatility2 %>%
    filter(taxon_id %in% top_taxa) %>%
    arrange(volatility) %>%
    pull(taxon_id)
  
  df <- df %>%
    mutate(
      taxon_id = factor(taxon_id, levels = row_order),
      col_id = factor(col_id, levels = unique(col_id[order(date, site_group, depth_bin)]))
    )
  
  p <- ggplot(df, aes(x = col_id, y = taxon_id, fill = z)) +
    geom_tile(color = "grey90", linewidth = 0.1) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-3, 3),
      name = "z-score\nlog1p RPM"
    ) +
    theme_nature_like(base_size = 8) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 5),
      panel.grid = element_blank()
    ) +
    labs(
      title = paste0(target_domain, " volatile genus heatmap"),
      subtitle = "Rows are top volatile taxa; columns are habitat-depth-time combinations",
      x = "Habitat group | depth | time",
      y = NULL
    )
  
  ggsave(
    file.path(outdir, paste0("Fig5_heatmap_", safe_filename(target_domain), "_top_volatile.pdf")),
    p,
    width = 13,
    height = 10
  )
  
  return(p)
}

plot_volatile_heatmap("Bacteria", n_top = 50)
plot_volatile_heatmap("Archaea", n_top = 50)
plot_volatile_heatmap("Fungi", n_top = 50)


# ------------------------------------------------------------
# 20. Figure 6: depth-time profile for candidate taxa
# ------------------------------------------------------------
# This function is useful for manually checking important taxa.
# Example:
#   plot_taxon_depth_time("Fungi", "Trichophyton")
#   plot_taxon_depth_time("Fungi", "Fungi_gen_Incertae_sedis")

plot_taxon_depth_time <- function(target_domain, target_genus) {
  
  message2("Plotting depth-time profile: ", target_domain, " | ", target_genus)
  
  df <- taxa_complete %>%
    filter(domain == target_domain, genus == target_genus) %>%
    group_by(date, year_month, site_group, site_id, depth_bin) %>%
    summarise(
      rpm = sum(rpm, na.rm = TRUE),
      .groups = "drop"
    )
  
  if (nrow(df) == 0) {
    warning("No data found for: ", target_domain, " | ", target_genus)
    return(NULL)
  }
  
  p <- ggplot(df, aes(x = date, y = depth_bin, fill = log1p(rpm))) +
    geom_tile(color = "grey85", linewidth = 0.25) +
    facet_grid(site_id ~ site_group, scales = "free_x", space = "free_x") +
    scale_y_discrete(limits = rev(depth_levels)) +
    scale_fill_gradientn(
      colors = c("#F7F7F7", "#92C5DE", "#2166AC", "#B2182B"),
      name = "log1p RPM"
    ) +
    theme_nature_like(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6)
    ) +
    labs(
      title = paste(target_domain, target_genus, sep = " | "),
      subtitle = "Depth-time RPM profile across fixed stations",
      x = NULL,
      y = "Depth"
    )
  
  outfile <- paste0(
    "Fig6_depth_time_",
    safe_filename(target_domain),
    "_",
    safe_filename(target_genus),
    ".pdf"
  )
  
  ggsave(
    file.path(outdir, outfile),
    p,
    width = 11,
    height = 8
  )
  
  return(p)
}


# Automatically plot top 3 candidate taxa from each domain
auto_taxa_to_plot <- candidate_taxa %>%
  group_by(domain) %>%
  slice_max(order_by = max_shock, n = 3, with_ties = FALSE) %>%
  ungroup() %>%
  select(domain, genus)

if (nrow(auto_taxa_to_plot) > 0) {
  for (i in seq_len(nrow(auto_taxa_to_plot))) {
    plot_taxon_depth_time(
      as.character(auto_taxa_to_plot$domain[i]),
      as.character(auto_taxa_to_plot$genus[i])
    )
  }
}


# ------------------------------------------------------------
# 21. Figure 7: top taxa in the strongest event interval
# ------------------------------------------------------------

top_event_taxa <- taxon_event_contribution %>%
  group_by(domain, interval) %>%
  slice_max(order_by = shocked_units, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    interval = factor(interval, levels = interval_order),
    label = paste(domain, genus, sep = " | "),
    label = fct_reorder(label, median_abs_log2fc)
  )

p_event_taxa <- ggplot(
  top_event_taxa,
  aes(x = median_abs_log2fc, y = label)
) +
  geom_point(aes(size = shocked_ratio, color = domain), alpha = 0.85) +
  facet_wrap(~ interval, scales = "free_y") +
  scale_color_npg() +
  scale_size_continuous(range = c(1.5, 6)) +
  theme_nature_like(base_size = 8) +
  theme(
    axis.text.y = element_text(size = 5),
    panel.grid.major.y = element_blank()
  ) +
  labs(
    title = "Top taxa contributing to each disturbance interval",
    subtitle = "Ranked by shocked units and median absolute log2FC",
    x = "Median |log2FC|",
    y = NULL,
    color = "Domain",
    size = "Shocked ratio"
  )

ggsave(
  file.path(outdir, "Fig7_top_taxa_by_event_interval.pdf"),
  p_event_taxa,
  width = 14,
  height = 10
)


# ------------------------------------------------------------
# 22. Final report summary
# ------------------------------------------------------------

summary_lines <- c(
  paste0("Input file: ", input_file),
  paste0("Output directory: ", outdir),
  paste0("Number of samples: ", n_samples),
  paste0("Number of taxa before filtering: ", nrow(taxa_pre_filter)),
  paste0("Number of taxa retained: ", nrow(taxa_keep)),
  paste0("Pseudo for log2FC: ", signif(pseudo, 4)),
  paste0("min_prevalence: ", min_prevalence),
  paste0("min_mean_rpm: ", min_mean_rpm),
  paste0("shock_log2fc: ", shock_log2fc, " ; fold-change = ", 2^shock_log2fc),
  paste0("strong_shock_log2fc: ", strong_shock_log2fc, " ; fold-change = ", 2^strong_shock_log2fc)
)

writeLines(summary_lines, con = file.path(outdir, "analysis_summary.txt"))

message2("Analysis finished.")
message2("Results written to: ", outdir)