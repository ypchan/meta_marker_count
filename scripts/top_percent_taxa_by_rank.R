#!/usr/bin/env Rscript

# ============================================================
# top99_taxa_by_rank.R
#
# Purpose:
#   For a specified taxonomic rank, identify the minimum number
#   of top taxa that cumulatively explain 99% relative abundance.
#
# Main outputs:
#   1. detailed taxa table with rpm, relative abundance and cumulative abundance
#   2. summary table with number of taxa required for 99%
#   3. taxa list table
#
# Input:
#   long-format table with at least:
#     sample_id, date, site_id, site_group, depth_bin, domain, lineage, rpm
#
# Author:
#   For domain-to-genus RPM dynamics
# ============================================================


suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
})


# ------------------------------------------------------------
# 1. User parameters
# ------------------------------------------------------------

input_file <- "domain_rank_rpm.tsv"

outdir <- "top99_taxa_results"

# 指定要分析的分类水平：
# 可选：
#   domain, phylum, class, order, family, genus, species
target_rank <- "genus"

# 累计相对丰度阈值
coverage_cutoff <- 0.99

# 是否在每个 domain 内部单独计算 99%
#
# TRUE:
#   对 Bacteria / Archaea / Fungi 分别计算。
#   推荐用于你的三域比较。
#
# FALSE:
#   把所有 domain 混在一起计算。
within_domain <- TRUE

# 分组变量。
# 也就是在什么尺度下分别计算 top99 taxa。
#
# 推荐：
#   每个时间 × 样点 × 深度 × domain 单独计算
#
# 如果想按地理组聚合，可以改成：
#   c("date", "site_group", "depth_bin", "domain")
#
# 如果想全局每个时间算：
#   c("date", "domain")
group_vars <- c("date", "site_group", "site_id", "depth_bin", "domain")

# 站点和深度顺序
site_group_levels <- c("MF", "MG", "TD")
site_levels <- c("MF1", "MF2", "MG1", "MG2", "TD1", "TD2")
depth_levels <- c("00-10", "10-20", "20-30", "30-40", "40-50", "50-60")
domain_levels <- c("Bacteria", "Archaea", "Fungi")


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
  x <- as.character(x)
  out <- suppressWarnings(as.Date(x))
  
  idx <- is.na(out)
  
  if (any(idx)) {
    x2 <- x[idx]
    
    # YYYY-MM
    id1 <- grepl("^\\d{4}-\\d{2}$", x2)
    if (any(id1)) {
      out[idx][id1] <- as.Date(paste0(x2[id1], "-01"))
    }
    
    # YYYYMM
    id2 <- grepl("^\\d{6}$", x2)
    if (any(id2)) {
      out[idx][id2] <- as.Date(
        paste0(substr(x2[id2], 1, 4), "-", substr(x2[id2], 5, 6), "-01")
      )
    }
    
    # YYYY/MM/DD
    id3 <- grepl("^\\d{4}/\\d{2}/\\d{2}$", x2)
    if (any(id3)) {
      out[idx][id3] <- as.Date(gsub("/", "-", x2[id3]))
    }
  }
  
  out
}


# ------------------------------------------------------------
# 3. Extract taxon name from lineage
# ------------------------------------------------------------
# 支持常见格式：
#
#   d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;...
#   k__Fungi;p__Ascomycota;c__Sordariomycetes;...
#   D_0__Bacteria;D_1__Proteobacteria;D_2__Gammaproteobacteria;...
#
# 如果 lineage 没有标准前缀，则尝试按顺序解析：
#   domain;phylum;class;order;family;genus;species

extract_rank_from_lineage <- function(lineage, rank) {
  
  rank <- tolower(rank)
  
  rank_prefix <- list(
    domain  = c("d__", "k__", "D_0__"),
    phylum  = c("p__", "D_1__"),
    class   = c("c__", "D_2__"),
    order   = c("o__", "D_3__"),
    family  = c("f__", "D_4__"),
    genus   = c("g__", "D_5__"),
    species = c("s__", "D_6__")
  )
  
  rank_index <- c(
    domain = 1,
    phylum = 2,
    class = 3,
    order = 4,
    family = 5,
    genus = 6,
    species = 7
  )
  
  if (!rank %in% names(rank_prefix)) {
    stop("Unsupported rank: ", rank, call. = FALSE)
  }
  
  purrr::map_chr(lineage, function(x) {
    
    if (is.na(x) || x == "") {
      return(paste0("Unclassified_", rank))
    }
    
    # 分割 lineage
    parts <- stringr::str_split(x, ";|\\|", simplify = FALSE)[[1]]
    parts <- stringr::str_trim(parts)
    parts <- parts[parts != ""]
    
    if (length(parts) == 0) {
      return(paste0("Unclassified_", rank))
    }
    
    # 先按 rank prefix 搜索
    prefixes <- rank_prefix[[rank]]
    
    hit <- NA_character_
    
    for (pf in prefixes) {
      idx <- which(stringr::str_starts(parts, fixed(pf)))
      if (length(idx) > 0) {
        hit <- parts[idx[1]]
        hit <- stringr::str_replace(hit, paste0("^", stringr::fixed(pf)), "")
        break
      }
    }
    
    # 如果没有 prefix，按 lineage 顺序取
    if (is.na(hit)) {
      idx <- rank_index[[rank]]
      if (length(parts) >= idx) {
        hit <- parts[idx]
        # 去掉可能存在的前缀
        hit <- stringr::str_replace(hit, "^[A-Za-z]__", "")
        hit <- stringr::str_replace(hit, "^D_[0-9]__", "")
      }
    }
    
    # 清理异常值
    if (is.na(hit) || hit == "" || hit %in% c("uncultured", "metagenome")) {
      hit <- paste0("Unclassified_", rank)
    }
    
    hit
  })
}


# ------------------------------------------------------------
# 4. Read input
# ------------------------------------------------------------

message2("Reading input: ", input_file)

df <- readr::read_tsv(input_file, show_col_types = FALSE)

required_cols <- c(
  "sample_id", "date", "site_id", "site_group",
  "depth_bin", "domain", "lineage", "rpm"
)

missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", "),
    call. = FALSE
  )
}

df <- df %>%
  mutate(
    sample_id = as.character(sample_id),
    date = parse_sample_date(date),
    year_month = format(date, "%Y-%m"),
    
    site_id = as.character(site_id),
    site_group = as.character(site_group),
    depth_bin = as.character(depth_bin),
    domain = as.character(domain),
    lineage = as.character(lineage),
    
    rpm = as.numeric(rpm),
    rpm = replace_na(rpm, 0),
    rpm = if_else(rpm < 0, 0, rpm),
    
    site_group = factor(site_group, levels = site_group_levels),
    site_id = factor(site_id, levels = site_levels),
    depth_bin = factor(depth_bin, levels = depth_levels),
    domain = factor(domain, levels = domain_levels)
  )

if (any(is.na(df$date))) {
  stop("Some date values cannot be parsed.", call. = FALSE)
}

message2("Rows: ", nrow(df))
message2("Samples: ", n_distinct(df$sample_id))
message2("Target rank: ", target_rank)
message2("Coverage cutoff: ", coverage_cutoff)


# ------------------------------------------------------------
# 5. Prepare target taxon
# ------------------------------------------------------------
# 如果输入表中已经有对应 rank 列，比如 genus/family/order，
# 优先使用已有列。
# 否则从 lineage 中解析。

if (target_rank %in% colnames(df)) {
  
  message2("Using existing column for rank: ", target_rank)
  
  df <- df %>%
    mutate(
      taxon = as.character(.data[[target_rank]]),
      taxon = if_else(
        is.na(taxon) | taxon == "" | taxon == "NA",
        paste0("Unclassified_", target_rank),
        taxon
      )
    )
  
} else {
  
  message2("Extracting rank from lineage: ", target_rank)
  
  df <- df %>%
    mutate(
      taxon = extract_rank_from_lineage(lineage, target_rank)
    )
}

df <- df %>%
  mutate(
    rank = target_rank,
    taxon_id = paste(rank, taxon, sep = "|")
  )


# ------------------------------------------------------------
# 6. Grouping design
# ------------------------------------------------------------

if (!within_domain) {
  group_vars <- setdiff(group_vars, "domain")
}

missing_group_vars <- setdiff(group_vars, colnames(df))

if (length(missing_group_vars) > 0) {
  stop(
    "Grouping variables not found in input: ",
    paste(missing_group_vars, collapse = ", "),
    call. = FALSE
  )
}

message2("Grouping variables: ", paste(group_vars, collapse = ", "))


# ------------------------------------------------------------
# 7. Aggregate RPM at target rank
# ------------------------------------------------------------
# 如果多个 lineage 映射到同一个 taxon，例如同一个 genus 下多个 species，
# 在指定 rank 下需要把 RPM 合并。
#
# lineage 这里保留该 taxon 下最常见或最高 RPM 的代表 lineage，
# 方便回溯。

taxa_rank <- df %>%
  group_by(across(all_of(group_vars)), rank, taxon) %>%
  summarise(
    rpm = sum(rpm, na.rm = TRUE),
    
    # 保留该 taxon 中 RPM 最大的 lineage 作为代表 lineage
    lineage = lineage[which.max(rpm)][1],
    
    .groups = "drop"
  )

readr::write_tsv(
  taxa_rank,
  file.path(outdir, paste0("rank_", target_rank, "_rpm_aggregated.tsv"))
)


# ------------------------------------------------------------
# 8. Calculate relative abundance and cumulative abundance
# ------------------------------------------------------------

top99_detail_all <- taxa_rank %>%
  group_by(across(all_of(group_vars))) %>%
  mutate(
    total_rpm = sum(rpm, na.rm = TRUE),
    
    relative_abundance = if_else(
      total_rpm > 0,
      rpm / total_rpm,
      0
    )
  ) %>%
  arrange(
    across(all_of(group_vars)),
    desc(relative_abundance),
    desc(rpm),
    taxon
  ) %>%
  mutate(
    rank_order = row_number(),
    cumulative_relative_abundance = cumsum(relative_abundance),
    
    # include_top99:
    #   保留累计丰度达到 99% 之前的所有 taxa，
    #   同时包括第一个使累计丰度超过 99% 的 taxon。
    previous_cumulative = lag(cumulative_relative_abundance, default = 0),
    include_top99 = previous_cumulative < coverage_cutoff
  ) %>%
  ungroup()

readr::write_tsv(
  top99_detail_all,
  file.path(outdir, paste0("top99_", target_rank, "_all_ranked_taxa.tsv"))
)

top99_detail <- top99_detail_all %>%
  filter(include_top99) %>%
  select(
    all_of(group_vars),
    rank,
    taxon,
    lineage,
    rpm,
    total_rpm,
    relative_abundance,
    cumulative_relative_abundance,
    rank_order
  )

readr::write_tsv(
  top99_detail,
  file.path(outdir, paste0("top99_", target_rank, "_detail.tsv"))
)


# ------------------------------------------------------------
# 9. Summary: number of taxa needed to cover 99%
# ------------------------------------------------------------

top99_summary <- top99_detail %>%
  group_by(across(all_of(group_vars)), rank) %>%
  summarise(
    top99_taxa_number = n(),
    top99_total_rpm = sum(rpm, na.rm = TRUE),
    total_rpm = max(total_rpm, na.rm = TRUE),
    top99_cumulative_relative_abundance = max(
      cumulative_relative_abundance,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  arrange(across(all_of(group_vars)))

readr::write_tsv(
  top99_summary,
  file.path(outdir, paste0("top99_", target_rank, "_summary.tsv"))
)


# ------------------------------------------------------------
# 10. Taxa list table
# ------------------------------------------------------------

top99_taxa_list <- top99_detail %>%
  arrange(across(all_of(group_vars)), rank_order) %>%
  group_by(across(all_of(group_vars)), rank) %>%
  summarise(
    top99_taxa_number = n(),
    top99_taxa_list = paste(taxon, collapse = "; "),
    top99_lineage_list = paste(lineage, collapse = " || "),
    .groups = "drop"
  )

readr::write_tsv(
  top99_taxa_list,
  file.path(outdir, paste0("top99_", target_rank, "_taxa_list.tsv"))
)


# ------------------------------------------------------------
# 11. Optional: summary across time, site group and depth
# ------------------------------------------------------------
# 这个表可以回答：
#   哪些样点/深度需要更多 taxa 才能覆盖 99%？
#   群落是否变得更分散？
#
# top99_taxa_number 越大：
#   表示优势类群不集中，多样性/均匀度更高。
#
# top99_taxa_number 越小：
#   表示少数类群占据绝对优势，群落更集中。

summary_by_group_depth <- top99_summary %>%
  group_by(site_group, depth_bin, domain, rank) %>%
  summarise(
    mean_top99_taxa_number = mean(top99_taxa_number, na.rm = TRUE),
    median_top99_taxa_number = median(top99_taxa_number, na.rm = TRUE),
    min_top99_taxa_number = min(top99_taxa_number, na.rm = TRUE),
    max_top99_taxa_number = max(top99_taxa_number, na.rm = TRUE),
    .groups = "drop"
  )

readr::write_tsv(
  summary_by_group_depth,
  file.path(outdir, paste0("top99_", target_rank, "_summary_by_group_depth.tsv"))
)


# ------------------------------------------------------------
# 12. Visualization
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
      plot.title = element_text(face = "bold", hjust = 0)
    )
}

# 图 1：每个时间点达到 99% 所需 taxa 数量
p_top99_number <- ggplot(
  top99_summary,
  aes(x = date, y = top99_taxa_number, color = domain, group = interaction(domain, site_id))
) +
  geom_line(alpha = 0.55, linewidth = 0.5) +
  geom_point(size = 1.6, alpha = 0.85) +
  facet_grid(depth_bin ~ site_group, scales = "free_y") +
  ggsci::scale_color_npg() +
  theme_nature_like(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7)
  ) +
  labs(
    title = paste0("Number of ", target_rank, " taxa required to cover 99% abundance"),
    subtitle = "Higher values indicate a more even or dispersed community",
    x = NULL,
    y = paste0("Top taxa number for 99% ", target_rank, " abundance"),
    color = "Domain"
  )

ggsave(
  file.path(outdir, paste0("Fig_top99_", target_rank, "_taxa_number.pdf")),
  p_top99_number,
  width = 12,
  height = 8
)


# 图 2：累计相对丰度曲线
# 每个 group-depth-domain-date 的曲线太多，这里按 site_group × depth_bin × domain 聚合显示。
cum_curve <- top99_detail_all %>%
  group_by(site_group, depth_bin, domain, rank_order) %>%
  summarise(
    mean_cumulative = mean(cumulative_relative_abundance, na.rm = TRUE),
    median_cumulative = median(cumulative_relative_abundance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(rank_order <= 100)

p_cum_curve <- ggplot(
  cum_curve,
  aes(x = rank_order, y = median_cumulative, color = domain)
) +
  geom_hline(yintercept = coverage_cutoff, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.8) +
  facet_grid(depth_bin ~ site_group) +
  ggsci::scale_color_npg() +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_nature_like(base_size = 10) +
  labs(
    title = paste0("Cumulative relative abundance curve at ", target_rank, " level"),
    subtitle = paste0("Dashed line indicates ", coverage_cutoff * 100, "% coverage"),
    x = paste0(target_rank, " rank order"),
    y = "Cumulative relative abundance",
    color = "Domain"
  )

ggsave(
  file.path(outdir, paste0("Fig_top99_", target_rank, "_cumulative_curve.pdf")),
  p_cum_curve,
  width = 12,
  height = 8
)


# 图 3：某个指定 rank 的 top99 taxa 组成
# 默认只画每个 domain 中总体 RPM 最高的前 20 个 taxa，其余合并为 Others。

top_taxa_for_bar <- top99_detail %>%
  group_by(domain, taxon) %>%
  summarise(total_rpm = sum(rpm, na.rm = TRUE), .groups = "drop") %>%
  group_by(domain) %>%
  slice_max(order_by = total_rpm, n = 20, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(taxon_keep = taxon)

bar_df <- top99_detail %>%
  left_join(
    top_taxa_for_bar %>% select(domain, taxon, taxon_keep),
    by = c("domain", "taxon")
  ) %>%
  mutate(
    taxon_plot = if_else(is.na(taxon_keep), "Others", taxon)
  ) %>%
  group_by(date, site_group, depth_bin, domain, taxon_plot) %>%
  summarise(
    rpm = sum(rpm, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(date, site_group, depth_bin, domain) %>%
  mutate(
    relative_abundance = rpm / sum(rpm, na.rm = TRUE)
  ) %>%
  ungroup()

p_bar <- ggplot(
  bar_df,
  aes(x = date, y = relative_abundance, fill = taxon_plot)
) +
  geom_col(width = 25) +
  facet_grid(depth_bin ~ site_group + domain, scales = "free_x", space = "free_x") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_nature_like(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
    legend.position = "right",
    legend.text = element_text(size = 6),
    panel.grid = element_blank()
  ) +
  labs(
    title = paste0("Composition of taxa covering 99% abundance at ", target_rank, " level"),
    subtitle = "Top taxa are shown individually; remaining taxa are grouped as Others",
    x = NULL,
    y = "Relative abundance within top99 set",
    fill = target_rank
  )

ggsave(
  file.path(outdir, paste0("Fig_top99_", target_rank, "_composition_bar.pdf")),
  p_bar,
  width = 16,
  height = 10
)


# ------------------------------------------------------------
# 13. Done
# ------------------------------------------------------------

message2("Finished.")
message2("Output directory: ", outdir)

message2("Main output files:")
message2("  ", file.path(outdir, paste0("top99_", target_rank, "_detail.tsv")))
message2("  ", file.path(outdir, paste0("top99_", target_rank, "_summary.tsv")))
message2("  ", file.path(outdir, paste0("top99_", target_rank, "_taxa_list.tsv")))