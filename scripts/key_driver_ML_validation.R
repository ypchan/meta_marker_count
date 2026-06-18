# Fig_key_driver_RF_XGBoost_validation_ASCII.R
#
# Purpose:
#   Screen key predictors with Random Forest and XGBoost, then validate
#   top predictors using linear mixed models or linear models.
#
# This script is ASCII-only to avoid garbled characters in RStudio,
# Windows terminals, Linux servers, and PDF output.
#
# Input 1: metadata_env.tsv
# Required columns:
#   sample_id
# Recommended covariates:
#   year, month, site_id, site_type, depth
# Candidate driver columns:
#   numeric or categorical environmental variables, for example:
#   pH, SOC, TN, TP, AP, salinity, moisture, year, depth, site_type
#
# Input 2: response_abundance.tsv
# Required columns:
#   sample_id
#   response_type
#   response_id
#   response_label
#   abundance
#
# Example response_type:
#   Bacteria_genus
#   Archaea_genus
#   Fungi_genus
#   ARG
#   Function

options(stringsAsFactors = FALSE)

packages <- c(
  "tidyverse", "data.table", "randomForest", "xgboost",
  "lme4", "lmerTest", "broom.mixed", "patchwork",
  "ggsci", "scales"
)

for (p in packages) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p)
  }
}

library(tidyverse)
library(data.table)
library(randomForest)
library(xgboost)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(patchwork)
library(ggsci)
library(scales)

# ============================================================
# Parameters
# ============================================================

metadata_file <- "metadata_env.tsv"
response_file <- "response_abundance.tsv"
output_dir <- "key_driver_RF_XGBoost_validation"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

response_types_keep <- NULL
# response_types_keep <- c("Bacteria_genus", "Archaea_genus", "Fungi_genus", "ARG")

top_n_response_per_type <- 20
min_detected_samples <- 8
min_samples_for_model <- 20

candidate_predictors <- NULL
# candidate_predictors <- c("year", "site_type", "depth", "pH", "SOC", "TN", "AP")

fixed_covariates <- c("year", "site_type", "depth")
random_effects <- c("site_id")
factor_variables <- c("year", "month", "site_id", "site_type", "depth")

log_transform_response <- TRUE
pseudo_count <- 1e-8

rf_ntree <- 1000
rf_seed <- 123

xgb_seed <- 123
xgb_nrounds <- 1000
xgb_nfold <- 5
xgb_early_stopping_rounds <- 30
xgb_params <- list(
  objective = "reg:squarederror",
  eval_metric = "rmse",
  eta = 0.03,
  max_depth = 3,
  min_child_weight = 3,
  subsample = 0.8,
  colsample_bytree = 0.8,
  lambda = 1,
  alpha = 0
)

top_n_predictors_to_validate <- 8
p_adjust_method <- "BH"

plot_top_n_predictors <- 20
plot_top_n_validation <- 40

# ============================================================
# Helper functions
# ============================================================

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

safe_log <- function(x) {
  if (log_transform_response) {
    log10(x + pseudo_count)
  } else {
    x
  }
}

safe_zscore <- function(x) {
  x[!is.finite(x)] <- NA_real_
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    rep(0, length(x))
  } else {
    as.numeric(scale(x))
  }
}

safe_name <- function(x) {
  make.names(x, unique = TRUE)
}

make_design_matrix <- function(df, predictors) {
  xdf <- df %>%
    select(all_of(predictors)) %>%
    mutate(across(where(is.character), as.factor))

  for (v in names(xdf)) {
    if (is.logical(xdf[[v]])) {
      xdf[[v]] <- as.factor(xdf[[v]])
    }
  }

  mm <- model.matrix(~ . - 1, data = xdf)
  mode(mm) <- "numeric"

  feature_map <- tibble(model_feature = colnames(mm)) %>%
    rowwise() %>%
    mutate(
      predictor = {
        hits <- predictors[
          purrr::map_lgl(
            predictors,
            function(p) model_feature == p || stringr::str_starts(model_feature, fixed(p))
          )
        ]
        if (length(hits) == 0) model_feature else hits[which.max(nchar(hits))]
      }
    ) %>%
    ungroup()

  list(matrix = mm, feature_map = feature_map)
}

get_rf_importance <- function(df, y_col, predictors) {
  mdf <- df %>% select(all_of(c(y_col, predictors))) %>% drop_na()
  if (nrow(mdf) < min_samples_for_model) return(NULL)
  if (length(unique(mdf[[y_col]])) < 2) return(NULL)

  x <- mdf %>%
    select(all_of(predictors)) %>%
    mutate(across(where(is.character), as.factor))

  for (v in names(x)) {
    if (is.logical(x[[v]])) x[[v]] <- as.factor(x[[v]])
  }

  y <- mdf[[y_col]]
  set.seed(rf_seed)

  fit <- randomForest(
    x = x,
    y = y,
    ntree = rf_ntree,
    importance = TRUE
  )

  imp <- randomForest::importance(fit, type = 1) %>%
    as.data.frame() %>%
    rownames_to_column("predictor")

  if (!"%IncMSE" %in% names(imp)) {
    value_col <- setdiff(names(imp), "predictor")[1]
    imp <- imp %>% rename(`%IncMSE` = all_of(value_col))
  }

tibble(
  predictor = imp$predictor,
  rf_importance = imp$`%IncMSE`,
  rf_oob_r2 = max(0, tail(fit$rsq, 1) * 100, na.rm = TRUE)
)
}

get_xgb_importance <- function(df, y_col, predictors) {
  mdf <- df %>% select(all_of(c(y_col, predictors))) %>% drop_na()
  if (nrow(mdf) < min_samples_for_model) return(NULL)
  if (length(unique(mdf[[y_col]])) < 2) return(NULL)
  
  design <- make_design_matrix(mdf, predictors)
  xmat <- design$matrix
  fmap <- design$feature_map
  y <- mdf[[y_col]]
  
  if (ncol(xmat) == 0) return(NULL)
  if (any(!is.finite(xmat))) return(NULL)
  
  dtrain <- xgboost::xgb.DMatrix(data = xmat, label = y)
  
  set.seed(xgb_seed)
  cv <- xgboost::xgb.cv(
    params = xgb_params,
    data = dtrain,
    nrounds = xgb_nrounds,
    nfold = min(xgb_nfold, nrow(mdf)),
    early_stopping_rounds = xgb_early_stopping_rounds,
    verbose = 0
  )
  
  best_n <- cv$best_iteration
  if (is.null(best_n) || is.na(best_n) || best_n < 1) best_n <- 100
  
  set.seed(xgb_seed)
  fit <- xgboost::xgb.train(
    params = xgb_params,
    data = dtrain,
    nrounds = best_n,
    verbose = 0
  )
  
  pred <- predict(fit, dtrain)
  ss_res <- sum((y - pred)^2, na.rm = TRUE)
  ss_tot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  train_r2 <- ifelse(ss_tot > 0, max(0, 1 - ss_res / ss_tot) * 100, NA_real_)
  
  imp <- xgboost::xgb.importance(
    feature_names = colnames(xmat),
    model = fit
  )
  
  if (nrow(imp) == 0) return(NULL)
  
  imp %>%
    as_tibble() %>%
    rename(model_feature = Feature) %>%
    left_join(fmap, by = "model_feature") %>%
    group_by(predictor) %>%
    summarise(
      xgb_gain = sum(Gain, na.rm = TRUE),
      xgb_cover = sum(Cover, na.rm = TRUE),
      xgb_frequency = sum(Frequency, na.rm = TRUE),
      xgb_best_nrounds = best_n,
      xgb_train_r2 = train_r2,
      .groups = "drop"
    )
}

build_formula <- function(response_col, predictor, fixed_covariates, random_effects, data_cols) {
  fixed_covariates <- fixed_covariates[fixed_covariates %in% data_cols]
  random_effects <- random_effects[random_effects %in% data_cols]
  fixed_covariates <- setdiff(fixed_covariates, predictor)
  random_effects <- setdiff(random_effects, predictor)
  
  fixed_terms <- unique(c(predictor, fixed_covariates))
  fixed_part <- paste(fixed_terms, collapse = " + ")
  
  if (length(random_effects) > 0) {
    random_part <- paste0("(1 | ", random_effects, ")", collapse = " + ")
    as.formula(paste(response_col, "~", fixed_part, "+", random_part))
  } else {
    as.formula(paste(response_col, "~", fixed_part))
  }
}

extract_p_from_anova <- function(aov_df, predictor) {
  p_col <- names(aov_df)[stringr::str_detect(names(aov_df), "Pr")][1]
  stat_col <- names(aov_df)[stringr::str_detect(names(aov_df), "F.value|F value|Chisq")][1]
  
  if (is.na(p_col) || !predictor %in% rownames(aov_df)) {
    return(tibble(test_statistic = NA_real_, p_value = NA_real_))
  }
  
  stat <- if (!is.na(stat_col)) aov_df[predictor, stat_col] else NA_real_
  
  tibble(
    test_statistic = as.numeric(stat),
    p_value = as.numeric(aov_df[predictor, p_col])
  )
}

validate_predictor <- function(df, predictor) {
  vars <- unique(c("response_value", predictor, fixed_covariates_safe, random_effects_safe))
  vars <- vars[vars %in% names(df)]
  mdf <- df %>% select(all_of(vars)) %>% drop_na()
  
  if (nrow(mdf) < min_samples_for_model) return(tibble())
  if (length(unique(mdf$response_value)) < 2) return(tibble())
  if (length(unique(mdf[[predictor]])) < 2) return(tibble())
  
  usable_random <- random_effects_safe[random_effects_safe %in% names(mdf)]
  usable_random <- usable_random[
    purrr::map_lgl(usable_random, function(v) length(unique(mdf[[v]])) >= 2)
  ]
  
  fml <- build_formula(
    response_col = "response_value",
    predictor = predictor,
    fixed_covariates = fixed_covariates_safe,
    random_effects = usable_random,
    data_cols = names(mdf)
  )
  
  if (length(usable_random) > 0) {
    fit <- lmerTest::lmer(fml, data = mdf, REML = FALSE)
    aov_df <- as.data.frame(lmerTest::anova(fit, type = 3))
    tidy_fit <- broom.mixed::tidy(fit, effects = "fixed")
    model_type <- "mixed_model"
  } else {
    fit <- lm(fml, data = mdf)
    aov_df <- as.data.frame(drop1(fit, test = "F"))
    tidy_fit <- broom.mixed::tidy(fit)
    model_type <- "linear_model"
  }
  
  p_res <- extract_p_from_anova(aov_df, predictor)
  
  if (is.numeric(mdf[[predictor]])) {
    erow <- tidy_fit %>% filter(term == predictor) %>% slice(1)
    est <- ifelse(nrow(erow) == 1, erow$estimate, NA_real_)
    se <- ifelse(nrow(erow) == 1, erow$std.error, NA_real_)
    pclass <- "numeric"
  } else {
    est <- NA_real_
    se <- NA_real_
    pclass <- "categorical"
  }
  
  tibble(
    predictor = predictor,
    predictor_class = pclass,
    estimate = est,
    std_error = se,
    test_statistic = p_res$test_statistic,
    p_value = p_res$p_value,
    n = nrow(mdf),
    model_type = model_type,
    model_formula = deparse(fml)
  )
}

# ============================================================
# Read input files
# ============================================================

if (!file.exists(metadata_file)) stop("Metadata file not found: ", metadata_file)
if (!file.exists(response_file)) stop("Response file not found: ", response_file)

metadata_raw <- data.table::fread(metadata_file, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)
response_raw <- data.table::fread(response_file, sep = "\t", header = TRUE, data.table = FALSE, check.names = FALSE)

names(metadata_raw) <- trimws(names(metadata_raw))
names(response_raw) <- trimws(names(response_raw))

if (!"sample_id" %in% names(metadata_raw)) stop("metadata_file must contain sample_id.")
required_response_cols <- c("sample_id", "response_type", "response_id", "response_label", "abundance")
missing_response_cols <- setdiff(required_response_cols, names(response_raw))
if (length(missing_response_cols) > 0) {
  stop("response_file missing columns: ", paste(missing_response_cols, collapse = ", "))
}

metadata_raw <- metadata_raw %>% mutate(sample_id = as.character(sample_id))
response_raw <- response_raw %>%
  mutate(
    sample_id = as.character(sample_id),
    response_type = as.character(response_type),
    response_id = as.character(response_id),
    response_label = as.character(response_label),
    abundance = as.numeric(abundance)
  ) %>%
  filter(!is.na(abundance), abundance >= 0)

if (!is.null(response_types_keep)) {
  response_raw <- response_raw %>% filter(response_type %in% response_types_keep)
}

common_samples <- intersect(metadata_raw$sample_id, response_raw$sample_id)
if (length(common_samples) < min_samples_for_model) {
  stop("Too few common samples between metadata and response tables.")
}

metadata_raw <- metadata_raw %>% filter(sample_id %in% common_samples)
response_raw <- response_raw %>% filter(sample_id %in% common_samples)

# ============================================================
# Prepare metadata
# ============================================================

old_names <- names(metadata_raw)
safe_names <- make.names(old_names, unique = TRUE)
names(metadata_raw) <- safe_names
name_map <- tibble(original = old_names, safe = safe_names)

sample_safe <- name_map$safe[name_map$original == "sample_id"][1]
names(metadata_raw)[names(metadata_raw) == sample_safe] <- "sample_id"
name_map$safe[name_map$original == "sample_id"] <- "sample_id"

metadata <- metadata_raw

map_vars <- function(v) {
  out <- name_map$safe[match(v, name_map$original)]
  out <- out[!is.na(out)]
  unique(out)
}

fixed_covariates_safe <- map_vars(fixed_covariates)
random_effects_safe <- map_vars(random_effects)
factor_variables_safe <- map_vars(factor_variables)

for (v in factor_variables_safe) {
  if (v %in% names(metadata)) metadata[[v]] <- factor(metadata[[v]])
}

if (is.null(candidate_predictors)) {
  candidate_predictors_safe <- setdiff(names(metadata), unique(c("sample_id", random_effects_safe)))
} else {
  candidate_predictors_safe <- map_vars(candidate_predictors)
}

candidate_predictors_safe <- candidate_predictors_safe[candidate_predictors_safe %in% names(metadata)]

if (length(candidate_predictors_safe) == 0) stop("No candidate predictors found.")

for (v in candidate_predictors_safe) {
  if (is.character(metadata[[v]])) metadata[[v]] <- factor(metadata[[v]])
}

predictor_name_table <- name_map %>%
  filter(safe %in% candidate_predictors_safe) %>%
  transmute(predictor = safe, predictor_original = original)

# ============================================================
# Select response features
# ============================================================

response_info <- response_raw %>%
  group_by(response_type, response_id, response_label) %>%
  summarise(
    total_abundance = sum(abundance, na.rm = TRUE),
    mean_abundance = mean(abundance, na.rm = TRUE),
    detected_samples = n_distinct(sample_id[abundance > 0]),
    .groups = "drop"
  ) %>%
  filter(detected_samples >= min_detected_samples) %>%
  group_by(response_type) %>%
  arrange(desc(total_abundance), .by_group = TRUE) %>%
  slice_head(n = top_n_response_per_type) %>%
  ungroup() %>%
  mutate(response_safe_id = safe_name(paste(response_type, response_id, sep = "__")))

write_tsv(response_info, file.path(output_dir, "selected_response_features.tsv"))

response_selected <- response_raw %>%
  inner_join(
    response_info %>% select(response_type, response_id, response_label, response_safe_id),
    by = c("response_type", "response_id", "response_label")
  )

# ============================================================
# Run ML and validation
# ============================================================

rf_list <- list()
xgb_list <- list()
val_list <- list()

feature_ids <- unique(response_selected$response_safe_id)
message("[INFO] Selected response features: ", length(feature_ids))

for (fid in feature_ids) {
  fmeta <- response_info %>% filter(response_safe_id == fid) %>% slice(1)
  
  ydf <- response_selected %>%
    filter(response_safe_id == fid) %>%
    group_by(sample_id) %>%
    summarise(abundance = sum(abundance, na.rm = TRUE), .groups = "drop")
  
  model_df <- metadata %>%
    select(sample_id, all_of(unique(c(fixed_covariates_safe, random_effects_safe, candidate_predictors_safe)))) %>%
    left_join(ydf, by = "sample_id") %>%
    mutate(
      abundance = ifelse(is.na(abundance), 0, abundance),
      response_value = safe_log(abundance)
    )
  
  rf_imp <- get_rf_importance(model_df, "response_value", candidate_predictors_safe)
  xgb_imp <- get_xgb_importance(model_df, "response_value", candidate_predictors_safe)
  
  if (!is.null(rf_imp) && nrow(rf_imp) > 0) {
    rf_list[[fid]] <- rf_imp %>%
      mutate(
        response_type = fmeta$response_type,
        response_id = fmeta$response_id,
        response_label = fmeta$response_label,
        response_safe_id = fid
      )
  }
  
  if (!is.null(xgb_imp) && nrow(xgb_imp) > 0) {
    xgb_list[[fid]] <- xgb_imp %>%
      mutate(
        response_type = fmeta$response_type,
        response_id = fmeta$response_id,
        response_label = fmeta$response_label,
        response_safe_id = fid
      )
  }
  
  top_rf <- character(0)
  top_xgb <- character(0)
  
  if (!is.null(rf_imp) && nrow(rf_imp) > 0) {
    top_rf <- rf_imp %>%
      arrange(desc(rf_importance)) %>%
      slice_head(n = top_n_predictors_to_validate) %>%
      pull(predictor)
  }
  
  if (!is.null(xgb_imp) && nrow(xgb_imp) > 0) {
    top_xgb <- xgb_imp %>%
      arrange(desc(xgb_gain)) %>%
      slice_head(n = top_n_predictors_to_validate) %>%
      pull(predictor)
  }
  
  top_predictors <- union(top_rf, top_xgb)
  
  val_one <- purrr::map_dfr(top_predictors, function(pred) validate_predictor(model_df, pred))
  
  if (nrow(val_one) > 0) {
    val_list[[fid]] <- val_one %>%
      mutate(
        response_type = fmeta$response_type,
        response_id = fmeta$response_id,
        response_label = fmeta$response_label,
        response_safe_id = fid
      )
  }
}

rf_importance <- bind_rows(rf_list)
xgb_importance <- bind_rows(xgb_list)
validation_results <- bind_rows(val_list)

if (nrow(rf_importance) == 0) stop("No Random Forest model was fitted.")
if (nrow(xgb_importance) == 0) stop("No XGBoost model was fitted.")

rf_importance <- rf_importance %>% left_join(predictor_name_table, by = "predictor")
xgb_importance <- xgb_importance %>% left_join(predictor_name_table, by = "predictor")

if (nrow(validation_results) > 0) {
  validation_results <- validation_results %>%
    left_join(predictor_name_table, by = "predictor") %>%
    group_by(response_type) %>%
    mutate(q_value = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(validated = !is.na(p_value) & p_value < 0.05)
}

write_tsv(rf_importance, file.path(output_dir, "ML_random_forest_predictor_importance.tsv"))
write_tsv(xgb_importance, file.path(output_dir, "ML_xgboost_predictor_importance.tsv"))
write_tsv(validation_results, file.path(output_dir, "LMM_validation_top_predictors.tsv"))

# ============================================================
# Summarize key drivers
# ============================================================

rf_summary <- rf_importance %>%
  group_by(response_type, predictor, predictor_original) %>%
  summarise(
    mean_rf_importance = mean(rf_importance, na.rm = TRUE),
    n_rf_features = n_distinct(response_safe_id),
    mean_rf_oob_r2 = mean(rf_oob_r2, na.rm = TRUE),
    .groups = "drop"
  )

xgb_summary <- xgb_importance %>%
  group_by(response_type, predictor, predictor_original) %>%
  summarise(
    mean_xgb_gain = mean(xgb_gain, na.rm = TRUE),
    n_xgb_features = n_distinct(response_safe_id),
    mean_xgb_train_r2 = mean(xgb_train_r2, na.rm = TRUE),
    .groups = "drop"
  )

driver_summary <- full_join(
  rf_summary,
  xgb_summary,
  by = c("response_type", "predictor", "predictor_original")
) %>%
  mutate(
    mean_rf_importance = replace_na(mean_rf_importance, 0),
    mean_xgb_gain = replace_na(mean_xgb_gain, 0),
    n_rf_features = replace_na(n_rf_features, 0),
    n_xgb_features = replace_na(n_xgb_features, 0)
  )

if (nrow(validation_results) > 0) {
  val_summary <- validation_results %>%
    group_by(response_type, predictor, predictor_original) %>%
    summarise(
      n_validated = sum(p_value < 0.05, na.rm = TRUE),
      n_tested = n(),
      validation_rate = n_validated / n_tested,
      median_p_value = median(p_value, na.rm = TRUE),
      median_abs_effect = median(abs(estimate), na.rm = TRUE),
      .groups = "drop"
    )
  
  driver_summary <- driver_summary %>%
    left_join(val_summary, by = c("response_type", "predictor", "predictor_original"))
}

driver_summary <- driver_summary %>%
  mutate(
    n_validated = replace_na(n_validated, 0),
    n_tested = replace_na(n_tested, 0),
    validation_rate = replace_na(validation_rate, 0),
    median_abs_effect = replace_na(median_abs_effect, 0)
  ) %>%
  group_by(response_type) %>%
  mutate(
    z_rf = safe_zscore(mean_rf_importance),
    z_xgb = safe_zscore(mean_xgb_gain),
    z_validation = safe_zscore(validation_rate),
    z_effect = safe_zscore(median_abs_effect),
    driver_score = z_rf + z_xgb + z_validation + z_effect
  ) %>%
  ungroup() %>%
  arrange(response_type, desc(driver_score))

write_tsv(driver_summary, file.path(output_dir, "key_driver_summary.tsv"))

# ============================================================
# Plot
# ============================================================

plot_rf <- driver_summary %>%
  group_by(response_type) %>%
  slice_max(order_by = mean_rf_importance, n = plot_top_n_predictors, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    predictor_label = ifelse(is.na(predictor_original), predictor, predictor_original),
    predictor_label = factor(predictor_label, levels = rev(unique(predictor_label)))
  )

p_rf <- ggplot(plot_rf, aes(x = response_type, y = predictor_label, fill = mean_rf_importance)) +
  geom_tile(color = "grey85", linewidth = 0.25) +
  scale_fill_gradient(low = "grey95", high = "#B2182B", name = "Mean RF importance") +
  labs(title = "Random forest screening", x = "Response type", y = "Predictor") +
  theme_nature(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.25))

plot_xgb <- driver_summary %>%
  group_by(response_type) %>%
  slice_max(order_by = mean_xgb_gain, n = plot_top_n_predictors, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    predictor_label = ifelse(is.na(predictor_original), predictor, predictor_original),
    predictor_label = factor(predictor_label, levels = rev(unique(predictor_label)))
  )

p_xgb <- ggplot(plot_xgb, aes(x = response_type, y = predictor_label, fill = mean_xgb_gain)) +
  geom_tile(color = "grey85", linewidth = 0.25) +
  scale_fill_gradient(low = "grey95", high = "#2166AC", name = "Mean XGBoost Gain") +
  labs(title = "XGBoost screening", x = "Response type", y = "Predictor") +
  theme_nature(base_size = 8) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.25))

plot_driver <- driver_summary %>%
  group_by(response_type) %>%
  slice_max(order_by = driver_score, n = 10, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    predictor_label = ifelse(is.na(predictor_original), predictor, predictor_original),
    driver_label = paste(response_type, predictor_label, sep = " | "),
    driver_label = factor(driver_label, levels = rev(driver_label))
  )

p_driver <- ggplot(plot_driver, aes(x = driver_score, y = driver_label, fill = response_type)) +
  geom_col(width = 0.70, color = "black", linewidth = 0.25) +
  scale_fill_npg(name = "Response type") +
  labs(title = "Integrated key driver score", x = "Driver score", y = NULL) +
  theme_nature(base_size = 8) +
  theme(legend.position = "bottom")

if (nrow(validation_results) > 0) {
  plot_val <- driver_summary %>%
    filter(n_tested > 0) %>%
    group_by(response_type) %>%
    slice_max(order_by = driver_score, n = plot_top_n_validation, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      predictor_label = ifelse(is.na(predictor_original), predictor, predictor_original),
      predictor_label = factor(predictor_label, levels = rev(unique(predictor_label))),
      neg_log10_p = -log10(pmax(median_p_value, 1e-300))
    )
  
  p_val <- ggplot(plot_val, aes(x = response_type, y = predictor_label)) +
    geom_point(aes(size = validation_rate, fill = neg_log10_p),
               shape = 21, color = "black", stroke = 0.35) +
    scale_size_continuous(range = c(1.0, 5.0), limits = c(0, 1),
                          labels = percent_format(accuracy = 1),
                          name = "Validation rate") +
    scale_fill_gradient(low = "grey95", high = "black", name = "Median P (-log10)") +
    labs(title = "Statistical validation", x = "Response type", y = "Predictor") +
    theme_nature(base_size = 8) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.25))
} else {
  p_val <- ggplot() +
    annotate("text", x = 0, y = 0, label = "No validation results") +
    theme_void()
}

final_fig <- (p_rf | p_xgb) / (p_driver | p_val) +
  plot_layout(heights = c(1, 1.15), guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    theme = theme(plot.tag = element_text(face = "bold", size = 12),
                  legend.position = "bottom")
  )

final_fig

ggsave(
  filename = file.path(output_dir, "key_driver_RF_XGBoost_validation.pdf"),
  plot = final_fig,
  width = 240,
  height = 200,
  units = "mm",
  device = "pdf"
)

ggsave(
  filename = file.path(output_dir, "key_driver_RF_XGBoost_validation.png"),
  plot = final_fig,
  width = 240,
  height = 200,
  units = "mm",
  dpi = 600
)

message("[DONE] Results written to: ", output_dir)

