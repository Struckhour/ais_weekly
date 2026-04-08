################################################################################
# Script:    AIS_eDNA_brt_thresholds.R
# Purpose:   Boosted Regression Trees & Temperature Thresholds
# Author:    Melissa K. Morrison
# Created:   20‑Jun‑2025
# Modified:  07‑Jul‑2025  (extensive annotation for long‑term storage)
#
# Overview ---------------------------------------------------------------------
# 1) loads cleaned data (source: AIS_eDNA_data_prep.R)
# 2) fits one GBM regression per region‑species combo
# 3) extracts temperature thresholds from each fitted model
# 4) saves tidy tables + gives a quick facetted plot by region
#
# WHY WE USE GBM ---------------------------------------------------------------
# • Handles non‑linear, interacting effects (temp × sal) without specifying them
# • Robust to unequal sampling effort (bagging + boosting)
# • Provides interpretable variable importance
###############################################################################

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
load_pkgs <- function() {
  pkgs <- c("tidyverse",    # dplyr, ggplot2, purrr, etc.
            "caret",        # unified model interface + CV tuning
            "gbm",          # original GBM back‑end
            "pracma")       # gradient() for inflection estimate
  invisible(lapply(pkgs, require, character.only = TRUE))
}

# ── 1. LOAD DATA (dfRawCleanCiona, etc.) ──────────────────────────────────────
load_data <- function(prep_script = "AIS_eDNA_data_prep.R") {
  source(prep_script)  # loads dfRawClean, etc. into env
  # Return whatever you need in a list
  list(df = dfRawClean)
}

# ── 2. PREP MODEL DATA ────────────────────────────────────────────────────────
prep_model_df <- function(df) {
  df %>%
    select(region, temp, sal, species, concentration, logConc) %>%
    mutate(
      region   = factor(region),
      reg.spp  = paste(region, species, sep = "-")
    ) %>%
    drop_na(concentration, temp, sal) %>%
    filter(region != "GOM")        # omit GOM due to no environmental data
}

# ── 3. FIT ONE GBM (with modest tuning) ───────────────────────────────────────
fit_gbm <- function(train_df, grid, ctrl) {
  train(logConc ~ temp + sal,
        data      = train_df,
        method    = "gbm",
        metric    = "RMSE",
        tuneGrid  = grid,
        trControl = ctrl,
        verbose   = FALSE,
        na.action = na.omit)
}

# ── 4. THRESHOLD EXTRACTOR  ---------------------------------------------------
extract_thresholds <- function(model, full_df, n.trees) {
  # Build prediction grid varying TEMPERATURE while salinity held at mean
  temp_seq <- seq(min(full_df$temp), max(full_df$temp), by = 0.1)
  newdata  <- tibble(
    temp = temp_seq,
    sal  = mean(full_df$sal, na.rm = TRUE)
  )
  newdata$pred <- predict(model, newdata = newdata, n.trees = n.trees)

  # 1. T_max  = peak predicted log‑conc
  t_max <- newdata$temp[which.max(newdata$pred)]

  # 2. T_inflect = highest first derivative
  d1    <- pracma::gradient(newdata$pred, temp_seq)
  t_inflect <- temp_seq[which.max(d1)]

  # 3. 10–90 % window
  max_val <- max(newdata$pred)
  inside  <- newdata$temp[newdata$pred >= 0.10 * max_val]
  t_low   <- min(inside)
  t_high  <- max(inside)

  tibble(t_max, t_inflect, t_low, t_high)
}

# ── 5. QUICK RESPONSE CURVE PLOT  --------------------------------------------
plot_response_curve <- function(pred_df, thresholds, group_label) {
  ggplot(pred_df, aes(temp, exp(pred))) +
    geom_line() +
    geom_vline(xintercept = thresholds$t_max,      linetype = "dashed", colour = "red") +
    geom_vline(xintercept = c(thresholds$t_low,
                              thresholds$t_high), linetype = "dotted") +
    labs(title = group_label,
         y = "Predicted eDNA (copies/L)",
         x = "Temperature (°C)") +
    theme_classic()
}

# ── 6. MAIN WORKFLOW ----------------------------------------------------------
run_brt_workflow <- function() {
  load_pkgs()
  dat  <- load_data()
  mod_df <- prep_model_df(dat$df)

  # caret control & modest grid (adjust if needed)
  ctrl <- trainControl(method = "cv", number = 5, verboseIter = FALSE)
  grid <- expand.grid(interaction.depth = 3,
                      n.trees           = 1500,
                      shrinkage         = 0.05,
                      n.minobsinnode    = 5)

  # Store results
  model_list      <- list()
  threshold_table <- tibble()
  curve_plots     <- list()

  # Loop over each region‑species combo
  for(rs in unique(mod_df$reg.spp)) {
    message("Fitting model for ", rs)
    df_sub <- filter(mod_df, reg.spp == rs)

    # If too few rows, skip
    if (nrow(df_sub) < 30 || dplyr::n_distinct(df_sub$logConc) < 5) {
      warning("Skipping ", rs, " (insufficient data)")
      next
    }

    model <- fit_gbm(df_sub, grid, ctrl)
    best_iter <- model$bestTune$n.trees

    # thresholds
    th <- extract_thresholds(model, df_sub, best_iter) %>%
      mutate(reg.spp = rs,
             region  = str_split(rs, "-", simplify = TRUE)[,1],
             species = str_split(rs, "-", simplify = TRUE)[,2])

    # prediction curve (for optional plotting)
    temp_seq <- seq(min(df_sub$temp), max(df_sub$temp), by = 0.1)
    pred_df  <- tibble(temp = temp_seq,
                       sal  = mean(df_sub$sal, na.rm = TRUE))
    pred_df$pred <- predict(model, newdata = pred_df, n.trees = best_iter)

    curve_plots[[rs]] <- plot_response_curve(pred_df, th, rs)
    model_list[[rs]]  <- model
    threshold_table   <- bind_rows(threshold_table, th)
  }

  # ── OUTPUTS ────────────────────────────────────────────────────────────────
  print(threshold_table |> select(-reg.spp))
  # e.g. save thresholds for later use
  write_csv(threshold_table, "./data/temperature_thresholds.csv")

  # Combine and print plots (optional)
  p <- patchwork::wrap_plots(curve_plots, ncol = 3) +
    patchwork::plot_layout(guides = "collect")

  print(p)

  invisible(list(models = model_list,
                 thresholds = threshold_table,
                 plots = curve_plots))
}

# --- uncomment to run ---------------------------------------------------------
run_brt_workflow()
