################################################################################
# Script: AIS_eDNA_GAM_seasonality.R
# Purpose: AIS eDNA Generalised Additive Models (seasonal signal & thresholds)
# Author: Melissa K. Morrison
# Created: 20 June 2025
# Updated: 07 July 2025  (extensive annotation for long‑term storage)
#
# Overview ---------------------------------------------------------------------
# 1. loads weekly‑aggregated eDNA + env data       (dfWeeks  from prep script)
# 2. fits a cyclic‑spline GAM:
#       logConc ~  s(sampWeek, bs = "cc")  +  s(replicate, bs = "re")
# 3. extracts useful seasonality info:
#       • peak week (max fitted value)
#       • low week  (min fitted value)
#       • SE bands for plotting
# 4. returns a tidy list (model, predictions, peaks) + optional diagnostic plot
###############################################################################

# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
load_pkgs <- function() {
  pkgs <- c("tidyverse",  # dplyr, ggplot2, purrr, tibble, etc.
            "mgcv",       # GAM engine
            "gratia")     # nice GAM helpers (draw, derivatives, etc.)
  invisible(lapply(pkgs, require, character.only = TRUE))
}

# ── 1. LOAD WEEK‑LEVEL DATA ───────────────────────────────────────────────────
load_data <- function(prep_script = "AIS_eDNA_data_prep.R") {
  source(prep_script)   # loads dfWeeks and others
  dfRawClean    # return object for use below
}

# ── 2. SINGLE‑GROUP MODEL FITTER  --------------------------------------------
fit_gam_one <- function(data = dfRawClean,
                        species_name,
                        region_name) {

  # Filter to target group
  df_sub <- dfRawClean |>
    filter(species != "Didemnum vexillum",
           species == species_name,
           region  == region_name) %>%
    select(region, species, date, replicate,
           logConc, month, temp, sal, month_reordered) %>%
    group_by(region, species) %>%
    # Week index relative to first sample in that reg.spp group
    mutate(
      sampWeek = floor(as.period(date - min(date)) / weeks()) + 1
    ) %>%
    ungroup() %>%
    mutate(replicate = factor(replicate))            # random‑effect term

  # If group is very small, bail early
  if (nrow(df_sub) < 30) {
    stop("Too few rows for ", species_name, " in ", region_name)
  }

  # ---- 2A  Fit cyclic spline on week  ---------------------------------------
  mod <- gam(logConc ~
               s(sampWeek, bs = "cc", k = 12) +   # cyclic cubic (k ≈ # knots)
               s(replicate, bs = "re"),       # random effect, replicate ID
             data   = df_sub,
             method = "REML")

  # ---- 2B  Predictions & SE --------------------------------------------------
  df_fit <- df_sub |>
    mutate(
      fit = predict(mod, type = "response", se.fit = TRUE)$fit,
      se  = predict(mod, type = "response", se.fit = TRUE)$se.fit
      ) |>
    select(region, species, date, sampWeek, fit, se) |> unique()

  # ---- 2C  Peak & low weeks --------------------------------------------------
  peaks <- df_fit |>
    group_by(region, species) %>%
    summarise(
      peak_date = date[which.max(fit)],
      peak_week = sampWeek[which.max(fit)],
      peak_val  = max(fit),
      low_week  = sampWeek[which.min(fit)],
      low_date  = date[which.min(fit)],
      low_val   = min(fit)
    )

  list(model = mod, preds = df_fit, peaks = peaks)
}

fit_gam_one_v2 <- function(df, species_name, region_name) {
  df_sub <- dfRawClean |>
    drop_na(sal) %>%
    filter(species != "Didemnum vexillum", !region %in% "GOM",
           species == species_name, region == region_name) %>%
    select(region, species, date, replicate,
           logConc, month, temp, sal, month_reordered) %>%
    group_by(region, species) %>%
    # Week index relative to first sample in that reg.spp group
    mutate(
      sampWeek = floor(as.period(date - min(date)) / weeks()) + 1
    ) %>%
    ungroup() %>%
    mutate(replicate = factor(replicate))        # random‑effect term


  if (nrow(df_sub) < 10) return(NULL)  # safety check

  mod <- gam(
    logConc ~ s(sampWeek, bs = "cc") +
      s(temp, bs = "tp") +
      s(sal, bs = "tp") +
      s(replicate, bs = "re"),
    data = df_sub,
    method = "REML", na.action = na.omit
  )

  df_fit <- df_sub |>
    mutate(fit = predict(mod, type = "response", se.fit = TRUE)$fit,
           se = predict(mod, type = "response", se.fit = TRUE)$se.fit) |>
    select(region, species, date, sampWeek, fit, se) |> unique()

  # Find peaks
  peaks <- df_fit |>
    group_by(region, species) |>
    summarise(
      peak_week = sampWeek[which.max(fit)],
      peak_date = date[which.max(fit)],
      peak_val = max(fit, na.rm = TRUE),
      low_week = sampWeek[which.min(fit)],
      low_val = min(fit, na.rm = TRUE),
      low_date  = date[which.min(fit)],
    )

  list(
    model = mod,
    pred = df_fit,
    peaks = peaks
  )
}


# ── 3. QUICK DIAGNOSTIC PLOT  -------------------------------------------------
plot_gam_curve <- function(pred_df, peaks, title_txt) {
  ggplot(pred_df, aes(sampWeek, fit)) +
    geom_line(colour = "steelblue", linewidth = 1) +
    geom_ribbon(aes(ymin = fit - se, ymax = fit + se),
                alpha = 0.25, fill = "steelblue") +
    geom_point(data = peaks,
               aes(x = peak_week, y = peak_val),
               colour = "firebrick", size = 3) +
    geom_point(data = peaks,
               aes(x = low_week,  y = low_val),
               colour = "black",   size = 3) +
    scale_x_continuous(breaks = seq(1, 52, 4), expand = c(0, 0)) +
    labs(x = "Sampling week",
         y = "Fitted log‑eDNA",
         title = title_txt,
         subtitle = paste0("Peak = ", month.abb[month(peaks$peak_date)],
                               ", Low = ", month.abb[month(peaks$low_date)])) +
    theme_classic()
}



# ── 4. DRIVER FOR ONE GROUP  --------------------------------------------------
run_gam_workflow <- function(species_name, region_name, show_plot = TRUE) {
  load_pkgs()
  wks <- load_data()

  res  <- fit_gam_one(wks, species_name, region_name)
  print(summary(res$model))          # param & smooth table
  if (show_plot) {
    plot_gam_curve(res$preds, res$peaks,
                   glue::glue("{species}  –  {region}")) |>
      print()
  }
  invisible(res)
}

# ── 5. LOOP ACROSS ALL REGION–SPECIES COMBOS  ---------------------------------
run_all_gams <- function(min_rows = 30, show_plot = FALSE) {
  load_pkgs()
  wks <- load_data()

  # ---- Identify valid species–region combinations ----------------------------
  combos <- wks |>
    drop_na(sal, temp) |>
    count(species, region, name = "n") |>
    filter(n >= min_rows, !species %in% "Didemnum vexillum", !region %in% "GOM") |>
  #  filter(!species %in% "Membranipora membranacea" | !region %in% "HAL") |>
    mutate(label = paste0(species, "-", region))

  # ---- Fit GAMs for each combo -----------------------------------------------
  results <- map2(combos$species, combos$region,
                  ~ fit_gam_one_v2(wks, .x, .y))
  names(results) <- combos$label

  # ---- Extract & bind peak results -------------------------------------------
  all_peaks <- map_dfr(results, ~ .x$peaks, .id = "group") |>
    separate(group, into = c("species", "region"), sep = "-")

  write_csv(all_peaks, "./data/GAM_peaks_by_group2.csv")
  message("✔ Saved peak values to: ./data/GAM_peaks_by_group2.csv")

  # ---- Optional: plot each smooth --------------------------------------------
  if (show_plot) {
    walk2(results, names(results), function(res, label) {
      if (!is.null(res)) {
        plot_gam_curve(res$preds, res$peaks, title_txt = label) |> print()
      }
    })
  }

  return(list(results = results, peaks = all_peaks))
}

# -----------------------------------------------------------------------------
# Uncomment to run
# -----------------------------------------------------------------------------
lapply(unique(dfRawClean$species), function(x)
  run_gam_workflow(x, "BOF"))# single demo

out <- run_all_gams()                                        # run every reg‑spp

# ── Extract smooth effects for temp and sal ───────────────────────────────────
extract_covariate_effects <- function(results_list, covariates = c("temp", "sal"), save_dir = NULL) {
#  require(gratia)
 # require(ggplot2)
  plots <- list()

 resultsComplete =  Filter(Negate(is.null), results_list)

  for (var in covariates) {
    plots[[var]] <- map2(resultsComplete, names(resultsComplete), function(res, label) {
      p <- draw(res$model, select = var, partial_match = TRUE) +
      #  ggtitle(glue::glue("Effect of {var}")) +
        labs(title = NULL, x = NULL, subtitle = str_split_fixed(label, "-", 2)[2]) +
        theme_minimal(base_size = 12)

      # Optionally save plots
      if (!is.null(save_dir)) {
        ggsave(
          filename = file.path(save_dir, paste0(label, "_", var, "_effect.png")),
          plot = p, width = 6, height = 4
        )
      }
      return(p)
    }) |> compact()  # remove NULLs
  }

  return(plots)
}
# Extract & view smooth plots
cov_plots <- extract_covariate_effects(out$results)#, save_dir = "./figures/gam_covariates")

cov_plotlistTemp <- lapply(unique(combos$species), function(spp)
  cowplot::plot_grid(plotlist = cov_plots$temp[grepl(spp, names(cov_plots$temp))],
                     ncol = 1) +
  plot_annotation(title = spp, subtitle = "Temperature effect"))

plot_grid(plotlist = cov_plotlistTemp, nrow = 1)

cov_plotlistSal <- lapply(unique(combos$species), function(spp)
  cowplot::plot_grid(plotlist = cov_plots$sal[grepl(spp, names(cov_plots$sal))],
                     ncol = 1) +
    plot_annotation(title = spp, subtitle = "Salinity"))

plot_grid(plotlist = cov_plotlistSal, nrow = 1)


# ── Summarise smooth significance (edf + p-values) ----------------------------
summarise_covariate_significance <- function(results_list, covariates = c("temp", "sal")) {
  map_dfr(names(results_list), function(label) {
    res <- results_list[[label]]
    if (is.null(res)) return(NULL)

    s_table <- summary(res$model)$s.table
    s_table_df <- as.data.frame(s_table) %>%
      rownames_to_column("term") %>%
      filter(str_detect(term, paste0("^s\\(", covariates, "\\)", collapse = "|"))) %>%
      mutate(label = label) %>%
      relocate(label, term, edf, `p-value`)
  }) %>%
    separate(label, into = c("species", "region"), sep = "-", remove = FALSE)
}

# Summarise smooth term significance
cov_summary <- summarise_covariate_significance(out$results)
print(cov_summary)

# Optionally save to CSV
write_csv(cov_summary, "./data/GAM_covariate_smooth_significance.csv")

# ── Check GAMs ----------------------------------------------------------------

model_list = lapply(resultsComplete, '[[', 1)

gam_checks = lapply(model_list, gam.check)
