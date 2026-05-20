################################################################################
# Script: AIS_eDNA_descriptive_stats.R
# Purpose: AIS eDNA Descriptive Stats & Visualization
# Author: Melissa K. Morrison
# Created: 20 June 2025
# Updated: 07 July 2025  (extensive annotation for long‑term storage)

# Overview ---------------------------------------------------------------------
# 1)  Loads cleaned data from AIS_eDNA_data_prep.R
# 2)  Produces:
#     • Boxplots of temp / sal / turb / pH by region
#     • Basic summaries (mean / sd / five‑number) per region
#     • Region‑wide correlation matrices and species‑specific correlations
#     • Month‑by‑month eDNA scatter + trend plots
#     • Weekly env‑variable bar plots
#     • Life‑stage timeline & plate‑abundance columns
# 3)  Writes correlation CSVs to ./data/pub2/
#
# HOW TO USE ------------------------------------------------------------------
#   source("AIS_descriptive_stats.R");  run_all()
#   – Everything is encapsulated in small functions so you can run parts
#     interactively (e.g., plot_env_boxplots(dfRaw))
###############################################################################

# ──────────────────────────────────────────────────────────────────────────────
# 0. Load libraries ------------------------------------------------------------
# ──────────────────────────────────────────────────────────────────────────────

load_pkgs <- function() {
  pkgs <- c("tidyverse",   # readr, dplyr, tidyr, ggplot2, stringr
            "patchwork",   # combine ggplots
            "correlation") # nice correlation tables
  invisible(lapply(pkgs, require, character.only = TRUE))
}

# ───────────────────────────────────────────────────────────────────────────────
# 1. Load pre-processed data ----------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
load_data <- function(prep_script = "AIS_eDNA_data_prep.R") {
  source(prep_script) # brings dfRaw, dfWeeks, dfStages, etc. into env
}

# ───────────────────────────────────────────────────────────────────────────────
# 2. Environmental boxplots (excl. GOM) -----------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
plot_env_boxplots <- function(df) {
  vars   <- c(sal = "Salinity (psu)",
              temp = "Temperature (°C)",
              turb = "Turbidity (NTU)",
              pH   = "pH")
  # imap(): .x = value (pretty label), .y = name (column in df)
  p_list <- imap(vars, \(label, varname) {
    ggplot(filter(df, region != "GOM"),
           aes(x = region, y = .data[[varname]])) +
      geom_boxplot() +
      theme_classic() +
      labs(x = "Region", y = label)
  })
  wrap_plots(p_list) + plot_layout(axis_titles = "collect")
}

# ───────────────────────────────────────────────────────────────────────────────
# 3. Quick summary stats by region ----------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
summary_by_region <- function(df, cols) {
  stats <- map(cols, ~ list(
    mean = tapply(df[[.x]], df$region, mean, na.rm = TRUE),
    sd       = tapply(df[[.x]], df$region, sd, na.rm = TRUE)
  ))
  set_names(stats, cols)
}

# ───────────────────────────────────────────────────────────────────────────────
# 4. Correlation helpers --------------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
write_region_cor <- function(data, region, vars, file_stub) {
  correlation(filter(data, region == !!region)[, vars],
              include_factors = TRUE, method = "pearson") %>%
    as.data.frame() %>%
    select(Parameter1, Parameter2, r, p) %>%
    write_csv(file = glue::glue("./data/{region}_{file_stub}.csv"))
}

write_species_cor <- function(data, region, vars, file_stub) {
  split(filter(data, region == !!region), ~species) |>
    map_df(~ correlation(.x[, vars], include_factors = TRUE) |>
             as.data.frame() |>
             select(Parameter1, Parameter2, r, p),
           .id = "species") |>
    write_csv(file = glue::glue("./data/{region}_{file_stub}.csv"))
}

# ───────────────────────────────────────────────────────────────────────────────
# 5. eDNA month-by-month plot ---------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
plot_edna_month <- function(df) {
  ggplot(df) +
    geom_jitter(aes(x = month_reordered, y = meanConc, fill = region),
                alpha = 0.8, shape = 21, size = 1) +
    geom_smooth(aes(x = month_reordered, y = meanConc, fill = region),
                method = "loess", colour = "grey10") +
    facet_grid(region ~ species, scales = "free_y") +
    theme_classic() +
    scale_fill_manual(values = mapCols) +
    scale_x_continuous(breaks = month_labels_split$HAL$month_reordered,
                       labels = month_labels_split$HAL$month_plot_label) +
    labs(x = "Month", y = "eDNA Copies/L", fill = "Region")
}

# ───────────────────────────────────────────────────────────────────────────────
# 6. Log mean eDNA by region ----------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
plot_log_mean_conc <- function(dfWeeks, region, month_labels) {
  ggplot(filter(dfWeeks, region == !!region),
         aes(x = month_reordered,
             y = exp(log(meanConc + 1)))) +
    geom_jitter() +
    facet_wrap(~species, nrow = 1) +
    theme_classic() +
    scale_x_continuous(name   = "Month",
                       breaks = month_labels_split[[region]]$month_reordered,
                       labels = month_labels_split[[region]]$month_plot_label) +
    scale_y_continuous(name   = "Log mean\n eDNA concentration",
                       trans  = "log",
                       breaks = c(1,10,100,1e3,1e4,1e6))
}

# ───────────────────────────────────────────────────────────────────────────────
# 7. Weekly bar plot (temp, sal) ------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
plot_env_week <- function(dfWeeks, region, var, y_lab) {
  ggplot(filter(dfWeeks, region == !!region),
         aes(x = sampWeek, y = .data[[var]])) +
    geom_col(width = 1, fill = "grey75", alpha = 0.9) +
    facet_wrap(~species, nrow = 1) +
    theme_classic() +
    scale_x_continuous("Sampling week", breaks = seq(1, 52, 9)) +
    labs(y = y_lab)
}

# ───────────────────────────────────────────────────────────────────────────────
# 8. Life‑stage observations & plate abundance ----------------------------------
# ───────────────────────────────────────────────────────────────────────────────
plot_life_stage <- function(dfStages, region) {
  ggplot(filter(dfStages, region == !!region, !is.na(state)),
         aes(x = sampWeek, y = state, colour = state)) +
    geom_point(size = 3, shape = 15) +
    facet_wrap(~species, nrow = 1) +
    theme_classic() +
    scale_x_continuous("Sampling week", breaks = seq(1, 52, 9)) +
    labs(y = NULL) +
    theme(legend.position = "none")
}

plot_plate_abundance <- function(plateDF, region) {
  ggplot(filter(plateDF, region == !!region),
         aes(x = sampWeek, y = Avg)) +
    geom_col(width = 1) +
    facet_wrap(~species, nrow = 1, scales = "free") +
    theme_classic() +
    scale_x_continuous("Sampling week", breaks = seq(1, 52, 9), limits = c(1, 52)) +
    labs(y = "Mean plate abundance")
}

# ───────────────────────────────────────────────────────────────────────────────
# 9. Run all functions ----------------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
run_all <- function() {
  load_pkgs()
  load_data()   # creates dfRaw, dfWeeks, dfCor, dfStages, plateAbundance, mapCols, etc.

  ## 2. Boxplots
  plot_env_boxplots(dfRawClean)

  ## 3. Quick summaries
  env_stats <- summary_by_region(dfRawClean, c("temp", "sal", "pH", "tds", "chl"))
  # (inspect via env_stats$sal$five_num etc.)

  ## 4. Correlation CSVs
  vars_cor <- c("concentration", "temp", "sal", "turb", "pH", "chl")
  lapply(dfCorRegion, function(x)
    write_region_cor(x, first(x$region), vars_cor, "RegionCorr"))
  lapply(dfCorRegion, function(x)
    write_species_cor(x, first(x$region), vars_cor, "SpeciesCorr"))

  ## 5. Month‑by‑month eDNA
  plot_edna_month(dfWeeks) |> print()

  ## 6. Log mean conc by region
  regions <- c("MAG", "PEI", "HAL", "BOF", "GOM")
  walk(regions, ~ plot_log_mean_conc(dfWeeks, .x, month_labels_split) |> print())

  ## 7. Weekly env bars - temp and sal
  lapply(split(dfWeeks, ~region), function(x)
      plot_env_week(x, first(x$region), "meanTemp", "Mean water temperature (°C)")
)
  lapply(unique(dfWeeks$region), function(x)
      plot_env_week(dfWeeks, x, "meanSal",  "Mean salinity (psu)")
      )

  ## 8. Life stage + plate abundance
  lapply(unique(dfStages$region), function(x)
    plot_life_stage(dfStages, x))
  lapply(unique(plateAbundance$region), function(x)
    plot_plate_abundance(plateAbundance, x))
}

# ─── Uncomment to run everything ───────────────────────────────────────────────
run_all()

# ───────────────────────────────── END OF FILE ─────────────────────────────────
