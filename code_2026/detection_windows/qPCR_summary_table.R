source('./AIS_eDNA_data_prep.R')


summary_table <- dfRawClean %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    # overall stats
    mean_log_conc = mean(logConc, na.rm = TRUE),
    sd_log_conc   = sd(logConc, na.rm = TRUE),

    # date range
    start_date = min(date, na.rm = TRUE),
    end_date   = max(date, na.rm = TRUE),

    # mean SD across months
    mean_monthly_sd = {
      df_month <- cur_data() %>%
        dplyr::group_by(month) %>%
        dplyr::summarise(sd_month = sd(logConc, na.rm = TRUE), .groups = "drop")

      mean(df_month$sd_month, na.rm = TRUE)
    },

    .groups = "drop"
  ) %>%
  dplyr::arrange(species, region)































####################
#SD TABLE
####################

summary_table <- dfRawClean %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    n_obs = dplyr::n(),
    n_samples = dplyr::n_distinct(materialSampleID),
    n_stations = dplyr::n_distinct(station),

    mean_log_conc = mean(logConc, na.rm = TRUE),
    sd_log_conc = sd(logConc, na.rm = TRUE),

    .groups = "drop"
  )


monthly_sd <- dfRawClean %>%
  dplyr::group_by(species, region, month) %>%
  dplyr::summarise(
    monthly_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    mean_monthly_sd = mean(monthly_sd, na.rm = TRUE),
    .groups = "drop"
  )


weekly_sd <- dfRawClean %>%
  dplyr::mutate(
    week = lubridate::floor_date(date, "week")
  ) %>%
  dplyr::group_by(species, region, week) %>%
  dplyr::summarise(
    weekly_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    mean_weekly_sd = mean(weekly_sd, na.rm = TRUE),
    .groups = "drop"
  )


weekly_station_sd <- dfRawClean %>%
  dplyr::mutate(
    week = lubridate::floor_date(date, "week")
  ) %>%
  dplyr::group_by(species, region, station, week) %>%
  dplyr::summarise(
    station_weekly_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    mean_weekly_sd_within_stations = mean(station_weekly_sd, na.rm = TRUE),
    .groups = "drop"
  )


sample_sd <- dfRawClean %>%
  dplyr::group_by(species, region, materialSampleID) %>%
  dplyr::summarise(
    sample_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    mean_sd_within_samples = mean(sample_sd, na.rm = TRUE),
    .groups = "drop"
  )


summary_table_deeper <- summary_table %>%
  dplyr::left_join(monthly_sd, by = c("species", "region")) %>%
  dplyr::left_join(weekly_sd, by = c("species", "region")) %>%
  dplyr::left_join(sample_sd, by = c("species", "region")) %>%
  dplyr::arrange(species, region)


library(gt)
summary_table_gt <- summary_table_deeper %>%
  dplyr::arrange(species, region) %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(
    species_display = dplyr::if_else(
      dplyr::row_number() == 1,
      as.character(species),
      ""
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(species_display, dplyr::everything(), -species)


summary_table_gt %>%
  gt() %>%
  fmt_number(
    columns = -c(species_display, region, n_obs, n_samples, n_stations),
    decimals = 2
  ) %>%
  cols_label(
    species_display = "Species",
    region = "Region",
    n_obs = "N obs",
    n_samples = "N samples",
    n_stations = "N stations",
    mean_log_conc = "Mean log conc",
    sd_log_conc = "SD log conc",
    mean_monthly_sd = "Mean SD (monthly)",
    mean_weekly_sd = "Mean SD (weekly)",
    mean_sd_within_samples = "Mean SD (within sample)"
  ) %>%
  tab_header(
    title = "Summary of eDNA variability by species and region"
  ) %>%
  tab_options(
    table.font.size = px(10),
    heading.title.font.size = px(12),
    column_labels.font.size = px(9),
    data_row.padding = px(2),
    column_labels.padding = px(3),
    heading.padding = px(4),
    table.width = pct(100)
  ) %>%
  cols_width(
    species_display ~ px(120),
    region ~ px(45),
    n_obs:n_stations ~ px(50),
    everything() ~ px(75)
  )


####################
#Pairwise difference TABLE
####################

mean_pairwise_diff <- function(x) {
  x <- x[!is.na(x)]

  if (length(x) < 2) return(NA_real_)

  mean(abs(stats::dist(x)))
}

summary_table <- dfRawClean %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    n_obs = dplyr::n(),
    n_samples = dplyr::n_distinct(materialSampleID),
    n_stations = dplyr::n_distinct(station),

    mean_log_conc = mean(logConc, na.rm = TRUE),
    overall_mean_pairwise_diff = mean_pairwise_diff(logConc),

    .groups = "drop"
  )


monthly_pairwise <- dfRawClean %>%
  dplyr::group_by(species, region, month) %>%
  dplyr::summarise(
    monthly_pairwise_diff = mean_pairwise_diff(logConc),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    mean_monthly_pairwise_diff = mean(monthly_pairwise_diff, na.rm = TRUE),
    .groups = "drop"
  )


weekly_pairwise <- dfRawClean %>%
  dplyr::mutate(
    week = lubridate::floor_date(date, "week")
  ) %>%
  dplyr::group_by(species, region, week) %>%
  dplyr::summarise(
    weekly_pairwise_diff = mean_pairwise_diff(logConc),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    mean_weekly_pairwise_diff = mean(weekly_pairwise_diff, na.rm = TRUE),
    .groups = "drop"
  )


weekly_station_pairwise <- dfRawClean %>%
  dplyr::mutate(
    week = lubridate::floor_date(date, "week")
  ) %>%
  dplyr::group_by(species, region, station, week) %>%
  dplyr::summarise(
    station_weekly_pairwise_diff = mean_pairwise_diff(logConc),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    mean_weekly_pairwise_diff_within_stations =
      mean(station_weekly_pairwise_diff, na.rm = TRUE),
    .groups = "drop"
  )


sample_pairwise <- dfRawClean %>%
  dplyr::group_by(species, region, materialSampleID) %>%
  dplyr::summarise(
    sample_pairwise_diff = mean_pairwise_diff(logConc),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    mean_pairwise_diff_within_samples =
      mean(sample_pairwise_diff, na.rm = TRUE),
    .groups = "drop"
  )


summary_table_deeper <- summary_table %>%
  dplyr::left_join(monthly_pairwise, by = c("species", "region")) %>%
  dplyr::left_join(weekly_pairwise, by = c("species", "region")) %>%
  dplyr::left_join(weekly_station_pairwise, by = c("species", "region")) %>%
  dplyr::left_join(sample_pairwise, by = c("species", "region")) %>%
  dplyr::arrange(species, region)


summary_table_gt <- summary_table_deeper %>%
  dplyr::arrange(species, region) %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(
    species_display = dplyr::if_else(
      dplyr::row_number() == 1,
      as.character(species),
      ""
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(species_display, dplyr::everything(), -species)


summary_table_gt %>%
  gt::gt() %>%
  gt::fmt_number(
    columns = -c(species_display, region, n_obs, n_samples, n_stations),
    decimals = 2
  ) %>%
  gt::cols_label(
    species_display = "Species",
    region = "Region",
    n_obs = "N obs",
    n_samples = "N samples",
    n_stations = "N stations",
    mean_log_conc = "Mean log conc",
    overall_mean_pairwise_diff = "Overall mean pairwise diff",
    mean_monthly_pairwise_diff = "Mean pairwise diff (monthly)",
    mean_weekly_pairwise_diff = "Mean pairwise diff (weekly)",
    mean_weekly_pairwise_diff_within_stations = "Mean pairwise diff (station × week)",
    mean_pairwise_diff_within_samples = "Mean pairwise diff (within sample)"
  ) %>%
  gt::tab_header(
    title = "Summary of eDNA variability by species and region"
  ) %>%
  gt::tab_options(
    table.font.size = gt::px(10),
    heading.title.font.size = gt::px(12),
    column_labels.font.size = gt::px(9),
    data_row.padding = gt::px(2),
    column_labels.padding = gt::px(3),
    heading.padding = gt::px(4),
    table.width = gt::pct(100)
  ) %>%
  gt::cols_width(
    species_display ~ gt::px(120),
    region ~ gt::px(45),
    n_obs:n_stations ~ gt::px(50),
    everything() ~ gt::px(75)
  )









###########
#BAR PLOTS
###########

plot_df <- summary_table_deeper %>%
  dplyr::select(
    species, region,
    sd_log_conc,
    mean_monthly_sd,
    mean_weekly_sd
  ) %>%
  tidyr::pivot_longer(
    cols = -c(species, region),
    names_to = "metric",
    values_to = "value"
  )

plot_df <- plot_df %>%
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      sd_log_conc = "Annual SD",
      mean_monthly_sd = "Monthly SD",
      mean_weekly_sd = "Weekly SD"
    )
  )

ggplot(plot_df, aes(x = metric, y = value, fill = metric)) +
  geom_col(position = "dodge") +
  facet_grid(region ~ species) +
  theme_classic() +
  labs(
    title = "Comparison of variability metrics by species and region",
    x = "",
    y = "SD (log10 concentration)",
    fill = "Metric"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.background = element_rect(color = "black", fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(0, "lines")
  )



############PAIRWISE DIFFERENCE
plot_df <- summary_table_deeper %>%
  dplyr::select(
    species, region,
    overall_mean_pairwise_diff,
    mean_monthly_pairwise_diff,
    mean_weekly_pairwise_diff,
    mean_weekly_pairwise_diff_within_stations,
    mean_pairwise_diff_within_samples
  ) %>%
  tidyr::pivot_longer(
    cols = -c(species, region),
    names_to = "metric",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    metric = dplyr::recode(
      metric,
      overall_mean_pairwise_diff = "Annual",
      mean_monthly_pairwise_diff = "Monthly",
      mean_weekly_pairwise_diff = "Weekly",
      mean_weekly_pairwise_diff_within_stations = "Station-week",
      mean_pairwise_diff_within_samples = "Within-sample"
    ),
    metric = factor(
      metric,
      levels = c("Annual", "Monthly", "Weekly", "Station-week", "Within-sample")
    )
  )

ggplot(plot_df, aes(x = metric, y = value, fill = metric)) +
  geom_col() +
  facet_grid(region ~ species) +
  theme_classic() +
  labs(
    title = "Comparison of pairwise variability metrics by species and region",
    x = "",
    y = "Mean absolute pairwise difference in log10 concentration",
    fill = "Metric"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "bottom",
    strip.background = element_rect(color = "black", fill = "white"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(0, "lines")
  )
































#################################
#REGION SUMMARIES PAIRWISE DIFF
#################################
summary_species <- dfRawClean %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    n_obs = dplyr::n(),
    n_samples = dplyr::n_distinct(materialSampleID),
    n_stations = dplyr::n_distinct(station),
    n_regions = dplyr::n_distinct(region),

    mean_log_conc = mean(logConc, na.rm = TRUE),
    min_log_conc  = min(logConc, na.rm = TRUE),
    max_log_conc  = max(logConc, na.rm = TRUE),

    overall_mean_pairwise_diff = mean_pairwise_diff(logConc),

    .groups = "drop"
  )

monthly_pairwise_species <- dfRawClean %>%
  dplyr::group_by(species, month) %>%
  dplyr::summarise(
    monthly_pairwise_diff = mean_pairwise_diff(logConc),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    mean_monthly_pairwise_diff = mean(monthly_pairwise_diff, na.rm = TRUE),
    .groups = "drop"
  )

weekly_pairwise_species <- dfRawClean %>%
  dplyr::mutate(
    week = lubridate::floor_date(date, "week")
  ) %>%
  dplyr::group_by(species, week) %>%
  dplyr::summarise(
    weekly_pairwise_diff = mean_pairwise_diff(logConc),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    mean_weekly_pairwise_diff = mean(weekly_pairwise_diff, na.rm = TRUE),
    .groups = "drop"
  )


weekly_station_pairwise_species <- dfRawClean %>%
  dplyr::mutate(
    week = lubridate::floor_date(date, "week")
  ) %>%
  dplyr::group_by(species, station, week) %>%
  dplyr::summarise(
    station_weekly_pairwise_diff = mean_pairwise_diff(logConc),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    mean_weekly_pairwise_diff_within_stations =
      mean(station_weekly_pairwise_diff, na.rm = TRUE),
    .groups = "drop"
  )

sample_pairwise_species <- dfRawClean %>%
  dplyr::group_by(species, materialSampleID) %>%
  dplyr::summarise(
    sample_pairwise_diff = mean_pairwise_diff(logConc),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    mean_pairwise_diff_within_samples =
      mean(sample_pairwise_diff, na.rm = TRUE),
    .groups = "drop"
  )

summary_species_deeper <- summary_species %>%
  dplyr::left_join(monthly_pairwise_species, by = "species") %>%
  dplyr::left_join(weekly_pairwise_species, by = "species") %>%
  dplyr::left_join(weekly_station_pairwise_species, by = "species") %>%
  dplyr::left_join(sample_pairwise_species, by = "species") %>%
  dplyr::arrange(species)

summary_species_deeper %>%
  gt::gt() %>%
  gt::fmt_number(
    columns = -c(n_obs, n_samples, n_stations, n_regions),
    decimals = 2
  ) %>%
  gt::cols_label(
    species = "Species",
    n_obs = "N obs",
    n_samples = "N samples",
    n_stations = "N stations",
    n_regions = "N regions",
    mean_log_conc = "Mean log conc",
    min_log_conc = "Min",
    max_log_conc = "Max",
    overall_mean_pairwise_diff = "Annual pairwise",
    mean_monthly_pairwise_diff = "Monthly pairwise",
    mean_weekly_pairwise_diff = "Weekly pairwise",
    mean_weekly_pairwise_diff_within_stations = "Station-week pairwise",
    mean_pairwise_diff_within_samples = "Within-sample pairwise"
  ) %>%
  gt::tab_header(
    title = "Summary of eDNA variability by species (all regions combined)"
  )



##################
#SPECIES SUMMARY PAIRWISE DIFF
##################

summary_region <- dfRawClean %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    n_obs = dplyr::n(),
    n_samples = dplyr::n_distinct(materialSampleID),
    n_stations = dplyr::n_distinct(station),
    n_species = dplyr::n_distinct(species),

    mean_log_conc = mean(logConc, na.rm = TRUE),
    min_log_conc  = min(logConc, na.rm = TRUE),
    max_log_conc  = max(logConc, na.rm = TRUE),

    overall_mean_pairwise_diff = mean_pairwise_diff(logConc),

    .groups = "drop"
  )

monthly_pairwise_region <- dfRawClean %>%
  dplyr::group_by(region, month) %>%
  dplyr::summarise(
    monthly_pairwise_diff = mean_pairwise_diff(logConc),
    .groups = "drop"
  ) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    mean_monthly_pairwise_diff = mean(monthly_pairwise_diff, na.rm = TRUE),
    .groups = "drop"
  )

weekly_pairwise_region <- dfRawClean %>%
  dplyr::mutate(
    week = lubridate::floor_date(date, "week")
  ) %>%
  dplyr::group_by(region, week) %>%
  dplyr::summarise(
    weekly_pairwise_diff = mean_pairwise_diff(logConc),
    .groups = "drop"
  ) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    mean_weekly_pairwise_diff = mean(weekly_pairwise_diff, na.rm = TRUE),
    .groups = "drop"
  )

weekly_station_pairwise_region <- dfRawClean %>%
  dplyr::mutate(
    week = lubridate::floor_date(date, "week")
  ) %>%
  dplyr::group_by(region, station, week) %>%
  dplyr::summarise(
    station_weekly_pairwise_diff = mean_pairwise_diff(logConc),
    .groups = "drop"
  ) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    mean_weekly_pairwise_diff_within_stations =
      mean(station_weekly_pairwise_diff, na.rm = TRUE),
    .groups = "drop"
  )

sample_pairwise_region <- dfRawClean %>%
  dplyr::group_by(region, materialSampleID) %>%
  dplyr::summarise(
    sample_pairwise_diff = mean_pairwise_diff(logConc),
    .groups = "drop"
  ) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    mean_pairwise_diff_within_samples =
      mean(sample_pairwise_diff, na.rm = TRUE),
    .groups = "drop"
  )

summary_region_deeper <- summary_region %>%
  dplyr::left_join(monthly_pairwise_region, by = "region") %>%
  dplyr::left_join(weekly_pairwise_region, by = "region") %>%
  dplyr::left_join(weekly_station_pairwise_region, by = "region") %>%
  dplyr::left_join(sample_pairwise_region, by = "region") %>%
  dplyr::arrange(region)

summary_region_deeper %>%
  gt::gt() %>%
  gt::fmt_integer(
    columns = c(n_obs, n_samples, n_stations, n_species)
  ) %>%
  gt::fmt_number(
    columns = -c(n_obs, n_samples, n_stations, n_species),
    decimals = 2
  ) %>%
  gt::cols_label(
    region = "Region",
    n_obs = "N obs",
    n_samples = "N samples",
    n_stations = "N stations",
    n_species = "N species",
    mean_log_conc = "Mean log conc",
    min_log_conc = "Min",
    max_log_conc = "Max",
    overall_mean_pairwise_diff = "Annual pairwise",
    mean_monthly_pairwise_diff = "Monthly pairwise",
    mean_weekly_pairwise_diff = "Weekly pairwise",
    mean_weekly_pairwise_diff_within_stations = "Station-week pairwise",
    mean_pairwise_diff_within_samples = "Within-sample pairwise"
  ) %>%
  gt::tab_header(
    title = "Summary of eDNA variability by region (all species combined)"
  )








































































############################
#SPECIES SD SUMMARY
############################
# --- restrict to species actually tested in each region ---
region_species_tested <- dfRawClean %>%
  dplyr::distinct(region, species)

all_samples <- dfRawClean %>%
  dplyr::distinct(materialSampleID, region, date, month)

sample_species_grid <- all_samples %>%
  dplyr::left_join(region_species_tested, by = "region")

# --- collapse to sample × species ---
df_full <- sample_species_grid %>%
  dplyr::left_join(
    dfRawClean %>%
      dplyr::group_by(materialSampleID, species) %>%
      dplyr::summarise(
        logConc = mean(logConc, na.rm = TRUE),
        sample_detected = as.integer(any(detected == 1)),
        .groups = "drop"
      ),
    by = c("materialSampleID", "species")
  ) %>%
  dplyr::mutate(
    sample_detected = tidyr::replace_na(sample_detected, 0),
    logConc = tidyr::replace_na(logConc, 0)
  )

# --- positive-only data ---
df_pos <- df_full %>%
  dplyr::filter(sample_detected == 1)

# --- sample-level detection by species ---
sample_detection_species <- df_full %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    n_obs = dplyr::n(),
    n_samples = dplyr::n_distinct(materialSampleID),
    n_regions = dplyr::n_distinct(region),
    n_sample_detected = sum(sample_detected == 1),
    sample_detection_rate = mean(sample_detected == 1),
    .groups = "drop"
  )

# --- replicate-level detection by species ---
rep_detection_species <- dfRawClean %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    n_reps = dplyr::n(),
    n_pos_reps = sum(detected == 1, na.rm = TRUE),
    rep_detection_rate = mean(detected == 1, na.rm = TRUE),
    .groups = "drop"
  )

# --- monthly SD by species, positive only ---
monthly_sd_species <- df_pos %>%
  dplyr::group_by(species, month) %>%
  dplyr::summarise(
    monthly_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    mean_monthly_sd = mean(monthly_sd, na.rm = TRUE),
    .groups = "drop"
  )

# --- weekly SD by species, positive only ---
weekly_sd_species <- df_pos %>%
  dplyr::mutate(
    week = lubridate::floor_date(date, "week")
  ) %>%
  dplyr::group_by(species, week) %>%
  dplyr::summarise(
    weekly_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    mean_weekly_sd = mean(weekly_sd, na.rm = TRUE),
    .groups = "drop"
  )

# --- within-sample SD by species, positive only ---
sample_sd_species <- dfRawClean %>%
  dplyr::filter(detected == 1) %>%
  dplyr::group_by(species, materialSampleID) %>%
  dplyr::summarise(
    sample_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    mean_sd_within_samples = mean(sample_sd, na.rm = TRUE),
    .groups = "drop"
  )

# --- concentration stats by species, positive only ---
summary_species_sd <- df_pos %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    mean_log_conc = mean(logConc, na.rm = TRUE),
    max_log_conc  = max(logConc, na.rm = TRUE),
    sd_log_conc   = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  )

# --- combine ---
summary_species_final <- sample_detection_species %>%
  dplyr::left_join(rep_detection_species, by = "species") %>%
  dplyr::left_join(summary_species_sd, by = "species") %>%
  dplyr::left_join(monthly_sd_species, by = "species") %>%
  dplyr::left_join(weekly_sd_species, by = "species") %>%
  dplyr::left_join(sample_sd_species, by = "species") %>%
  dplyr::arrange(species)

summary_species_final %>%
  dplyr::select(-n_samples) %>%
  dplyr::rename(
    n_samples = n_obs
  ) %>%
  dplyr::select(
    species,
    n_regions,
    n_samples,
    n_sample_detected,
    sample_detection_rate,
    n_reps,
    n_pos_reps,
    rep_detection_rate,
    mean_log_conc,
    max_log_conc,
    sd_log_conc,
    mean_monthly_sd,
    mean_weekly_sd,
    mean_sd_within_samples
  ) %>%
  gt::gt() %>%
  gt::fmt_integer(
    columns = c(n_regions, n_samples, n_sample_detected, n_reps, n_pos_reps)
  ) %>%
  gt::fmt_percent(
    columns = c(sample_detection_rate, rep_detection_rate),
    decimals = 1
  ) %>%
  gt::fmt_number(
    columns = c(
      mean_log_conc, max_log_conc, sd_log_conc,
      mean_monthly_sd, mean_weekly_sd, mean_sd_within_samples
    ),
    decimals = 2
  ) %>%
  gt::cols_label(
    species = "Species",
    n_regions = "N regions",
    n_samples = "N samples",
    n_sample_detected = "N sample detections",
    sample_detection_rate = "Sample detection rate",
    n_reps = "N replicates",
    n_pos_reps = "N replicate detections",
    rep_detection_rate = "Replicate detection rate",
    mean_log_conc = "Mean log conc. (positive only)",
    max_log_conc = "Max log conc.",
    sd_log_conc = "Annual SD",
    mean_monthly_sd = "Monthly SD",
    mean_weekly_sd = "Weekly SD",
    mean_sd_within_samples = "Within-sample SD"
  ) %>%
  gt::tab_header(
    title = "Summary of eDNA detection and concentration variability by species"
  ) %>%
  gt::tab_options(
    table.font.size = gt::px(13),
    heading.title.font.size = gt::px(16),
    column_labels.font.size = gt::px(14),
    data_row.padding = gt::px(6),
    column_labels.padding = gt::px(2),
    heading.padding = gt::px(4),
    table.width = gt::pct(90)
  )

#####################
#new version of region with samples x species
#####################


# --- restrict to species actually tested in each region ---
region_species_tested <- dfRawClean %>%
  distinct(region, species)

all_samples <- dfRawClean %>%
  distinct(materialSampleID, region, date, month)

sample_species_grid <- all_samples %>%
  left_join(region_species_tested, by = "region")

# --- collapse to sample × species ---
df_full <- sample_species_grid %>%
  left_join(
    dfRawClean %>%
      group_by(materialSampleID, species) %>%
      summarise(
        logConc = mean(logConc, na.rm = TRUE),
        sample_detected = as.integer(any(detected == 1)),
        .groups = "drop"
      ),
    by = c("materialSampleID", "species")
  ) %>%
  mutate(
    sample_detected = replace_na(sample_detected, 0),
    logConc = replace_na(logConc, 0)
  )

# --- positive-only data ---
df_pos <- df_full %>%
  filter(sample_detected == 1)

# --- sample-level detection ---
sample_detection_region <- df_full %>%
  group_by(region) %>%
  summarise(
    n_obs = n(),
    n_samples = n_distinct(materialSampleID),
    n_species = n_distinct(species),
    n_sample_detected = sum(sample_detected == 1),
    sample_detection_rate = mean(sample_detected == 1),
    .groups = "drop"
  )

# --- replicate-level detection ---
rep_detection_region <- dfRawClean %>%
  group_by(region) %>%
  summarise(
    n_reps = n(),
    n_pos_reps = sum(detected == 1, na.rm = TRUE),
    rep_detection_rate = mean(detected == 1, na.rm = TRUE),
    .groups = "drop"
  )

# --- monthly SD (positive only) ---
monthly_sd_region <- df_pos %>%
  dplyr::group_by(region, month) %>%
  dplyr::summarise(
    monthly_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    mean_monthly_sd = mean(monthly_sd, na.rm = TRUE),
    .groups = "drop"
  )

# --- weekly SD (positive only) ---
weekly_sd_region <- df_pos %>%
  dplyr::mutate(
    week = lubridate::floor_date(date, "week")
  ) %>%
  dplyr::group_by(region, week) %>%
  dplyr::summarise(
    weekly_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    mean_weekly_sd = mean(weekly_sd, na.rm = TRUE),
    .groups = "drop"
  )

# --- within-sample SD (positive only) ---
sample_sd_region <- dfRawClean %>%
  dplyr::filter(detected == 1) %>%
  dplyr::group_by(region, materialSampleID, species) %>%
  dplyr::summarise(
    sample_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    mean_sd_within_samples = mean(sample_sd, na.rm = TRUE),
    .groups = "drop"
  )
# --- concentration stats (positive only) ---
summary_region_sd <- df_pos %>%
  group_by(region) %>%
  summarise(
    mean_log_conc = mean(logConc, na.rm = TRUE),
    max_log_conc  = max(logConc, na.rm = TRUE),
    sd_log_conc   = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  )

# --- combine ---
summary_region_final <- sample_detection_region %>%
  dplyr::left_join(rep_detection_region, by = "region") %>%
  dplyr::left_join(summary_region_sd, by = "region") %>%
  dplyr::left_join(monthly_sd_region, by = "region") %>%
  dplyr::left_join(weekly_sd_region, by = "region") %>%
  dplyr::left_join(sample_sd_region, by = "region") %>%
  dplyr::arrange(region)


summary_region_final %>%
  dplyr::select(-n_samples) %>%
  dplyr::rename(
    n_samples = n_obs
  ) %>%
  dplyr::select(
    region,
    n_species,
    n_samples,
    n_sample_detected,
    sample_detection_rate,
    n_reps,
    n_pos_reps,
    rep_detection_rate,
    mean_log_conc,
    max_log_conc,
    sd_log_conc,
    mean_monthly_sd,
    mean_weekly_sd,
    mean_sd_within_samples
  ) %>%
  gt::gt() %>%
  gt::fmt_integer(
    columns = c(n_species, n_samples, n_sample_detected, n_reps, n_pos_reps)
  ) %>%
  gt::fmt_percent(
    columns = c(sample_detection_rate, rep_detection_rate),
    decimals = 1
  ) %>%
  gt::fmt_number(
    columns = c(
      mean_log_conc, max_log_conc, sd_log_conc,
      mean_monthly_sd, mean_weekly_sd, mean_sd_within_samples
    ),
    decimals = 2
  ) %>%
  gt::cols_label(
    region = "Region",
    n_species = "N species",
    n_samples = "N samples",
    n_sample_detected = "N sample detections",
    sample_detection_rate = "Sample detection rate",
    n_reps = "N replicates",
    n_pos_reps = "N replicate detections",
    rep_detection_rate = "Replicate detection rate",
    mean_log_conc = "Mean log conc. (positive only)",
    max_log_conc = "Max log conc.",
    sd_log_conc = "Annual SD",
    mean_monthly_sd = "Monthly SD",
    mean_weekly_sd = "Weekly SD",
    mean_sd_within_samples = "Within-sample SD"
  ) %>%
  gt::tab_header(
    title = "Summary of eDNA detection and concentration variability by region"
  ) %>%
  gt::tab_options(
    table.font.size = gt::px(13),
    heading.title.font.size = gt::px(16),
    column_labels.font.size = gt::px(14),
    data_row.padding = gt::px(6),
    column_labels.padding = gt::px(2),
    heading.padding = gt::px(4),
    table.width = gt::pct(90)
  )



###########################
#all species x region
###########################

# --- sample-level detection (species × region) ---
sample_detection_sr <- df_full %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    n_obs = dplyr::n(),
    n_samples = dplyr::n_distinct(materialSampleID),
    n_sample_detected = sum(sample_detected == 1),
    sample_detection_rate = mean(sample_detected == 1),
    .groups = "drop"
  )

# --- replicate-level detection ---
rep_detection_sr <- dfRawClean %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    n_reps = dplyr::n(),
    n_pos_reps = sum(detected == 1, na.rm = TRUE),
    rep_detection_rate = mean(detected == 1, na.rm = TRUE),
    .groups = "drop"
  )

# --- monthly SD ---
monthly_sd_sr <- df_pos %>%
  dplyr::group_by(species, region, month) %>%
  dplyr::summarise(
    monthly_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    mean_monthly_sd = mean(monthly_sd, na.rm = TRUE),
    .groups = "drop"
  )

# --- weekly SD ---
weekly_sd_sr <- df_pos %>%
  dplyr::mutate(
    week = lubridate::floor_date(date, "week")
  ) %>%
  dplyr::group_by(species, region, week) %>%
  dplyr::summarise(
    weekly_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    mean_weekly_sd = mean(weekly_sd, na.rm = TRUE),
    .groups = "drop"
  )

# --- within-sample SD ---
sample_sd_sr <- dfRawClean %>%
  dplyr::filter(detected == 1) %>%
  dplyr::group_by(species, region, materialSampleID) %>%
  dplyr::summarise(
    sample_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    mean_sd_within_samples = mean(sample_sd, na.rm = TRUE),
    .groups = "drop"
  )

# --- concentration stats ---
summary_sd_sr <- df_pos %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    mean_log_conc = mean(logConc, na.rm = TRUE),
    max_log_conc  = max(logConc, na.rm = TRUE),
    sd_log_conc   = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  )

# --- combine ---
summary_sr_final <- sample_detection_sr %>%
  dplyr::left_join(rep_detection_sr, by = c("species", "region")) %>%
  dplyr::left_join(summary_sd_sr, by = c("species", "region")) %>%
  dplyr::left_join(monthly_sd_sr, by = c("species", "region")) %>%
  dplyr::left_join(weekly_sd_sr, by = c("species", "region")) %>%
  dplyr::left_join(sample_sd_sr, by = c("species", "region")) %>%
  dplyr::arrange(species, region)

summary_sr_final_clean <- summary_sr_final %>%
  dplyr::arrange(species, region) %>%
  dplyr::mutate(
    species = as.character(species),
    species_display = dplyr::if_else(
      species == dplyr::lag(species, default = ""),
      "",
      species
    )
  )

summary_sr_final_clean %>%
  dplyr::select(-n_samples) %>%
  dplyr::rename(n_samples = n_obs) %>%
  dplyr::select(
    species_display, region,
    n_samples,
    n_sample_detected,
    sample_detection_rate,
    n_reps,
    n_pos_reps,
    rep_detection_rate,
    mean_log_conc,
    max_log_conc,
    sd_log_conc,
    mean_monthly_sd,
    mean_weekly_sd,
    mean_sd_within_samples
  ) %>%
  gt::gt() %>%
  gt::fmt_integer(
    columns = c(n_samples, n_sample_detected, n_reps, n_pos_reps)
  ) %>%
  gt::cols_label(
    species_display = "Species",
    region = "Region",
    n_samples = "N samples",
    n_sample_detected = "N sample detections",
    sample_detection_rate = "Sample detection rate",
    n_reps = "N replicates",
    n_pos_reps = "N replicate detections",
    rep_detection_rate = "Replicate detection rate",
    mean_log_conc = "Mean log conc. (positive only)",
    max_log_conc = "Max log conc.",
    sd_log_conc = "Annual SD",
    mean_monthly_sd = "Monthly SD",
    mean_weekly_sd = "Weekly SD",
    mean_sd_within_samples = "Within-sample SD"
  ) %>%
  gt::cols_align(
    align = "center",
    everything()
  ) %>%
  gt::cols_width(
    everything() ~ gt::px(90)   # key line
  ) %>%
  gt::fmt_percent(
    columns = c(sample_detection_rate, rep_detection_rate),
    decimals = 1
  ) %>%
  gt::fmt_number(
    columns = c(
      mean_log_conc, max_log_conc, sd_log_conc,
      mean_monthly_sd, mean_weekly_sd, mean_sd_within_samples
    ),
    decimals = 2
  ) %>%
  gt::tab_header(
    title = "Summary of eDNA detection and concentration variability by species and region"
  ) %>%
  gt::cols_width(
    species_display ~ gt::px(95),
    region ~ gt::px(45),
    n_samples ~ gt::px(55),
    n_sample_detected ~ gt::px(65),
    sample_detection_rate ~ gt::px(70),
    n_reps ~ gt::px(55),
    n_pos_reps ~ gt::px(65),
    rep_detection_rate ~ gt::px(70),
    mean_log_conc ~ gt::px(75),
    max_log_conc ~ gt::px(55),
    sd_log_conc ~ gt::px(55),
    mean_monthly_sd ~ gt::px(55),
    mean_weekly_sd ~ gt::px(55),
    mean_sd_within_samples ~ gt::px(65)
  ) %>%
  gt::tab_options(
    table.font.size = gt::px(11),
    heading.title.font.size = gt::px(13),
    column_labels.font.size = gt::px(11),
    data_row.padding = gt::px(3),
    column_labels.padding = gt::px(1),
    heading.padding = gt::px(2),
    table.width = gt::pct(100)
  )
