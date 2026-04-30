source('./AIS_eDNA_data_prep.R')

str(dfRawClean)

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
