source('./AIS_eDNA_data_prep.R')

bof_monthly_station <- dfRawClean %>%
  dplyr::filter(region == "BOF") %>%
  dplyr::group_by(station, species, month) %>%
  dplyr::summarise(
    decimalLatitude  = dplyr::first(decimalLatitude),
    decimalLongitude = dplyr::first(decimalLongitude),
    mean_log_conc = mean(logConc, na.rm = TRUE),
    sd_log_conc   = sd(logConc, na.rm = TRUE),
    n             = dplyr::n(),
    se_log_conc   = sd_log_conc / sqrt(n),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    species = factor(species, levels = species_order),
    month = factor(month, levels = 1:12, labels = month.abb)
  )


ggplot(
  bof_monthly_station,
  aes(x = month, y = mean_log_conc, color = station, group = station)
) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = mean_log_conc - se_log_conc,
      ymax = mean_log_conc + se_log_conc
    ),
    width = 0.2,
    alpha = 0.6
  ) +
  facet_wrap(~ species, ncol = 1, scales = "free_y") +
  theme_classic() +
  labs(
    title = "Monthly mean qPCR concentration across BOF stations",
    x = "Month",
    y = "Mean log10(concentration + 1)",
    color = "Station"
  )

































make_region_monthly_station <- function(df, region_name) {

  df %>%
    dplyr::filter(region == region_name) %>%
    dplyr::group_by(station, species, month) %>%
    dplyr::summarise(
      decimalLatitude  = dplyr::first(decimalLatitude),
      decimalLongitude = dplyr::first(decimalLongitude),
      mean_log_conc = mean(logConc, na.rm = TRUE),
      sd_log_conc   = sd(logConc, na.rm = TRUE),
      n             = dplyr::n(),
      se_log_conc   = sd_log_conc / sqrt(n),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      species = factor(species, levels = species_order),
      month = factor(month, levels = 1:12, labels = month.abb),
      region = region_name
    )
}

plot_region_monthly <- function(df_region) {

  ggplot(
    df_region,
    aes(x = month, y = mean_log_conc, color = station, group = station)
  ) +
    geom_line(linewidth = 1, alpha = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(
      aes(
        ymin = mean_log_conc - se_log_conc,
        ymax = mean_log_conc + se_log_conc
      ),
      width = 0.2,
      alpha = 0.6
    ) +
    facet_wrap(~ species, ncol = 1, scales = "free_y") +
    theme_classic() +
    labs(
      title = paste("Monthly mean qPCR concentration -", unique(df_region$region)),
      x = "Month",
      y = "Mean log10(concentration + 1)",
      color = "Station"
    )
}


mag_df <- make_region_monthly_station(dfRawClean, "MAG")
plot_region_monthly(mag_df)

pei_df <- make_region_monthly_station(dfRawClean, "PEI")
plot_region_monthly(pei_df)



hal_df <- make_region_monthly_station(dfRawClean, "HAL")
plot_region_monthly(hal_df)

bof_df <- make_region_monthly_station(dfRawClean, "BOF")
plot_region_monthly(bof_df)

gom_df <- make_region_monthly_station(dfRawClean, "GOM")
plot_region_monthly(gom_df)




























###################
#big summary plots with se bars
###################
region_monthly_df <- dfRawClean %>%
  dplyr::group_by(region, species, month) %>%
  dplyr::summarise(
    mean_log_conc = mean(logConc, na.rm = TRUE),
    sd_log_conc   = sd(logConc, na.rm = TRUE),
    n             = dplyr::n(),
    se_log_conc   = sd_log_conc / sqrt(n),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    species = factor(species, levels = species_order),
    region = factor(region, levels = region_order),
    month = factor(month, levels = 1:12, labels = month.abb)
  )

ggplot(
  region_monthly_df,
  aes(x = month, y = mean_log_conc, color = region, group = region)
) +
  geom_line(linewidth = 1, alpha = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(
    aes(
      ymin = mean_log_conc - se_log_conc,
      ymax = mean_log_conc + se_log_conc
    ),
    width = 0.2,
    alpha = 0.6
  ) +
  facet_wrap(~ species, ncol = 1, scales = "free_y") +
  scale_color_manual(values = region_colors) +
  theme_classic() +
  labs(
    title = "Monthly mean qPCR concentration by region",
    x = "Month",
    y = "Mean log10(concentration + 1)",
    color = "Region"
  )




