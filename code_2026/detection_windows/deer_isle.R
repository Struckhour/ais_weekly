source('./AIS_eDNA_data_prep.R')

gom_monthly_station <- dfRawClean %>%
  dplyr::filter(region == "GOM") %>%
  dplyr::group_by(station, decimalLatitude, decimalLongitude, species, month) %>%
  dplyr::summarise(
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
  gom_monthly_station,
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
    title = "Monthly mean qPCR concentration across GOM stations",
    x = "Month",
    y = "Mean log10(concentration + 1)",
    color = "Station"
  )

