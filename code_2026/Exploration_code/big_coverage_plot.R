source('./AIS_eDNA_data_prep.R')

meta_cov <- metadata %>%
  dplyr::select(region, date, waterTemp_C, salinity_ppt, pH) %>%
  dplyr::mutate(
    dplyr::across(c(waterTemp_C, salinity_ppt, pH), as.character)
  ) %>%
  tidyr::pivot_longer(
    cols = c(waterTemp_C, salinity_ppt, pH),
    names_to = "property",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    value = dplyr::na_if(value, "N/A"),
    value = dplyr::na_if(value, ""),
    property = dplyr::recode(
      property,
      waterTemp_C  = "Temperature",
      salinity_ppt = "Salinity",
      pH           = "pH"
    )
  ) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::distinct(region, date, property) %>%
  dplyr::mutate(
    track = paste("ENV", property, sep = " - "),
    source = "Environment"
  ) %>%
  dplyr::select(region, date, track, source)

life_cov <- plateAbundance %>%
  dplyr::filter(date >= as.Date("2022-06-01") & date <= as.Date("2028-05-31")) %>%
  tidyr::separate_rows(State, sep = ",\\s*") %>%
  dplyr::mutate(State = trimws(State)) %>%
  dplyr::filter(State != "") %>%
  dplyr::distinct(region, date, species, State) %>%
  dplyr::mutate(
    track = paste("LIFE", species, State, sep = " - "),
    source = "Life stage"
  ) %>%
  dplyr::select(region, date, track, source)

dfRawClean <- dfRawClean %>%
  filter(date >= as.Date("2013-06-01") & date <= as.Date("2024-10-31"))
qpcr_cov <- dfRawClean %>%
  dplyr::filter(!is.na(concentration)) %>%
  dplyr::distinct(region, date, species) %>%
  dplyr::mutate(
    track = paste("qPCR", species, sep = " - "),
    source = "qPCR"
  ) %>%
  dplyr::select(region, date, track, source)


super_coverage <- dplyr::bind_rows(
  meta_cov,
  life_cov,
  qpcr_cov
)

super_coverage <- super_coverage %>%
  dplyr::mutate(
    region = factor(region, levels = c("MAG", "PEI", "HAL", "BOF", "GOM"))
  )

super_coverage <- super_coverage %>%
  dplyr::mutate(
    track = factor(track, levels = unique(track))
  )


ggplot(super_coverage, aes(x = date, y = track, color = source)) +
  geom_point(size = 1.8, alpha = 0.8) +
  facet_grid(region ~ ., scales = "free_y", space = "free_y") +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b\n%Y"
  ) +
  scale_color_manual(
    values = c(
      "Environment" = "#1f78b4",
      "Life stage"  = "#e66101",
      "qPCR"        = "#984ea3"
    )
  ) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
    panel.spacing.y = unit(0.6, "lines"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank()
  ) +
  labs(
    title = "Data coverage by region across environmental, life-stage, and qPCR sources",
    x = "Date",
    y = NULL,
    color = NULL
  )
