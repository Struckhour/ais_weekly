source('./AIS_eDNA_data_prep.R')


library(dplyr)
library(ggplot2)
library(forcats)
library(lubridate)

# Make sure "GOM" is included as a level
plateAbundance <- plateAbundance %>%
  mutate(region = fct_expand(region, "GOM"),
         region = replace_na(region, "GOM")) %>%
  filter(date >= as.Date("2023-06-01") & date <= as.Date("2024-05-31")) %>%
  mutate(month_letter = substr(month.abb[month(date)], 1, 1))


hybrid_color <- c("#00A08A", "#446455", "#CCAA4F", "#5BBCD6", "#fb8072")


#Log
ggplot(plateAbundance, aes(x = date, y = Avg + 0.01, color = region)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  scale_color_manual(values = hybrid_color) +
  scale_y_log10() +  # log scale if abundances vary widely
  facet_grid(region ~ species) +  # rows = regions, columns = species
  scale_x_date(
    date_breaks = "1 month",
    labels = function(d) substr(month.abb[month(d)], 1, 1)
  ) +
  theme_classic() +
  labs(
    title = "Plate Abundance Over Time by Species and Region",
    x = "Month",
    y = "Log10(Avg + 0.01)"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  )

#NOT LOG
ggplot(plateAbundance, aes(x = date, y = Avg, color = region)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  scale_color_manual(values = hybrid_color) +
  facet_grid(region ~ species) +
  scale_x_date(
    date_breaks = "1 month",
    labels = function(d) substr(month.abb[month(d)], 1, 1)
  ) +
  theme_classic() +
  labs(
    title = "Plate Abundance Over Time by Species and Region",
    x = "Month",
    y = "Avg Abundance"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  )













library(dplyr)
library(tidyr)
library(ggplot2)
library(lubridate)
library(gt)

# overlapping species and regions
keep_species <- intersect(unique(dfRawClean$species), unique(plateAbundance$species))
keep_regions <- intersect(as.character(unique(dfRawClean$region)),
                          as.character(unique(plateAbundance$region)))

# qPCR summary
qpcr_plot_df <- dfRawClean %>%
  filter(
    species %in% keep_species,
    as.character(region) %in% keep_regions
  ) %>%
  group_by(region, species, date) %>%
  summarise(
    value = mean(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(week = floor_date(date, "week"))

# Plate summary
plate_plot_df <- plateAbundance %>%
  filter(
    species %in% keep_species,
    as.character(region) %in% keep_regions
  ) %>%
  group_by(region, species, date) %>%
  summarise(
    value = mean(Avg, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(week = floor_date(date, "week"))

# weekly aggregated
qpcr_weekly <- qpcr_plot_df %>%
  group_by(region, species, week) %>%
  summarise(
    value = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

plate_weekly <- plate_plot_df %>%
  group_by(region, species, week) %>%
  summarise(
    value = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

# normalize by overall max only, within each source
qpcr_max <- max(qpcr_weekly$value, na.rm = TRUE)
plate_max <- max(plate_weekly$value, na.rm = TRUE)

qpcr_weekly <- qpcr_weekly %>%
  mutate(
    value_norm = value / qpcr_max,
    source = "qPCR logConc"
  )

plate_weekly <- plate_weekly %>%
  mutate(
    value_norm = value / plate_max,
    source = "Plate Avg"
  )

# for correlations
corr_df_weekly <- full_join(
  qpcr_weekly,
  plate_weekly,
  by = c("region", "species", "week"),
  suffix = c("_qpcr", "_plate")
)

corr_results_weekly <- corr_df_weekly %>%
  group_by(region, species) %>%
  summarise(
    n = sum(!is.na(value_norm_qpcr) & !is.na(value_norm_plate)),
    cor = cor(value_norm_qpcr, value_norm_plate, use = "complete.obs"),
    p_value = if (n > 2) {
      cor.test(value_norm_qpcr, value_norm_plate)$p.value
    } else {
      NA_real_
    },
    .groups = "drop"
  )

# gt table
corr_results_weekly %>%
  arrange(species, region) %>%
  group_by(species) %>%
  mutate(
    species_display = ifelse(row_number() == 1, species, "")
  ) %>%
  ungroup() %>%
  select(
    species = species_display,
    region,
    n,
    cor,
    p_value
  ) %>%
  gt() %>%
  tab_header(
    title = "Correlation between Plate and qPCR Abundance"
  ) %>%
  fmt_number(
    columns = c(cor, p_value),
    decimals = 3
  ) %>%
  cols_label(
    species = "Species",
    region = "Region",
    n = "N",
    cor = "Correlation",
    p_value = "p-value"
  ) %>%
  tab_style(
    style = cell_text(color = "black"),
    locations = cells_body()
  ) %>%
  tab_options(
    table.background.color = "white",
    heading.background.color = "white",
    column_labels.background.color = "white"
  )

# recreate plot_df for plotting
plot_df <- bind_rows(
  qpcr_weekly %>%
    transmute(region, species, date = week, value_norm, source),
  plate_weekly %>%
    transmute(region, species, date = week, value_norm, source)
)

# plot
ggplot(plot_df, aes(x = date, y = value_norm, color = source)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_line(alpha = 0.8) +
  # geom_smooth(se = FALSE, method = "loess", span = 0.75) +
  facet_grid(region ~ species) +
  scale_x_date(
    date_breaks = "1 month",
    labels = function(d) substr(month.abb[month(d)], 1, 1)
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  theme_classic() +
  labs(
    title = "Normalized qPCR and Plate Abundance Over Time",
    x = "Month",
    y = "Normalized abundance (0–1)",
    color = "Source"
  )
