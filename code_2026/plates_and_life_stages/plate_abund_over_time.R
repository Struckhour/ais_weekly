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

ggplot(plateAbundance, aes(x = date, y = Avg + 0.01, color = region)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
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
