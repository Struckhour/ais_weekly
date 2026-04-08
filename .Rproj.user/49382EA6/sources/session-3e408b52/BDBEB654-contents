source('./AIS_eDNA_data_prep.R')







library(ggplot2)
library(dplyr)
library(forcats)
library(lubridate)

# Prepare data: add "GOM", filter dates, and add first-letter month column
dfCor_withGOM <- dfRawClean %>%
  filter(date >= as.Date("2013-06-01") & date <= as.Date("2024-05-31")) %>%
  mutate(month_letter = substr(month.abb[month(date)], 1, 1))

# Plot with actual dates on x-axis but labeled with month_letter
ggplot(dfCor_withGOM, aes(x = date, y = concentration + 1, color = region)) +
  geom_point(size = 2, alpha = 0.7) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  scale_y_log10() +
  facet_grid(region ~ species) +  # rows = regions, columns = species
  scale_x_date(
    date_breaks = "1 month",      # keep spacing by actual date
    labels = function(d) substr(month.abb[month(d)], 1, 1)  # first-letter labels
  ) +
  theme_classic() +
  labs(
    title = "eDNA Concentration Over Time by Species and Region",
    x = "Month",
    y = "Log10(Concentration + 1)"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  )
