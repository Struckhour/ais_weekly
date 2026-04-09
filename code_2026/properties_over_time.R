View(metadata)


#SALINITY

metadata_clean <- metadata %>%
  mutate(
    salinity_ppt = na_if(salinity_ppt, "N/A"),
    salinity_ppt = as.numeric(salinity_ppt),
    region = factor(region)              # converts to factor with only present levels
  ) %>%
  filter(!is.na(salinity_ppt))          # remove rows with no salinity data


ggplot(metadata_clean, aes(x = date, y = salinity_ppt)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  facet_wrap(~ region, ncol = 1) +
  theme_classic() +
  labs(
    title = "Salinity Over Time by Region",
    x = "Date",
    y = "Salinity (ppt)"
  )

ggplot(metadata_clean, aes(x = date, y = salinity_ppt, color = region)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_classic() +
  labs(
    title = "Salinity Over Time",
    x = "Date",
    y = "Salinity (ppt)",
    color = "Region"
  )



















#TEMPERATURE

library(tidyverse)
library(lubridate)

# Clean and prepare temperature data
metadata_temp <- metadata %>%
  mutate(
    waterTemp_C = as.numeric(waterTemp_C),  # convert to numeric
    region = factor(region)                 # only levels that exist
  ) %>%
  filter(!is.na(waterTemp_C))              # remove missing temperatures

# Plot temperature over time by region
ggplot(metadata_temp, aes(x = date, y = waterTemp_C, color = region)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_classic() +
  labs(
    title = "Water Temperature Over Time by Region",
    x = "Date",
    y = "Temperature (°C)",
    color = "Region"
  )

