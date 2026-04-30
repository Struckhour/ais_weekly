source('./AIS_eDNA_data_prep.R')

View(metadata)


#SALINITY

metadata_clean <- metadata %>%
  mutate(
    region = factor(region)
  ) %>%
  filter(!is.na(salinity_ppt))         # remove rows with no salinity data


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
  geom_smooth(method = "loess", se = FALSE, span = 0.1) +
  scale_color_manual(values = hybrid_color) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b"
  ) +
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
  geom_smooth(method = "loess", se = FALSE, span = 0.1) +
  scale_color_manual(values = hybrid_color) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b"
  ) +
  theme_classic() +
  labs(
    title = "Temperature Over Time",
    x = "Date",
    y = "Temp (C)",
    color = "Region"
  )





#pH

library(tidyverse)
library(lubridate)

# Clean and prepare temperature data
metadata_pH <- metadata %>%
  mutate(
    pH = as.numeric(pH),  # convert to numeric
    region = factor(region)                 # only levels that exist
  ) %>%
  filter(!is.na(pH))              # remove missing temperatures


# metadata_pH <- metadata_pH %>% filter(region == "BOF")

# Plot temperature over time by region
ggplot(metadata_pH, aes(x = date, y = pH, color = region)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE) +
  theme_classic() +
  labs(
    title = "pH Over Time by Region",
    x = "Date",
    y = "pH",
    color = "Region"
  )







##################
#DATA COVERAGE
##################

library(dplyr)
library(tidyr)
library(ggplot2)

coverage_df <- metadata %>%
  dplyr::select(
    region,
    date,
    waterTemp_C,
    salinity_ppt,
    pH,
    turbidity_NTU,
    chlorophyll
  ) %>%
  dplyr::mutate(
    dplyr::across(
      c(waterTemp_C, salinity_ppt, pH, turbidity_NTU, chlorophyll),
      as.character
    )
  ) %>%
  tidyr::pivot_longer(
    cols = c(waterTemp_C, salinity_ppt, pH, turbidity_NTU, chlorophyll),
    names_to = "property",
    values_to = "value"
  ) %>%
  dplyr::mutate(
    value = na_if(value, "N/A"),
    value = na_if(value, ""),
    property = dplyr::recode(
      property,
      waterTemp_C   = "Temperature",
      salinity_ppt  = "Salinity",
      pH            = "pH",
      turbidity_NTU = "Turbidity",
      chlorophyll   = "Chlorophyll"
    ),
    region_property = paste(region, property, sep = " - ")
  ) %>%
  dplyr::filter(!is.na(value))


ggplot(coverage_df, aes(x = date, y = property, color = property)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(
    values = c(
      "Temperature" = "#e66101",
      "Salinity" = "#1f78b4",
      "pH" = "#984ea3",
      "Turbidity" = "#4daf4a",
      "Chlorophyll" = "#a65628"
    )
  ) +
  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b\n%Y"
  ) +
  facet_grid(region ~ ., scales = "free_y", space = "free_y") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "grey70", fill = NA, linewidth = 0.5),
    panel.spacing.y = unit(0.6, "lines"),
    panel.grid.minor = element_blank()
  ) +
  labs(
    x = "Date",
    y = NULL,
    title = "Environmental data coverage by region and property"
  )
