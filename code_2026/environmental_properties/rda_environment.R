library(dplyr)
library(tidyr)
library(lubridate)

qpcr_weekly <- dfRawClean %>%
  group_by(region, species, week = floor_date(date, "week")) %>%
  summarise(value = mean(logConc, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = species,
    values_from = value
  )

env_weekly <- dfRawClean %>%
  group_by(region, week = floor_date(date, "week")) %>%
  summarise(
    temp = mean(temp, na.rm = TRUE),
    sal  = mean(sal, na.rm = TRUE),
    .groups = "drop"
  )

plate_weekly <- plateAbundance %>%
  group_by(region, week = floor_date(date, "week")) %>%
  summarise(
    plate_total = mean(Avg, na.rm = TRUE),
    .groups = "drop"
  )

rda_df <- qpcr_weekly %>%
  left_join(env_weekly, by = c("region", "week")) %>%
  left_join(plate_weekly, by = c("region", "week")) %>%
  drop_na()

Y <- rda_df %>%
  select(all_of(species_order))   # your 5 species



library(vegan)
Y <- decostand(Y, method = "hellinger")
X <- scale(X)

rda(Y ~ temp + sal + plate_total, data = as.data.frame(X))
summary(rda_model)
plot(rda_model)
