new_data <- read_xlsx(
  "./data/BIO_Stn24_SBE_CT_2023_2024_dateorderd.xlsx",
  na = ""
)



hal_met <- metadata %>% filter(region == "HAL")

hal_met <- hal_met %>% select(date, year, month, region, salinity_ppt, waterTemp_C)

new_data <- new_data %>% select(date, temp, PSU, ts) %>%
  filter(ts >= as.Date("2023-04-27") & ts <= as.Date("2024-06-3"))



hal_met_clean <- hal_met %>%
  mutate(
    date = as.Date(date),
    waterTemp_C = as.numeric(waterTemp_C),
    salinity_ppt = as.numeric(salinity_ppt)
  ) %>%
  filter(!is.na(waterTemp_C) | !is.na(salinity_ppt))

new_data_clean <- new_data %>%
  mutate(
    date = as.Date(date)
  )



sbe_daily <- new_data_clean %>%
  group_by(date) %>%
  summarise(
    temp = mean(temp, na.rm = TRUE),
    PSU = mean(PSU, na.rm = TRUE),
    .groups = "drop"
  )



library(tidyr)

hal_long <- hal_met_clean %>%
  select(date, waterTemp_C, salinity_ppt) %>%
  pivot_longer(
    cols = c(waterTemp_C, salinity_ppt),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(source = "old")

sbe_long <- sbe_daily %>%
  rename(
    waterTemp_C = temp,
    salinity_ppt = PSU
  ) %>%
  pivot_longer(
    cols = c(waterTemp_C, salinity_ppt),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(source = "new")

plot_df <- bind_rows(hal_long, sbe_long)

ggplot(plot_df, aes(x = date, y = value,
                    color = source,
                    shape = source,
                    linetype = source,
                    group = interaction(variable, source))) +

  geom_line() +
  geom_point(size = 1.3) +

  facet_wrap(~ variable, scales = "free_y", ncol = 1,
             labeller = as_labeller(c(
               waterTemp_C = "Temperature (°C)",
               salinity_ppt = "Salinity"
             ))) +

  scale_color_manual(
    values = c(
      new = "#66c2a5",
      old = "#fc8d62"
    )
  ) +

  scale_x_date(
    date_breaks = "1 month",
    date_labels = "%b"
  ) +

  labs(
    x = "Date",
    y = NULL,
    color = "Source",
    shape = "Source",
    linetype = "Source",
    title = "Temperature and Salinity Comparison",
    subtitle = "Old data vs new data"
  ) +

  theme_minimal()
