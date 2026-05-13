source('./AIS_eDNA_data_prep.R')

library(ggplot2)
library(dplyr)
library(forcats)
library(lubridate)

species_order <- c(
  "Membranipora membranacea",
  "Botrylloides violaceus",
  "Didemnum vexillum",
  "Ciona intestinalis",
  "Carcinus maenas"
)



color_vec <- c("#00A08A", "#446455", "#Fdd262", "#5BBCD6", "#046c9a", "#ABDDDE", "#d3dddc")

other_color_vec <- c("#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3")

hybrid_color <- c(
  MAG = "#00A08A",
  PEI = "#446455",
  HAL = "#CCAA4F",
  BOF = "#5BBCD6",
  GOM = "#fb8072"
)


################
# ALL REGIONS — raw concentrations by actual sampling date
################

df_plot <- dfRawClean %>%
  dplyr::mutate(
    species = factor(species, levels = species_order),
    region = factor(region, levels = c("MAG", "PEI", "HAL", "BOF", "GOM"))
  ) %>%
  filter(
    date >= as.Date("2023-01-01"),
    date <= as.Date("2024-08-30")
  )

ggplot(
  df_plot,
  aes(
    x = date,
    y = concentration + 1,
    color = region
  )
) +
  geom_point(size = 0.8, alpha = 0.75) +
  geom_smooth(
    aes(group = 1),
    method = "loess",
    se = TRUE,
    color = "black",
    fill = "black",     # controls ribbon color
    alpha = 0.35,        # controls transparency
    linewidth = 0.6,
    span = 0.35
  ) +
  scale_color_manual(values = hybrid_color) +
  scale_y_log10() +
  facet_grid(region ~ species) +
  scale_x_date(
    date_breaks = "1 month",
    labels = function(d) substr(month.abb[lubridate::month(d)], 1, 1)
  ) +
  theme_classic() +
  labs(
    title = "eDNA Concentration Over Time by Species and Region",
    x = "Month",
    y = "Log10(Concentration + 1)",
    color = "Region"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 9),
    strip.background = element_rect(color = "black", fill = "white"),
    strip.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(0, "lines")
  )


################
# ALL REGIONS — raw concentrations over full sampling year
################

df_plot <- dfRawClean %>%
  dplyr::mutate(
    species = factor(species, levels = species_order),
    region = factor(region, levels = c("MAG", "PEI", "HAL", "BOF", "GOM")),
    month_plot_label = substr(month.abb[month], 1, 1)
  )

ggplot(
  df_plot,
  aes(
    x = month_reordered,
    y = concentration + 1,
    color = region
  )
) +
  geom_point(size = 0.8, alpha = 0.75) +
  geom_smooth(
    aes(group = 1),
    method = "loess",
    se = FALSE,
    color = "black",
    linewidth = 0.8,
    span = 0.45
  ) +
  scale_color_manual(values = hybrid_color) +
  scale_y_log10() +
  facet_grid(region ~ species) +
  scale_x_continuous(
    breaks = 1:12,
    labels = function(x) substr(month.abb[((x - 1) %% 12) + 1], 1, 1)
  ) +
  theme_classic() +
  labs(
    title = "eDNA Concentration Over Time by Species and Region",
    x = "Month",
    y = "Log10(Concentration + 1)",
    color = "Region"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 7),
    strip.background = element_rect(color = "black", fill = "white"),
    strip.text = element_text(size = 7),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(0, "lines")
  )




################
#ALL REGIONS zoomed in on Halifax rainfall event
################
df_weekly <- dfRawClean %>%
  filter(
    species != "Didemnum vexillum",
    date >= as.Date("2023-06-01"),
    date <= as.Date("2023-12-30")
  ) %>%
  mutate(
    species = factor(species, levels = species_order),
    week = lubridate::floor_date(date, "week"),
    alpha_flag = ifelse(region == "HAL", "HAL", "other")
  ) %>%
  group_by(region, species, week, alpha_flag) %>%
  summarise(
    mean_conc = mean(concentration, na.rm = TRUE),
    .groups = "drop"
  )

ggplot(df_weekly, aes(x = week, y = mean_conc + 1,
                      color = region, group = region, alpha = alpha_flag)) +
  geom_point(size = 2) +
  geom_line() +
  geom_vline(
    xintercept = as.Date("2023-07-20"),
    color = "red",
    linewidth = 1,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = as.Date("2023-10-22"),
    color = "red",
    linewidth = 1,
    linetype = "dashed"
  ) +
  scale_alpha_manual(values = c(HAL = 1, other = 0.2)) +
  scale_color_manual(values = hybrid_color) +

  scale_y_log10() +
  facet_grid(. ~ species) +

  scale_x_date(
    date_breaks = "1 month",
    labels = function(d) substr(month.abb[lubridate::month(d)], 1, 1)
  ) +

  theme_classic() +
  labs(
    title = "Weekly Mean eDNA Concentration by Region",
    x = "Month",
    y = "Log10(Mean Concentration + 1)",
    color = "Region"
  ) +

  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 0, vjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
    panel.spacing = unit(0, "lines")
  )



















# Prepare data: add "GOM", filter dates, and add first-letter month column
dfRaw_window <- dfRawClean %>%
  # filter(date >= as.Date("2013-06-01") & date <= as.Date("2024-10-31")) %>%
  filter(date >= as.Date("2023-07-01") & date <= as.Date("2023-09-30")) %>%
  mutate(month_letter = substr(month.abb[month(date)], 1, 1))

# Plot with actual dates on x-axis but labeled with month_letter
ggplot(dfRaw_window, aes(x = date, y = concentration + 1, color = region)) +
  geom_point(size = 1, alpha = 0.5) +
  geom_smooth(method = "loess", se = TRUE, color = "black") +
  scale_color_manual(values = hybrid_color) +
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


###########################
#MONTHLY COMPARISON OF NORMALIZED CONCENTRATION AND DETECTION RATES
###########################
library(dplyr)
library(tidyr)
library(ggplot2)

df_plot <- dfMonths %>%
  select(region, species, month_reordered, normMeanLogConc, det_rate) %>%
  pivot_longer(
    cols = c(normMeanLogConc, det_rate),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    normMeanLogConc = "Normalized log concentration",
                    det_rate = "Detection rate"),
    month_label = substr(month.abb[month_reordered], 1, 1),
    # Set factor levels for regions to control facet row order
    region = factor(region, levels = c("MAG", "PEI", "HAL", "BOF", "GOM"))
  )

ggplot(df_plot, aes(x = month_reordered, y = value, group = metric, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +

  facet_grid(region ~ species) +

  scale_x_continuous(
    breaks = 1:12,
    labels = substr(month.abb, 1, 1)
  ) +

  scale_color_manual(
    values = c(
      "Normalized log concentration" = "steelblue",
      "Detection rate" = "darkorange"
    )
  ) +

  scale_y_continuous(limits = c(0, 1)) +

  theme_classic() +
  labs(
    title = "Detection Rate vs Normalized Log Concentration",
    x = "Month",
    y = "Scaled Value",
    color = NULL
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  )




####################
#JUST THE POSITIVE ABUNDANCES IN NORMALIZED MEAN CONCENTRATION
####################
df_plot <- dfMonths %>%
  select(region, species, month_reordered, normPosMeanLogConc, det_rate) %>%
  pivot_longer(
    cols = c(normPosMeanLogConc, det_rate),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    normPosMeanLogConc = "Normalized Positive log concentration",
                    det_rate = "Detection rate"),
    month_label = substr(month.abb[month_reordered], 1, 1),
    region = factor(region, levels = c("MAG", "PEI", "HAL", "BOF", "GOM")),
    species = factor(species, levels = species_order)
  )

ggplot(df_plot, aes(x = month_reordered, y = value, group = metric, color = metric)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +

  facet_grid(region ~ species) +

  scale_x_continuous(
    breaks = 1:12,
    labels = substr(month.abb, 1, 1)
  ) +

  scale_color_manual(
    values = c(
      "Normalized Positive log concentration" = "steelblue",
      "Detection rate" = "darkorange"
    )
  ) +

  scale_y_continuous(limits = c(0, 1)) +

  theme_classic() +
  labs(
    title = "Detection Rate vs Normalized Positive Log Concentration",
    x = "Month",
    y = "Scaled Value",
    color = NULL
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  )


########################
#ALL THREE
########################

df_plot <- dfMonths %>%
  select(region, species, month_original,
         normMeanLogConc, normPosMeanLogConc, det_rate) %>%
  pivot_longer(
    cols = c(normMeanLogConc, normPosMeanLogConc, det_rate),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    normMeanLogConc     = "Normalized log concentration",
                    normPosMeanLogConc  = "Normalized Positive log concentration",
                    det_rate            = "Detection rate"
    ),
    month_plot_label = factor(
      month_original,
      levels = 1:12,
      labels = month.abb
    ),
    region = factor(region, levels = c("MAG", "PEI", "HAL", "BOF", "GOM")),
    species = factor(species, levels = species_order)
  )
# Plot
ggplot(df_plot, aes(x = month_plot_label, y = value, color = metric)) +

  geom_line(
    data = df_plot %>% filter(metric != "Normalized Positive log concentration"),
    aes(group = interaction(metric, region, species)),
    linewidth = 1,
    alpha = 1
  ) +

  geom_line(
    data = df_plot %>% filter(metric == "Normalized Positive log concentration"),
    aes(group = interaction(metric, region, species)),
    linewidth = 1,
    alpha = 0.5
  ) +

  geom_point(aes(group = interaction(metric, region, species)), size = 2) +

  facet_grid(region ~ species) +

  scale_color_manual(
    values = c(
      "Normalized log concentration"           = "steelblue",
      "Normalized Positive log concentration" = "lightsteelblue",
      "Detection rate"                        = "darkorange"
    )
  ) +

  scale_y_continuous(limits = c(0, 1)) +

  scale_x_discrete(labels = c(
    "Jan" = "J",
    "Feb" = "F",
    "Mar" = "M",
    "Apr" = "A",
    "May" = "M",
    "Jun" = "J",
    "Jul" = "J",
    "Aug" = "A",
    "Sep" = "S",
    "Oct" = "O",
    "Nov" = "N",
    "Dec" = "D"
  )) +

  theme_classic() +
  labs(
    title = "Detection Rate vs Normalized Log Concentrations",
    x = "Month",
    y = "Scaled Value",
    color = NULL
  ) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 0, vjust = 0.5)
  )
