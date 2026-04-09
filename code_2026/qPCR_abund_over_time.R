source('./AIS_eDNA_data_prep.R')

library(ggplot2)
library(dplyr)
library(forcats)
library(lubridate)



# Prepare data: add "GOM", filter dates, and add first-letter month column
dfRaw_window <- dfRawClean %>%
  filter(date >= as.Date("2013-06-01") & date <= as.Date("2024-05-31")) %>%
  mutate(month_letter = substr(month.abb[month(date)], 1, 1))

# Plot with actual dates on x-axis but labeled with month_letter
ggplot(dfRaw_window, aes(x = date, y = concentration + 1, color = region)) +
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


# Prepare the data
df_plot <- dfMonths %>%
  select(region, species, month_reordered, normMeanLogConc, normPosMeanLogConc, det_rate) %>%
  pivot_longer(
    cols = c(normMeanLogConc, normPosMeanLogConc, det_rate),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(metric,
                    normMeanLogConc    = "Normalized log concentration",
                    normPosMeanLogConc = "Normalized Positive log concentration",
                    det_rate           = "Detection rate"),
    month_label = substr(month.abb[month_reordered], 1, 1),
    region = factor(region, levels = c("MAG", "PEI", "HAL", "BOF", "GOM"))
  )

# Plot
ggplot(df_plot, aes(x = month_reordered, y = value, group = metric, color = metric)) +
  geom_line(aes(alpha = if_else(metric == "Normalized Positive log concentration", 0.5, 1)),
            linewidth = 1) +
  geom_point(size = 2) +

  facet_grid(region ~ species) +

  scale_x_continuous(
    breaks = 1:12,
    labels = substr(month.abb, 1, 1)
  ) +

  scale_color_manual(
    values = c(
      "Normalized log concentration"           = "steelblue",
      "Normalized Positive log concentration" = "lightsteelblue",
      "Detection rate"                        = "darkorange"
    )
  ) +

  scale_alpha_identity() +  # uses the alpha values directly without a legend

  scale_y_continuous(limits = c(0, 1)) +

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
