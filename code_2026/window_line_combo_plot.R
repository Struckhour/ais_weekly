library(dplyr)
library(ggplot2)
library(patchwork)

###################################
# BUILD COMBINED DATA
###################################

window_base <- collect_window_plot_data(
  df_monthly = df,
  threshold = 0.8
) %>%
  dplyr::filter(
    !(species == "Didemnum vexillum" & region %in% c("MAG", "PEI", "HAL"))
  ) %>%
  dplyr::mutate(
    species = factor(species, levels = species_order),
    region  = factor(region, levels = rev(region_order))
  )

plot_df_combined <- window_base %>%
  left_join(
    wilcox_results_df %>%
      dplyr::select(species, region, prob_superiority),
    by = c("species", "region")
  ) %>%
  prep_window_segments() %>%
  dplyr::mutate(
    species = factor(species, levels = species_order),
    region = factor(region, levels = rev(region_order)),
    prob_label = sprintf("%d", round(prob_superiority * 100))
  ) %>%
  dplyr::arrange(species, region)

###################################
# AXIS LIMITS
###################################

x_min <- 3
x_max <- 15

label_x <- 16.2
plot_x_max <- 17.1

###################################
# ONE-TIME HEADER DATA
###################################

ps_header_df <- tibble::tibble(
  species = factor(species_order[1], levels = species_order),
  region  = factor(rev(region_order)[1], levels = rev(region_order))
)

###################################
# PLOT
###################################

ggplot(plot_df_combined) +

  geom_segment(
    aes(
      x = start_plot,
      xend = end_plot,
      y = region,
      yend = region,
      color = region
    ),
    linewidth = 4,
    lineend = "round"
  ) +

  geom_point(
    aes(
      x = peak_plot,
      y = region
    ),
    color = "black",
    size = 2.5
  ) +

  geom_text(
    aes(
      x = label_x,
      y = region,
      label = prob_label
    ),
    hjust = 0,
    size = 3.2
  ) +



  facet_wrap(
    ~ factor(species, levels = species_order),
    ncol = 1,
    scales = "free_y",
    strip.position = "top"
  ) +

  scale_color_manual(values = region_colors) +

  scale_x_continuous(
    limits = c(x_min, plot_x_max),
    breaks = seq(x_min, x_max, by = 1),
    labels = shifted_month_labels,
    expand = expansion(mult = c(0.02, 0.02))
  ) +

  coord_cartesian(clip = "off") +

  labs(
    x = "Month",
    y = NULL,
    title = "Optimal timing windows by species and region",
    subtitle = "Colored line = optimal window, black point = peak month, PS = probability of superiority (%)"
  ) +


  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    strip.text = element_text(face = "italic"),
    axis.text.y = element_text(face = "bold"),
    plot.margin = margin(10, 40, 10, 10)
  )
