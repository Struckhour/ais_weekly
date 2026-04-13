source('./AIS_eDNA_data_prep.R')

df <- dfWeeks %>%
  dplyr::select(region, species, week_of_year, scaleLogConc)


prep_weekly_signal <- function(df) {
  df %>%
    dplyr::group_by(week_of_year) %>%
    dplyr::summarise(
      value = mean(scaleLogConc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(week = week_of_year) %>%
    dplyr::arrange(week)
}


interp_weekly_circular <- function(df, n_weeks = 52) {

  full <- data.frame(week = 1:n_weeks)

  x <- full %>%
    dplyr::left_join(df, by = "week") %>%
    dplyr::arrange(week)

  x$is_observed <- !is.na(x$value)

  if (all(is.na(x$value))) {
    x$value <- NA
    x$is_imputed <- FALSE
    return(x)
  }

  if (sum(x$is_observed) == 1) {
    val <- x$value[x$is_observed]
    x$value <- rep(val, n_weeks)
    x$is_imputed <- !x$is_observed
    return(x)
  }

  idx <- which(x$is_observed)
  y <- x$value[idx]

  idx_ext <- c(idx - n_weeks, idx, idx + n_weeks)
  y_ext   <- c(y, y, y)

  interp <- approx(
    x = idx_ext,
    y = y_ext,
    xout = 1:n_weeks,
    method = "linear",
    rule = 2
  )$y

  x$value <- interp
  x$is_imputed <- !x$is_observed

  x
}


classify_weeks <- function(df, threshold = 0.75) {
  df %>%
    dplyr::mutate(
      status = dplyr::case_when(
        value >= threshold ~ "above",
        TRUE ~ "below"
      )
    )
}



plot_weekly_wheel <- function(df, title = NULL) {

  week_breaks <- seq(1, 52, by = 4)

  ggplot2::ggplot(df) +

    # radial lines
    ggplot2::geom_vline(
      xintercept = 1:52,
      color = "black",
      linewidth = 0.2,
      alpha = 0.25
    ) +

    # circular grid lines
    ggplot2::geom_hline(
      yintercept = c(0.25, 0.5, 0.75, 1),
      color = "black",
      linewidth = 0.3,
      alpha = 0.5
    ) +

    # base bars
    ggplot2::geom_col(
      ggplot2::aes(
        x = week,
        y = value,
        fill = status
      ),
      width = 0.95
    ) +

    ggplot2::scale_fill_manual(
      values = c(
        above = "#5B2C83",
        below = "grey70"
      ),
      labels = c(
        above = "Above threshold",
        below = "Below threshold"
      )
    ) +

    # striped overlay for imputed weeks
    ggpattern::geom_col_pattern(
      data = dplyr::filter(df, is_imputed == TRUE),
      ggplot2::aes(x = week, y = value),
      width = 0.95,
      fill = NA,
      pattern = "stripe",
      pattern_fill = "white",
      pattern_color = "white",
      pattern_density = 0.1,
      pattern_spacing = 0.015
    ) +

    ggplot2::coord_polar(start = 0) +

    ggplot2::scale_x_continuous(
      limits = c(0.5, 52.5),
      breaks = week_breaks,
      labels = week_breaks
    ) +

    ggplot2::ylim(0, 1) +

    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 9),
      legend.title = ggplot2::element_blank(),
      legend.position = "none"
    ) +

    ggplot2::labs(title = title)
}


# all_species <- unique(df$species)
# all_regions <- unique(df$region)
# sp <- all_species[1]
# reg <- all_regions[1]
# df_sub <- df %>%
#   dplyr::filter(
#     species == sp,
#     region == reg
#   )
#
# plot_df <- df_sub %>%
#   prep_weekly_signal() %>%
#   interp_weekly_circular(n_weeks = 52) %>%
#   classify_weeks(threshold = 0.5)
#
# plot_weekly_wheel(plot_df)



df <- dfWeeks %>%
  dplyr::select(region, species, week_of_year, scaleLogConc)

if (!dir.exists("abundance_wheels_weekly")) {
  dir.create("abundance_wheels_weekly")
}

all_species <- unique(df$species)
all_regions <- unique(df$region)

for (sp in all_species) {
  for (reg in all_regions) {

    df_sub <- df %>%
      dplyr::filter(
        species == sp,
        region == reg
      )

    if (nrow(df_sub) == 0) next

    plot_df <- df_sub %>%
      prep_weekly_signal() %>%
      interp_weekly_circular(n_weeks = 52) %>%
      classify_weeks(threshold = 0.5)

    p <- plot_weekly_wheel(plot_df)

    file_name <- paste0(
      "abundance_wheels_weekly/",
      gsub(" ", "_", sp),
      "__",
      gsub(" ", "_", reg),
      ".png"
    )

    ggplot2::ggsave(
      filename = file_name,
      plot = p,
      width = 6,
      height = 6,
      dpi = 300
    )
  }
}



























plot_df <- dfWeeks %>%
  select(region, species, week_of_year, scaleLogConc, scaleTemp) %>%
  pivot_longer(
    cols = c(scaleLogConc, scaleTemp),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = recode(
      variable,
      scaleLogConc = "Log concentration",
      scaleTemp    = "Temperature",
      scaleSal     = "Salinity"
    )
  )

ggplot(plot_df, aes(x = week_of_year, y = value, color = variable)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.2) +
  facet_grid(region ~ species) +
  scale_x_continuous(breaks = seq(1, 53, by = 4)) +
  labs(
    x = "Week of year",
    y = "Scaled value",
    color = NULL
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )






plot_df <- dfWeeks %>%
  select(region, species, week_of_year, scaleLogConc, scaleSal) %>%
  pivot_longer(
    cols = c(scaleLogConc, scaleSal),
    names_to = "variable",
    values_to = "value"
  ) %>%
  mutate(
    variable = recode(
      variable,
      scaleLogConc = "Log concentration",
      scaleTemp    = "Temperature",
      scaleSal     = "Salinity"
    )
  )

ggplot(plot_df, aes(x = week_of_year, y = value, color = variable)) +
  geom_line(linewidth = 0.9) +
  geom_point(size = 1.2) +
  facet_grid(region ~ species) +
  scale_x_continuous(breaks = seq(1, 53, by = 4)) +
  labs(
    x = "Week of year",
    y = "Scaled value",
    color = NULL
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
