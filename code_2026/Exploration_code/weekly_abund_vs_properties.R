source('./AIS_eDNA_data_prep.R')


plot_scaled_compare <- function(df, env_var, env_label, colors = NULL) {

  plot_df <- df %>%
    dplyr::select(region, species, week_of_year, scaleLogConc, {{ env_var }}) %>%
    tidyr::pivot_longer(
      cols = c(scaleLogConc, {{ env_var }}),
      names_to = "variable",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      variable = dplyr::recode(
        variable,
        scaleLogConc = "Log concentration",
        !!rlang::as_name(rlang::ensym(env_var)) := env_label
      )
    )

  p <- ggplot(plot_df, aes(x = week_of_year, y = value, color = variable)) +
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

  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }

  p
}
plot_scaled_compare(dfWeeks, scaleTemp, "Temperature", colors = c("Log concentration" = "black", "Temperature" = "#d7301f"))
plot_scaled_compare(dfWeeks, scaleSal, "Salinity", colors = c("Log concentration" = "black", "Salinity" = "blue"))
plot_scaled_compare(dfWeeks, scalePH, "pH", colors = c("Log concentration" = "black", "pH" = "#984ea3"))




plot_region_species_env <- function(
    df,
    region_name = "HAL",
    temp_var = meanTemp,
    sal_var = meanSal,
    temp_label = "Temperature",
    sal_label = "Salinity",
    species_colors = NULL
) {

  plot_df <- df %>%
    dplyr::filter(region == region_name) %>%
    dplyr::select(
      region,
      species,
      week_of_year,
      scaleLogConc,
      {{ temp_var }},
      {{ sal_var }}
    ) %>%
    dplyr::rename(
      temp_value = {{ temp_var }},
      sal_value  = {{ sal_var }}
    ) %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(
      qPCR_norm = scaleLogConc / max(scaleLogConc, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(species, week_of_year, qPCR_norm, temp_value, sal_value) %>%
    tidyr::pivot_longer(
      cols = c(qPCR_norm, temp_value, sal_value),
      names_to = "variable",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      variable = dplyr::case_when(
        variable == "qPCR_norm" ~ as.character(species),
        variable == "temp_value" ~ temp_label,
        variable == "sal_value" ~ sal_label,
        TRUE ~ variable
      )
    ) %>%
    dplyr::group_by(variable) %>%
    dplyr::mutate(
      value_norm = (value - min(value, na.rm = TRUE)) /
        (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      is_env = variable %in% c(temp_label, sal_label)
    )

  species_levels <- unique(plot_df$species)
  species_levels <- species_levels[!is.na(species_levels)]

  if (is.null(species_colors)) {
    species_colors <- setNames(
      scales::hue_pal()(length(species_levels)),
      species_levels
    )
  }

  color_values <- c(
    species_colors,
    setNames(c("grey50", "black"), c(temp_label, sal_label))
  )

  p <- ggplot(plot_df, aes(x = week_of_year, y = value_norm)) +

    # Species lines (all dotted)
    geom_line(
      data = subset(plot_df, !is_env),
      aes(color = variable),
      linewidth = 0.9,
      linetype = "solid"
    ) +

    # Species points (all same marker)
    geom_point(
      data = subset(plot_df, !is_env),
      aes(color = variable),
      size = 1.8,
      shape = 16
    ) +

    # Temp line
    geom_line(
      data = subset(plot_df, variable == temp_label),
      aes(color = variable),
      linewidth = 1.1
    ) +

    # Salinity line
    geom_line(
      data = subset(plot_df, variable == sal_label),
      aes(color = variable),
      linewidth = 1.1
    ) +

    scale_color_manual(values = color_values) +

    scale_x_continuous(breaks = seq(1, 53, by = 4)) +
    labs(
      x = "Week of year",
      y = "Normalized value",
      color = NULL
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )

  p
}
plot_region_species_env(dfWeeks)

plot_region_species_sal_loess <- function(
    df,
    region_name = "HAL",
    sal_var = meanSal,
    sal_label = "Salinity",
    species_colors = NULL,
    loess_span = 0.25
) {

  plot_df <- df %>%
    dplyr::filter(region == region_name) %>%
    dplyr::select(
      region,
      species,
      week_of_year,
      scaleLogConc,
      {{ sal_var }}
    ) %>%
    dplyr::rename(
      sal_value = {{ sal_var }}
    ) %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(
      qPCR_norm = (scaleLogConc - min(scaleLogConc, na.rm = TRUE)) /
        (max(scaleLogConc, na.rm = TRUE) - min(scaleLogConc, na.rm = TRUE))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(species, week_of_year, qPCR_norm, sal_value) %>%
    tidyr::pivot_longer(
      cols = c(qPCR_norm, sal_value),
      names_to = "variable",
      values_to = "value"
    ) %>%
    dplyr::mutate(
      variable = dplyr::case_when(
        variable == "qPCR_norm" ~ as.character(species),
        variable == "sal_value" ~ sal_label,
        TRUE ~ variable
      )
    ) %>%
    dplyr::group_by(variable) %>%
    dplyr::mutate(
      value_norm = (value - min(value, na.rm = TRUE)) /
        (max(value, na.rm = TRUE) - min(value, na.rm = TRUE))
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      is_salinity = variable == sal_label
    )

  species_levels <- unique(plot_df$species)
  species_levels <- species_levels[!is.na(species_levels)]

  if (is.null(species_colors)) {
    species_colors <- setNames(
      scales::hue_pal()(length(species_levels)),
      species_levels
    )
  }

  color_values <- c(
    species_colors,
    setNames("black", sal_label)
  )

  linetypes <- c(
    setNames(
      rep(c("solid", "dashed", "dotted", "dotdash", "longdash"),
          length.out = length(species_levels)),
      species_levels
    ),
    setNames("solid", sal_label)
  )

  ggplot(plot_df, aes(x = week_of_year, y = value_norm)) +

    geom_point(
      data = subset(plot_df, !is_salinity),
      aes(color = variable, shape = variable),
      size = 1.5,
      alpha = 0.6
    ) +

    geom_smooth(
      data = subset(plot_df, !is_salinity),
      aes(color = variable, fill = variable, linetype = variable),
      method = "loess",
      formula = y ~ x,
      span = loess_span,
      se = TRUE,
      linewidth = 1,
      alpha = 0.18
    ) +

    geom_point(
      data = subset(plot_df, is_salinity),
      aes(color = variable),
      size = 1.5,
      alpha = 0.5
    ) +

    geom_smooth(
      data = subset(plot_df, is_salinity),
      aes(color = variable, fill = variable),
      method = "loess",
      formula = y ~ x,
      span = loess_span,
      se = TRUE,
      linewidth = 1.2,
      alpha = 0.18
    ) +

    scale_color_manual(values = color_values) +
    scale_fill_manual(values = color_values) +
    scale_linetype_manual(values = linetypes) +

    scale_x_continuous(breaks = seq(1, 53, by = 4)) +
    labs(
      x = "Week of year",
      y = "Normalized value",
      color = NULL,
      fill = NULL,
      linetype = NULL,
      shape = NULL
    ) +
    theme_minimal() +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
}

plot_region_species_sal_loess(dfWeeks, loess_span = 0.15)




























df2 <- df_model %>%
  filter(region == "HAL") %>%
  arrange(species, week_of_year) %>%
  group_by(species) %>%
  mutate(
    d_qpcr = scaleLogConc - lag(scaleLogConc),
    d_sal  = meanSal - lag(meanSal)
  ) %>%
  ungroup()

threshold_q <- quantile(df2$d_qpcr, 0.1, na.rm = TRUE)
threshold_s <- quantile(df2$d_sal, 0.1, na.rm = TRUE)

print(threshold_q)
print(threshold_s)

df2 <- df2 %>%
  mutate(
    qpcr_drop = d_qpcr <= threshold_q,
    sal_drop  = d_sal  <= threshold_s
  ) %>%
  group_by(species) %>%
  arrange(week_of_year) %>%
  mutate(
    qpcr_drop_window =
      qpcr_drop |
      dplyr::lag(qpcr_drop, 1, default = FALSE) |
      dplyr::lead(qpcr_drop, 1, default = FALSE)
  ) %>%
  ungroup()

table(df2$sal_drop, df2$qpcr_drop_window)

mean(df2$qpcr_drop_window[df2$sal_drop], na.rm = TRUE)
