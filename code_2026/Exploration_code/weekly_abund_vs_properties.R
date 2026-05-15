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

































plot_halifax_raw_qpcr_env <- function(
    df,
    region_name = "HAL",
    temp_var = temp,
    sal_var = sal,
    temp_label = "Temperature",
    sal_label = "Salinity",
    species_colors = c(
      "#e41a1c",
      "#377eb8",
      "#4daf4a",
      "#984ea3",
      "#ff7f00"
    ),
    date_min = as.Date("2023-01-01"),
    date_max = as.Date("2024-08-30"),
    loess_span = 0.03,
    salinity_valley_threshold = 20
) {

  df_qpcr <- df %>%
    dplyr::filter(
      region == region_name,
      date >= date_min,
      date <= date_max
    ) %>%
    dplyr::mutate(
      species = factor(species, levels = species_order)
    )

  species_levels <- levels(droplevels(df_qpcr$species))

  species_colors <- setNames(
    species_colors[seq_along(species_levels)],
    species_levels
  )

  df_env <- df_qpcr %>%
    dplyr::select(date, {{ temp_var }}, {{ sal_var }}) %>%
    dplyr::distinct() %>%
    dplyr::rename(
      temp_value = {{ temp_var }},
      sal_value  = {{ sal_var }}
    ) %>%
    dplyr::filter(
      is.finite(temp_value),
      is.finite(sal_value)
    ) %>%
    dplyr::arrange(date)

  temp_min <- min(df_env$temp_value, na.rm = TRUE)
  temp_max <- max(df_env$temp_value, na.rm = TRUE)

  sal_min <- min(df_env$sal_value, na.rm = TRUE)
  sal_max <- max(df_env$sal_value, na.rm = TRUE)

  df_env <- df_env %>%
    dplyr::mutate(
      temp_scaled_to_sal =
        ((temp_value - temp_min) / (temp_max - temp_min)) *
        (sal_max - sal_min) + sal_min
    )

  df_vlines <- df_env %>%
    dplyr::mutate(
      below_threshold = sal_value < salinity_valley_threshold,
      valley_group = cumsum(below_threshold != dplyr::lag(below_threshold, default = FALSE))
    ) %>%
    dplyr::filter(below_threshold) %>%
    dplyr::group_by(valley_group) %>%
    dplyr::slice_min(sal_value, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(date)

  p_qpcr <- ggplot(
    df_qpcr,
    aes(
      x = date,
      y = logConc,
      color = species,
      linetype = species
    )
  ) +

    geom_vline(
      data = df_vlines,
      aes(xintercept = date),
      inherit.aes = FALSE,
      linetype = "dotted",
      color = "grey35",
      linewidth = 0.5
    ) +

    geom_smooth(
      method = "loess",
      se = FALSE,
      linewidth = 0.9,
      span = loess_span
    ) +

    scale_color_manual(values = species_colors) +

    scale_linetype_manual(
      values = c(
        "longdash",
        "solid",
        "twodash",
        "dotdash",
        "twodash"
      )[seq_along(species_levels)]
    ) +

    scale_x_date(
      date_breaks = "1 month",
      labels = function(d) substr(month.abb[lubridate::month(d)], 1, 1)
    ) +

    labs(
      title = paste("(b) Log10 qPCR Concentrations —", region_name),
      x = NULL,
      y = "log10(copies/L + 1)",
      color = NULL,
      linetype = NULL
    ) +

    guides(
      color = guide_legend(order = 1),
      linetype = guide_legend(order = 1)
    ) +

    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.key.width = unit(1.5, "lines"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )

  p_env <- ggplot(df_env, aes(x = date)) +

    geom_vline(
      data = df_vlines,
      aes(xintercept = date),
      inherit.aes = FALSE,
      linetype = "dotted",
      color = "grey35",
      linewidth = 0.5
    ) +

    geom_line(
      aes(y = sal_value, linetype = sal_label),
      color = "black",
      linewidth = 0.7
    ) +
    labs(
      title = paste("(a) Salinity and Temperature —", region_name),
      x = "Month",
      linetype = NULL
    ) +
    geom_line(
      aes(y = temp_scaled_to_sal, linetype = temp_label),
      color = "black",
      linewidth = 0.7
    ) +

    scale_linetype_manual(
      values = c(
        "Temperature" = "dashed",
        "Salinity" = "solid"
      )
    ) +

    scale_y_continuous(
      name = "Salinity (PSU)",
      sec.axis = sec_axis(
        trans = ~ ((. - sal_min) / (sal_max - sal_min)) *
          (temp_max - temp_min) + temp_min,
        name = "Temperature (°C)"
      )
    ) +

    scale_x_date(
      date_breaks = "1 month",
      labels = function(d) substr(month.abb[lubridate::month(d)], 1, 1)
    ) +

    labs(
      x = "Month",
      linetype = NULL
    ) +

    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.text = element_text(size = 10),
      legend.key.width = unit(1.5, "lines"),
      legend.key.height = unit(0.6, "lines"),   # tighter rows
      legend.spacing.y = unit(0.2, "lines"),    # less space between entries
      legend.box.spacing = unit(0.2, "lines"),  # less space above legend
      axis.text.x = element_text(size = 9),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )

  list(
    env = p_env,
    qpcr = p_qpcr
  )
}

plot_halifax_raw_qpcr_env(dfRawClean, region_name = "HAL", loess_span = 0.03)








































plot_halifax_valley_bars <- function(
    df,
    wilcox_results = NULL,
    region_name = "HAL",
    temp_var = temp,
    sal_var = sal,
    species_colors = c(
      "#e41a1c",
      "#377eb8",
      "#4daf4a",
      "#984ea3",
      "#ff7f00"
    ),
    date_min = as.Date("2023-01-01"),
    date_max = as.Date("2024-08-30"),
    salinity_valley_threshold = 20
) {

  df_hal <- df %>%
    dplyr::filter(
      region == region_name,
      date >= date_min,
      date <= date_max
    ) %>%
    dplyr::mutate(
      species = factor(species, levels = species_order)
    )

  species_levels <- levels(droplevels(df_hal$species))

  species_colors <- setNames(
    species_colors[seq_along(species_levels)],
    species_levels
  )

  df_env <- df_hal %>%
    dplyr::select(date, {{ temp_var }}, {{ sal_var }}) %>%
    dplyr::distinct() %>%
    dplyr::rename(
      temp_value = {{ temp_var }},
      sal_value  = {{ sal_var }}
    ) %>%
    dplyr::filter(is.finite(sal_value)) %>%
    dplyr::arrange(date)

  valley_dates <- df_env %>%
    dplyr::mutate(
      below_threshold = sal_value < salinity_valley_threshold,
      valley_group = cumsum(
        below_threshold != dplyr::lag(below_threshold, default = FALSE)
      )
    ) %>%
    dplyr::filter(below_threshold) %>%
    dplyr::group_by(valley_group) %>%
    dplyr::slice_min(sal_value, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(
      event = paste0("Event ", dplyr::row_number()),
      event_label = paste0(
        "Event ", dplyr::row_number(),
        "\n", format(date, "%Y-%m-%d")
      )
    ) %>%
    dplyr::select(event, event_label, event_date = date, sal_value)

  sampling_dates <- sort(unique(df_hal$date))

  event_windows <- purrr::map_dfr(seq_len(nrow(valley_dates)), function(i) {

    event_date <- valley_dates$event_date[i]
    event_index <- match(event_date, sampling_dates)

    selected_dates <- sampling_dates[
      (event_index - 2):(event_index + 2)
    ]

    tibble::tibble(
      event = valley_dates$event[i],
      event_label = valley_dates$event_label[i],
      timing = factor(
        c(
          "2 weeks before",
          "1 week before",
          "Same week",
          "1 week after",
          "2 weeks after"
        ),
        levels = c(
          "2 weeks before",
          "1 week before",
          "Same week",
          "1 week after",
          "2 weeks after"
        )
      ),
      sample_date = selected_dates
    )
  }) %>%
    dplyr::filter(!is.na(sample_date))

  plot_df <- df_hal %>%
    dplyr::inner_join(
      event_windows,
      by = c("date" = "sample_date")
    ) %>%
    dplyr::group_by(event_label, species, timing) %>%
    dplyr::summarise(
      mean_logConc = mean(logConc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    tidyr::complete(
      event_label,
      species = factor(species_levels, levels = species_levels),
      timing = factor(
        c(
          "2 weeks before",
          "1 week before",
          "Same week",
          "1 week after",
          "2 weeks after"
        ),
        levels = c(
          "2 weeks before",
          "1 week before",
          "Same week",
          "1 week after",
          "2 weeks after"
        )
      )
    )

  if (is.null(wilcox_results)) {
    wilcox_results <- run_halifax_valley_wilcox(
      df = df,
      region_name = region_name,
      temp_var = {{ temp_var }},
      sal_var = {{ sal_var }},
      date_min = date_min,
      date_max = date_max,
      salinity_valley_threshold = salinity_valley_threshold
    )
  }

  sig_df <- wilcox_results %>%
    dplyr::filter(p_value < 0.1) %>%
    dplyr::mutate(
      species = factor(species, levels = species_levels),
      sig_label = dplyr::case_when(
        p_value < 0.001 ~ "****",
        p_value < 0.01  ~ "***",
        p_value < 0.05  ~ "**",
        p_value < 0.1 ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::left_join(
      plot_df %>%
        dplyr::group_by(event_label, species) %>%
        dplyr::summarise(
          y_pos = max(mean_logConc, na.rm = TRUE) + 0.35,
          .groups = "drop"
        ),
      by = c("event_label", "species")
    )

  ggplot(
    plot_df,
    aes(
      x = species,
      y = mean_logConc,
      fill = species,
      alpha = timing
    )
  ) +
    geom_col(
      position = position_dodge(width = 0.85),
      width = 0.75,
      color = "black",
      linewidth = 0.25
    ) +

    geom_text(
      data = sig_df,
      aes(
        x = species,
        y = y_pos,
        label = sig_label
      ),
      inherit.aes = FALSE,
      size = 6,
      fontface = "bold"
    ) +

    facet_wrap(~ event_label, nrow = 1) +

    scale_alpha_manual(
      values = c(
        "2 weeks before" = 0.30,
        "1 week before" = 0.55,
        "Same week" = 1,
        "1 week after" = 0.70,
        "2 weeks after" = 0.45
      ),
      guide = guide_legend(order = 1, nrow = 1)
    ) +

    scale_fill_manual(
      values = species_colors,
      guide = guide_legend(order = 2, nrow = 1)
    ) +

    scale_y_continuous(
      expand = expansion(mult = c(0, 0.12))
    ) +

    labs(
      title = paste("(c) qPCR Concentration Around Salinity Valleys —", region_name),
      x = NULL,
      y = "Mean log10(copies/L + 1)",
      fill = NULL,
      alpha = NULL
    ) +

    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 10),

      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),

      strip.background = element_rect(color = "black", fill = "white"),
      strip.text = element_text(size = 10),

      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )
}

plot_halifax_valley_boxplots <- function(
    df,
    wilcox_results = NULL,
    region_name = "HAL",
    temp_var = temp,
    sal_var = sal,
    species_colors = c(
      "#e41a1c",
      "#377eb8",
      "#4daf4a",
      "#984ea3",
      "#ff7f00"
    ),
    date_min = as.Date("2023-01-01"),
    date_max = as.Date("2024-08-30"),
    salinity_valley_threshold = 20
) {

  timing_levels <- c(
    "2 weeks before",
    "1 week before",
    "Same week",
    "1 week after",
    "2 weeks after"
  )

  timing_alpha <- c(
    "2 weeks before" = 0.30,
    "1 week before" = 0.55,
    "Same week" = 1.00,
    "1 week after" = 0.70,
    "2 weeks after" = 0.45
  )

  df_hal <- df %>%
    dplyr::filter(
      region == region_name,
      date >= date_min,
      date <= date_max
    ) %>%
    dplyr::mutate(
      species = factor(species, levels = species_order)
    )

  species_levels <- levels(droplevels(df_hal$species))

  species_colors <- setNames(
    species_colors[seq_along(species_levels)],
    species_levels
  )

  df_env <- df_hal %>%
    dplyr::select(date, {{ temp_var }}, {{ sal_var }}) %>%
    dplyr::distinct() %>%
    dplyr::rename(
      temp_value = {{ temp_var }},
      sal_value  = {{ sal_var }}
    ) %>%
    dplyr::filter(is.finite(sal_value)) %>%
    dplyr::arrange(date)

  valley_dates <- df_env %>%
    dplyr::mutate(
      below_threshold = sal_value < salinity_valley_threshold,
      valley_group = cumsum(
        below_threshold != dplyr::lag(below_threshold, default = FALSE)
      )
    ) %>%
    dplyr::filter(below_threshold) %>%
    dplyr::group_by(valley_group) %>%
    dplyr::slice_min(sal_value, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(
      event = paste0("Event ", dplyr::row_number()),
      event_label = paste0(
        "Event ", dplyr::row_number(),
        "\n", format(date, "%Y-%m-%d")
      )
    ) %>%
    dplyr::select(event, event_label, event_date = date)

  sampling_dates <- sort(unique(df_hal$date))

  event_windows <- purrr::map_dfr(seq_len(nrow(valley_dates)), function(i) {

    event_date <- valley_dates$event_date[i]
    event_index <- match(event_date, sampling_dates)

    selected_dates <- sampling_dates[
      (event_index - 2):(event_index + 2)
    ]

    tibble::tibble(
      event = valley_dates$event[i],
      event_label = valley_dates$event_label[i],
      timing = factor(
        timing_levels,
        levels = timing_levels
      ),
      sample_date = selected_dates
    )
  }) %>%
    dplyr::filter(!is.na(sample_date))

  plot_df <- df_hal %>%
    dplyr::inner_join(
      event_windows,
      by = c("date" = "sample_date")
    ) %>%
    dplyr::mutate(
      timing = factor(timing, levels = timing_levels)
    )

  if (is.null(wilcox_results)) {
    wilcox_results <- run_halifax_valley_wilcox(
      df = df,
      region_name = region_name,
      temp_var = {{ temp_var }},
      sal_var = {{ sal_var }},
      date_min = date_min,
      date_max = date_max,
      salinity_valley_threshold = salinity_valley_threshold
    )
  }

  y_lookup <- plot_df %>%
    dplyr::group_by(event_label, species) %>%
    dplyr::summarise(
      y_pos = max(logConc, na.rm = TRUE) + 0.4,
      .groups = "drop"
    )

  sig_df <- wilcox_results %>%
    dplyr::filter(p_value < 0.05) %>%
    dplyr::mutate(
      species = factor(species, levels = species_levels),
      sig_label = dplyr::case_when(
        p_value < 0.001 ~ "***",
        p_value < 0.01  ~ "**",
        p_value < 0.05  ~ "*",
        TRUE ~ ""
      )
    ) %>%
    dplyr::left_join(
      y_lookup,
      by = c("event_label", "species")
    )

  ggplot(
    plot_df,
    aes(
      x = species,
      y = logConc,
      fill = species,
      alpha = timing
    )
  ) +
    geom_boxplot(
      position = position_dodge(width = 0.85),
      width = 0.7,
      outlier.size = 0.8,
      color = "black",
      linewidth = 0.25
    ) +

    geom_text(
      data = sig_df,
      aes(
        x = species,
        y = y_pos,
        label = sig_label
      ),
      inherit.aes = FALSE,
      size = 6,
      fontface = "bold"
    ) +

    facet_wrap(~ event_label, nrow = 1) +

    scale_fill_manual(
      values = species_colors,
      guide = guide_legend(
        order = 2,
        nrow = 1,
        override.aes = list(alpha = 1)
      )
    ) +

    scale_alpha_manual(
      values = timing_alpha,
      guide = guide_legend(
        order = 1,
        nrow = 1,
        override.aes = list(
          fill = "grey40",
          color = "black"
        )
      )
    ) +

    scale_y_continuous(
      expand = expansion(mult = c(0, 0.15))
    ) +

    labs(
      title = paste("(c) qPCR Distributions Around Salinity Valleys —", region_name),
      x = NULL,
      y = "log10(copies/L + 1)",
      fill = NULL,
      alpha = NULL
    ) +

    theme_classic() +
    theme(
      legend.position = "bottom",
      legend.box = "vertical",
      legend.text = element_text(size = 10),

      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),

      strip.background = element_rect(color = "black", fill = "white"),
      strip.text = element_text(size = 10),

      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    )
}

plot_halifax_valley_bars(dfRawClean)
plot_halifax_valley_boxplots(dfRawClean)

plots_top <- plot_halifax_raw_qpcr_env(dfRawClean)
plot_bottom <- plot_halifax_valley_bars(dfRawClean)

final_plot <-
  (plots_top$env / plots_top$qpcr) /
  plot_bottom +
  patchwork::plot_layout(
    heights = c(2, 2, 2),
    guides = "keep"
  ) &
  theme(
    legend.position = "bottom"
  )

final_plot









inspect_halifax_valley_timing <- function(
    df,
    region_name = "HAL",
    temp_var = temp,
    sal_var = sal,
    date_min = as.Date("2023-01-01"),
    date_max = as.Date("2024-08-30"),
    salinity_valley_threshold = 20
) {

  df_hal <- df %>%
    dplyr::filter(
      region == region_name,
      date >= date_min,
      date <= date_max
    ) %>%
    dplyr::mutate(
      species = factor(species, levels = species_order)
    )

  df_env <- df_hal %>%
    dplyr::select(date, {{ temp_var }}, {{ sal_var }}) %>%
    dplyr::distinct() %>%
    dplyr::rename(
      temp_value = {{ temp_var }},
      sal_value  = {{ sal_var }}
    ) %>%
    dplyr::filter(is.finite(sal_value)) %>%
    dplyr::arrange(date)

  valley_dates <- df_env %>%
    dplyr::mutate(
      below_threshold = sal_value < salinity_valley_threshold,
      valley_group = cumsum(
        below_threshold != dplyr::lag(below_threshold, default = FALSE)
      )
    ) %>%
    dplyr::filter(below_threshold) %>%
    dplyr::group_by(valley_group) %>%
    dplyr::slice_min(sal_value, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(date) %>%
    dplyr::mutate(
      event = paste0("Event ", dplyr::row_number()),
      event_label = paste0(
        "Event ", dplyr::row_number(),
        "\n", format(date, "%Y-%m-%d")
      )
    ) %>%
    dplyr::select(event, event_label, event_date = date, sal_value)

  sampling_dates <- sort(unique(df_hal$date))

  event_windows <- purrr::map_dfr(seq_len(nrow(valley_dates)), function(i) {

    event_date <- valley_dates$event_date[i]
    event_index <- match(event_date, sampling_dates)

    selected_dates <- sampling_dates[
      (event_index - 2):(event_index + 2)
    ]

    tibble::tibble(
      event = valley_dates$event[i],
      event_label = valley_dates$event_label[i],
      event_date = event_date,
      valley_salinity = valley_dates$sal_value[i],
      event_index = event_index,
      timing = factor(
        c(
          "2 weeks before",
          "1 week before",
          "Same week",
          "1 week after",
          "2 weeks after"
        ),
        levels = c(
          "2 weeks before",
          "1 week before",
          "Same week",
          "1 week after",
          "2 weeks after"
        )
      ),
      sample_date = selected_dates
    )
  }) %>%
    dplyr::filter(!is.na(sample_date))

  timing_check <- event_windows %>%
    dplyr::left_join(
      df_env,
      by = c("sample_date" = "date")
    ) %>%
    dplyr::arrange(event, sample_date)

  value_check <- df_hal %>%
    dplyr::inner_join(
      event_windows,
      by = c("date" = "sample_date")
    ) %>%
    dplyr::group_by(
      event,
      event_label,
      event_date,
      timing,
      date,
      species
    ) %>%
    dplyr::summarise(
      n_rows = dplyr::n(),
      n_detected = sum(detected, na.rm = TRUE),
      mean_concentration = mean(concentration, na.rm = TRUE),
      mean_logConc = mean(logConc, na.rm = TRUE),
      min_logConc = min(logConc, na.rm = TRUE),
      max_logConc = max(logConc, na.rm = TRUE),
      raw_logConc_values = paste(round(logConc, 3), collapse = ", "),
      .groups = "drop"
    ) %>%
    dplyr::arrange(event, species, timing)

  list(
    valley_dates = valley_dates,
    event_windows = event_windows,
    timing_check = timing_check,
    value_check = value_check
  )
}
check <- inspect_halifax_valley_timing(dfRawClean)
View(check$valley_dates)
View(check$timing_check)


run_halifax_valley_ttests <- function(
    df,
    region_name = "HAL",
    temp_var = temp,
    sal_var = sal,
    date_min = as.Date("2023-01-01"),
    date_max = as.Date("2024-08-30"),
    salinity_valley_threshold = 20
) {

  check <- inspect_halifax_valley_timing(
    df = df,
    region_name = region_name,
    temp_var = {{ temp_var }},
    sal_var = {{ sal_var }},
    date_min = date_min,
    date_max = date_max,
    salinity_valley_threshold = salinity_valley_threshold
  )

  event_windows <- check$event_windows

  df_test <- df %>%
    dplyr::filter(
      region == region_name,
      date >= date_min,
      date <= date_max
    ) %>%
    dplyr::mutate(
      species = factor(species, levels = species_order)
    ) %>%
    dplyr::inner_join(
      event_windows,
      by = c("date" = "sample_date")
    ) %>%
    dplyr::mutate(
      group = dplyr::case_when(
        timing %in% c("2 weeks before", "1 week before") ~ "before",
        timing %in% c("Same week", "1 week after", "2 weeks after") ~ "after",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(group), is.finite(logConc))

  results <- df_test %>%
    dplyr::group_by(event, event_label, species) %>%
    dplyr::group_modify(~ {

      before_vals <- .x$logConc[.x$group == "before"]
      after_vals  <- .x$logConc[.x$group == "after"]

      can_test <- length(before_vals) >= 2 &&
        length(after_vals) >= 2 &&
        stats::sd(before_vals, na.rm = TRUE) > 0 &&
        stats::sd(after_vals, na.rm = TRUE) > 0

      if (can_test) {
        tt <- stats::t.test(before_vals, after_vals)

        tibble::tibble(
          n_before = length(before_vals),
          n_after = length(after_vals),
          mean_before = mean(before_vals, na.rm = TRUE),
          mean_after = mean(after_vals, na.rm = TRUE),
          difference_after_minus_before = mean(after_vals, na.rm = TRUE) -
            mean(before_vals, na.rm = TRUE),
          statistic = unname(tt$statistic),
          p_value = tt$p.value
        )
      } else {
        tibble::tibble(
          n_before = length(before_vals),
          n_after = length(after_vals),
          mean_before = mean(before_vals, na.rm = TRUE),
          mean_after = mean(after_vals, na.rm = TRUE),
          difference_after_minus_before = mean(after_vals, na.rm = TRUE) -
            mean(before_vals, na.rm = TRUE),
          statistic = NA_real_,
          p_value = NA_real_
        )
      }
    }) %>%
    dplyr::ungroup()

  results
}
ttest_results <- run_halifax_valley_ttests(dfRawClean)

View(ttest_results)



##############################
#Wilcoxon
############
run_halifax_valley_wilcox <- function(
    df,
    region_name = "HAL",
    temp_var = temp,
    sal_var = sal,
    date_min = as.Date("2023-01-01"),
    date_max = as.Date("2024-08-30"),
    salinity_valley_threshold = 20
) {

  check <- inspect_halifax_valley_timing(
    df = df,
    region_name = region_name,
    temp_var = {{ temp_var }},
    sal_var = {{ sal_var }},
    date_min = date_min,
    date_max = date_max,
    salinity_valley_threshold = salinity_valley_threshold
  )

  event_windows <- check$event_windows

  df_test <- df %>%
    dplyr::filter(
      region == region_name,
      date >= date_min,
      date <= date_max
    ) %>%
    dplyr::mutate(
      species = factor(species, levels = species_order)
    ) %>%
    dplyr::inner_join(
      event_windows,
      by = c("date" = "sample_date")
    ) %>%
    dplyr::mutate(
      group = dplyr::case_when(
        timing == "1 week before" ~ "before",
        timing == "1 week after"  ~ "after",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(group), is.finite(logConc))

  results <- df_test %>%
    dplyr::group_by(event, event_label, species) %>%
    dplyr::group_modify(~ {

      before_vals <- .x$logConc[.x$group == "before"]
      after_vals  <- .x$logConc[.x$group == "after"]

      can_test <- length(before_vals) >= 1 &&
        length(after_vals) >= 1

      if (can_test) {
        wt <- stats::wilcox.test(
          before_vals,
          after_vals,
          exact = FALSE
        )

        tibble::tibble(
          n_before = length(before_vals),
          n_after = length(after_vals),
          median_before = median(before_vals, na.rm = TRUE),
          median_after = median(after_vals, na.rm = TRUE),
          mean_before = mean(before_vals, na.rm = TRUE),
          mean_after = mean(after_vals, na.rm = TRUE),
          difference_median_after_minus_before =
            median(after_vals, na.rm = TRUE) - median(before_vals, na.rm = TRUE),
          statistic = unname(wt$statistic),
          p_value = wt$p.value
        )
      } else {
        tibble::tibble(
          n_before = length(before_vals),
          n_after = length(after_vals),
          median_before = median(before_vals, na.rm = TRUE),
          median_after = median(after_vals, na.rm = TRUE),
          mean_before = mean(before_vals, na.rm = TRUE),
          mean_after = mean(after_vals, na.rm = TRUE),
          difference_median_after_minus_before =
            median(after_vals, na.rm = TRUE) - median(before_vals, na.rm = TRUE),
          statistic = NA_real_,
          p_value = NA_real_
        )
      }
    }) %>%
    dplyr::ungroup()

  overall_paired <- df_test %>%
    dplyr::group_by(event, event_label, species, group) %>%
    dplyr::summarise(
      value = mean(logConc, na.rm = TRUE),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = group,
      values_from = c(value, n)
    ) %>%
    dplyr::filter(
      is.finite(value_before),
      is.finite(value_after)
    )

  if (nrow(overall_paired) >= 2) {
    overall_test <- stats::wilcox.test(
      overall_paired$value_after,
      overall_paired$value_before,
      paired = TRUE,
      exact = FALSE
    )

    overall_results <- tibble::tibble(
      comparison = "overall paired: 1 week after vs 1 week before (means)",
      n_pairs = nrow(overall_paired),
      mean_before = mean(overall_paired$value_before, na.rm = TRUE),
      mean_after = mean(overall_paired$value_after, na.rm = TRUE),
      mean_difference_after_minus_before =
        mean(overall_paired$value_after - overall_paired$value_before, na.rm = TRUE),
      statistic = unname(overall_test$statistic),
      p_value = overall_test$p.value
    )
  } else {
    overall_results <- tibble::tibble(
      comparison = "overall paired: 1 week after vs 1 week before (means)",
      n_pairs = nrow(overall_paired),
      mean_before = NA_real_,
      mean_after = NA_real_,
      mean_difference_after_minus_before = NA_real_,
      statistic = NA_real_,
      p_value = NA_real_
    )
  }
  list(
    by_event_species = results,
    overall_paired = overall_results,
    overall_paired_data = overall_paired
  )
}
wilcox_results <- run_halifax_valley_wilcox(dfRawClean)

View(wilcox_results$overall_paired)




plot_valley_distributions <- function(
    df,
    event_to_plot = "Event 1"
) {

  check <- inspect_halifax_valley_timing(df)
  event_windows <- check$event_windows

  df_plot <- df %>%
    dplyr::filter(region == "HAL") %>%
    dplyr::inner_join(
      event_windows,
      by = c("date" = "sample_date")
    ) %>%
    dplyr::filter(event == event_to_plot) %>%
    dplyr::mutate(
      species = factor(species, levels = species_order),
      group = dplyr::case_when(
        timing %in% c("2 weeks before", "1 week before") ~ "Before",
        timing %in% c("Same week", "1 week after", "2 weeks after") ~ "After"
      )
    )

  ggplot(df_plot, aes(x = logConc, fill = group)) +

    geom_histogram(
      bins = 20,
      alpha = 0.6,
      position = "identity"
    ) +

    facet_wrap(~ species, scales = "free_y") +

    labs(
      title = paste("Distribution of log10 eDNA concentrations —", event_to_plot),
      x = "log10(copies/L + 1)",
      y = "Count",
      fill = NULL
    ) +

    theme_classic() +
    theme(
      legend.position = "top",
      strip.background = element_rect(color = "black", fill = "white"),
      strip.text = element_text(size = 10)
    )
}
plot_valley_distributions(dfRawClean, "Event 1")

plot_valley_distributions(dfRawClean, "Event 2")
plot_valley_distributions(dfRawClean, "Event 3")


run_halifax_valley_zero_prop <- function(
    df,
    region_name = "HAL",
    temp_var = temp,
    sal_var = sal,
    date_min = as.Date("2023-01-01"),
    date_max = as.Date("2024-08-30"),
    salinity_valley_threshold = 20
) {

  check <- inspect_halifax_valley_timing(
    df = df,
    region_name = region_name,
    temp_var = {{ temp_var }},
    sal_var = {{ sal_var }},
    date_min = date_min,
    date_max = date_max,
    salinity_valley_threshold = salinity_valley_threshold
  )

  event_windows <- check$event_windows

  df_test <- df %>%
    dplyr::filter(
      region == region_name,
      date >= date_min,
      date <= date_max
    ) %>%
    dplyr::mutate(
      species = factor(species, levels = species_order)
    ) %>%
    dplyr::inner_join(
      event_windows,
      by = c("date" = "sample_date")
    ) %>%
    dplyr::mutate(
      group = dplyr::case_when(
        timing == "1 week before" ~ "before",
        timing == "1 week after"  ~ "after",
        TRUE ~ NA_character_
      )
    ) %>%
    dplyr::filter(!is.na(group), is.finite(logConc))

  results <- df_test %>%
    dplyr::group_by(event, event_label, species, group) %>%
    dplyr::summarise(
      n = dplyr::n(),
      n_zero = sum(logConc == 0, na.rm = TRUE),
      prop_zero = n_zero / n,
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = group,
      values_from = c(n, n_zero, prop_zero)
    ) %>%
    dplyr::mutate(
      diff_prop_zero_after_minus_before =
        prop_zero_after - prop_zero_before
    )

  results
}
res <- run_halifax_valley_zero_prop(dfRawClean)
View(res)
