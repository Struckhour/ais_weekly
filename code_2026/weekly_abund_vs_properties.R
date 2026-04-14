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
