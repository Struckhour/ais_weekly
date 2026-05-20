source('./AIS_eDNA_data_prep.R')

run_env_qpcr_circular_lag <- function(
    df,
    env_var,
    env_name,
    qpcr_var = logConc
) {

  env_var  <- rlang::ensym(env_var)
  qpcr_var <- rlang::ensym(qpcr_var)

  qpcr_weekly <- df %>%
    dplyr::group_by(region, species, date) %>%
    dplyr::summarise(
      qpcr = mean(!!qpcr_var, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(week = lubridate::isoweek(date)) %>%
    dplyr::group_by(region, species, week) %>%
    dplyr::summarise(
      qpcr = mean(qpcr, na.rm = TRUE),
      .groups = "drop"
    )

  env_weekly <- df %>%
    dplyr::group_by(region, species, date) %>%
    dplyr::summarise(
      env = mean(!!env_var, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(week = lubridate::isoweek(date)) %>%
    dplyr::group_by(region, species, week) %>%
    dplyr::summarise(
      env = mean(env, na.rm = TRUE),
      .groups = "drop"
    )

  season_df <- dplyr::full_join(
    qpcr_weekly,
    env_weekly,
    by = c("region", "species", "week")
  ) %>%
    dplyr::group_by(region, species) %>%
    tidyr::complete(week = 1:52) %>%
    dplyr::arrange(region, species, week) %>%
    dplyr::ungroup()

  season_df %>%
    dplyr::group_by(region, species) %>%
    dplyr::group_modify(~ {
      prof <- circular_lag_profile(.x$env, .x$qpcr)

      if (nrow(prof) == 0 || all(is.na(prof$cor))) {
        return(tibble::tibble(
          best_lag = NA_real_,
          max_cor = NA_real_
        ))
      }

      best_i <- which.max(prof$cor)

      tibble::tibble(
        best_lag = prof$lag[best_i],
        max_cor = prof$cor[best_i]
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(env = env_name)
}

temp_results <- run_env_qpcr_circular_lag(
  df = dfRawClean,
  env_var = temp,
  env_name = "Temperature"
)

sal_results <- run_env_qpcr_circular_lag(
  df = dfRawClean,
  env_var = sal,
  env_name = "Salinity"
)

combined_circular_results <- dplyr::bind_rows(
  temp_results,
  sal_results
) %>%
  dplyr::filter(!is.na(best_lag), !is.na(max_cor)) %>%
  dplyr::mutate(
    region = factor(region, levels = region_order),
    species = factor(species, levels = species_order)
  ) %>%
  dplyr::arrange(species, region) %>%
  tidyr::pivot_wider(
    names_from = env,
    values_from = c(best_lag, max_cor),
    names_glue = "{env}_{.value}"
  )


combined_circular_results %>%
  dplyr::arrange(species, region) %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(
    species_display = dplyr::if_else(
      dplyr::row_number() == 1,
      as.character(species),
      ""
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::select(
    species = species_display,
    region,
    Temperature_best_lag,
    Temperature_max_cor,
    Salinity_best_lag,
    Salinity_max_cor
  ) %>%
  gt::gt() %>%
  gt::tab_header(
    title = "Cross-Correlation Lags Between Environmental Conditions and qPCR Abundance"
  ) %>%
  gt::tab_spanner(
    label = "Temperature",
    columns = c(Temperature_best_lag, Temperature_max_cor)
  ) %>%
  gt::tab_spanner(
    label = "Salinity",
    columns = c(Salinity_best_lag, Salinity_max_cor)
  ) %>%
  gt::fmt_number(
    columns = c(
      Temperature_best_lag,
      Temperature_max_cor,
      Salinity_best_lag,
      Salinity_max_cor
    ),
    decimals = 2
  ) %>%
  gt::cols_label(
    species = "Species",
    region = "Region",
    Temperature_best_lag = "Lag",
    Temperature_max_cor = "r",
    Salinity_best_lag = "Lag",
    Salinity_max_cor = "r"
  ) %>%
  gt::tab_style(
    style = gt::cell_text(color = "black"),
    locations = gt::cells_body()
  ) %>%
  gt::tab_options(
    table.background.color = "white",
    heading.background.color = "white",
    column_labels.background.color = "white",
    data_row.padding = gt::px(2),
    column_labels.padding = gt::px(2),
    heading.padding = gt::px(2)
  )


library(dplyr)
library(writexl)

combined_circular_results %>%
  arrange(species, region) %>%
  group_by(species) %>%
  mutate(
    species_display = if_else(
      row_number() == 1,
      as.character(species),
      ""
    )
  ) %>%
  ungroup() %>%
  select(
    species = species_display,
    region,
    Temperature_best_lag,
    Temperature_max_cor,
    Salinity_best_lag,
    Salinity_max_cor
  ) %>%
  rename(
    Species = species,
    Region = region,
    `Temperature lag` = Temperature_best_lag,
    `Temperature r` = Temperature_max_cor,
    `Salinity lag` = Salinity_best_lag,
    `Salinity r` = Salinity_max_cor
  ) %>%
  writexl::write_xlsx("combined_circular_results.xlsx")
