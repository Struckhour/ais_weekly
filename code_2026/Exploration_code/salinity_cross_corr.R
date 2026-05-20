# ----------------------------
# 1. Build weekly profiles
# ----------------------------
qpcr_weekly <- dfRawClean %>%
  group_by(region, species, date) %>%
  summarise(
    qpcr = mean(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(week = isoweek(date)) %>%
  group_by(region, species, week) %>%
  summarise(
    qpcr = mean(qpcr, na.rm = TRUE),
    .groups = "drop"
  )

sal_weekly <- dfRawClean %>%
  group_by(region, species, date) %>%
  summarise(
    sal = mean(sal, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(week = isoweek(date)) %>%
  group_by(region, species, week) %>%
  summarise(
    sal = mean(sal, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------
# 2. Join and complete 52-week grid
# ----------------------------
season_df <- full_join(
  qpcr_weekly,
  sal_weekly,
  by = c("region", "species", "week")
) %>%
  group_by(region, species) %>%
  complete(week = 1:52) %>%
  arrange(region, species, week) %>%
  ungroup()

# ----------------------------
# 3. Run by region × species
# ----------------------------
circular_results <- season_df %>%
  group_by(region, species) %>%
  group_modify(~{
    prof <- circular_lag_profile(.x$sal, .x$qpcr)

    if (nrow(prof) == 0 || all(is.na(prof$cor))) {
      return(tibble(
        best_lag = NA_real_,
        max_cor = NA_real_
      ))
    }

    best_i <- which.max(prof$cor)

    tibble(
      best_lag = prof$lag[best_i],
      max_cor = prof$cor[best_i]
    )
  }) %>%
  ungroup()

# ----------------------------
# 4. Table
# ----------------------------
circular_results %>%
  filter(!is.na(best_lag), !is.na(max_cor)) %>%
  mutate(
    region = factor(region, levels = region_order),
    species = factor(species, levels = species_order)
  ) %>%
  arrange(species, region) %>%
  group_by(species) %>%
  mutate(
    species_display = ifelse(row_number() == 1, as.character(species), "")
  ) %>%
  ungroup() %>%
  select(
    species = species_display,
    region,
    best_lag,
    max_cor
  ) %>%
  gt() %>%
  tab_header(
    title = "Cross Correlation Lag Between Salinity and qPCR Abundance"
  ) %>%
  fmt_number(
    columns = c(best_lag, max_cor),
    decimals = 2
  ) %>%
  cols_label(
    species = "Species",
    region = "Region",
    best_lag = "Lag (weeks)",
    max_cor = "Correlation"
  ) %>%
  tab_style(
    style = cell_text(color = "black"),
    locations = cells_body()
  ) %>%
  tab_options(
    table.background.color = "white",
    heading.background.color = "white",
    column_labels.background.color = "white",
    data_row.padding = gt::px(2),
    column_labels.padding = gt::px(2),
    heading.padding = gt::px(2)
  )
