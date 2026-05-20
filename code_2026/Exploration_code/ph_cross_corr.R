

library(dplyr)
library(lubridate)
library(tidyr)

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

pH_weekly <- dfRawClean %>%
  group_by(region, species, date) %>%
  summarise(
    pH = mean(pH, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(week = isoweek(date)) %>%
  group_by(region, species, week) %>%
  summarise(
    pH = mean(pH, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------
# 2. Join and complete 52-week grid
# ----------------------------
season_df <- full_join(
  qpcr_weekly,
  pH_weekly,
  by = c("region", "species", "week")
) %>%
  group_by(region, species) %>%
  complete(week = 1:52) %>%
  arrange(region, species, week) %>%
  ungroup()

# ----------------------------
# 3. Helper functions
# ----------------------------

# linear interpolation, leaving ends extended
fill_circular_series <- function(x) {
  idx <- which(is.finite(x))

  if (length(idx) == 0) return(rep(NA_real_, length(x)))
  if (length(idx) == 1) return(rep(x[idx], length(x)))

  approx(
    x = idx,
    y = x[idx],
    xout = seq_along(x),
    method = "linear",
    rule = 2
  )$y
}

# z-score helper
zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(s) || s == 0) return(rep(NA_real_, length(x)))
  (x - m) / s
}

# circular shift to the right by k
circ_shift <- function(x, k) {
  n <- length(x)
  k <- k %% n
  if (k == 0) return(x)
  c(tail(x, k), head(x, n - k))
}

# circular correlation profile
# positive lag means pH leads qPCR by that many weeks
circular_lag_profile <- function(pH, qpcr) {
  n <- length(pH)

  pH_filled <- fill_circular_series(pH)
  qpcr_filled <- fill_circular_series(qpcr)

  pH_z <- zscore(pH_filled)
  qpcr_z <- zscore(qpcr_filled)

  if (all(is.na(pH_z)) || all(is.na(qpcr_z))) {
    return(tibble(lag = integer(), cor = numeric()))
  }

  lags <- 0:(n - 1)
  cors <- sapply(lags, function(k) {
    # shift pH forward by k weeks relative to qpcr
    # so positive lag means pH occurs earlier / leads
    cor(circ_shift(pH_z, k), qpcr_z, use = "complete.obs")
  })

  tibble(
    lag_raw = lags,
    cor = cors
  ) %>%
    mutate(
      lag = ifelse(lag_raw <= n / 2, lag_raw, lag_raw - n)
    ) %>%
    select(lag, cor) %>%
    arrange(lag)
}

# ----------------------------
# 4. Run by region × species
# ----------------------------
circular_results <- season_df %>%
  group_by(region, species) %>%
  group_modify(~{
    prof <- circular_lag_profile(.x$pH, .x$qpcr)

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

circular_results



bound_circular_lag_profile <- function(pH, qpcr, max_lag = 12) {
  n <- length(pH)

  pH_filled <- fill_circular_series(pH)
  qpcr_filled <- fill_circular_series(qpcr)

  pH_z <- zscore(pH_filled)
  qpcr_z <- zscore(qpcr_filled)

  # bail out early if either series is unusable
  if (all(is.na(pH_z)) || all(is.na(qpcr_z))) {
    return(tibble(
      lag = (-max_lag):(max_lag),
      cor = NA_real_
    ))
  }

  lags <- (-max_lag):(max_lag)

  cors <- sapply(lags, function(k) {
    shifted_pH <- circ_shift(pH_z, k)

    keep <- is.finite(shifted_pH) & is.finite(qpcr_z)

    if (sum(keep) < 3) {
      return(NA_real_)
    }

    cor(shifted_pH[keep], qpcr_z[keep])
  })

  tibble(
    lag = lags,
    cor = cors
  )
}


bound_circular_results <- season_df %>%
  group_by(region, species) %>%
  group_modify(~{
    prof <- bound_circular_lag_profile(.x$pH, .x$qpcr)

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



library(dplyr)
library(gt)

# your preferred orders (adjust if needed)
region_order <- c("MAG", "PEI", "HAL", "BOF")
species_order <- c(
  "Membranipora membranacea",
  "Botrylloides violaceus",
  "Didemnum vexillum",
  "Ciona intestinalis",
  "Carcinus maenas"
)




circular_results %>%
  # remove NA rows
  filter(!is.na(best_lag), !is.na(max_cor)) %>%

  # enforce ordering
  mutate(
    region = factor(region, levels = region_order),
    species = factor(species, levels = species_order)
  ) %>%
  arrange(species, region) %>%

  # show species once per block
  group_by(species) %>%
  mutate(
    species_display = ifelse(row_number() == 1, as.character(species), "")
  ) %>%
  ungroup() %>%

  # select + rename columns
  select(
    species = species_display,
    region,
    best_lag,
    max_cor
  ) %>%

  gt() %>%

  tab_header(
    title = "Cross Correlation Lag Between pH and qPCR Abundance"
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

  # black text everywhere
  tab_style(
    style = cell_text(color = "black"),
    locations = cells_body()
  ) %>%

  # clean white styling
  tab_options(
    table.background.color = "white",
    heading.background.color = "white",
    column_labels.background.color = "white"
  )
