source('./AIS_eDNA_data_prep.R')


library(dplyr)
library(lubridate)

# ----------------------------
# 1. Daily summaries from dfRawClean
# ----------------------------
qpcr_daily <- dfRawClean %>%
  group_by(region, species, date) %>%
  summarise(
    qpcr = mean(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(week = floor_date(date, "week"))

temp_daily <- dfRawClean %>%
  group_by(region, species, date) %>%
  summarise(
    temp = mean(temp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(week = floor_date(date, "week"))

# ----------------------------
# 2. Weekly summaries
# ----------------------------
qpcr_weekly <- qpcr_daily %>%
  group_by(region, species, week) %>%
  summarise(
    qpcr = mean(qpcr, na.rm = TRUE),
    .groups = "drop"
  )

temp_weekly <- temp_daily %>%
  group_by(region, species, week) %>%
  summarise(
    temp = mean(temp, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------
# 3. Join qPCR and temperature
# ----------------------------
lag_df <- full_join(
  qpcr_weekly,
  temp_weekly,
  by = c("region", "species", "week")
) %>%
  filter(
    !is.na(qpcr),
    !is.na(temp)
  ) %>%
  arrange(region, species, week)

# ----------------------------
# 4. Lag analysis
#    Positive lag = temp leads qPCR
# ----------------------------
lag_results <- lag_df %>%
  group_by(region, species) %>%
  group_modify(~{
    df <- arrange(.x, week)

    # need enough points after differencing
    if (nrow(df) < 4) {
      return(tibble(
        n = nrow(df),
        best_lag = NA_real_,
        max_cor = NA_real_
      ))
    }

    temp_diff <- diff(df$temp)
    qpcr_diff <- diff(df$qpcr)

    # remove any NA pairs just in case
    keep <- is.finite(temp_diff) & is.finite(qpcr_diff)
    temp_diff <- temp_diff[keep]
    qpcr_diff <- qpcr_diff[keep]

    if (length(temp_diff) < 3 || length(qpcr_diff) < 3) {
      return(tibble(
        n = nrow(df),
        best_lag = NA_real_,
        max_cor = NA_real_
      ))
    }

    ccf_res <- ccf(temp_diff, qpcr_diff, plot = FALSE, na.action = na.pass)

    tibble(
      n = nrow(df),
      best_lag = ccf_res$lag[which.max(ccf_res$acf)],
      max_cor = max(ccf_res$acf, na.rm = TRUE)
    )
  }) %>%
  ungroup()

lag_results


















































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

temp_weekly <- dfRawClean %>%
  group_by(region, species, date) %>%
  summarise(
    temp = mean(temp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(week = isoweek(date)) %>%
  group_by(region, species, week) %>%
  summarise(
    temp = mean(temp, na.rm = TRUE),
    .groups = "drop"
  )

# ----------------------------
# 2. Join and complete 52-week grid
# ----------------------------
season_df <- full_join(
  qpcr_weekly,
  temp_weekly,
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
  n <- length(x)
  idx <- which(is.finite(x))

  if (length(idx) == 0) return(rep(NA_real_, n))
  if (length(idx) == 1) return(rep(x[idx], n))

  # extend one full cycle on both sides
  idx_ext <- c(idx - n, idx, idx + n)
  y_ext   <- rep(x[idx], 3)

  out <- approx(
    x = idx_ext,
    y = y_ext,
    xout = seq_len(n),
    method = "linear",
    rule = 2
  )$y

  out
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
# positive lag means temp leads qPCR by that many weeks
circular_lag_profile <- function(temp, qpcr) {
  n <- length(temp)

  temp_filled <- fill_circular_series(temp)
  qpcr_filled <- fill_circular_series(qpcr)

  temp_z <- zscore(temp_filled)
  qpcr_z <- zscore(qpcr_filled)

  if (all(is.na(temp_z)) || all(is.na(qpcr_z))) {
    return(tibble(lag = integer(), cor = numeric()))
  }

  lags <- 0:(n - 1)
  cors <- sapply(lags, function(k) {
    # shift temp forward by k weeks relative to qpcr
    # so positive lag means temp occurs earlier / leads
    cor(circ_shift(temp_z, k), qpcr_z, use = "complete.obs")
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
    prof <- circular_lag_profile(.x$temp, .x$qpcr)

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
    title = "Cross Correlation Lag Between Temperature and qPCR Abundance"
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
