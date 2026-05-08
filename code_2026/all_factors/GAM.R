library(dplyr)
library(tidyr)
library(stringr)
library(mgcv)
library(ggplot2)

# =========================================================
# GAM starter for weekly qPCR + environment + plate data
# =========================================================
# Assumes these objects already exist in your workspace:
#   - dfWeeks
#   - plateAbundance
#   - state_levels
#
# This script:
#   1) rebuilds the joined weekly modelling dataframe
#   2) prepares GAM-ready species-specific data
#   3) fits a small set of starting GAMs
#   4) provides quick diagnostics and plotting helpers
#
# Recommended first response variable:
#   scaleLogConc
# because that is what you used in the random forest workflow.
#
# Later, if useful, you can split this into:
#   - detection GAM (binomial)
#   - positive-only abundance GAM (Gaussian/Gamma)
# =========================================================


# --------------------------
# 1. Prep plate abundance weekly
# --------------------------
plate_abund_weekly <- plateAbundance %>%
  mutate(
    region = as.character(region),
    species = as.character(species),
    sampWeek = as.numeric(sampWeek),
    Avg = as.numeric(Avg)
  ) %>%
  filter(
    !is.na(species),
    !is.na(region),
    !is.na(sampWeek)
  ) %>%
  group_by(species, region, sampWeek) %>%
  summarise(
    mean_plate_abund = mean(Avg, na.rm = TRUE),
    n_plate_abund_obs = sum(!is.na(Avg)),
    .groups = "drop"
  )


# --------------------------
# 2. Prep weekly stage indicators
# --------------------------
plate_stage_weekly <- plateAbundance %>%
  mutate(
    region = as.character(region),
    species = as.character(species),
    sampWeek = as.numeric(sampWeek),
    State = as.character(State)
  ) %>%
  filter(
    !is.na(species),
    !is.na(region),
    !is.na(sampWeek),
    !is.na(State),
    State != ""
  ) %>%
  separate_rows(State, sep = ",\\s*") %>%
  mutate(
    State = str_trim(State),
    State = if_else(State %in% state_levels, State, NA_character_)
  ) %>%
  filter(!is.na(State)) %>%
  distinct(species, region, sampWeek, State) %>%
  mutate(value = 1L) %>%
  pivot_wider(
    names_from = State,
    values_from = value,
    values_fill = 0
  )


# --------------------------
# 3. Combine weekly plate abundance + stage indicators
# --------------------------
plate_weekly <- plate_abund_weekly %>%
  full_join(
    plate_stage_weekly,
    by = c("species", "region", "sampWeek")
  )

for (st in state_levels) {
  if (!st %in% names(plate_weekly)) {
    plate_weekly[[st]] <- 0L
  }
}

plate_weekly <- plate_weekly %>%
  mutate(
    across(all_of(state_levels), ~ replace_na(.x, 0L))
  ) %>%
  select(
    species,
    region,
    sampWeek,
    mean_plate_abund,
    n_plate_abund_obs,
    all_of(state_levels)
  )


# --------------------------
# 4. Prep qPCR weekly keys for join
# --------------------------
qpcr_weekly <- dfWeeks %>%
  mutate(
    region = as.character(region),
    species = as.character(species),
    week_of_year = as.numeric(week_of_year)
  )


# --------------------------
# 5. Join everything together
# --------------------------
df_model <- qpcr_weekly %>%
  left_join(
    plate_weekly %>%
      rename(week_of_year = sampWeek),
    by = c("species", "region", "week_of_year")
  ) %>%
  mutate(
    across(all_of(state_levels), ~ replace_na(.x, 0L))
  )


# =========================================================
# Helpers
# =========================================================

prep_species_gam_data <- function(df, species_name, response = "scaleLogConc") {

  keep_vars <- c(
    response,
    "species",
    "region",
    "week_of_year",
    "meanTemp",
    "meanSal",
    "meanPH",
    "meanLat",
    "mean_plate_abund",
    "n_plate_abund_obs",
    state_levels
  )

  keep_vars <- keep_vars[keep_vars %in% names(df)]

  out <- df %>%
    filter(species == species_name) %>%
    select(all_of(keep_vars)) %>%
    mutate(
      region = factor(region),
      week_of_year = as.numeric(week_of_year)
    )

  out <- out %>%
    filter(
      is.finite(.data[[response]]),
      is.finite(week_of_year),
      is.finite(meanTemp),
      is.finite(meanSal),
      is.finite(meanLat)
    )

  out <- out %>%
    group_by(region) %>%
    filter(sum(.data[[response]] > 0, na.rm = TRUE) >= 5) %>%
    ungroup()

  out
}


fit_gam_set <- function(df_sp, response = "scaleLogConc") {

  week_knots <- list(week_of_year = c(0.5, 52.5))

  n_week_unique <- dplyr::n_distinct(df_sp$week_of_year)
  n_region <- dplyr::n_distinct(df_sp$region)

  if (nrow(df_sp) < 20) {
    stop("Too few rows for GAM.")
  }

  if (n_week_unique < 8) {
    stop("Too few unique weeks for a seasonal GAM.")
  }

  if (n_region < 1) {
    stop("No regions available.")
  }

  k_week <- min(8, n_week_unique - 1)
  k_temp <- min(6, sum(is.finite(df_sp$meanTemp)) - 1)
  k_sal  <- min(6, sum(is.finite(df_sp$meanSal)) - 1)

  m_season <- gam(
    formula = as.formula(
      paste0(response, " ~ s(week_of_year, bs = 'cc', k = ", k_week, ") + s(region, bs = 're')")
    ),
    data = df_sp,
    method = "REML",
    knots = week_knots
  )

  m_env <- gam(
    formula = as.formula(
      paste0(
        response,
        " ~ s(week_of_year, bs = 'cc', k = ", k_week, ")",
        " + s(meanTemp, k = ", k_temp, ")",
        " + s(meanSal, k = ", k_sal, ")",
        " + s(region, bs = 're')"
      )
    ),
    data = df_sp,
    method = "REML",
    knots = week_knots
  )

  df_plate <- df_sp %>% filter(is.finite(mean_plate_abund))

  m_env_plate <- NULL
  if (nrow(df_plate) >= 15 && dplyr::n_distinct(df_plate$mean_plate_abund) >= 5) {
    k_plate <- min(5, dplyr::n_distinct(df_plate$mean_plate_abund) - 1)

    m_env_plate <- gam(
      formula = as.formula(
        paste0(
          response,
          " ~ s(week_of_year, bs = 'cc', k = ", k_week, ")",
          " + s(meanTemp, k = ", k_temp, ")",
          " + s(meanSal, k = ", k_sal, ")",
          " + s(mean_plate_abund, k = ", k_plate, ")",
          " + s(region, bs = 're')"
        )
      ),
      data = df_plate,
      method = "REML",
      knots = week_knots
    )
  }

  list(
    season = m_season,
    env = m_env,
    env_plate = m_env_plate
  )
}


extract_gam_summary <- function(fit, species_name, model_name) {
  s <- summary(fit)

  tibble(
    species = species_name,
    model = model_name,
    n = nrow(fit$model),
    adj_r2 = unname(s$r.sq),
    dev_expl = unname(s$dev.expl),
    scale_est = unname(s$scale),
    aic = AIC(fit)
  )
}


plot_gam_smooths <- function(fit, pages = 1) {
  par(mfrow = c(2, 2))
  plot(fit, pages = pages, shade = TRUE, scale = 0)
}


make_prediction_grid <- function(df_sp) {
  df_sp %>%
    group_by(region, week_of_year) %>%
    summarise(
      meanTemp = mean(meanTemp, na.rm = TRUE),
      meanSal = mean(meanSal, na.rm = TRUE),
      mean_plate_abund = mean(mean_plate_abund, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    complete(region, week_of_year = 1:52) %>%
    group_by(region) %>%
    arrange(week_of_year) %>%
    mutate(
      meanTemp = approx(week_of_year[!is.na(meanTemp)], meanTemp[!is.na(meanTemp)], xout = week_of_year, rule = 2)$y,
      meanSal = approx(week_of_year[!is.na(meanSal)], meanSal[!is.na(meanSal)], xout = week_of_year, rule = 2)$y,
      mean_plate_abund = if (all(is.na(mean_plate_abund))) NA_real_ else approx(
        week_of_year[!is.na(mean_plate_abund)],
        mean_plate_abund[!is.na(mean_plate_abund)],
        xout = week_of_year,
        rule = 2
      )$y
    ) %>%
    ungroup() %>%
    mutate(region = factor(region, levels = levels(df_sp$region)))
}


plot_predicted_seasonal_curve <- function(fit, df_sp, title_text = NULL) {
  newdat <- make_prediction_grid(df_sp)
  newdat$pred <- predict(fit, newdata = newdat, type = "response")

  ggplot(newdat, aes(x = week_of_year, y = pred, color = region)) +
    geom_line(linewidth = 1.1) +
    scale_x_continuous(breaks = seq(1, 52, by = 4)) +
    labs(
      title = title_text,
      x = "Week of year",
      y = "Predicted response",
      color = "Region"
    ) +
    theme_classic()
}


# =========================================================
# Example: fit one species
# =========================================================

sp <- "Didemnum vexillum"

# other examples:
# sp <- "Membranipora membranacea"
sp <- "Botrylloides violaceus"
# sp <- "Ciona intestinalis"
# sp <- "Carcinus maenas"


df_sp_gam <- prep_species_gam_data(
  df = df_model,
  species_name = sp,
  response = "scaleLogConc"
)

models <- fit_gam_set(df_sp_gam, response = "scaleLogConc")

summary_table_one_species <- bind_rows(
  extract_gam_summary(models$season, sp, "season"),
  extract_gam_summary(models$env, sp, "season + temp + sal"),
  extract_gam_summary(models$env_plate, sp, "season + temp + sal + plate")
)

summary_table_one_species

gam.check(models$season)
gam.check(models$env)
gam.check(models$env_plate)

# plot_gam_smooths(models$env)
# plot_predicted_seasonal_curve(models$env, df_sp_gam, title_text = sp)


# =========================================================
# Loop through all species
# =========================================================

all_species <- unique(df_model$species)
all_gam_results <- list()
all_gam_models <- list()

for (sp in all_species) {

  df_sp_gam <- prep_species_gam_data(
    df = df_model,
    species_name = sp,
    response = "scaleLogConc"
  )

  if (nrow(df_sp_gam) < 20) next
  if (nlevels(droplevels(df_sp_gam$region)) < 2) next

  fitted_models <- fit_gam_set(df_sp_gam, response = "scaleLogConc")

  all_gam_models[[sp]] <- fitted_models

  all_gam_results[[sp]] <- bind_rows(
    extract_gam_summary(fitted_models$season, sp, "season"),
    extract_gam_summary(fitted_models$env, sp, "season + temp + sal"),
    extract_gam_summary(fitted_models$env_plate, sp, "season + temp + sal + plate")
  )
}

gam_summary_table <- bind_rows(all_gam_results) %>%
  mutate(
    across(where(is.numeric), ~ round(.x, 3))
  )

gam_summary_table


# =========================================================
# Next likely upgrades
# =========================================================
# 1) Add stage variables one at a time, not all at once
#    Example: + R + G + B or as factors if that makes more biological sense
#
# 2) Add tensor terms if needed:
#    te(meanTemp, meanSal)
#
# 3) Add region-specific seasonal shapes:
#    s(week_of_year, region, bs = 'fs', k = 8)
#    or separate smooths by region
#
# 4) Build a delta GAM:
#    - detected ~ ... family = binomial
#    - positive-only scale/log concentration ~ ...
#
# 5) Try lagged covariates once the basic GAMs are stable
# =========================================================


















df_sp <- df_model %>%
  dplyr::filter(species == "Membranipora membranacea")


df_sp <- df_model %>%
  dplyr::filter(species == "Botrylloides violaceus")

unique(df_sp$species)
#starting over

nrow(df_sp)

df_sp %>%
  dplyr::count(region)

df_sp <- df_sp %>%
  dplyr::filter(
    is.finite(scaleLogConc),
    is.finite(week_of_year)
  )


fit <- mgcv::gam(
  scaleLogConc ~ s(week_of_year, bs = "cc"),
  data = df_sp,
  method = "REML",
  knots = list(week_of_year = c(0.5, 52.5))
)

plot(fit)

View(df_sp_gam)

df_sp$pred <- predict(fit)

ggplot(df_sp, aes(x = week_of_year, y = scaleLogConc)) +
  geom_point(alpha = 0.3) +
  geom_line(aes(y = pred), color = "red", linewidth = 1.2) +
  theme_classic()

ggplot(df_sp, aes(x = week_of_year, y = scaleLogConc)) +
  geom_point(alpha = 0.3) +
  geom_line(aes(y = pred), color = "red", linewidth = 1.2) +
  facet_wrap(~ region) +
  theme_classic()

####ADD TEMP

summary(df_sp$meanTemp)
df_sp2 <- df_sp %>%
  dplyr::filter(is.finite(meanTemp))

fit_temp <- mgcv::gam(
  scaleLogConc ~
    s(week_of_year, bs = "cc") +
    s(meanTemp),
  data = df_sp2,
  method = "REML",
  knots = list(week_of_year = c(0.5, 52.5))
)
plot(fit_temp)

####Lag
# df_sp_lag <- df_sp %>%
#   arrange(region, week_of_year) %>%
#   group_by(region) %>%
#   mutate(
#     temp_lag4 = dplyr::lag(meanTemp, 4)
#   ) %>%
#   ungroup()
#
# df_sp_lag <- df_sp_lag %>%
#   dplyr::filter(is.finite(temp_lag4))
#
# fit_temp_lag4 <- mgcv::gam(
#   scaleLogConc ~
#     s(week_of_year, bs = "cc") +
#     s(temp_lag4),
#   data = df_sp_lag,
#   method = "REML",
#   knots = list(week_of_year = c(0.5, 52.5))
# )
#
# summary(fit_temp)
# summary(fit_temp_lag4)
#
# plot(fit_temp_lag4)
#
# df_sp_lag$pred <- predict(fit_temp_lag4)
#
# ggplot(df_sp_lag, aes(x = week_of_year, y = scaleLogConc)) +
#   geom_point(alpha = 0.3) +
#   geom_line(aes(y = pred), color = "red", linewidth = 1.2) +
#   facet_wrap(~ region) +
#   theme_classic()

############################
#adding sal
############################

summary(df_sp$meanSal)
df_sp3 <- df_sp %>%
  dplyr::filter(
    is.finite(meanTemp),
    is.finite(meanSal)
  )

fit_temp_sal <- mgcv::gam(
  scaleLogConc ~
    s(week_of_year, bs = "cc") +
    s(meanTemp) +
    s(meanSal),
  data = df_sp3,
  method = "REML",
  knots = list(week_of_year = c(0.5, 52.5))
)

summary(fit_temp)
summary(fit_temp_sal)
plot(fit_temp_sal)


##########
#adding plate
##########
df_sp_plate <- df_sp %>%
  dplyr::filter(
    is.finite(meanTemp),
    is.finite(mean_plate_abund)
  )
nrow(df_sp_plate)

df_sp_plate %>%
  dplyr::summarise(
    n = n(),
    n_plate = sum(is.finite(mean_plate_abund))
  )
fit_temp_plate <- mgcv::gam(
  scaleLogConc ~
    s(week_of_year, bs = "cc") +
    s(meanTemp) +
    s(mean_plate_abund),
  data = df_sp_plate,
  method = "REML",
  knots = list(week_of_year = c(0.5, 52.5))
)
fit_temp_only <- mgcv::gam(
  scaleLogConc ~
    s(week_of_year, bs = "cc") +
    s(meanTemp),
  data = df_sp_plate,
  method = "REML",
  knots = list(week_of_year = c(0.5, 52.5))
)
summary(fit_temp_only)
summary(fit_temp_plate)










































########
#big gam analysis for temp and salinity
########
extract_gam_summary <- function(fit, species_name, model_name) {
  s <- summary(fit)

  tibble::tibble(
    species = species_name,
    model = model_name,
    n = nrow(fit$model),
    dev_expl = s$dev.expl,
    adj_r2 = s$r.sq,
    edf_week = if ("s(week_of_year)" %in% rownames(s$s.table)) s$s.table["s(week_of_year)", "edf"] else NA,
    edf_temp = if ("s(meanTemp)" %in% rownames(s$s.table)) s$s.table["s(meanTemp)", "edf"] else NA,
    edf_sal  = if ("s(meanSal)" %in% rownames(s$s.table)) s$s.table["s(meanSal)", "edf"] else NA,
    p_week = if ("s(week_of_year)" %in% rownames(s$s.table)) s$s.table["s(week_of_year)", "p-value"] else NA,
    p_temp = if ("s(meanTemp)" %in% rownames(s$s.table)) s$s.table["s(meanTemp)", "p-value"] else NA,
    p_sal  = if ("s(meanSal)" %in% rownames(s$s.table)) s$s.table["s(meanSal)", "p-value"] else NA
  )
}


library(mgcv)
library(dplyr)
library(tidyr)

all_species <- unique(df_model$species)

all_results <- list()

for (sp in all_species) {

  # --------------------------
  # Prepare data (ONE dataset per species)
  # --------------------------
  df_sp <- df_model %>%
    filter(species == sp) %>%
    filter(
      is.finite(scaleLogConc),
      is.finite(week_of_year),
      is.finite(meanTemp),
      is.finite(meanSal)
    ) %>%
    mutate(region = factor(region))

  # Skip if too small
  if (nrow(df_sp) < 30) next

  # --------------------------
  # Model 1 — season
  # --------------------------
  m_season <- gam(
    scaleLogConc ~
      s(week_of_year, bs = "cc") +
      s(region, bs = "re"),
    data = df_sp,
    method = "REML",
    knots = list(week_of_year = c(0.5, 52.5))
  )

  # --------------------------
  # Model 2 — season + temp
  # --------------------------
  m_temp <- gam(
    scaleLogConc ~
      s(week_of_year, bs = "cc") +
      s(meanTemp) +
      s(region, bs = "re"),
    data = df_sp,
    method = "REML",
    knots = list(week_of_year = c(0.5, 52.5))
  )

  # --------------------------
  # Model 3 — season + temp + sal
  # --------------------------
  m_temp_sal <- gam(
    scaleLogConc ~
      s(week_of_year, bs = "cc") +
      s(meanTemp) +
      s(meanSal) +
      s(region, bs = "re"),
    data = df_sp,
    method = "REML",
    knots = list(week_of_year = c(0.5, 52.5))
  )

  # --------------------------
  # Collect results
  # --------------------------
  res <- bind_rows(
    extract_gam_summary(m_season, sp, "season"),
    extract_gam_summary(m_temp, sp, "season + temp"),
    extract_gam_summary(m_temp_sal, sp, "season + temp + sal")
  )

  all_results[[sp]] <- res
}


gam_summary_table <- bind_rows(all_results) %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

gam_summary_table


gam_summary_wide <- gam_summary_table %>%
  select(species, model, dev_expl) %>%
  pivot_wider(names_from = model, values_from = dev_expl) %>%
  mutate(
    temp_gain = `season + temp` - season,
    sal_gain  = `season + temp + sal` - `season + temp`
  )

gam_summary_wide














#####inspect ciona

df_sp <- df_model %>%
  filter(species == "Ciona intestinalis") %>%
  filter(
    is.finite(scaleLogConc),
    is.finite(week_of_year),
    is.finite(meanTemp),
    is.finite(meanSal)
  ) %>%
  mutate(region = factor(region))

m_temp_sal <- mgcv::gam(
  scaleLogConc ~
    s(week_of_year, bs = "cc") +
    s(meanTemp) +
    s(meanSal) +
    s(region, bs = "re"),
  data = df_sp,
  method = "REML",
  knots = list(week_of_year = c(0.5, 52.5))
)

plot(m_temp_sal)
