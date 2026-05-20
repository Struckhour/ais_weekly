library(dplyr)
library(tidyr)
library(stringr)
library(mgcv)
library(ggplot2)

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




######################
#GAM summary
######################

library(dplyr)
library(purrr)
library(tibble)
library(mgcv)

all_species <- unique(df_model$species)

fit_gam_models <- function(sp) {

  df_sp <- df_model %>%
    filter(species == sp) %>%
    filter(
      is.finite(scaleLogConc),
      is.finite(week_of_year),
      is.finite(meanTemp),
      is.finite(meanSal)
    ) %>%
    mutate(region = factor(region, levels = c("MAG","PEI","HAL","BOF","GOM")))

  if (nrow(df_sp) < 30) return(NULL)

  m_season <- gam(
    scaleLogConc ~ s(week_of_year, bs = "cc"),
    data = df_sp,
    method = "REML",
    knots = list(week_of_year = c(0.5, 52.5))
  )

  m_temp <- gam(
    scaleLogConc ~
      s(week_of_year, by = region, bs = "cc") +
      s(meanTemp) +
      region,
    data = df_sp,
    method = "REML",
    knots = list(week_of_year = c(0.5, 52.5))
  )

  m_temp_sal <- gam(
    scaleLogConc ~
      s(week_of_year, by = region, bs = "cc") +
      s(meanTemp) +
      s(meanSal) +
      region,
    data = df_sp,
    method = "REML",
    knots = list(week_of_year = c(0.5, 52.5))
  )

  list(
    species = sp,
    data = df_sp,
    season = m_season,
    temp = m_temp,
    temp_sal = m_temp_sal
  )
}

gam_fits <- all_species %>%
  setNames(all_species) %>%
  map(fit_gam_models) %>%
  compact()

gam_model_summary <- imap_dfr(gam_fits, function(x, sp) {

  tibble(
    species = sp,
    model = c("season", "season_temp", "season_temp_sal"),
    n = nrow(x$data),
    deviance_explained = c(
      summary(x$season)$dev.expl,
      summary(x$temp)$dev.expl,
      summary(x$temp_sal)$dev.expl
    ),
    AIC = c(
      AIC(x$season),
      AIC(x$temp),
      AIC(x$temp_sal)
    )
  )
})

gam_model_summary


extract_smooth_table <- function(fit, species_name, model_name) {

  s <- summary(fit)

  as.data.frame(s$s.table) %>%
    rownames_to_column("term") %>%
    as_tibble() %>%
    transmute(
      species = species_name,
      model = model_name,
      term,
      edf,
      F = `F`,
      p_value = `p-value`
    )
}

gam_smooth_summary <- imap_dfr(gam_fits, function(x, sp) {
  bind_rows(
    extract_smooth_table(x$season, sp, "season"),
    extract_smooth_table(x$temp, sp, "season_temp"),
    extract_smooth_table(x$temp_sal, sp, "season_temp_sal")
  )
})

gam_model_summary %>%
  arrange(species, AIC)

printable_summary_smooth <- gam_smooth_summary %>%
  filter(model == "season_temp_sal") %>%
  arrange(species, term)
print(printable_summary_smooth, n=Inf, width=Inf)


season_smooth_df <- purrr::imap_dfr(gam_fits, function(x, sp) {

  fit <- x$temp_sal
  dat <- x$data

  pred_grid <- tidyr::expand_grid(
    week_of_year = seq(1, 52, length.out = 300),
    region = unique(dat$region)
  ) %>%
    dplyr::mutate(
      meanTemp = mean(dat$meanTemp, na.rm = TRUE),
      meanSal = mean(dat$meanSal, na.rm = TRUE)
    )

  pred <- predict(
    fit,
    newdata = pred_grid,
    se.fit = TRUE,
    type = "response"
  )

  pred_grid %>%
    dplyr::mutate(
      species = sp,
      fit = pred$fit,
      se = pred$se.fit,
      lower = fit - 1.96 * se,
      upper = fit + 1.96 * se
    )
})

species_order <- c(
  "Membranipora membranacea",
  "Botrylloides violaceus",
  "Didemnum vexillum",
  "Ciona intestinalis",
  "Carcinus maenas"
)

season_smooth_df <- season_smooth_df %>%
  dplyr::mutate(
    species = factor(species, levels = species_order)
  )

p_season_smooths <- ggplot(
  season_smooth_df,
  aes(x = week_of_year, y = fit, color = region, fill = region)
) +
  geom_ribbon(
    aes(ymin = lower, ymax = upper),
    alpha = 0.15,
    color = NA
  ) +
  geom_line(linewidth = 1) +
  facet_wrap(~ species, scales = "free_y") +
  scale_color_manual(values = hybrid_color) +
  scale_fill_manual(values = hybrid_color) +
  scale_x_continuous(
    breaks = c(1, 13, 26, 39, 52),
    labels = c("Jan", "Apr", "Jul", "Oct", "Dec")
  ) +
  labs(
    x = "Week of year",
    y = "Predicted scaled log qPCR concentration",
    color = "Region",
    fill = "Region"
  ) +
  theme_classic() +
  theme(
    legend.position = "bottom",
    strip.text = element_text(face = "italic"),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_season_smooths
ggsave("manuscript_figures/figure_8.png", p_season_smooths, width = 8, height = 6, dpi = 300)


temp_effect_df <- purrr::imap_dfr(gam_fits, function(x, sp) {

  fit <- x$temp_sal
  dat <- x$data

  temp_seq <- seq(
    min(dat$meanTemp, na.rm = TRUE),
    max(dat$meanTemp, na.rm = TRUE),
    length.out = 200
  )

  pred_grid <- tibble::tibble(
    meanTemp = temp_seq,
    meanSal = mean(dat$meanSal, na.rm = TRUE),
    week_of_year = 26,  # mid-year
    region = unique(dat$region)[1]  # arbitrary, doesn't matter here
  )

  pred <- predict(
    fit,
    newdata = pred_grid,
    se.fit = TRUE,
    type = "terms"
  )

  # extract ONLY temperature term
  temp_term <- pred$fit[, "s(meanTemp)"]
  temp_se <- pred$se.fit[, "s(meanTemp)"]

  pred_grid %>%
    mutate(
      species = sp,
      fit = temp_term,
      lower = fit - 1.96 * temp_se,
      upper = fit + 1.96 * temp_se
    )
})

species_order <- c(
  "Membranipora membranacea",
  "Botrylloides violaceus",
  "Didemnum vexillum",
  "Ciona intestinalis",
  "Carcinus maenas"
)

temp_effect_df <- temp_effect_df %>%
  dplyr::mutate(
    species = factor(species, levels = species_order)
  )

p_temp <- ggplot(
  temp_effect_df,
  aes(x = meanTemp, y = fit)
) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2) +
  geom_line(linewidth = 1) +
  facet_wrap(~ species, scales = "free_y") +
  labs(
    x = "Temperature (°C)",
    y = "Partial effect on qPCR concentration"
  ) +
  theme_classic()

p_temp
ggsave("manuscript_figures/figure_9.png", p_temp, width = 8, height = 6, dpi = 300)
