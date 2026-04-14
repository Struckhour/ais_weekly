library(randomForest)
library(zoo)

source('./AIS_eDNA_data_prep.R')

df_sp <- dfWeeks %>%
  filter(species == sp)


df_sp_clean <- df_sp %>%
  as.data.frame() %>%
  dplyr::filter(
    is.finite(scaleLogConc),
    is.finite(week_of_year),
    is.finite(meanTemp),
    is.finite(meanSal),
    is.finite(meanPH),
    is.finite(meanLat)
  )



extract_rf_summary <- function(model_result, species_name, model_name) {
  fit <- model_result$fit

  imp <- randomForest::importance(fit)

  if ("%IncMSE" %in% colnames(imp)) {
    imp_vals <- imp[, "%IncMSE"]
  } else {
    imp_vals <- rep(NA_real_, nrow(imp))
  }

  imp_df <- tibble::tibble(
    variable = rownames(imp),
    incMSE = imp_vals
  ) %>%
    tidyr::pivot_wider(
      names_from = variable,
      values_from = incMSE,
      names_prefix = "IncMSE_"
    )

  dplyr::bind_cols(
    tibble::tibble(
      species = species_name,
      model = model_name
    ),
    model_result$metrics,
    imp_df
  )
}

run_rf_model <- function(df_sp_clean,
                         model_formula,
                         profile_vars,
                         importance = TRUE,
                         ntree = 500) {

  df_model <- df_sp_clean

  fit <- randomForest::randomForest(
    formula = model_formula,
    data = df_model,
    importance = importance,
    ntree = ntree
  )

  profiles <- df_model %>%
    dplyr::group_by(region, week_of_year) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(profile_vars), ~ mean(.x, na.rm = TRUE)),
      .groups = "drop"
    ) %>%
    tidyr::complete(region, week_of_year = 1:52) %>%
    dplyr::group_by(region) %>%
    dplyr::arrange(week_of_year) %>%
    dplyr::mutate(
      dplyr::across(
        dplyr::all_of(profile_vars),
        ~ zoo::na.approx(.x, week_of_year, na.rm = FALSE, rule = 2)
      )
    ) %>%
    dplyr::ungroup()

  profiles$pred <- predict(fit, newdata = profiles)
  df_model$pred <- predict(fit, newdata = df_model)

  r2_train <- cor(df_model$pred, df_model$scaleLogConc)^2
  rmse <- sqrt(mean((df_model$scaleLogConc - df_model$pred)^2))
  var_explained <- fit$rsq[length(fit$rsq)] * 100

  list(
    fit = fit,
    profiles = profiles,
    data = df_model,
    metrics = tibble::tibble(
      r2_train = r2_train,
      rmse = rmse,
      var_explained = var_explained,
      n = nrow(df_model)
    )
  )
}

prep_species_data <- function(df, species_name, required_vars) {
  df %>%
    dplyr::filter(species == species_name) %>%
    as.data.frame() %>%
    dplyr::filter(dplyr::if_all(dplyr::all_of(required_vars), is.finite))
}








#####################
#LOOP THROUGH SPECIES AND BUILD MODELS
#####################

all_species <- unique(dfWeeks$species)
all_results <- list()

for (sp in all_species) {

  df_sp_clean <- prep_species_data(
    df = dfWeeks,
    species_name = sp,
    required_vars = c("scaleLogConc", "week_of_year", "meanTemp", "meanSal", "meanPH", "meanLat")
  )

  m1 <- run_rf_model(
    df_sp_clean = df_sp_clean,
    model_formula = scaleLogConc ~ week_of_year + meanTemp + meanSal + meanPH + meanLat,
    profile_vars = c("meanTemp", "meanSal", "meanPH", "meanLat")
  )

  m2 <- run_rf_model(
    df_sp_clean = df_sp_clean,
    model_formula = scaleLogConc ~ meanTemp + meanSal + meanPH,
    profile_vars = c("meanTemp", "meanSal", "meanPH")
  )

  m3 <- run_rf_model(
    df_sp_clean = df_sp_clean,
    model_formula = scaleLogConc ~ meanTemp + meanSal + meanPH + meanLat,
    profile_vars = c("meanTemp", "meanSal", "meanPH", "meanLat")
  )


  df_sp_typical <- df_sp_clean %>%
    dplyr::filter(meanSal > 25, meanPH > 7.8)
  m4 <- run_rf_model(
    df_sp_clean = df_sp_typical,
    model_formula = scaleLogConc ~ meanTemp + meanSal + meanPH + meanLat,
    profile_vars = c("meanTemp", "meanSal", "meanPH", "meanLat")
  )

  summary_table <- dplyr::bind_rows(
    extract_rf_summary(m1, sp, "week + temp + sal + pH + lat"),
    extract_rf_summary(m2, sp, "temp + sal + pH"),
    extract_rf_summary(m3, sp, "temp + sal + pH + lat"),
    extract_rf_summary(m4, sp, "typical sal/pH + temp + sal + pH + lat")
  ) %>%
    dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 3)))

  all_results[[sp]] <- summary_table
}

summary_table <- dplyr::bind_rows(all_results)
























sp <- "Membranipora membranacea"
sp <- "Botrylloides violaceus"
sp <- "Ciona intestinalis"
sp <- "Carcinus maenas"



fit <- randomForest(
  scaleLogConc ~ week_of_year + meanTemp + meanSal + meanPH + meanLat,
  data = df_sp_clean,
  importance = TRUE
)


profiles <- df_sp_clean %>%
  group_by(region, week_of_year) %>%
  summarise(
    meanTemp = mean(meanTemp, na.rm = TRUE),
    meanSal  = mean(meanSal, na.rm = TRUE),
    meanPH  = mean(meanPH, na.rm = TRUE),
    meanLat = mean(meanLat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(region, week_of_year = 1:52)


profiles <- profiles %>%
  group_by(region) %>%
  arrange(week_of_year) %>%
  mutate(
    meanTemp = zoo::na.approx(meanTemp, week_of_year, na.rm = FALSE, rule = 2),
    meanSal  = zoo::na.approx(meanSal,  week_of_year, na.rm = FALSE, rule = 2),
    meanPH  = zoo::na.approx(meanPH,  week_of_year, na.rm = FALSE, rule = 2),
    meanLat  = zoo::na.approx(meanLat,  week_of_year, na.rm = FALSE, rule = 2)
  ) %>%
  ungroup()


profiles$pred <- predict(fit, newdata = profiles)




ggplot(profiles, aes(x = week_of_year, y = pred, color = region)) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(breaks = seq(1, 52, by = 4)) +
  labs(
    title = sp,
    x = "Week of year",
    y = "Predicted scaled log concentration",
    color = "Region"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank()
  )


peaks <- profiles %>%
  group_by(region) %>%
  summarise(
    peak_week = week_of_year[which.max(pred)],
    peak_val  = max(pred),
    .groups = "drop"
  )

ggplot(profiles, aes(x = week_of_year, y = pred, color = region)) +
  geom_line(linewidth = 1.2) +
  geom_point(data = peaks, aes(x = peak_week, y = peak_val), size = 3) +
  labs(
    title = paste0(sp, " (predicted seasonal curves)"),
    x = "Week of year",
    y = "Predicted scaled log concentration"
  ) +
  theme_minimal()


df_sp_clean$pred <- predict(fit, newdata = df_sp_clean)



ggplot(df_sp_clean, aes(x = pred, y = scaleLogConc)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(
    x = "Predicted",
    y = "Observed",
    title = paste("Predicted vs Observed —", sp)
  ) +
  theme_minimal()



#VALUES
importance(fit)
varImpPlot(fit)
cor(df_sp_clean$pred, df_sp_clean$scaleLogConc)^2
rmse <- sqrt(mean((df_sp_clean$scaleLogConc - df_sp_clean$pred)^2))
rmse

fit

# df_sp_base <- as.data.frame(df_sp_clean)

partialPlot(fit, df_sp_clean, x.var = "week_of_year")
partialPlot(fit, df_sp_clean, x.var = "meanTemp")
partialPlot(fit, df_sp_clean, x.var = "meanSal")
partialPlot(fit, df_sp_clean, x.var = "meanPH")
partialPlot(fit, df_sp_clean, x.var = "meanLat")












###############################################################
###############################################################
###############################################################
#no week_of_year or Lat as input

fit <- randomForest(
  scaleLogConc ~ meanTemp + meanSal + meanPH,
  data = df_sp_clean,
  importance = TRUE
)


profiles <- df_sp_clean %>%
  group_by(region, week_of_year) %>%
  summarise(
    meanTemp = mean(meanTemp, na.rm = TRUE),
    meanSal  = mean(meanSal, na.rm = TRUE),
    meanPH = mean(meanPH, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(region, week_of_year = 1:52)


profiles <- profiles %>%
  group_by(region) %>%
  arrange(week_of_year) %>%
  mutate(
    meanTemp = zoo::na.approx(meanTemp, week_of_year, na.rm = FALSE, rule = 2),
    meanSal  = zoo::na.approx(meanSal,  week_of_year, na.rm = FALSE, rule = 2),
    meanPH = zoo::na.approx(meanPH, week_of_year, na.rm = FALSE, rule = 2)
  ) %>%
  ungroup()


profiles$pred <- predict(fit, newdata = profiles)




ggplot(profiles, aes(x = week_of_year, y = pred, color = region)) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(breaks = seq(1, 52, by = 4)) +
  labs(
    title = sp,
    x = "Week of year",
    y = "Predicted scaled log concentration",
    color = "Region"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank()
  )

df_sp_clean$pred <- predict(fit, newdata = df_sp_clean)

ggplot(df_sp_clean, aes(x = pred, y = scaleLogConc)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(
    x = "Predicted",
    y = "Observed",
    title = paste("Predicted vs Observed —", sp)
  ) +
  theme_minimal()


#VALUES
importance(fit)
varImpPlot(fit)

partialPlot(fit, df_sp_clean, x.var = "meanTemp")
partialPlot(fit, df_sp_clean, x.var = "meanSal")
partialPlot(fit, df_sp_clean, x.var = "meanPH")


fit



###############################################################
###############################################################
###############################################################
#no week_of_year, but latitude added as input

fit <- randomForest(
  scaleLogConc ~ meanTemp + meanSal + meanPH + meanLat,
  data = df_sp_clean,
  importance = TRUE
)


profiles <- df_sp_clean %>%
  group_by(region, week_of_year) %>%
  summarise(
    meanTemp = mean(meanTemp, na.rm = TRUE),
    meanSal  = mean(meanSal, na.rm = TRUE),
    meanPH = mean(meanPH, na.rm = TRUE),
    meanLat = mean(meanLat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(region, week_of_year = 1:52)


profiles <- profiles %>%
  group_by(region) %>%
  arrange(week_of_year) %>%
  mutate(
    meanTemp = zoo::na.approx(meanTemp, week_of_year, na.rm = FALSE, rule = 2),
    meanSal  = zoo::na.approx(meanSal,  week_of_year, na.rm = FALSE, rule = 2),
    meanPH = zoo::na.approx(meanPH, week_of_year, na.rm = FALSE, rule = 2),
    meanLat = zoo::na.approx(meanLat, week_of_year, na.rm = FALSE, rule = 2)
  ) %>%
  ungroup()


profiles$pred <- predict(fit, newdata = profiles)




ggplot(profiles, aes(x = week_of_year, y = pred, color = region)) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(breaks = seq(1, 52, by = 4)) +
  labs(
    title = sp,
    x = "Week of year",
    y = "Predicted scaled log concentration",
    color = "Region"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank()
  )

df_sp_clean$pred <- predict(fit, newdata = df_sp_clean)

ggplot(df_sp_clean, aes(x = pred, y = scaleLogConc)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(
    x = "Predicted",
    y = "Observed",
    title = paste("Predicted vs Observed —", sp)
  ) +
  theme_minimal()


#VALUES
importance(fit)
varImpPlot(fit)

partialPlot(fit, df_sp_clean, x.var = "meanTemp")
partialPlot(fit, df_sp_clean, x.var = "meanSal")
partialPlot(fit, df_sp_clean, x.var = "meanPH")
partialPlot(fit, df_sp_clean, x.var = "meanLat")

fit







#################################################
#TEST THRESHOLD???
#################################################
obs_peaks <- df_sp_clean %>%
  group_by(region, week_of_year) %>%
  summarise(
    obs = mean(scaleLogConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(region) %>%
  summarise(
    peak_week = week_of_year[which.max(obs)],
    .groups = "drop"
  )

ph_profiles <- df_sp_clean %>%
  group_by(region, week_of_year) %>%
  summarise(
    meanPH = mean(meanPH, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(region, week_of_year = 1:52) %>%
  group_by(region) %>%
  arrange(week_of_year) %>%
  mutate(
    meanPH = zoo::na.approx(meanPH, week_of_year, na.rm = FALSE, rule = 2)
  ) %>%
  ungroup()

ph_threshold <- 8.0
ph_cross <- ph_profiles %>%
  group_by(region) %>%
  summarise(
    cross_week = min(week_of_year[meanPH >= ph_threshold], na.rm = TRUE),
    .groups = "drop"
  )

comparison <- left_join(obs_peaks, ph_cross, by = "region")

comparison

comparison <- comparison %>%
  mutate(diff = peak_week - cross_week)

ggplot(comparison, aes(x = cross_week, y = peak_week)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(
    x = "Week pH crosses threshold",
    y = "Peak week",
    title = "Does pH threshold predict peak timing?"
  ) +
  theme_minimal()

cor.test(comparison$cross_week, comparison$peak_week)
lm(peak_week ~ cross_week, data = comparison)
###############################################################
###############################################################
###############################################################
















###############################################################
###############################################################
###############################################################
#TYPICAL SALINITY and PH
df_sp_marine <- df_sp_clean %>%
  filter(meanSal > 25, meanPH < 8.5 & meanPH > 7.5)

fit_marine <- randomForest(
  scaleLogConc ~ meanTemp + meanSal + meanPH + meanLat,
  data = df_sp_marine,
  importance = TRUE
)



profiles <- df_sp_marine %>%
  group_by(region, week_of_year) %>%
  summarise(
    meanTemp = mean(meanTemp, na.rm = TRUE),
    meanSal  = mean(meanSal, na.rm = TRUE),
    meanPH = mean(meanPH, na.rm = TRUE),
    meanLat = mean(meanLat, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  complete(region, week_of_year = 1:52)


profiles <- profiles %>%
  group_by(region) %>%
  arrange(week_of_year) %>%
  mutate(
    meanTemp = zoo::na.approx(meanTemp, week_of_year, na.rm = FALSE, rule = 2),
    meanSal  = zoo::na.approx(meanSal,  week_of_year, na.rm = FALSE, rule = 2),
    meanPH = zoo::na.approx(meanPH, week_of_year, na.rm = FALSE, rule = 2),
    meanLat = zoo::na.approx(meanLat, week_of_year, na.rm = FALSE, rule = 2)
  ) %>%
  ungroup()


profiles$pred <- predict(fit_marine, newdata = profiles)




ggplot(profiles, aes(x = week_of_year, y = pred, color = region)) +
  geom_line(linewidth = 1.2) +
  scale_x_continuous(breaks = seq(1, 52, by = 4)) +
  labs(
    title = sp,
    x = "Week of year",
    y = "Predicted scaled log concentration",
    color = "Region"
  ) +
  theme_minimal() +
  theme(
    panel.grid.minor = element_blank()
  )

df_sp_marine$pred <- predict(fit_marine, newdata = df_sp_marine)

ggplot(df_sp_marine, aes(x = pred, y = scaleLogConc)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(
    x = "Predicted",
    y = "Observed",
    title = paste("Predicted vs Observed —", sp)
  ) +
  theme_minimal()


cor(df_sp_marine$pred, df_sp_marine$scaleLogConc)^2
rmse <- sqrt(mean((df_sp_marine$scaleLogConc - df_sp_marine$pred)^2))
rmse

fit_marine
#VALUES
importance(fit_marine)
varImpPlot(fit_marine)


partialPlot(fit_marine, df_sp_marine, x.var = "week_of_year")
partialPlot(fit_marine, df_sp_marine, x.var = "meanTemp")
partialPlot(fit_marine, df_sp_marine, x.var = "meanSal")
partialPlot(fit_marine, df_sp_marine, x.var = "meanPH")

###############################################################
###############################################################
###############################################################



