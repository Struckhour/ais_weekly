library(randomForest)
library(zoo)

source('./AIS_eDNA_data_prep.R')

sp <- "Membranipora membranacea"
sp <- "Botrylloides violaceus"
sp <- "Ciona intestinalis"
sp <- "Carcinus maenas"

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


fit <- randomForest(
  scaleLogConc ~ week_of_year + meanTemp + meanSal + meanPH + meanLat,
  data = df_sp_clean
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
#no week_of_year as input

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
  data = df_sp_marine
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



