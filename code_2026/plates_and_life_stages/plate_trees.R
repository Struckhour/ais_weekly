source('./AIS_eDNA_data_prep.R')

library(dplyr)
library(tidyr)
library(stringr)

state_levels <- c("R", "G", "B", "S", "D")

#--------------------------
# 1. Prep plate abundance weekly
#--------------------------
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

#--------------------------
# 2. Prep weekly stage indicators
#--------------------------
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

#--------------------------
# 3. Combine weekly plate abundance + stage indicators
#--------------------------
plate_weekly <- plate_abund_weekly %>%
  full_join(
    plate_stage_weekly,
    by = c("species", "region", "sampWeek")
  )

# make sure all stage columns exist
for (st in state_levels) {
  if (!st %in% names(plate_weekly)) {
    plate_weekly[[st]] <- 0L
  }
}

plate_weekly <- plate_weekly %>%
  mutate(
    across(all_of(state_levels), ~ tidyr::replace_na(.x, 0L))
  ) %>%
  select(
    species,
    region,
    sampWeek,
    mean_plate_abund,
    n_plate_abund_obs,
    all_of(state_levels)
  )

#--------------------------
# 4. Prep dfWeeks keys for join
#--------------------------
# Assumes dfWeeks uses week_of_year
qpcr_weekly <- dfWeeks %>%
  mutate(
    region = as.character(region),
    species = as.character(species),
    week_of_year = as.numeric(week_of_year)
  )

#--------------------------
# 5. Join everything together
#--------------------------
df_model <- qpcr_weekly %>%
  left_join(
    plate_weekly %>%
      rename(week_of_year = sampWeek),
    by = c("species", "region", "week_of_year")
  ) %>%
  mutate(
    across(all_of(state_levels), ~ tidyr::replace_na(.x, 0L))
  )

# inspect
# glimpse(df_model)
# View(df_model)
#
# summary(df_model$mean_plate_abund)
#
# colSums(
#   dplyr::select(df_model, all_of(state_levels)),
#   na.rm = TRUE
# )
#
# df_model %>%
#   count(species, region, week_of_year) %>%
#   filter(n > 1)



df_tree <- df_model %>%
  filter(!is.na(mean_plate_abund)) %>%
  mutate(
    region = factor(region),
    species = factor(species)
  )

names(df_tree)
summary(df_tree$mean_plate_abund)

df_tree %>%
  select(species, region, week_of_year, mean_plate_abund, R, G, B, S, D)

#######################################


library(dplyr)
library(rpart)
library(rpart.plot)

#--------------------------
# 1. Keep only modeling columns we need
#--------------------------
df_tree_model <- df_tree %>%
  select(
    species,
    region,
    week_of_year,
    scaleLogConc,
    mean_plate_abund,
    R, G, B, S, D
  ) %>%
  filter(
    !is.na(scaleLogConc),
    !is.na(mean_plate_abund)
  ) %>%
  mutate(
    region = factor(region),
    species = factor(species)
  )

# quick check
df_tree_model %>%
  group_by(species) %>%
  summarise(
    n_rows = n(),
    n_regions = n_distinct(region),
    .groups = "drop"
  )





species_list <- levels(df_tree_model$species)

tree_models <- lapply(species_list, function(sp) {
  df_sp <- df_tree_model %>%
    filter(species == sp)

  fit <- rpart(
    scaleLogConc ~ mean_plate_abund + R + G + B + S + D + region,
    data = df_sp,
    method = "anova",
    control = rpart.control(
      cp = 0.01,
      minsplit = 10,
      minbucket = 5,
      maxdepth = 4
    )
  )

  list(
    species = sp,
    data = df_sp,
    fit = fit
  )
})

names(tree_models) <- species_list



for (sp in names(tree_models)) {
  cat("\n=============================\n")
  cat("Species:", sp, "\n")
  cat("=============================\n")
  print(tree_models[[sp]]$fit)
  cat("\n")
  printcp(tree_models[[sp]]$fit)
}


for (sp in names(tree_models)) {
  rpart.plot(
    tree_models[[sp]]$fit,
    main = sp,
    type = 2,
    extra = 101,
    fallen.leaves = TRUE,
    tweak = 1.1
  )
}




##########################################
#NO REGION
##########################################
tree_models_noregion <- lapply(species_list, function(sp) {
  df_sp <- df_tree_model %>%
    filter(species == sp)

  fit <- rpart(
    scaleLogConc ~ mean_plate_abund + R + G + B + S + D,
    data = df_sp,
    method = "anova",
    model = TRUE,
    control = rpart.control(
      cp = 0.01,
      minsplit = 10,
      minbucket = 5,
      maxdepth = 4
    )
  )

  list(
    species = sp,
    data = df_sp,
    fit = fit
  )
})

names(tree_models_noregion) <- species_list


for (sp in names(tree_models_noregion)) {
  cat("\n=============================\n")
  cat("Species:", sp, "\n")
  cat("=============================\n")
  print(tree_models_noregion[[sp]]$fit)
  cat("\n")
  printcp(tree_models_noregion[[sp]]$fit)
}





































fit_and_eval <- function(df, formula) {
  fit <- lm(formula, data = df)

  pred <- predict(fit, newdata = df)

  tibble(
    rmse = sqrt(mean((df$scaleLogConc - pred)^2, na.rm = TRUE)),
    rsq = cor(df$scaleLogConc, pred, use = "complete.obs")^2
  )
}



species_list <- levels(df_tree_model$species)

model_results <- lapply(species_list, function(sp) {

  df_sp <- df_tree_model %>%
    filter(species == sp)

  res_A <- fit_and_eval(df_sp, scaleLogConc ~ mean_plate_abund)
  res_B <- fit_and_eval(df_sp, scaleLogConc ~ mean_plate_abund + R + G + B + S + D)
  res_C <- fit_and_eval(df_sp, scaleLogConc ~ mean_plate_abund + week_of_year)

  bind_rows(
    res_A %>% mutate(model = "A: plate only"),
    res_B %>% mutate(model = "B: plate + stage"),
    res_C %>% mutate(model = "C: plate + time")
  ) %>%
    mutate(species = sp)
})

model_results <- bind_rows(model_results) %>%
  select(species, model, rmse, rsq)



model_results %>%
  arrange(species, rmse)































#################
#Cyclic Time
#################


df_tree_model <- df_tree_model %>%
  mutate(
    week_sin = sin(2 * pi * week_of_year / 52),
    week_cos = cos(2 * pi * week_of_year / 52)
  )


library(rpart)

model_formulas <- list(
  "1_plate"              = scaleLogConc ~ mean_plate_abund,
  "2_plate_stage"        = scaleLogConc ~ mean_plate_abund + R + G + B + S + D,
  "3_plate_stage_cycle"  = scaleLogConc ~ mean_plate_abund + R + G + B + S + D + week_sin + week_cos,
  "4_cycle_only"         = scaleLogConc ~ week_sin + week_cos,
  "5_plate_cycle"        = scaleLogConc ~ mean_plate_abund + week_sin + week_cos
)


species_list <- levels(df_tree_model$species)

tree_models_5 <- lapply(species_list, function(sp) {
  df_sp <- df_tree_model %>%
    filter(
      species == sp,
      !is.na(scaleLogConc),
      !is.na(mean_plate_abund)
    )

  fits <- lapply(model_formulas, function(fm) {
    rpart(
      formula = fm,
      data = df_sp,
      method = "anova",
      model = TRUE,
      control = rpart.control(
        cp = 0.01,
        minsplit = 10,
        minbucket = 5,
        maxdepth = 4
      )
    )
  })

  list(
    species = sp,
    data = df_sp,
    fits = fits
  )
})

names(tree_models_5) <- species_list


for (sp in names(tree_models_5)) {
  cat("\n====================================\n")
  cat("Species:", sp, "\n")
  cat("====================================\n")

  for (nm in names(tree_models_5[[sp]]$fits)) {
    cat("\n-----------------------------\n")
    cat("Model:", nm, "\n")
    cat("-----------------------------\n")
    print(tree_models_5[[sp]]$fits[[nm]])
    cat("\n")
    printcp(tree_models_5[[sp]]$fits[[nm]])
  }
}


tree_results_5 <- lapply(names(tree_models_5), function(sp) {
  obj <- tree_models_5[[sp]]
  df_sp <- obj$data

  bind_rows(lapply(names(obj$fits), function(nm) {
    fit <- obj$fits[[nm]]
    pred <- predict(fit, newdata = df_sp)

    tibble(
      species = sp,
      model = nm,
      rmse = sqrt(mean((df_sp$scaleLogConc - pred)^2, na.rm = TRUE)),
      rsq = cor(df_sp$scaleLogConc, pred, use = "complete.obs")^2
    )
  }))
}) %>%
  bind_rows() %>%
  arrange(species, rmse)

tree_results_5



library(gt)


tree_results_5 %>%
  arrange(species, rmse) %>%
  group_by(species) %>%
  mutate(
    species_display = ifelse(row_number() == 1, species, ""),
    is_best = rmse == min(rmse)
  ) %>%
  ungroup() %>%
  select(
    species = species_display,
    model,
    rmse,
    rsq,
    is_best
  ) %>%
  gt() %>%
  tab_header(
    title = "Regression Tree Model Performance",
    subtitle = "Scaled eDNA concentration (within region)"
  ) %>%
  fmt_number(
    columns = c(rmse, rsq),
    decimals = 3
  ) %>%
  cols_label(
    species = "Species",
    model   = "Model",
    rmse    = "RMSE",
    rsq     = "R²"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = model,
      rows = is_best
    )
  ) %>%
  cols_hide(is_best) %>%
  tab_options(
    table.background.color = "white",
    heading.background.color = "white",
    column_labels.background.color = "white"
  )


fit_ciona_m3 <- tree_models_5[["Ciona intestinalis"]]$fits[["3_plate_stage_cycle"]]
library(rpart.plot)

rpart.plot(
  fit_ciona_m3,
  main = "Ciona intestinalis — Plate + Stage + Cyclic Time",
  type = 2,
  extra = 101,
  fallen.leaves = TRUE,
  tweak = 1.1
)












































###############################
#TIME LAGS
###############################

df_tree_lagged <- df_tree_model %>%
  arrange(species, region, week_of_year) %>%
  group_by(species, region) %>%
  mutate(
    y_lag0 = scaleLogConc,
    y_lag1 = lead(scaleLogConc, 1),
    y_lag2 = lead(scaleLogConc, 2),
    y_lag3 = lead(scaleLogConc, 3)
  ) %>%
  ungroup()

run_tree_models_for_response <- function(df, response_var) {

  model_formulas <- list(
    "1_plate"             = as.formula(paste(response_var, "~ mean_plate_abund")),
    "2_plate_stage"       = as.formula(paste(response_var, "~ mean_plate_abund + R + G + B + S + D")),
    "3_plate_stage_cycle" = as.formula(paste(response_var, "~ mean_plate_abund + R + G + B + S + D + week_sin + week_cos")),
    "4_cycle_only"        = as.formula(paste(response_var, "~ week_sin + week_cos")),
    "5_plate_cycle"       = as.formula(paste(response_var, "~ mean_plate_abund + week_sin + week_cos"))
  )

  species_list <- levels(df$species)

  tree_models <- lapply(species_list, function(sp) {
    df_sp <- df %>%
      filter(
        species == sp,
        !is.na(.data[[response_var]]),
        !is.na(mean_plate_abund)
      )

    fits <- lapply(model_formulas, function(fm) {
      rpart(
        formula = fm,
        data = df_sp,
        method = "anova",
        model = TRUE,
        control = rpart.control(
          cp = 0.01,
          minsplit = 10,
          minbucket = 5,
          maxdepth = 4
        )
      )
    })

    list(
      species = sp,
      data = df_sp,
      fits = fits
    )
  })

  names(tree_models) <- species_list

  tree_results <- lapply(names(tree_models), function(sp) {
    obj <- tree_models[[sp]]
    df_sp <- obj$data

    bind_rows(lapply(names(obj$fits), function(nm) {
      fit <- obj$fits[[nm]]
      pred <- predict(fit, newdata = df_sp)
      obs  <- df_sp[[response_var]]

      tibble(
        species = sp,
        model = nm,
        rmse = sqrt(mean((obs - pred)^2, na.rm = TRUE)),
        rsq = cor(obs, pred, use = "complete.obs")^2,
        n = nrow(df_sp)
      )
    }))
  }) %>%
    bind_rows() %>%
    arrange(species, rmse)

  list(
    models = tree_models,
    results = tree_results
  )
}



lag_runs <- list(
  lag0 = run_tree_models_for_response(df_tree_lagged, "y_lag0"),
  lag1 = run_tree_models_for_response(df_tree_lagged, "y_lag1"),
  lag2 = run_tree_models_for_response(df_tree_lagged, "y_lag2"),
  lag3 = run_tree_models_for_response(df_tree_lagged, "y_lag3")
)


lag_results_all <- bind_rows(
  lag_runs$lag0$results %>% mutate(lag = "lag0"),
  lag_runs$lag1$results %>% mutate(lag = "lag1"),
  lag_runs$lag2$results %>% mutate(lag = "lag2"),
  lag_runs$lag3$results %>% mutate(lag = "lag3")
) %>%
  mutate(
    lag = factor(lag, levels = c("lag0", "lag1", "lag2", "lag3")),
    model = recode(
      model,
      "1_plate" = "Plate",
      "2_plate_stage" = "Plate + Stage",
      "3_plate_stage_cycle" = "Plate + Stage + Time",
      "4_cycle_only" = "Time",
      "5_plate_cycle" = "Plate + Time"
    )
  ) %>%
  arrange(species, lag, rmse)

lag_results_all

best_by_lag <- lag_results_all %>%
  group_by(species, lag) %>%
  slice_min(order_by = rmse, n = 1, with_ties = FALSE) %>%
  ungroup()

best_by_lag


best_by_lag %>%
  group_by(species) %>%
  mutate(
    species_display = ifelse(row_number() == 1, as.character(species), "")
  ) %>%
  ungroup() %>%
  select(
    species = species_display,
    lag,
    model,
    rmse,
    rsq,
    n
  ) %>%
  gt() %>%
  tab_header(
    title = "Best Regression Tree Model by Time Lag",
    subtitle = "Predictors at week t, response at week t + lag"
  ) %>%
  fmt_number(
    columns = c(rmse, rsq),
    decimals = 3
  ) %>%
  cols_label(
    species = "Species",
    lag = "Lag",
    model = "Best model",
    rmse = "RMSE",
    rsq = "R²",
    n = "N"
  ) %>%
  tab_options(
    table.background.color = "white",
    heading.background.color = "white",
    column_labels.background.color = "white"
  )


lag_results_all %>%
  group_by(species, lag) %>%
  mutate(is_best = rmse == min(rmse)) %>%
  ungroup() %>%
  group_by(species) %>%
  mutate(
    species_display = ifelse(row_number() == 1, as.character(species), "")
  ) %>%
  ungroup() %>%
  select(
    species = species_display,
    lag,
    model,
    rmse,
    rsq,
    is_best
  ) %>%
  gt() %>%
  tab_header(
    title = "Regression Tree Model Performance Across Time Lags",
    subtitle = "Predictors at week t, response at week t + lag"
  ) %>%
  fmt_number(
    columns = c(rmse, rsq),
    decimals = 3
  ) %>%
  cols_label(
    species = "Species",
    lag = "Lag",
    model = "Model",
    rmse = "RMSE",
    rsq = "R²"
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = model,
      rows = is_best
    )
  ) %>%
  cols_hide(is_best) %>%
  tab_options(
    table.background.color = "white",
    heading.background.color = "white",
    column_labels.background.color = "white"
  )
