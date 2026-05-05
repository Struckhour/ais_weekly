source('./AIS_eDNA_data_prep.R')


##########################
#NO MEDIAN REF, NOW WE USE PAIRWISE LEAST SQUARES
##########################


# ----------------------------
# 1. Weekly qPCR curves
# ----------------------------
qpcr_weekly <- dfRawClean %>%
  group_by(region, species, week = isoweek(date)) %>%
  summarise(
    value = mean(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(region, species) %>%
  complete(week = 1:52) %>%
  arrange(region, species, week) %>%
  ungroup()

qpcr_weekly_flagged <- qpcr_weekly %>%
  mutate(
    is_observed = !is.na(value)
  )

interp_summary_species <- qpcr_weekly_flagged %>%
  group_by(species) %>%
  summarise(
    total_weeks = n(),
    observed_weeks = sum(is_observed),
    interpolated_weeks = sum(!is_observed),
    pct_interpolated = 100 * interpolated_weeks / total_weeks,
    .groups = "drop"
  )

interp_summary_species

interp_summary_region <- qpcr_weekly_flagged %>%
  group_by(species, region) %>%
  summarise(
    total_weeks = n(),
    observed_weeks = sum(is_observed),
    interpolated_weeks = sum(!is_observed),
    pct_interpolated = 100 * interpolated_weeks / total_weeks,
    .groups = "drop"
  )
interp_summary_region

# ----------------------------
# 2. Pairwise region cross-correlations within species
# ----------------------------
pairwise_lags <- qpcr_weekly %>%
  group_by(species) %>%
  group_modify(~{

    dat <- .x
    regs <- sort(unique(as.character(dat$region)))

    region_pairs <- combn(regs, 2, simplify = FALSE)

    map_dfr(region_pairs, function(pair) {

      r1 <- pair[1]
      r2 <- pair[2]

      curve1 <- dat %>%
        filter(as.character(region) == r1) %>%
        arrange(week) %>%
        pull(value)

      curve2 <- dat %>%
        filter(as.character(region) == r2) %>%
        arrange(week) %>%
        pull(value)

      prof <- circular_lag_profile(curve1, curve2)

      if (nrow(prof) == 0 || all(is.na(prof$cor))) {
        return(tibble(
          region1 = r1,
          region2 = r2,
          lag = NA_real_,
          cor = NA_real_
        ))
      }

      best_i <- which.max(prof$cor)

      tibble(
        region1 = r1,
        region2 = r2,
        lag = prof$lag[best_i],
        cor = prof$cor[best_i]
      )
    })
  }) %>%
  ungroup()

pairwise_lags

compute_lag_fit_error <- function(lag_df, pos_df, this_species) {

  lag_df %>%
    mutate(species = this_species) %>%
    inner_join(pos_df, by = c("species", "region1" = "region")) %>%
    rename(pos1 = position) %>%
    inner_join(pos_df, by = c("species", "region2" = "region")) %>%
    rename(pos2 = position) %>%
    mutate(
      lag_pred = pos1 - pos2,
      residual = lag - lag_pred,
      weight = pmax(cor, 0)^2
    )
}

lag_fit_summary <- pairwise_lags %>%
  group_by(species) %>%
  group_modify(~{

    this_species <- .y$species[[1]]

    pos_df <- region_positions_by_species %>%
      filter(species == this_species)

    res <- compute_lag_fit_error(.x, pos_df, this_species)

    tibble(
      rmse = sqrt(weighted.mean(res$residual^2, res$weight, na.rm = TRUE)),
      mean_cor = mean(res$cor, na.rm = TRUE),
      n_pairs = nrow(res)
    )
  }) %>%
  ungroup()

library(ggrepel)

ggplot(lag_fit_summary, aes(x = mean_cor, y = rmse, label = species)) +
  geom_point(size = 3) +
  geom_text_repel(
    size = 5,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  theme_classic() +
  labs(
    x = "Mean pairwise correlation",
    y = "Lag model RMSE (weeks)",
    title = "Consistency of regional timing estimates across species"
  )

#####

# ----------------------------
# 3. Solve global regional alignment per species
# ----------------------------
solve_region_positions <- function(lag_df) {

  lag_df <- lag_df %>%
    filter(!is.na(lag), !is.na(cor)) %>%
    mutate(
      weight = pmax(cor, 0)^2
    ) %>%
    filter(weight > 0)

  regs <- sort(unique(c(lag_df$region1, lag_df$region2)))

  if (nrow(lag_df) == 0 || length(regs) < 2) {
    return(tibble(region = regs, position = NA_real_))
  }

  # use last region as temporary reference
  ref_region <- tail(regs, 1)
  fit_regions <- setdiff(regs, ref_region)

  X <- matrix(0, nrow = nrow(lag_df), ncol = length(fit_regions))
  colnames(X) <- fit_regions

  for (i in seq_len(nrow(lag_df))) {
    r1 <- lag_df$region1[i]
    r2 <- lag_df$region2[i]

    if (r1 != ref_region) X[i, r1] <-  1
    if (r2 != ref_region) X[i, r2] <- -1
  }

  y <- lag_df$lag
  w <- lag_df$weight

  fit <- lm(y ~ X - 1, weights = w)

  beta <- coef(fit)
  names(beta) <- fit_regions

  # add reference region back as 0
  positions <- c(beta, setNames(0, ref_region))

  # center so species mean position = 0
  positions <- positions - mean(positions, na.rm = TRUE)

  tibble(
    region = names(positions),
    position = as.numeric(positions)
  )
}

region_positions_by_species <- pairwise_lags %>%
  group_by(species) %>%
  group_modify(~ solve_region_positions(.x)) %>%
  ungroup()

region_positions_by_species


#####
# ----------------------------
# 4. Regional timing summary across species
# ----------------------------
region_position_summary <- region_positions_by_species %>%
  group_by(region) %>%
  summarise(
    mean_position = mean(position, na.rm = TRUE),
    sd_position = sd(position, na.rm = TRUE),
    n_species = sum(!is.na(position)),
    .groups = "drop"
  ) %>%
  arrange(mean_position)

region_position_summary


#############
#with threshold
#############

solve_region_positions_thresh <- function(lag_df, cor_min = 0.6) {

  lag_df <- lag_df %>%
    filter(!is.na(lag), !is.na(cor)) %>%
    filter(cor >= cor_min) %>%
    mutate(
      weight = pmax(cor, 0)^2
    ) %>%
    filter(weight > 0)

  regs <- sort(unique(c(lag_df$region1, lag_df$region2)))

  if (nrow(lag_df) == 0 || length(regs) < 2) {
    return(tibble(region = regs, position = NA_real_))
  }

  ref_region <- tail(regs, 1)
  fit_regions <- setdiff(regs, ref_region)

  X <- matrix(0, nrow = nrow(lag_df), ncol = length(fit_regions))
  colnames(X) <- fit_regions

  for (i in seq_len(nrow(lag_df))) {
    r1 <- lag_df$region1[i]
    r2 <- lag_df$region2[i]

    if (r1 != ref_region) X[i, r1] <-  1
    if (r2 != ref_region) X[i, r2] <- -1
  }

  y <- lag_df$lag
  w <- lag_df$weight

  fit <- lm(y ~ X - 1, weights = w)

  beta <- coef(fit)
  names(beta) <- fit_regions

  positions <- c(beta, setNames(0, ref_region))
  positions <- positions - mean(positions, na.rm = TRUE)

  tibble(
    region = names(positions),
    position = as.numeric(positions),
    n_pairs_used = nrow(lag_df)
  )
}

region_positions_by_species_thresh <- pairwise_lags %>%
  group_by(species) %>%
  group_modify(~ solve_region_positions_thresh(.x, cor_min = 0.6)) %>%
  ungroup()

region_position_summary_thresh <- region_positions_by_species_thresh %>%
  group_by(region) %>%
  summarise(
    mean_position = mean(position, na.rm = TRUE),
    sd_position = sd(position, na.rm = TRUE),
    n_species = sum(!is.na(position)),
    .groups = "drop"
  ) %>%
  arrange(mean_position)

print(region_positions_by_species_thresh, n = Inf, width = Inf)
region_position_summary_thresh



###############
#FILTER FOR JUST STRONG SIGNAL SPECIES
###############

species_strength <- pairwise_lags %>%
  group_by(species) %>%
  summarise(
    mean_cor = mean(cor, na.rm = TRUE),
    n_strong_pairs = sum(cor > 0.6, na.rm = TRUE),
    .groups = "drop"
  )

species_strength

strong_species <- species_strength %>%
  filter(mean_cor > 0.6 & n_strong_pairs >= 5) %>%
  pull(species)

pairwise_lags_strong <- pairwise_lags %>%
  filter(species %in% strong_species)
pairwise_lags_strong


region_positions_by_species_strong <- pairwise_lags_strong %>%
  group_by(species) %>%
  group_modify(~ solve_region_positions(.x)) %>%
  ungroup()

print(region_positions_by_species_strong, n = Inf, width = Inf)

region_position_summary_strong <- region_positions_by_species_strong %>%
  group_by(region) %>%
  summarise(
    mean_position = mean(position, na.rm = TRUE),
    sd_position = sd(position, na.rm = TRUE),
    n_species = n(),
    .groups = "drop"
  ) %>%
  arrange(mean_position)

region_position_summary_strong


################
#UPDATED PART 5 WITH RESTRICTIONS
################
# ----------------------------
# 5. Restrict to species with coherent regional timing
# ----------------------------

valid_species <- lag_fit_summary %>%
  filter(rmse < 3) %>%   # interpretable cutoff: lag model error < ~3 weeks
  pull(species)

valid_species

pairwise_lags_valid <- pairwise_lags %>%
  filter(species %in% valid_species)

region_positions_by_species_valid <- pairwise_lags_valid %>%
  group_by(species) %>%
  group_modify(~ solve_region_positions(.x)) %>%
  ungroup()

print(region_positions_by_species_valid, n = Inf, width = Inf)

region_position_summary_valid <- region_positions_by_species_valid %>%
  group_by(region) %>%
  summarise(
    mean_position = mean(position, na.rm = TRUE),
    sd_position = sd(position, na.rm = TRUE),
    n_species = sum(!is.na(position)),
    .groups = "drop"
  ) %>%
  arrange(mean_position)

region_position_summary_valid

library(dplyr)
library(gt)

region_position_summary_valid %>%
  mutate(
    rank = rank(-mean_position, ties.method = "first")  # early = 1
  ) %>%
  arrange(rank) %>%
  rename(
    `Mean adjusted position (weeks)` = mean_position,
    `SD (weeks)` = sd_position,
    `Number of species` = n_species,
    `Region` = region,
    `Rank (early → late)` = rank
  ) %>%
  gt() %>%
  fmt_number(
    columns = c(`Mean adjusted position (weeks)`, `SD (weeks)`),
    decimals = 2
  ) %>%
  cols_align(
    align = "center",
    -Region
  ) %>%
  tab_header(
    title = "Regional ordering of seasonal timing"
  )


#PERMUTATION TEST#####################################
obs_stat <- region_position_summary_valid %>%
  summarise(stat = var(mean_position, na.rm = TRUE)) %>%
  pull(stat)

permute_once <- function(df) {
  df %>%
    group_by(species) %>%
    mutate(region = sample(region)) %>%   # shuffle labels
    ungroup()
}

compute_stat <- function(df) {

  pairwise_lags_perm <- df %>%
    group_by(species) %>%
    group_modify(~{

      dat <- .x
      regs <- sort(unique(as.character(dat$region)))
      region_pairs <- combn(regs, 2, simplify = FALSE)

      purrr::map_dfr(region_pairs, function(pair) {

        r1 <- pair[1]
        r2 <- pair[2]

        curve1 <- dat %>%
          filter(region == r1) %>%
          arrange(week) %>%
          pull(value)

        curve2 <- dat %>%
          filter(region == r2) %>%
          arrange(week) %>%
          pull(value)

        prof <- circular_lag_profile(curve1, curve2)

        if (nrow(prof) == 0 || all(is.na(prof$cor))) {
          return(tibble(region1 = r1, region2 = r2, lag = NA, cor = NA))
        }

        best_i <- which.max(prof$cor)

        tibble(
          region1 = r1,
          region2 = r2,
          lag = prof$lag[best_i],
          cor = prof$cor[best_i]
        )
      })
    }) %>%
    ungroup()

  region_positions_perm <- pairwise_lags_perm %>%
    group_by(species) %>%
    group_modify(~ solve_region_positions(.x)) %>%
    ungroup()

  region_summary_perm <- region_positions_perm %>%
    group_by(region) %>%
    summarise(mean_position = mean(position, na.rm = TRUE), .groups = "drop")

  var(region_summary_perm$mean_position, na.rm = TRUE)
}

set.seed(123)

n_perm <- 1000
perm_stats <- numeric(n_perm)

start_time <- Sys.time()

for (i in seq_len(n_perm)) {

  df_perm <- permute_once(
    qpcr_weekly %>%
      filter(species %in% valid_species)
  )

  perm_stats[i] <- compute_stat(df_perm)

  if (i %% 50 == 0) {
    elapsed <- as.numeric(Sys.time() - start_time, units = "mins")
    rate <- elapsed / i
    eta <- rate * (n_perm - i)

    cat(sprintf(
      "Completed %d / %d | elapsed: %.1f mins | ETA: %.1f mins\n",
      i, n_perm, elapsed, eta
    ))
  }
}

p_value <- mean(perm_stats >= obs_stat)
p_value





























###############
#PLOTS
###############

pairwise_lags %>%
  mutate(
    pair = paste(region1, region2, sep = " vs "),
    species = factor(species, levels = rev(species_order))
  ) %>%
  ggplot(aes(x = pair, y = species, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(
    aes(
      label = paste0("lag=", lag, "\nr=", round(cor, 2)),
      color = cor > 0.65   # threshold tweakable
    ),
    size = 3
  ) +
  scale_color_manual(
    values = c("black", "white"),
    guide = "none"
  ) +
  scale_fill_gradient(low = "grey90", high = "black") +
  theme_classic() +
  labs(
    x = "Region pair",
    y = "Species",
    fill = "Max correlation",
    title = "Pairwise circular cross-correlation between regional qPCR curves"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

library(dplyr)
library(ggplot2)
library(patchwork)

# ----------------------------
# Compare pairwise-optimal lags with least-squares reconciled lags
# ----------------------------

pairwise_lags_fit <- pairwise_lags %>%
  inner_join(
    region_positions_by_species %>%
      select(species, region, position) %>%
      rename(pos1 = position),
    by = c("species", "region1" = "region")
  ) %>%
  inner_join(
    region_positions_by_species %>%
      select(species, region, position) %>%
      rename(pos2 = position),
    by = c("species", "region2" = "region")
  ) %>%
  mutate(
    lag_pairwise = lag,
    lag_fitted = pos1 - pos2,

    # circular disagreement between pairwise and fitted lag
    lag_error = abs(((lag_pairwise - lag_fitted + 26) %% 52) - 26),

    # 1 = perfect agreement, 0 = maximum circular disagreement
    fit_score = 1 - lag_error / 26,

    pair = paste(region1, region2, sep = " vs "),
    species = factor(species, levels = rev(species_order))
  )

pairwise_lags_fit <- pairwise_lags_fit %>%
  mutate(
    # rescale error based on observed range
    error_min = min(lag_error, na.rm = TRUE),
    error_max = max(lag_error, na.rm = TRUE),

    fit_score_scaled = 1 - (lag_error - error_min) / (error_max - error_min)
  )

p_pairwise <- pairwise_lags_fit %>%
  ggplot(aes(x = pair, y = species, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(
    aes(
      label = paste0(
        "lag=", round(lag_pairwise, 1),
        "\nr=", round(cor, 2)
      ),
      color = cor > 0.65
    ),
    size = 3
  ) +
  scale_fill_gradient(
    low = "grey90",
    high = "black"
  ) +
  scale_color_manual(
    values = c("black", "white"),
    guide = "none"
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = "Species",
    fill = "Max r",
    title = "Pairwise optimal lags"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

p_fitted <- pairwise_lags_fit %>%
  ggplot(aes(x = pair, y = species, fill = fit_score_scaled)) +
  geom_tile(color = "white") +
  geom_text(
    aes(
      label = paste0(
        "lag=", round(lag_fitted, 1),
        "\nfit=", round(fit_score_scaled, 2)
      ),
      color = fit_score < 0.75
    ),
    size = 3
  ) +
  scale_fill_gradient(
    low = "grey90",
    high = "black",
    limits = c(0, 1)
  ) +
  scale_color_manual(
    values = c("white", "black"),
    guide = "none"
  ) +
  theme_classic() +
  labs(
    x = "Region pair",
    y = "Species",
    fill = "Fit score",
    title = "Least-squares reconciled lags"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_pairwise / p_fitted
#
# plot_alignment_for_species <- function(target_species) {
#
#   shifts <- region_positions_by_species %>%
#     filter(species == target_species) %>%
#     select(region, position)
#
#   plot_dat <- qpcr_weekly %>%
#     filter(species == target_species) %>%
#     left_join(shifts, by = "region") %>%
#     group_by(region) %>%
#     mutate(
#       value_norm = value / max(value, na.rm = TRUE),
#       week_shifted = ((week + round(position) - 1) %% 52) + 1
#     ) %>%
#     ungroup()
#
#   p_before <- ggplot(plot_dat, aes(x = week, y = value_norm, color = region)) +
#     geom_line(linewidth = 1, alpha = 0.9) +
#     scale_color_manual(values = hybrid_color) +
#     theme_classic() +
#     labs(
#       title = paste(target_species, "- before alignment"),
#       x = "Week",
#       y = "Normalized mean log concentration"
#     )
#
#   p_after <- ggplot(plot_dat, aes(x = week_shifted, y = value_norm, color = region)) +
#     geom_line(linewidth = 1, alpha = 0.9) +
#     scale_color_manual(values = hybrid_color) +
#     theme_classic() +
#     labs(
#       title = paste(target_species, "- after alignment"),
#       x = "Shifted week",
#       y = "Normalized mean log concentration"
#     ) +
#     theme(legend.position = "none")
#
#   p_before / p_after
# }

plot_alignment_for_species <- function(target_species) {

  shifts <- region_positions_by_species %>%
    filter(species == target_species) %>%
    select(region, position)

  plot_dat <- qpcr_weekly %>%
    filter(species == target_species) %>%
    left_join(shifts, by = "region") %>%
    group_by(region) %>%
    mutate(
      is_obs = !is.na(value),
      value_filled = fill_circular_series(value),
      value_filled_norm = value_filled / max(value_filled, na.rm = TRUE),
      value_norm = value / max(value, na.rm = TRUE),
      week_shifted = ((week + round(position) - 1) %% 52) + 1
    ) %>%
    ungroup()

  # BEFORE
  p_before <- ggplot(plot_dat, aes(x = week, color = region)) +

    # thin continuous line (interpolated curve)
    geom_line(
      aes(y = value_filled_norm),
      linewidth = 0.7,
      alpha = 0.9
    ) +

    # dots = real observations only
    geom_point(
      data = ~ dplyr::filter(.x, is_obs),
      aes(y = value_norm),
      size = 1.8,
      alpha = 0.9
    ) +

    scale_color_manual(values = hybrid_color) +
    theme_classic() +
    labs(
      title = paste(target_species, "- before alignment"),
      x = "Week",
      y = "Normalized mean log concentration"
    )

  # AFTER
  p_after <- ggplot(plot_dat, aes(x = week_shifted, color = region)) +

    geom_line(
      aes(y = value_filled_norm),
      linewidth = 0.7,
      alpha = 0.9
    ) +

    geom_point(
      data = ~ dplyr::filter(.x, is_obs),
      aes(y = value_norm),
      size = 1.8,
      alpha = 0.9
    ) +

    scale_color_manual(values = hybrid_color) +
    theme_classic() +
    labs(
      title = paste(target_species, "- after alignment"),
      x = "Shifted week",
      y = "Normalized mean log concentration"
    ) +
    theme(legend.position = "none")

  p_before / p_after
}



plot_alignment_for_species("Membranipora membranacea")
plot_alignment_for_species("Botrylloides violaceus")
plot_alignment_for_species("Ciona intestinalis")
plot_alignment_for_species("Didemnum vexillum")
plot_alignment_for_species("Carcinus maenas")

plot_alignment_all_species <- function() {

  plot_dat <- qpcr_weekly %>%
    left_join(
      region_positions_by_species %>%
        select(species, region, position),
      by = c("species", "region")
    ) %>%
    group_by(species, region) %>%
    mutate(
      is_obs = !is.na(value),
      value_filled = fill_circular_series(value),

      value_min = min(value_filled, na.rm = TRUE),
      value_max = max(value_filled, na.rm = TRUE),

      value_filled_norm = dplyr::if_else(
        value_max > value_min,
        (value_filled - value_min) / (value_max - value_min),
        NA_real_
      ),

      value_norm = dplyr::if_else(
        value_max > value_min,
        (value - value_min) / (value_max - value_min),
        NA_real_
      ),

      week_shifted = ((week + round(position) - 1) %% 52) + 1
    ) %>%
    ungroup() %>%
    mutate(
      species = factor(species, levels = species_order)
    )

  plot_long_line <- plot_dat %>%
    select(species, region, week, week_shifted, value_filled_norm) %>%
    tidyr::pivot_longer(
      cols = c(week, week_shifted),
      names_to = "alignment",
      values_to = "plot_week"
    ) %>%
    mutate(
      alignment = recode(
        alignment,
        week = "Before alignment",
        week_shifted = "After alignment"
      ),
      alignment = factor(alignment, levels = c("Before alignment", "After alignment"))
    )

  plot_long_points <- plot_dat %>%
    filter(is_obs) %>%
    select(species, region, week, week_shifted, value_norm) %>%
    tidyr::pivot_longer(
      cols = c(week, week_shifted),
      names_to = "alignment",
      values_to = "plot_week"
    ) %>%
    mutate(
      alignment = recode(
        alignment,
        week = "Before alignment",
        week_shifted = "After alignment"
      ),
      alignment = factor(alignment, levels = c("Before alignment", "After alignment"))
    )

  ggplot(plot_long_line, aes(x = plot_week, color = region)) +
    geom_line(
      aes(y = value_filled_norm, group = region),
      linewidth = 0.7,
      alpha = 0.9
    ) +
    geom_point(
      data = plot_long_points,
      aes(y = value_norm),
      size = 1.8,
      alpha = 0.9
    ) +
    facet_grid(alignment ~ species) +
    scale_color_manual(values = hybrid_color) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
    ) +
    labs(
      x = "Week",
      y = "Normalized mean log concentration",
      color = "Region"
    )
}

plot_alignment_all_species()
