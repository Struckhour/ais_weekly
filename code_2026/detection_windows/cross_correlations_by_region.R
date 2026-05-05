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


###############
#PLOTS
###############

pairwise_lags %>%
  mutate(
    pair = paste(region1, region2, sep = " vs "),
    species = factor(species, levels = species_order)
  ) %>%
  ggplot(aes(x = pair, y = species, fill = cor)) +
  geom_tile(color = "white") +
  geom_text(aes(label = paste0("lag=", lag, "\nr=", round(cor, 2))), size = 3) +
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
