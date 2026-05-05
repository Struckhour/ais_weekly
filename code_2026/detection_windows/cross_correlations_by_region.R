source('./AIS_eDNA_data_prep.R')

qpcr_weekly <- dfRawClean %>%
  group_by(region, species, week = isoweek(date)) %>%
  summarise(
    value = mean(logConc, na.rm = TRUE),
    .groups = "drop"
  )

ref_curve <- qpcr_weekly %>%
  group_by(species, week) %>%
  summarise(
    value = median(value, na.rm = TRUE),
    .groups = "drop"
  )


qpcr_weekly <- qpcr_weekly %>%
  group_by(region, species) %>%
  complete(week = 1:52) %>%
  arrange(region, species, week) %>%
  ungroup()

ref_curve <- ref_curve %>%
  group_by(species) %>%
  complete(week = 1:52) %>%
  arrange(species, week) %>%
  ungroup()

qpcr_weekly %>% count(region, species)
ref_curve %>% count(species)

lag_vs_ref <- qpcr_weekly %>%
  group_by(region, species) %>%
  group_modify(~{

    sp <- .y$species[[1]]

    ref <- ref_curve %>%
      filter(species == sp) %>%
      arrange(week)

    dat <- .x %>%
      arrange(week)

    prof <- circular_lag_profile(ref$value, dat$value)

    if (nrow(prof) == 0 || all(is.na(prof$cor))) {
      return(tibble(
        lag = NA_real_,
        cor = NA_real_
      ))
    }

    best_i <- which.max(prof$cor)

    tibble(
      lag = prof$lag[best_i],
      cor = prof$cor[best_i]
    )
  }) %>%
  ungroup()


lag_summary <- lag_vs_ref %>%
  group_by(region) %>%
  summarise(
    mean_lag = mean(lag, na.rm = TRUE),
    sd_lag = sd(lag, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_lag)


lag_filtered <- lag_vs_ref %>%
  filter(cor > 0.7)
lag_summary_filtered <- lag_filtered %>%
  group_by(region) %>%
  summarise(
    mean_lag = mean(lag, na.rm = TRUE),
    sd_lag = sd(lag, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(mean_lag)



lag_summary_weighted <- lag_vs_ref %>%
  filter(!is.na(cor), cor > 0) %>%
  group_by(region) %>%
  summarise(
    mean_lag = weighted.mean(lag, w = cor, na.rm = TRUE),
    sd_lag = sqrt(weighted.mean((lag - mean_lag)^2, w = cor, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  arrange(mean_lag)



###################
#PERMUTATION TEST
###################

obs_stat <- lag_summary_weighted %>%
  summarise(var = var(mean_lag)) %>%
  pull(var)

set.seed(1)

n_perm <- 2000

perm_stats <- replicate(n_perm, {

  permuted <- lag_vs_ref %>%
    group_by(species) %>%
    mutate(region_perm = sample(region)) %>%
    ungroup()

  perm_summary <- permuted %>%
    group_by(region_perm) %>%
    summarise(
      mean_lag = weighted.mean(lag, w = cor, na.rm = TRUE),
      .groups = "drop"
    )

  var(perm_summary$mean_lag, na.rm = TRUE)
})

p_value <- mean(perm_stats >= obs_stat)
p_value

ggplot() +
  geom_line(data = qpcr_weekly,
            aes(x = week, y = value, color = region),
            alpha = 0.6) +
  geom_line(data = ref_curve,
            aes(x = week, y = value),
            color = "black",
            linewidth = 1.2) +
  facet_wrap(~ species, scales = "free_y") +
  theme_classic()










##########################
#NO MEDIAN REF, NO WE USE PAIRWISE LEAST SQUARES
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
