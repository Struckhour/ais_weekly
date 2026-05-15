source('./AIS_eDNA_data_prep.R')


###################
#NEW CODE -- THIS SECTION (TILL LINE #760) SHOULD BE WORKABLE
#FOR ALL COMBINATIONS OF SPECIES AND REGION
#EVERYTHING BELOW THIS SECTION (#760) SHOULD NOW BE OBSOLETE
###################
library(dplyr)
library(tidyr)
library(purrr)
library(tibble)
library(ggplot2)
library(gt)

###############################
# USER CONTROLS
###############################

# species to include
species_include <- species_order
# species_include <- c(
#   "Membranipora membranacea",
#   "Botrylloides violaceus",
#   "Didemnum vexillum",
#   "Ciona intestinalis",
#   "Carcinus maenas"
# )

# optional species to exclude
species_exclude <- c()
# species_exclude <- c("Didemnum vexillum")

# regions to include in the pipeline
regions_include <- c("MAG", "PEI", "HAL", "BOF", "GOM")
# regions_include <- c("BOF", "GOM", "HAL")

# exact region x species combinations to remove
region_species_exclude <- tibble::tribble(
  ~region, ~species,
  # "MAG", "Carcinus maenas",

)

species_keep <- setdiff(species_include, species_exclude)

###############################
# 1. FILTER RAW DATA FIRST
###############################

df_for_pipeline <- dfRawClean %>%
  filter(
    species %in% species_keep,
    region %in% regions_include
  ) %>%
  anti_join(region_species_exclude, by = c("region", "species"))

###############################
# 2. WEEKLY QPCR CURVES
###############################

qpcr_weekly_selected <- df_for_pipeline %>%
  group_by(region, species, week = isoweek(date)) %>%
  summarise(
    value = mean(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(region, species) %>%
  complete(week = 1:52) %>%
  arrange(region, species, week) %>%
  ungroup()

###############################
# 3. PAIRWISE REGION CROSS-CORRELATIONS
###############################

pairwise_lags_selected <- qpcr_weekly_selected %>%
  group_by(species) %>%
  group_modify(~{

    dat <- .x
    regs <- sort(unique(as.character(dat$region)))

    if (length(regs) < 2) {
      return(tibble(
        region1 = character(),
        region2 = character(),
        lag = numeric(),
        cor = numeric()
      ))
    }

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

print(pairwise_lags_selected, n = Inf)

###############################
# 4. SOLVE RECONCILED REGIONAL POSITIONS
###############################

solve_region_positions <- function(lag_df) {

  lag_df <- lag_df %>%
    filter(!is.na(lag), !is.na(cor)) %>%
    mutate(
      weight = pmax(cor, 0)^2
    ) %>%
    filter(weight > 0)

  regs <- sort(unique(c(lag_df$region1, lag_df$region2)))

  if (nrow(lag_df) == 0 || length(regs) < 2) {
    return(tibble(
      region = regs,
      position = NA_real_
    ))
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
    position = as.numeric(positions)
  )
}

region_positions_selected <- pairwise_lags_selected %>%
  group_by(species) %>%
  group_modify(~ solve_region_positions(.x)) %>%
  ungroup()

print(region_positions_selected, n = Inf)

library(ggrepel)


species_fit_selected <- pairwise_lags_selected %>%
  left_join(
    region_positions_selected %>%
      rename(region1 = region, pos1 = position),
    by = c("species", "region1")
  ) %>%
  left_join(
    region_positions_selected %>%
      rename(region2 = region, pos2 = position),
    by = c("species", "region2")
  ) %>%
  mutate(
    lag_pred = pos1 - pos2,
    residual = lag - lag_pred
  ) %>%
  group_by(species) %>%
  summarise(
    mean_cor = mean(cor, na.rm = TRUE),
    rmse = sqrt(mean(residual^2, na.rm = TRUE)),
    n_pairs = sum(!is.na(residual)),
    .groups = "drop"
  ) %>%
  mutate(
    species = factor(species, levels = species_order)
  )

manual_point <- tibble::tibble(
  species = "Carcinus maenas (MAG removed)",
  mean_cor = 0.462,
  rmse = .892
)

ggplot(species_fit_selected, aes(x = mean_cor, y = rmse, label = species)) +
  geom_point(size = 3) +
  geom_text_repel(
    size = 5,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf
  ) +
  geom_point(
    data = manual_point,
    size = 3
  ) +
  geom_text_repel(
    data = manual_point,
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
#
# ggplot(species_fit_selected, aes(x = mean_cor, y = rmse, label = species)) +
#   geom_point(size = 3) +
#   geom_text_repel(
#     size = 5,
#     box.padding = 0.4,
#     point.padding = 0.3,
#     segment.color = "grey50",
#     max.overlaps = Inf
#   ) +
#   theme_classic() +
#   labs(
#     x = "Mean pairwise correlation",
#     y = "Lag model RMSE (weeks)",
#     title = "Consistency of regional timing estimates across species"
#   )

library(patchwork)

pairwise_lags_fit_selected <- pairwise_lags_selected %>%
  inner_join(
    region_positions_selected %>%
      select(species, region, position) %>%
      rename(pos1 = position),
    by = c("species", "region1" = "region")
  ) %>%
  inner_join(
    region_positions_selected %>%
      select(species, region, position) %>%
      rename(pos2 = position),
    by = c("species", "region2" = "region")
  ) %>%
  mutate(
    lag_pairwise = lag,
    lag_fitted = pos1 - pos2,

    # signed circular difference: pairwise - reconciled
    lag_error = ((lag_pairwise - lag_fitted + 26) %% 52) - 26,

    # absolute difference in weeks
    abs_error = abs(lag_error),

    pair = paste(region1, region2, sep = " vs "),
    species = factor(species, levels = rev(species_order))
  )

p_pairwise <- pairwise_lags_fit_selected %>%
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

max_abs_error <- max(pairwise_lags_fit_selected$abs_error, na.rm = TRUE)

text_cutoff <- max_abs_error / 2

p_fitted <- pairwise_lags_fit_selected %>%
  ggplot(aes(x = pair, y = species, fill = abs_error)) +
  geom_tile(color = "white") +
  geom_text(
    aes(
      label = paste0(
        "lag=", round(lag_fitted, 1),
        "\nΔ=", round(abs_error, 1)
      ),
      color = abs_error <= text_cutoff
    ),
    size = 3
  ) +
  scale_fill_gradient(
    low = "black",
    high = "grey90",
    limits = c(0, max_abs_error)
  ) +
  scale_color_manual(
    values = c("black", "white"),
    guide = "none"
  ) +
  theme_classic() +
  labs(
    x = "Region pair",
    y = "Species",
    fill = "Δ lag (weeks)",
    title = "Difference between pairwise and reconciled lags"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

p_pairwise / p_fitted


#############################
#COMBINED PLOT
#############################
library(patchwork)

p_consistency <- ggplot(species_fit_selected, aes(x = mean_cor, y = rmse, label = species)) +
  geom_point(size = 3) +
  geom_text_repel(
    size = 5,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf,
    nudge_x = ifelse(species_fit_selected$species == "Carcinus maenas", 0.02, 0),
    nudge_y = ifelse(species_fit_selected$species == "Carcinus maenas", -0.8, 0)
  ) +
  geom_point(
    data = manual_point,
    size = 3
  ) +
  geom_text_repel(
    data = manual_point,
    size = 5,
    box.padding = 0.4,
    point.padding = 0.3,
    segment.color = "grey50",
    max.overlaps = Inf,
    nudge_x = 0.02,
    nudge_y = 0.8
  ) +
  theme_classic() +
  labs(
    x = "Mean pairwise correlation",
    y = "Lag model RMSE (weeks)",
    title = "(a) Consistency of regional timing estimates across species"
  )

p_pairwise <- p_pairwise +
  labs(title = "(b) Pairwise optimal lags")

p_fitted <- p_fitted +
  labs(title = "(c) Difference between pairwise and reconciled lags")

combined_plot <- p_consistency /
  (p_pairwise + theme(legend.position = "right")) /
  (p_fitted + theme(legend.position = "right")) +
  plot_layout(
    heights = c(1.8, 1.2, 1.2),
    guides = "keep"
  ) &
  theme(
    plot.title = element_text(face = "bold", hjust = 0)
  )

combined_plot

###############################
# 5. REGIONAL SUMMARY + RANKS
###############################

region_position_summary_selected <- region_positions_selected %>%
  group_by(region) %>%
  summarise(
    mean_position = mean(position, na.rm = TRUE),
    sd_position = sd(position, na.rm = TRUE),
    n_species = n_distinct(species),
    .groups = "drop"
  ) %>%
  mutate(
    rank = rank(-mean_position, ties.method = "first")
  ) %>%
  arrange(rank)

print(region_position_summary_selected)

region_position_summary_selected %>%
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

###############################
# 6. OPTIONAL LAG TIMELINE PLOT
###############################

region_offsets <- c(
  MAG = -0.1,
  PEI = -0.05,
  HAL =  0.0,
  BOF =  0.05,
  GOM =  0.1
)

lag_timeline_df <- region_positions_selected %>%
  filter(!is.na(position)) %>%
  mutate(
    position_plot = -position,
    species = factor(species, levels = rev(species_order)),
    region = factor(region, levels = names(region_offsets))
  ) %>%
  filter(!is.na(species), !is.na(region)) %>%
  mutate(
    species = droplevels(species),
    species_num = as.numeric(species),
    y_pos = species_num + region_offsets[as.character(region)]
  )

species_breaks <- seq_along(levels(lag_timeline_df$species))
between_species_lines <- species_breaks[-length(species_breaks)] + 0.5

ggplot(lag_timeline_df, aes(x = position_plot, y = y_pos, color = region, shape = region)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(
    yintercept = between_species_lines,
    color = "grey25",
    linewidth = 0.6
  ) +
  geom_point(size = 3, alpha = 1) +
  scale_y_continuous(
    breaks = species_breaks,
    labels = levels(lag_timeline_df$species)
  ) +
  scale_color_manual(values = hybrid_color, drop = FALSE) +
  scale_shape_manual(values = c(
    MAG = 16,
    PEI = 17,
    HAL = 15,
    BOF = 18,
    GOM = 8
  ), drop = FALSE) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  labs(
    x = "Relative seasonal timing (weeks)",
    y = "Species",
    color = "Region",
    shape = "Region"
  )

plot_alignment_all_species <- function() {

  library(patchwork)

  region_levels <- regions_include

  region_shapes <- c(
    MAG = 16,
    PEI = 17,
    HAL = 15,
    BOF = 18,
    GOM = 8
  )

  legend_guides <- guides(
    color = guide_legend(
      override.aes = list(
        linetype = 1,
        shape = region_shapes[region_levels],
        linewidth = 0.8,
        size = 3
      )
    ),
    shape = "none"
  )

  plot_dat <- qpcr_weekly_selected %>%
    left_join(
      region_positions_selected %>%
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

      week_shifted = ((week + round(position) - 1) %% 52) + 1
    ) %>%
    ungroup() %>%
    mutate(
      species = factor(species, levels = species_order),
      region = factor(region, levels = rev(region_levels))
    ) %>%
    filter(!is.na(species), !is.na(region), !is.na(position))

  plot_long_line <- plot_dat %>%
    select(species, region, week, week_shifted, value_filled_norm) %>%
    tidyr::pivot_longer(
      cols = c(week, week_shifted),
      names_to = "alignment",
      values_to = "plot_week"
    ) %>%
    mutate(
      alignment = dplyr::recode(
        alignment,
        week = "Before alignment",
        week_shifted = "After alignment"
      ),
      alignment = factor(
        alignment,
        levels = c("Before alignment", "After alignment")
      )
    )

  make_curve_plot <- function(alignment_name, y_label, show_strips = TRUE) {

    p <- plot_long_line %>%
      filter(alignment == alignment_name) %>%
      ggplot(aes(x = plot_week, color = region, shape = region)) +
      geom_smooth(
        aes(
          y = value_filled_norm,
          group = region
        ),
        method = "loess",
        formula = y ~ x,
        span = 0.3,
        se = FALSE,
        linewidth = 0.8,
        alpha = 0.9
      ) +
      facet_grid(. ~ species) +
      scale_color_manual(
        values = hybrid_color[region_levels],
        drop = FALSE
      ) +
      scale_shape_manual(
        values = region_shapes[region_levels],
        drop = FALSE
      ) +
      scale_x_continuous(
        breaks = seq(0, 52, by = 13),
        limits = c(1, 52)
      ) +
      legend_guides +
      theme_classic() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
      ) +
      labs(
        x = NULL,
        y = y_label,
        color = "Region",
        shape = "Region"
      )

    if (!show_strips) {
      p <- p +
        theme(
          strip.text = element_blank(),
          strip.background = element_blank()
        )
    }

    p
  }

  lag_panel_df <- region_positions_selected %>%
    filter(!is.na(position)) %>%
    mutate(
      position_plot = -position,
      species = factor(species, levels = species_order),
      region = factor(region, levels = rev(region_levels))
    ) %>%
    filter(!is.na(species), !is.na(region))

  p_before <- make_curve_plot(
    alignment_name = "Before alignment",
    y_label = "Before\nalignment",
    show_strips = TRUE
  )

  p_lags <- ggplot(
    lag_panel_df,
    aes(x = position_plot, y = region, color = region, shape = region)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(size = 3, alpha = 0.95) +
    facet_grid(. ~ species) +
    scale_color_manual(
      values = hybrid_color[region_levels],
      drop = FALSE
    ) +
    scale_shape_manual(
      values = region_shapes[region_levels],
      drop = FALSE
    ) +
    legend_guides +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      strip.text = element_blank(),
      strip.background = element_blank()
    ) +
    labs(
      x = NULL,
      y = "Regional\nlag",
      color = "Region",
      shape = "Region"
    )

  p_after <- make_curve_plot(
    alignment_name = "After alignment",
    y_label = "After\nalignment",
    show_strips = FALSE
  ) +
    labs(x = "Week")

  (p_before / p_lags / p_after) +
    plot_layout(
      heights = c(1, 0.55, 1),
      guides = "collect"
    ) &
    theme(
      legend.position = "bottom"
    )
}

plot_alignment_all_species()


##############################
#PERMUTATION
###############
###############################
# 7. PERMUTATION TEST
# uses the same selected weekly workflow
###############################

obs_stat <- region_positions_selected %>%
  group_by(region) %>%
  summarise(
    mean_position = mean(position, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  summarise(
    stat = var(mean_position, na.rm = TRUE)
  ) %>%
  pull(stat)

obs_stat

permute_once <- function(df) {
  df %>%
    group_by(species) %>%
    mutate(region = sample(region)) %>%
    ungroup()
}

compute_perm_stat <- function(df_perm) {

  pairwise_lags_perm <- df_perm %>%
    group_by(species) %>%
    group_modify(~{

      dat <- .x
      regs <- sort(unique(as.character(dat$region)))

      if (length(regs) < 2) {
        return(tibble(
          region1 = character(),
          region2 = character(),
          lag = numeric(),
          cor = numeric()
        ))
      }

      region_pairs <- combn(regs, 2, simplify = FALSE)

      purrr::map_dfr(region_pairs, function(pair) {

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

  region_positions_perm <- pairwise_lags_perm %>%
    group_by(species) %>%
    group_modify(~ solve_region_positions(.x)) %>%
    ungroup()

  region_summary_perm <- region_positions_perm %>%
    group_by(region) %>%
    summarise(
      mean_position = mean(position, na.rm = TRUE),
      .groups = "drop"
    )

  var(region_summary_perm$mean_position, na.rm = TRUE)
}

set.seed(123)

n_perm <- 1000
perm_stats <- numeric(n_perm)

start_time <- Sys.time()

for (i in seq_len(n_perm)) {

  df_perm <- permute_once(qpcr_weekly_selected)

  perm_stats[i] <- compute_perm_stat(df_perm)

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

p_value <- mean(perm_stats >= obs_stat, na.rm = TRUE)

p_value



















































































































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

#####
#Figure showing ranks

region_offsets <- c(
  MAG = -0.1,
  PEI = -0.05,
  HAL =  0.0,
  BOF =  0.05,
  GOM =  0.1
)

lag_timeline_df <- region_positions_by_species_valid %>%
  mutate(
    position_plot = -position,
    species = factor(species, levels = rev(species_order)),
    region = factor(region, levels = names(region_offsets)),
    species_num = as.numeric(species),
    y_pos = species_num + region_offsets[as.character(region)]
  )

lag_timeline_df <- region_positions_by_species_valid %>%
  filter(!is.na(position)) %>%
  mutate(
    position_plot = -position,
    species = factor(species, levels = rev(species_order)),
    region = factor(region, levels = c("MAG", "PEI", "HAL", "BOF", "GOM"))
  ) %>%
  filter(!is.na(species)) %>%
  mutate(
    species = droplevels(species),
    species_num = as.numeric(species),
    y_pos = species_num + region_offsets[as.character(region)]
  )

species_breaks <- seq_along(levels(lag_timeline_df$species))
between_species_lines <- species_breaks[-length(species_breaks)] + 0.5


###############
#JUST FULL SETS
###############

complete_species <- region_positions_by_species_valid %>%
  group_by(species) %>%
  filter(
    all(c("MAG", "PEI", "HAL", "BOF", "GOM") %in% region),
    all(!is.na(position))
  ) %>%
  ungroup()

region_position_summary_balanced <- complete_species %>%
  group_by(region) %>%
  summarise(
    mean_position = mean(position, na.rm = TRUE),
    sd_position = sd(position, na.rm = TRUE),
    n_species = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_position)) %>%
  mutate(
    rank = row_number()
  )

region_position_summary_balanced %>%
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


##################
#select species to keep FOR SUMMARY TABLE
##################
regions_keep <- c("MAG", "PEI", "HAL", "BOF", "GOM")


# -------------------------------
# 1. choose species to include/exclude
# -------------------------------
species_include <- species_order
# species_include <- c(
#   "Membranipora membranacea",
#   "Botrylloides violaceus",
#   "Didemnum vexillum"
# )

species_exclude <- c()
# species_exclude <- c("Didemnum vexillum")

# -------------------------------
# 2. choose specific region x species combos to exclude
# -------------------------------
region_species_exclude <- tibble::tribble(
  ~region, ~species,
  "MAG", "Carcinus maenas"
  # ,"GOM", "Membranipora membranacea"
)

# -------------------------------
# 3. define species to summarize
# -------------------------------
species_for_summary <- species_include %>%
  setdiff(species_exclude)

# -------------------------------
# 4. start from the FULL unfiltered table
#    replace this object name with your real upstream object
# -------------------------------
positions_full <- region_positions_by_species

# -------------------------------
# 5. apply only the exclusions you actually want
# -------------------------------
positions_filtered <- positions_full %>%
  filter(species %in% species_for_summary) %>%
  anti_join(region_species_exclude, by = c("region", "species"))

# -------------------------------
# 6. keep only species that still have all required regions
#    and non-missing positions for those remaining rows
# -------------------------------
complete_species <- positions_filtered %>%
  group_by(species) %>%
  filter(
    all(regions_keep %in% region),
    all(!is.na(position[match(regions_keep, region)]))
  ) %>%
  ungroup()

region_position_summary_balanced <- complete_species %>%
  group_by(region) %>%
  summarise(
    mean_position = mean(position, na.rm = TRUE),
    sd_position = sd(position, na.rm = TRUE),
    n_species = n(),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_position)) %>%
  mutate(rank = row_number())

region_position_summary_balanced %>%
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

#PERMUTATION TEST - BASIC#####################################
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

######################
#permutation with species selected
######################
regions_keep <- c("MAG", "PEI", "HAL", "BOF", "GOM")



# ----------------------------
# Choose species to include/exclude
# ----------------------------

species_include <- species_order
# or manually:
# species_include <- c(
#   "Membranipora membranacea",
#   "Botrylloides violaceus",
#   "Didemnum vexillum"
#
# )

species_exclude <- c()
# example:
# species_exclude <- c("Didemnum vexillum")



species_for_test <- valid_species %>%
  intersect(species_include) %>%
  setdiff(species_exclude)
species_for_test


obs_stat <- region_positions_by_species_valid %>%
  filter(species %in% species_for_test) %>%
  group_by(region) %>%
  summarise(
    mean_position = mean(position, na.rm = TRUE),
    .groups = "drop"
  ) %>%
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
    filter(species %in% species_for_test) %>%
    group_by(region) %>%
    summarise(
      mean_position = mean(position, na.rm = TRUE),
      .groups = "drop"
    )

  var(region_summary_perm$mean_position, na.rm = TRUE)
}

set.seed(123)

n_perm <- 1000
perm_stats <- numeric(n_perm)

start_time <- Sys.time()

for (i in seq_len(n_perm)) {

  df_perm <- permute_once(
    qpcr_weekly %>%
      filter(species %in% species_for_test)
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

  # ggplot(plot_long_line, aes(x = plot_week, color = region)) +
  #   geom_line(
  #     aes(y = value_filled_norm, group = region),
  #     linewidth = 0.7,
  #     alpha = 0.9
  #   ) +
  #   geom_point(
  #     data = plot_long_points,
  #     aes(y = value_norm),
  #     size = 1.8,
  #     alpha = 0.9
  #   ) +
  #   facet_grid(alignment ~ species) +
  #   scale_color_manual(values = hybrid_color) +
  #   theme_classic() +
  #   theme(
  #     panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  #   ) +
  #   labs(
  #     x = "Week",
  #     y = "Normalized mean log concentration",
  #     color = "Region"
  #   )

  ggplot(plot_long_line, aes(x = plot_week, color = region)) +
    geom_smooth(
      aes(
        y = value_filled_norm,
        group = region
      ),
      method = "loess",
      formula = y ~ x,
      span = 0.3,
      se = FALSE,
      linewidth = 0.8,
      alpha = 0.9
    ) +
    # geom_point(
    #   data = plot_long_points,
    #   aes(y = value_norm),
    #   size = 1.8,
    #   alpha = 0.5
    # ) +
    facet_grid(alignment ~ species) +
    scale_color_manual(values = hybrid_color) +
    scale_x_continuous(
      breaks = seq(0, 52, by = 13),
      limits = c(1, 52)
    ) +
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



##################
#NEW TRIPLE PLOT
##################
plot_alignment_all_species <- function() {

  library(patchwork)

  region_levels <- c("MAG", "PEI", "HAL", "BOF", "GOM")

  region_shapes <- c(
    MAG = 16,
    PEI = 17,
    HAL = 15,
    BOF = 18,
    GOM = 8
  )

  legend_guides <- guides(
    color = guide_legend(
      override.aes = list(
        linetype = 1,
        shape = region_shapes[rev(region_levels)],
        linewidth = 0.8,
        size = 3
      )
    ),
    shape = "none"
  )

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

      week_shifted = ((week + round(position) - 1) %% 52) + 1
    ) %>%
    ungroup() %>%
    mutate(
      species = factor(species, levels = species_order),
      region = factor(region, levels = rev(region_levels))
    ) %>%
    filter(!is.na(species), !is.na(region))

  plot_long_line <- plot_dat %>%
    select(species, region, week, week_shifted, value_filled_norm) %>%
    tidyr::pivot_longer(
      cols = c(week, week_shifted),
      names_to = "alignment",
      values_to = "plot_week"
    ) %>%
    mutate(
      alignment = dplyr::recode(
        alignment,
        week = "Before alignment",
        week_shifted = "After alignment"
      ),
      alignment = factor(alignment, levels = c("Before alignment", "After alignment"))
    )

  make_curve_plot <- function(alignment_name, y_label, show_strips = TRUE) {

    p <- plot_long_line %>%
      filter(alignment == alignment_name) %>%
      ggplot(aes(x = plot_week, color = region, shape = region)) +
      geom_smooth(
        aes(
          y = value_filled_norm,
          group = region
        ),
        method = "loess",
        formula = y ~ x,
        span = 0.3,
        se = FALSE,
        linewidth = 0.8,
        alpha = 0.9
      ) +
      facet_grid(. ~ species) +
      scale_color_manual(values = hybrid_color, drop = FALSE) +
      scale_shape_manual(values = region_shapes, drop = FALSE) +
      scale_x_continuous(
        breaks = seq(0, 52, by = 13),
        limits = c(1, 52)
      ) +
      legend_guides +
      theme_classic() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
      ) +
      labs(
        x = NULL,
        y = y_label,
        color = "Region",
        shape = "Region"
      )

    if (!show_strips) {
      p <- p +
        theme(
          strip.text = element_blank(),
          strip.background = element_blank()
        )
    }

    p
  }

  lag_panel_df <- region_positions_by_species %>%
    filter(!is.na(position)) %>%
    mutate(
      position_plot = -position,
      species = factor(species, levels = species_order),
      region = factor(region, levels = rev(region_levels))
    ) %>%
    filter(!is.na(species), !is.na(region))

  p_before <- make_curve_plot(
    alignment_name = "Before alignment",
    y_label = "Before\nalignment",
    show_strips = TRUE
  )

  p_lags <- ggplot(
    lag_panel_df,
    aes(x = position_plot, y = region, color = region, shape = region)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(size = 3, alpha = 0.95) +
    facet_grid(. ~ species) +
    scale_color_manual(values = hybrid_color, drop = FALSE) +
    scale_shape_manual(values = region_shapes, drop = FALSE) +
    legend_guides +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      strip.text = element_blank(),
      strip.background = element_blank()
    ) +
    labs(
      x = NULL,
      y = "Regional\nlag",
      color = "Region",
      shape = "Region"
    )

  p_after <- make_curve_plot(
    alignment_name = "After alignment",
    y_label = "After\nalignment",
    show_strips = FALSE
  ) +
    labs(x = "Week")

  (p_before / p_lags / p_after) +
    plot_layout(
      heights = c(1, 0.55, 1),
      guides = "collect"
    ) &
    theme(
      legend.position = "bottom"
    )
}

plot_alignment_all_species()







































































#######################
#ALL SPECIES ALL REGIONS
#######################
# ----------------------------
# Libraries
# ----------------------------
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(gt)
library(patchwork)

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

# ----------------------------
# 3. Solve reconciled regional positions per species
# ----------------------------
solve_region_positions <- function(lag_df) {

  lag_df <- lag_df %>%
    filter(!is.na(lag), !is.na(cor)) %>%
    mutate(weight = pmax(cor, 0)^2) %>%
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

  fit <- lm(lag_df$lag ~ X - 1, weights = lag_df$weight)

  beta <- coef(fit)
  names(beta) <- fit_regions

  positions <- c(beta, setNames(0, ref_region))
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

# ----------------------------
# 4. Regional timing summary table
# Higher mean_position = earlier
# ----------------------------
region_position_summary <- region_positions_by_species %>%
  group_by(region) %>%
  summarise(
    mean_position = mean(position, na.rm = TRUE),
    sd_position = sd(position, na.rm = TRUE),
    n_species = sum(!is.na(position)),
    .groups = "drop"
  ) %>%
  mutate(
    rank = rank(-mean_position, ties.method = "first")
  ) %>%
  arrange(rank)

region_position_summary %>%
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

# ----------------------------
# 5. Lag timeline figure
# ----------------------------
region_offsets <- c(
  MAG = -0.1,
  PEI = -0.05,
  HAL =  0.0,
  BOF =  0.05,
  GOM =  0.1
)

lag_timeline_df <- region_positions_by_species %>%
  filter(!is.na(position)) %>%
  mutate(
    position_plot = -position,
    species = factor(species, levels = rev(species_order)),
    region = factor(region, levels = names(region_offsets))
  ) %>%
  filter(!is.na(species), !is.na(region)) %>%
  mutate(
    species = droplevels(species),
    species_num = as.numeric(species),
    y_pos = species_num + region_offsets[as.character(region)]
  )

species_breaks <- seq_along(levels(lag_timeline_df$species))
between_species_lines <- species_breaks[-length(species_breaks)] + 0.5

ggplot(lag_timeline_df, aes(x = position_plot, y = y_pos, color = region, shape = region)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_hline(
    yintercept = between_species_lines,
    color = "grey25",
    linewidth = 0.6
  ) +
  geom_point(size = 3, alpha = 1) +
  scale_y_continuous(
    breaks = species_breaks,
    labels = levels(lag_timeline_df$species)
  ) +
  scale_color_manual(values = hybrid_color) +
  scale_shape_manual(values = c(
    MAG = 16,
    PEI = 17,
    HAL = 15,
    BOF = 18,
    GOM = 8
  )) +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  ) +
  labs(
    x = "Relative seasonal timing (weeks)",
    y = "Species",
    color = "Region",
    shape = "Region"
  )

# ----------------------------
# 6. Before / lag / after alignment figure
# ----------------------------
plot_alignment_all_species <- function() {

  region_levels <- c("MAG", "PEI", "HAL", "BOF", "GOM")

  region_shapes <- c(
    MAG = 16,
    PEI = 17,
    HAL = 15,
    BOF = 18,
    GOM = 8
  )

  plot_dat <- qpcr_weekly %>%
    left_join(
      region_positions_by_species %>%
        select(species, region, position),
      by = c("species", "region")
    ) %>%
    group_by(species, region) %>%
    mutate(
      value_filled = fill_circular_series(value),

      value_min = min(value_filled, na.rm = TRUE),
      value_max = max(value_filled, na.rm = TRUE),

      value_filled_norm = dplyr::if_else(
        value_max > value_min,
        (value_filled - value_min) / (value_max - value_min),
        NA_real_
      ),

      week_shifted = ((week + round(position) - 1) %% 52) + 1
    ) %>%
    ungroup() %>%
    mutate(
      species = factor(species, levels = species_order),
      region = factor(region, levels = rev(region_levels))
    ) %>%
    filter(!is.na(species), !is.na(region), !is.na(position))

  plot_long_line <- plot_dat %>%
    select(species, region, week, week_shifted, value_filled_norm) %>%
    pivot_longer(
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
      alignment = factor(
        alignment,
        levels = c("Before alignment", "After alignment")
      )
    )

  make_curve_plot <- function(alignment_name, y_label, show_strips = TRUE) {

    p <- plot_long_line %>%
      filter(alignment == alignment_name) %>%
      ggplot(aes(x = plot_week, color = region)) +
      geom_smooth(
        aes(
          y = value_filled_norm,
          group = region
        ),
        method = "loess",
        formula = y ~ x,
        span = 0.3,
        se = FALSE,
        linewidth = 0.8,
        alpha = 0.9
      ) +
      facet_grid(. ~ species) +
      scale_color_manual(values = hybrid_color, drop = FALSE) +
      scale_x_continuous(
        breaks = seq(0, 52, by = 13),
        limits = c(1, 52)
      ) +
      theme_classic() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        legend.position = "none"
      ) +
      labs(
        x = NULL,
        y = y_label
      )

    if (!show_strips) {
      p <- p +
        theme(
          strip.text = element_blank(),
          strip.background = element_blank()
        )
    }

    p
  }

  lag_panel_df <- region_positions_by_species %>%
    filter(!is.na(position)) %>%
    mutate(
      position_plot = -position,
      species = factor(species, levels = species_order),
      region = factor(region, levels = rev(region_levels))
    ) %>%
    filter(!is.na(species), !is.na(region))

  p_before <- make_curve_plot(
    alignment_name = "Before alignment",
    y_label = "Before\nalignment",
    show_strips = TRUE
  )

  p_lags <- ggplot(
    lag_panel_df,
    aes(x = position_plot, y = region, color = region, shape = region)
  ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(size = 3, alpha = 0.95) +
    facet_grid(. ~ species) +
    scale_color_manual(values = hybrid_color, drop = FALSE) +
    scale_shape_manual(values = region_shapes, drop = FALSE) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      strip.text = element_blank(),
      strip.background = element_blank()
    ) +
    labs(
      x = NULL,
      y = "Regional\nlag",
      color = "Region",
      shape = "Region"
    )

  p_after <- make_curve_plot(
    alignment_name = "After alignment",
    y_label = "After\nalignment",
    show_strips = FALSE
  ) +
    labs(x = "Week")

  p_before / p_lags / p_after +
    plot_layout(
      heights = c(1, 0.55, 1)
    )
}

plot_alignment_all_species()

# ----------------------------
# 7. Permutation test — unrestricted
# ----------------------------
obs_stat <- region_position_summary %>%
  summarise(stat = var(mean_position, na.rm = TRUE)) %>%
  pull(stat)

permute_once <- function(df) {
  df %>%
    group_by(species) %>%
    mutate(region = sample(region)) %>%
    ungroup()
}

compute_stat <- function(df) {

  pairwise_lags_perm <- df %>%
    group_by(species) %>%
    group_modify(~{

      dat <- .x
      regs <- sort(unique(as.character(dat$region)))
      region_pairs <- combn(regs, 2, simplify = FALSE)

      map_dfr(region_pairs, function(pair) {

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
          return(tibble(region1 = r1, region2 = r2, lag = NA_real_, cor = NA_real_))
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
    summarise(
      mean_position = mean(position, na.rm = TRUE),
      .groups = "drop"
    )

  var(region_summary_perm$mean_position, na.rm = TRUE)
}

set.seed(123)

n_perm <- 1000
perm_stats <- numeric(n_perm)

start_time <- Sys.time()

for (i in seq_len(n_perm)) {

  df_perm <- permute_once(qpcr_weekly)

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




