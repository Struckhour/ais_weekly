source('./AIS_eDNA_data_prep.R')


###################
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

abbr_species <- function(x) {
  sapply(strsplit(as.character(x), " "), function(parts) {
    if (length(parts) >= 2) {
      paste0(substr(parts[1], 1, 1), ". ", parts[2])
    } else {
      x
    }
  })
}

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

# linear interpolation, leaving ends extended
fill_circular_series <- function(x) {
  idx <- which(is.finite(x))

  if (length(idx) == 0) return(rep(NA_real_, length(x)))
  if (length(idx) == 1) return(rep(x[idx], length(x)))

  approx(
    x = idx,
    y = x[idx],
    xout = seq_along(x),
    method = "linear",
    rule = 2
  )$y
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
  ungroup() %>%
  mutate(
    species_abbr = abbr_species(species)
  )

print(pairwise_lags_selected, n = Inf)

pairwise_lags_fit_selected <- pairwise_lags_selected %>%
  mutate(
    species_abbr = abbr_species(species)
  ) %>%
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
  mutate(
    species_abbr = abbr_species(species)
  ) %>%
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
    species = factor(
      abbr_species(species),
      levels = abbr_species(species_order)
    )
  )


manual_point <- tibble::tibble(
  species = "C. maenas (MAG removed)",
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
    species = factor(
      abbr_species(species),
      levels = rev(abbr_species(species_order))
    )
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
    point.padding = 0.6,
    segment.color = "grey50",
    max.overlaps = Inf,
    nudge_x = dplyr::case_when(
      species_fit_selected$species == "C. maenas" ~ 0.02,
      species_fit_selected$species == "M. membranacea" ~ 0.05,
      TRUE ~ 0
    ),
    nudge_y = dplyr::case_when(
      species_fit_selected$species == "C. intestinalis" ~ 0.1,
      species_fit_selected$species == "C. maenas" ~ -0.8,
      species_fit_selected$species == "M. membranacea" ~ -0.55,
      TRUE ~ 0
    )
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
ggsave("manuscript_figures/figure_5.png", combined_plot, width = 10, height = 10, dpi = 300)




























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
  scale_fill_manual(
    values = hybrid_color,
    drop = FALSE
  ) +
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

hybrid_color <- c(
  MAG = "#00A08A",
  PEI = "#446455",
  HAL = "#CCAA4F",
  BOF = "#5BBCD6",
  GOM = "#fb8072"
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
    filter(
      !is.na(species),
      !is.na(region),
      !is.na(position),
      !is.na(value_filled_norm)
    )

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
      ggplot(aes(x = plot_week, y = value_filled_norm, color = region)) +
      geom_smooth(
        aes(group = region),
        method = "loess",
        formula = y ~ x,
        span = 0.3,
        se = FALSE,
        linewidth = 0.8,
        alpha = 0.9
      ) +
      facet_grid(. ~ species) +
      scale_color_manual(
        values = hybrid_color,
        breaks = rev(region_levels),
        drop = FALSE,
        name = "Region"
      ) +
      scale_x_continuous(
        breaks = seq(0, 52, by = 13),
        limits = c(1, 52)
      ) +
      theme_classic() +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
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
      values = hybrid_color,
      breaks = rev(region_levels),
      drop = FALSE,
      name = "Region"
    ) +
    scale_shape_manual(
      values = region_shapes,
      breaks = rev(region_levels),
      drop = FALSE,
      name = "Region"
    ) +
    theme_classic() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      strip.text = element_blank(),
      strip.background = element_blank()
    ) +
    labs(
      x = NULL,
      y = "Regional\nlag"
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

p <- plot_alignment_all_species()
p
ggsave("manuscript_figures/figure_4.png", p, width = 10, height = 10, dpi = 300)

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




