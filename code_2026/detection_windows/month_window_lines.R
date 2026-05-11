library(dplyr)
library(purrr)
library(ggplot2)
library(tidyr)



###################################
# COLLECT WINDOW + PEAK FOR PLOTTING
###################################

collect_window_plot_data <- function(df_monthly, threshold = 0.75) {

  all_species <- sort(unique(df_monthly$species))
  all_regions <- c("MAG", "PEI", "HAL", "BOF", "GOM")

  combos <- expand.grid(
    species = all_species,
    region = all_regions,
    stringsAsFactors = FALSE
  )

  purrr::pmap_dfr(combos, function(species, region) {

    df_sub_monthly <- df_monthly %>%
      dplyr::filter(
        species == !!species,
        region == !!region
      )

    if (nrow(df_sub_monthly) == 0) {
      return(NULL)
    }

    plot_df <- df_sub_monthly %>%
      prep_monthly_signal() %>%
      interp_monthly_circular() %>%
      classify_months(threshold = threshold)

    window_res <- calc_window_simple(plot_df, threshold = threshold)

    if (is.null(window_res)) {
      return(NULL)
    }

    peak_month <- plot_df$month[which.max(plot_df$value)]

    inside_months <- if (!window_res$wrap_around) {
      window_res$start_month:window_res$end_month
    } else {
      c(window_res$start_month:12, 1:window_res$end_month)
    }

    above_outside_months <- plot_df %>%
      dplyr::filter(value >= threshold, !month %in% inside_months) %>%
      dplyr::pull(month)

    tibble::tibble(
      species = species,
      region = region,
      start_month = window_res$start_month,
      end_month = window_res$end_month,
      wrap_around = window_res$wrap_around,
      peak_month = peak_month,
      above_outside_months = list(above_outside_months)
    )
  })
}

###################################
# FIXED MONTH SHIFT: MARCH -> NEXT MARCH
###################################

shift_month_may <- function(m) {
  ifelse(m < 5, m + 12, m)
}

###################################
# PREP WINDOW SEGMENTS FOR PLOTTING
###################################

prep_window_segments <- function(window_df) {

  window_df %>%
    dplyr::mutate(
      start_plot = shift_month_may(start_month),
      end_plot   = shift_month_may(end_month),
      peak_plot  = shift_month_may(peak_month)
    ) %>%
    dplyr::mutate(
      end_plot = ifelse(end_plot < start_plot, end_plot + 12, end_plot),
      peak_plot = ifelse(peak_plot < start_plot, peak_plot + 12, peak_plot),

      # expand windows to cover full month bins
      start_plot = start_plot - 0.5,
      end_plot   = end_plot + 0.5
    )
}

###################################
# AXIS LABELS
###################################

shifted_month_labels <- function(x) {
  m <- ((round(x) - 1) %% 12) + 1
  c("J","F","M","A","M","J","J","A","S","O","N","D")[m]
}

###################################
# BUILD PLOT DATA
###################################

species_order <- c(
  "Membranipora membranacea",
  "Botrylloides violaceus",
  "Didemnum vexillum",
  "Ciona intestinalis",
  "Carcinus maenas"
)

region_order <- c("MAG", "PEI", "HAL", "BOF", "GOM")

region_colors <- c(
  MAG = "#00A08A",
  PEI = "#446455",
  HAL = "#CCAA4F",
  BOF = "#5BBCD6",
  GOM = "#fb8072"
)

region_colors <- c(
  MAG = "#008C78",
  PEI = "#3A5549",
  HAL = "#B89945",
  BOF = "#4FA7C0",
  GOM = "#E57266"
)

region_colors <- c(
  MAG = "#007566",
  PEI = "#31493E",
  HAL = "#A8893D",
  BOF = "#3E8FA6",
  GOM = "#D96459"
)



window_plot_df <- collect_window_plot_data(
  df_monthly = dfMonths,
  threshold = 0.9
) %>%
  dplyr::left_join(
    wilcox_results_df %>%
      dplyr::select(species, region, prob_superiority),
    by = c("species", "region")
  ) %>%
  dplyr::mutate(
    species = factor(species, levels = species_order),
    region = factor(region, levels = rev(region_order)),
    alpha_pos = pmax(0, pmin(1, (prob_superiority - 0.5) / 0.5))
  )

window_plot_df2 <- prep_window_segments(window_plot_df)

x_min <- 4.5
x_max <- 16.5


above_outside_df <- window_plot_df %>%
  tidyr::unnest_longer(above_outside_months, values_to = "month") %>%
  dplyr::filter(!is.na(month)) %>%
  dplyr::mutate(
    month_plot = shift_month_may(month),
    x_start = month_plot - 0.5,
    x_end   = month_plot + 0.5,
    species = factor(species, levels = species_order),
    region = factor(region, levels = rev(region_order))
  )
###################################
# PLOT
###################################

ggplot(window_plot_df2) +
  geom_segment(
    data = above_outside_df,
    aes(
      x = x_start,
      xend = x_end,
      y = region,
      yend = region
    ),
    color = "grey75",
    linewidth = 4,
    lineend = "round"
  ) +
  geom_segment(
    aes(
      x = start_plot,
      xend = end_plot,
      y = region,
      yend = region,
      color = region,
      alpha = alpha_pos
    ),
    linewidth = 4,
    lineend = "round"
  ) +
  scale_alpha_continuous(range = c(0, 1), limits = c(0, 1), guide = "none") +
  geom_point(
    aes(
      x = peak_plot,
      y = region
    ),
    color = "black",
    size = 2.5
  ) +
  facet_wrap(~ species, ncol = 1) +
  scale_color_manual(values = region_colors) +
  scale_x_continuous(
    limits = c(x_min, x_max + 1),
    breaks = seq(5, 16, by = 1),
    labels = shifted_month_labels,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
  scale_y_discrete(limits = rev(region_order)) +
  labs(
    x = "Month",
    y = NULL,
    title = "Optimal timing windows by species and region",
    subtitle = "Colored line = optimal window, black point = peak month"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "italic"),
    axis.text.y = element_text(face = "bold")
  ) +
  geom_text(
    aes(
      x = x_max + 0.2,   # moved left (was +0.5)
      y = region,
      label = sprintf("%.0f%%", prob_superiority * 100)
    ),
    hjust = 0,
    size = 3,
    color = "grey10"
  )


avg_ps_species <- window_plot_df2 %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    mean_ps = mean(prob_superiority, na.rm = TRUE),
    sd_ps = sd(prob_superiority, na.rm = TRUE),
    n = dplyr::n(),
    .groups = "drop"
  )

avg_ps_species


























##################################
#FRIEDMAN
##################################

circular_midpoint <- function(start, end, n = 12) {
  d <- ((end - start + n) %% n)
  ((start + d / 2 - 1) %% n) + 1
}

min_n_regions <- 2

window_centers_df <- window_plot_df %>%
  dplyr::mutate(
    center_month = circular_midpoint(start_month, end_month)
  )

# keep only species with full regional coverage
window_centers_df_complete <- window_centers_df %>%
  dplyr::group_by(species) %>%
  dplyr::filter(dplyr::n_distinct(region) >= min_n_regions) %>%
  dplyr::ungroup()

center_matrix <- window_centers_df_complete %>%
  dplyr::select(species, region, center_month) %>%
  tidyr::pivot_wider(names_from = region, values_from = center_month) %>%
  dplyr::arrange(species)

center_matrix_complete <- center_matrix %>%
  tidyr::drop_na()

friedman_res <- friedman.test(as.matrix(center_matrix_complete[, -1]))

k <- ncol(center_matrix_complete[, -1])  # regions
n <- nrow(center_matrix_complete)        # species

W <- friedman_res$statistic / (n * (k - 1))
W

rank_df <- window_centers_df_complete %>%
  dplyr::mutate(
    center_shifted = shift_month_may(center_month)
  ) %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(
    rank = rank(center_shifted, ties.method = "average")
  ) %>%
  dplyr::ungroup()

region_rank_summary <- rank_df %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    mean_rank = mean(rank, na.rm = TRUE),
    sd_rank = sd(rank, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(mean_rank)


library(gt)
region_rank_summary %>%
  gt() %>%
  fmt_number(columns = c(mean_rank, sd_rank), decimals = 2) %>%
  tab_header(
    title = "Regional timing ranks"
  ) %>%
  cols_align(align = "center", -region)



###############
#NEW TIMING ANALYSIS
###############
month_offset <- function(x, ref, n = 12) {
  ((x - ref + n / 2) %% n) - n / 2
}

window_center_offsets <- window_centers_df_complete %>%
  dplyr::mutate(
    center_shifted = shift_month_may(center_month)
  ) %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(
    relative_month = center_shifted - mean(center_shifted, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()
region_center_summary <- window_center_offsets %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    mean_relative_month = mean(relative_month, na.rm = TRUE),
    sd_relative_month = sd(relative_month, na.rm = TRUE),
    n_species = dplyr::n_distinct(species),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    rank = rank(mean_relative_month, ties.method = "first")
  ) %>%
  dplyr::arrange(rank)

region_center_summary %>%
  gt() %>%
  fmt_number(columns = c(rank, sd_relative_month), decimals = 2) %>%
  tab_header(
    title = "Regional timing ranks"
  ) %>%
  cols_align(align = "center", -region)

obs_stat_centers <- region_center_summary %>%
  dplyr::summarise(
    stat = var(mean_relative_month, na.rm = TRUE)
  ) %>%
  dplyr::pull(stat)

obs_stat_centers


permute_once <- function(df) {
  df %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(region = sample(region)) %>%
    dplyr::ungroup()
}
compute_perm_stat_centers <- function(df_perm) {

  offsets <- df_perm %>%
    dplyr::mutate(
      center_shifted = shift_month_may(center_month)
    ) %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(
      relative_month = center_shifted - mean(center_shifted, na.rm = TRUE)
    ) %>%
    dplyr::ungroup()

  region_summary <- offsets %>%
    dplyr::group_by(region) %>%
    dplyr::summarise(
      mean_relative_month = mean(relative_month, na.rm = TRUE),
      .groups = "drop"
    )

  var(region_summary$mean_relative_month, na.rm = TRUE)
}
set.seed(123)

perm_stats <- replicate(
  1000,
  compute_perm_stat_centers(
    permute_once(window_centers_df_complete)
  )
)

p_value <- mean(perm_stats >= obs_stat_centers)





###############
# WINDOW CENTER TIMING ANALYSIS
# pairwise circular month differences + reconciled regional positions
###############
###############
# WINDOW CENTER TIMING ANALYSIS
# coherent per-species timing positions
###############

circular_midpoint <- function(start, end, n = 12) {
  d <- ((end - start + n) %% n)
  ((start + d / 2 - 1) %% n) + 1
}

shift_month_may <- function(m) {
  ifelse(m < 5, m + 12, m)
}

min_n_regions <- 2

window_centers_df <- window_plot_df %>%
  dplyr::mutate(
    center_month = circular_midpoint(start_month, end_month)
  )

window_centers_df_use <- window_centers_df %>%
  dplyr::group_by(species) %>%
  dplyr::filter(dplyr::n_distinct(region) >= min_n_regions) %>%
  dplyr::ungroup()

window_center_offsets <- window_centers_df_use %>%
  dplyr::mutate(
    center_unwrapped = shift_month_may(center_month)
  ) %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(
    relative_month = center_unwrapped - mean(center_unwrapped, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

region_center_summary <- window_center_offsets %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    mean_relative_month = mean(relative_month, na.rm = TRUE),
    sd_relative_month = sd(relative_month, na.rm = TRUE),
    n_species = dplyr::n_distinct(species),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    rank = rank(mean_relative_month, ties.method = "first")
  ) %>%
  dplyr::arrange(rank)

region_center_summary

obs_stat_centers <- region_center_summary %>%
  dplyr::summarise(
    stat = var(mean_relative_month, na.rm = TRUE)
  ) %>%
  dplyr::pull(stat)

obs_stat_centers


###############
# PERMUTATION TEST FOR WINDOW CENTERS
###############

permute_once_centers <- function(df) {
  df %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(region = sample(region)) %>%
    dplyr::ungroup()
}

compute_perm_stat_centers <- function(df_perm) {

  offsets_perm <- df_perm %>%
    dplyr::mutate(
      center_unwrapped = shift_month_may(center_month)
    ) %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(
      relative_month = center_unwrapped - mean(center_unwrapped, na.rm = TRUE)
    ) %>%
    dplyr::ungroup()

  summary_perm <- offsets_perm %>%
    dplyr::group_by(region) %>%
    dplyr::summarise(
      mean_relative_month = mean(relative_month, na.rm = TRUE),
      .groups = "drop"
    )

  var(summary_perm$mean_relative_month, na.rm = TRUE)
}

set.seed(123)

n_perm <- 1000
perm_stats_centers <- numeric(n_perm)

start_time <- Sys.time()

for (i in seq_len(n_perm)) {

  perm_stats_centers[i] <- compute_perm_stat_centers(
    permute_once_centers(window_centers_df_use)
  )

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

p_value_centers <- mean(perm_stats_centers >= obs_stat_centers, na.rm = TRUE)

p_value_centers

#peak month rather than window center

###############
# PEAK MONTH TIMING ANALYSIS
# coherent per-species timing positions
###############

min_n_regions <- 2

peak_month_df_use <- window_plot_df %>%
  dplyr::group_by(species) %>%
  dplyr::filter(dplyr::n_distinct(region) >= min_n_regions) %>%
  dplyr::ungroup()

peak_month_offsets <- peak_month_df_use %>%
  dplyr::mutate(
    peak_unwrapped = shift_month_may(peak_month)
  ) %>%
  dplyr::group_by(species) %>%
  dplyr::mutate(
    relative_month = peak_unwrapped - mean(peak_unwrapped, na.rm = TRUE)
  ) %>%
  dplyr::ungroup()

region_peak_summary <- peak_month_offsets %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    mean_relative_month = mean(relative_month, na.rm = TRUE),
    sd_relative_month = sd(relative_month, na.rm = TRUE),
    n_species = dplyr::n_distinct(species),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    rank = rank(mean_relative_month, ties.method = "first")
  ) %>%
  dplyr::arrange(rank)

region_peak_summary

obs_stat_peaks <- region_peak_summary %>%
  dplyr::summarise(
    stat = var(mean_relative_month, na.rm = TRUE)
  ) %>%
  dplyr::pull(stat)

obs_stat_peaks


###############
# PERMUTATION TEST FOR PEAK MONTHS
###############

permute_once_peaks <- function(df) {
  df %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(region = sample(region)) %>%
    dplyr::ungroup()
}

compute_perm_stat_peaks <- function(df_perm) {

  offsets_perm <- df_perm %>%
    dplyr::mutate(
      peak_unwrapped = shift_month_may(peak_month)
    ) %>%
    dplyr::group_by(species) %>%
    dplyr::mutate(
      relative_month = peak_unwrapped - mean(peak_unwrapped, na.rm = TRUE)
    ) %>%
    dplyr::ungroup()

  summary_perm <- offsets_perm %>%
    dplyr::group_by(region) %>%
    dplyr::summarise(
      mean_relative_month = mean(relative_month, na.rm = TRUE),
      .groups = "drop"
    )

  var(summary_perm$mean_relative_month, na.rm = TRUE)
}

set.seed(123)

n_perm <- 1000
perm_stats_peaks <- numeric(n_perm)

start_time <- Sys.time()

for (i in seq_len(n_perm)) {

  perm_stats_peaks[i] <- compute_perm_stat_peaks(
    permute_once_peaks(peak_month_df_use)
  )

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

p_value_peaks <- mean(perm_stats_peaks >= obs_stat_peaks, na.rm = TRUE)

p_value_peaks





































######################
#HOW CLUSTERED ARE THE WINDOWS?
######################

library(circular)
compute_species_clustering <- function(df, sp) {

  dat <- df %>%
    dplyr::filter(species == sp)

  n_regions <- nrow(dat)

  if (n_regions < 2) {
    return(tibble::tibble(
      species = sp,
      n_regions = n_regions,
      clustering_R = NA_real_,
      p_value = NA_real_
    ))
  }

  theta <- 2 * pi * dat$center_month / 12
  circ <- circular::circular(theta)

  R <- circular::rho.circular(circ)

  p_val <- if (n_regions >= 3) {
    circular::rayleigh.test(circ)$p.value
  } else {
    NA_real_
  }

  tibble::tibble(
    species = sp,
    n_regions = n_regions,
    clustering_R = R,
    p_value = p_val
  )
}

circular_midpoint <- function(start, end, n = 12) {
  d <- ((end - start + n) %% n)
  ((start + d / 2 - 1) %% n) + 1
}

window_plot_df_with_centers <- window_plot_df %>%
  dplyr::mutate(
    center_month = circular_midpoint(start_month, end_month)
  )

species_clustering_table <- window_plot_df_with_centers %>%
  dplyr::distinct(species) %>%
  dplyr::pull(species) %>%
  purrr::map_dfr(~ compute_species_clustering(window_plot_df_with_centers, .x)) %>%
  dplyr::mutate(
    species = factor(species, levels = species_order)
  ) %>%
  dplyr::arrange(species)

species_clustering_table %>%
  gt::gt() %>%
  gt::fmt_number(columns = clustering_R, decimals = 2) %>%
  gt::tab_header(title = "Clustering of Window Centers Across Regions")














################################
#JACCARD
################################

jaccard_region_consistency <- function(
    window_df,
    species_name,
    weight_col = NULL
) {

  # subset to one species
  df <- window_df %>%
    dplyr::filter(species == species_name)

  if (nrow(df) < 2) {
    return("Not enough regions")
  }

  # convert window to binary
  window_to_binary <- function(start_month, end_month, wrap_around, n = 12) {
    out <- rep(0, n)

    if (!wrap_around) {
      out[start_month:end_month] <- 1
    } else {
      out[c(start_month:n, 1:end_month)] <- 1
    }

    out
  }

  # build binary matrix
  region_windows <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      bin = list(window_to_binary(start_month, end_month, wrap_around))
    ) %>%
    tidyr::unnest_wider(bin, names_sep = "_") %>%
    dplyr::ungroup()

  mat <- region_windows %>%
    dplyr::select(dplyr::starts_with("bin_")) %>%
    as.matrix()

  rownames(mat) <- region_windows$region

  # optional weights
  if (!is.null(weight_col) && weight_col %in% names(region_windows)) {
    weights <- region_windows[[weight_col]]
  } else {
    weights <- rep(1, nrow(mat))
  }

  # pairwise comparisons
  pair_idx <- utils::combn(seq_len(nrow(mat)), 2, simplify = FALSE)

  jacc_df <- purrr::map_dfr(pair_idx, function(idx) {

    i <- idx[1]
    j <- idx[2]

    v1 <- mat[i, ]
    v2 <- mat[j, ]

    n_int <- sum(v1 == 1 & v2 == 1)

    # same "union" logic as before
    n_union <- max(sum(v1), sum(v2))

    n_wt <- weights[i] + weights[j]

    tibble::tibble(
      region_pair = paste(rownames(mat)[i], rownames(mat)[j], sep = "_"),
      n_int = n_int,
      n_union = n_union,
      J_index = ifelse(n_union == 0, NA_real_, n_int / n_union),
      n_wt = n_wt,
      J_wt = J_index * n_wt
    )
  })

  jacc_wt_mean <- jacc_df %>%
    dplyr::summarise(
      wt_mean = sum(J_wt, na.rm = TRUE) / sum(n_wt) * 100
    ) %>%
    dplyr::mutate(
      wt_text = dplyr::case_when(
        dplyr::between(wt_mean, 0, 9.99) ~ "Very low",
        dplyr::between(wt_mean, 10, 29.99) ~ "Low",
        dplyr::between(wt_mean, 30, 69.99) ~ "Medium",
        dplyr::between(wt_mean, 70, 89.99) ~ "High",
        dplyr::between(wt_mean, 90, 100) ~ "Very high"
      )
    )

  list(
    summary = jacc_wt_mean,
    pairwise = jacc_df,
    region_windows = region_windows
  )
}

###################################################################
jaccard_true_region_consistency <- function(
    window_df,
    species_name,
    weight_col = NULL
) {

  # subset to one species
  df <- window_df %>%
    dplyr::filter(species == species_name)

  if (nrow(df) < 2) {
    return("Not enough regions")
  }

  # convert window to binary
  window_to_binary <- function(start_month, end_month, wrap_around, n = 12) {
    out <- rep(0, n)

    if (!wrap_around) {
      out[start_month:end_month] <- 1
    } else {
      out[c(start_month:n, 1:end_month)] <- 1
    }

    out
  }

  # build binary matrix
  region_windows <- df %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      bin = list(window_to_binary(start_month, end_month, wrap_around))
    ) %>%
    tidyr::unnest_wider(bin, names_sep = "_") %>%
    dplyr::ungroup()

  mat <- region_windows %>%
    dplyr::select(dplyr::starts_with("bin_")) %>%
    as.matrix()

  rownames(mat) <- region_windows$region

  # optional weights
  if (!is.null(weight_col) && weight_col %in% names(region_windows)) {
    weights <- region_windows[[weight_col]]
  } else {
    weights <- rep(1, nrow(mat))
  }

  # pairwise comparisons
  pair_idx <- utils::combn(seq_len(nrow(mat)), 2, simplify = FALSE)

  jacc_df <- purrr::map_dfr(pair_idx, function(idx) {

    i <- idx[1]
    j <- idx[2]

    v1 <- mat[i, ]
    v2 <- mat[j, ]

    n_int <- sum(v1 == 1 & v2 == 1)

    # true jaccard logic?
    n_union <- min(sum(v1), sum(v2))

    n_wt <- weights[i] + weights[j]

    tibble::tibble(
      region_pair = paste(rownames(mat)[i], rownames(mat)[j], sep = "_"),
      n_int = n_int,
      n_union = n_union,
      J_index = ifelse(n_union == 0, NA_real_, n_int / n_union),
      n_wt = n_wt,
      J_wt = J_index * n_wt
    )
  })

  jacc_wt_mean <- jacc_df %>%
    dplyr::summarise(
      wt_mean = sum(J_wt, na.rm = TRUE) / sum(n_wt) * 100
    ) %>%
    dplyr::mutate(
      wt_text = dplyr::case_when(
        dplyr::between(wt_mean, 0, 9.99) ~ "Very low",
        dplyr::between(wt_mean, 10, 29.99) ~ "Low",
        dplyr::between(wt_mean, 30, 69.99) ~ "Medium",
        dplyr::between(wt_mean, 70, 89.99) ~ "High",
        dplyr::between(wt_mean, 90, 100) ~ "Very high"
      )
    )

  list(
    summary = jacc_wt_mean,
    pairwise = jacc_df,
    region_windows = region_windows
  )
}



res <- jaccard_true_region_consistency(
  window_df = window_plot_df,
  species_name = "Carcinus maenas"
)

res$summary
res$pairwise


compute_species_overlap <- function(window_df, sp) {

  res <- tryCatch(
    jaccard_true_region_consistency(
      window_df = window_df,
      species_name = sp
    ),
    error = function(e) NULL
  )

  if (is.null(res) || is.character(res)) {
    return(tibble::tibble(
      species = sp,
      n_regions = NA,
      mean_overlap = NA,
      overlap_class = NA
    ))
  }

  tibble::tibble(
    species = sp,
    n_regions = nrow(res$region_windows),
    mean_overlap = res$summary$wt_mean,
    overlap_class = res$summary$wt_text
  )
}


species_overlap_table <- window_plot_df %>%
  dplyr::distinct(species) %>%
  dplyr::pull(species) %>%
  purrr::map_dfr(~ compute_species_overlap(window_plot_df, .x))

species_overlap_table <- species_overlap_table %>%
  dplyr::mutate(
    species = factor(species, levels = species_order)
  ) %>%
  dplyr::arrange(species)


species_overlap_table %>%
  gt() %>%
  fmt_number(columns = mean_overlap, decimals = 2) %>%
  tab_header(title = "Overlap Between Windows - Across Regions & Within Species") %>%
  cols_align(
    align = "center",
    columns = everything()
  ) %>%
  tab_options(
    table.font.color = "black"
  )

#####################################
#WINDOW LENGTHS
#####################################
calc_window_length <- function(start_month, end_month, wrap_around, n = 12) {
  dplyr::if_else(
    !wrap_around,
    end_month - start_month + 1,
    (n - start_month + 1) + end_month
  )
}


window_length_df <- window_plot_df %>%
  dplyr::mutate(
    window_length = calc_window_length(
      start_month = start_month,
      end_month = end_month,
      wrap_around = wrap_around
    )
  )

####################################
#WHICH REGIONS ARE SIMILAR ACROSS SPECIES
####################################

calc_window_months <- function(start_month, end_month, wrap_around) {
  if (!wrap_around) {
    start_month:end_month
  } else {
    c(start_month:12, 1:end_month)
  }
}

region_pair_overlap_df <- window_plot_df %>%
  dplyr::mutate(
    window_months = purrr::pmap(
      list(start_month, end_month, wrap_around),
      calc_window_months
    )
  ) %>%
  dplyr::select(species, region, window_months) %>%
  dplyr::group_by(species) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    pairwise = purrr::map(data, function(dat) {
      region_pairs <- combn(dat$region, 2, simplify = FALSE)

      purrr::map_dfr(region_pairs, function(pair) {
        m1 <- dat$window_months[[which(dat$region == pair[1])]]
        m2 <- dat$window_months[[which(dat$region == pair[2])]]

        tibble::tibble(
          region_1 = pair[1],
          region_2 = pair[2],
          overlap = length(intersect(m1, m2)) / length(union(m1, m2))
        )
      })
    })
  ) %>%
  dplyr::select(species, pairwise) %>%
  tidyr::unnest(pairwise) %>%
  dplyr::ungroup()


region_similarity_table <- region_pair_overlap_df %>%
  dplyr::group_by(region_1, region_2) %>%
  dplyr::summarise(
    n_species = dplyr::n(),
    mean_overlap = mean(overlap, na.rm = TRUE),
    sd_overlap = sd(overlap, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(mean_overlap))


region_similarity_table %>%
  gt::gt() %>%
  gt::fmt_number(columns = c(mean_overlap, sd_overlap), decimals = 2) %>%
  gt::tab_header(
    title = "Similarity Between Regions Across Species"
  ) %>%
  gt::cols_align(
    align = "center",
    columns = everything()
  )






















calc_window_months <- function(start_month, end_month, wrap_around) {
  if (!wrap_around) {
    start_month:end_month
  } else {
    c(start_month:12, 1:end_month)
  }
}

region_pair_overlap_df <- window_plot_df %>%
  dplyr::mutate(
    window_months = purrr::pmap(
      list(start_month, end_month, wrap_around),
      calc_window_months
    )
  ) %>%
  dplyr::select(species, region, window_months) %>%
  dplyr::group_by(species) %>%
  tidyr::nest() %>%
  dplyr::mutate(
    pairwise = purrr::map(data, function(dat) {
      region_pairs <- combn(dat$region, 2, simplify = FALSE)

      purrr::map_dfr(region_pairs, function(pair) {
        m1 <- dat$window_months[[which(dat$region == pair[1])]]
        m2 <- dat$window_months[[which(dat$region == pair[2])]]

        tibble::tibble(
          region_1 = pair[1],
          region_2 = pair[2],
          overlap = length(intersect(m1, m2)) / min(length(m1), length(m2))
        )
      })
    })
  ) %>%
  dplyr::select(species, pairwise) %>%
  tidyr::unnest(pairwise) %>%
  dplyr::ungroup()


region_similarity_table <- region_pair_overlap_df %>%
  dplyr::group_by(region_1, region_2) %>%
  dplyr::summarise(
    n_species = dplyr::n(),
    mean_overlap = mean(overlap, na.rm = TRUE),
    sd_overlap = sd(overlap, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(mean_overlap))


region_similarity_table %>%
  gt::gt() %>%
  gt::fmt_number(columns = c(mean_overlap, sd_overlap), decimals = 2) %>%
  gt::tab_header(
    title = "Similarity Between Regions Across Species"
  ) %>%
  gt::cols_align(
    align = "center",
    columns = everything()
  )














####################
#Jaccard Apr 29
####################

window_months <- function(start_month, end_month, wrap_around) {
  if (!wrap_around) {
    start_month:end_month
  } else {
    c(start_month:12, 1:end_month)
  }
}

jaccard_overlap <- function(a, b) {
  length(intersect(a, b)) / length(union(a, b))
}

compute_species_overlap_jaccard <- function(window_df, sp) {

  dat <- window_df %>%
    dplyr::filter(species == sp) %>%
    dplyr::mutate(
      months = purrr::pmap(
        list(start_month, end_month, wrap_around),
        window_months
      )
    )

  if (nrow(dat) < 2) {
    return(tibble::tibble(
      species = sp,
      n_regions = nrow(dat),
      mean_overlap = NA_real_,
      overlap_class = NA_character_
    ))
  }

  pairs <- combn(seq_len(nrow(dat)), 2, simplify = FALSE)

  pairwise <- purrr::map_dfr(pairs, function(idx) {
    tibble::tibble(
      region_1 = as.character(dat$region[idx[1]]),
      region_2 = as.character(dat$region[idx[2]]),
      jaccard = jaccard_overlap(
        dat$months[[idx[1]]],
        dat$months[[idx[2]]]
      )
    )
  })

  tibble::tibble(
    species = sp,
    n_regions = nrow(dat),
    mean_overlap = mean(pairwise$jaccard, na.rm = TRUE),
    overlap_class = dplyr::case_when(
      mean_overlap >= 0.75 ~ "High",
      mean_overlap >= 0.50 ~ "Moderate",
      mean_overlap >= 0.25 ~ "Low",
      TRUE ~ "Very low"
    )
  )
}


species_overlap_table <- window_plot_df %>%
  dplyr::distinct(species) %>%
  dplyr::pull(species) %>%
  purrr::map_dfr(~ compute_species_overlap_jaccard(window_plot_df, .x)) %>%
  dplyr::mutate(
    species = factor(species, levels = species_order)
  ) %>%
  dplyr::arrange(species)

species_overlap_table %>%
  gt::gt() %>%
  gt::fmt_number(columns = mean_overlap, decimals = 2) %>%
  gt::tab_header(title = "Jaccard Overlap Between Windows - Across Regions & Within Species") %>%
  gt::cols_align(
    align = "center",
    columns = everything()
  ) %>%
  gt::tab_options(
    table.font.color = "black"
  )



##########################
#JACCARD ALL MONTHS ABOVE THRESHOLD
##########################

collect_above_threshold_months <- function(df_monthly, threshold = 0.9) {

  df_monthly %>%
    dplyr::group_by(species, region) %>%
    dplyr::group_modify(~ {
      plot_df <- .x %>%
        prep_monthly_signal() %>%
        interp_monthly_circular() %>%
        classify_months(threshold = threshold)

      tibble::tibble(
        months = list(plot_df$month[plot_df$value >= threshold])
      )
    }) %>%
    dplyr::ungroup()
}

above_threshold_months_df <- collect_above_threshold_months(
  df_monthly = dfMonths,
  threshold = 0.9
)

compute_species_above_threshold_jaccard <- function(month_df, sp) {

  dat <- month_df %>%
    dplyr::filter(species == sp)

  if (nrow(dat) < 2) {
    return(tibble::tibble(
      species = sp,
      n_regions = nrow(dat),
      mean_overlap = NA_real_
    ))
  }

  pairs <- combn(seq_len(nrow(dat)), 2, simplify = FALSE)

  pairwise <- purrr::map_dfr(pairs, function(idx) {
    tibble::tibble(
      region_1 = as.character(dat$region[idx[1]]),
      region_2 = as.character(dat$region[idx[2]]),
      jaccard = jaccard_overlap(
        dat$months[[idx[1]]],
        dat$months[[idx[2]]]
      )
    )
  })

  tibble::tibble(
    species = sp,
    n_regions = nrow(dat),
    mean_overlap = mean(pairwise$jaccard, na.rm = TRUE)
  )
}


species_above_threshold_overlap_table <- above_threshold_months_df %>%
  dplyr::distinct(species) %>%
  dplyr::pull(species) %>%
  purrr::map_dfr(~ compute_species_above_threshold_jaccard(
    above_threshold_months_df, .x
  )) %>%
  dplyr::mutate(
    species = factor(species, levels = species_order)
  ) %>%
  dplyr::arrange(species)
