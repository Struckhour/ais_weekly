
##############################
#WILCOXON TEST BETWEEN OPTIMAL WINDOW OF MONTHS AND SUBOPTIMAL MONTHS
##############################

calc_window_simple <- function(plot_df, threshold = 0.75) {

  df <- plot_df %>%
    dplyr::arrange(month) %>%
    dplyr::mutate(above = value >= threshold)

  x <- df$above
  n <- length(x)

  # handle full-year window
  if (all(x)) {
    return(list(
      start_month = 1,
      end_month = 12,
      wrap_around = FALSE
    ))
  }

  peak_month <- which.max(df$value)

  # expand right (+1 direction)
  right <- peak_month
  repeat {
    next_r <- (right %% n) + 1
    if (!x[next_r] || next_r == peak_month) break
    right <- next_r
  }

  # expand left (-1 direction)
  left <- peak_month
  repeat {
    next_l <- ((left - 2 + n) %% n) + 1
    if (!x[next_l] || next_l == peak_month) break
    left <- next_l
  }

  list(
    start_month = left,
    end_month = right,
    wrap_around = left > right
  )
}


test_window_wilcox <- function(df_raw, species, region, window_res) {

  if (is.null(window_res)) {
    warning("No window defined")
    return(NULL)
  }

  start <- window_res$start_month
  end   <- window_res$end_month
  wrap  <- window_res$wrap_around

  df_sub <- df_raw %>%
    dplyr::filter(
      species == !!species,
      region == !!region
    )

  # ---- NEW: drop all-zero or no-variation cases ----
  if (nrow(df_sub) == 0 || all(df_sub$logConc == 0, na.rm = TRUE)) {
    warning("All zero values — skipping")
    return(NULL)
  }

  if (sd(df_sub$logConc, na.rm = TRUE) == 0) {
    warning("No variation — skipping")
    return(NULL)
  }
  if (nrow(df_sub) == 0) {
    warning("No data for this species/region")
    return(NULL)
  }

  # define window membership
  if (!wrap) {
    in_window <- df_sub$month >= start & df_sub$month <= end
  } else {
    in_window <- df_sub$month >= start | df_sub$month <= end
  }

  df_sub <- df_sub %>%
    dplyr::mutate(window = ifelse(in_window, "in", "out"))

  if (length(unique(df_sub$window)) < 2) {
    warning("Only one group present (all in or all out)")
    return(NULL)
  }

  # Wilcoxon test
  w <- wilcox.test(logConc ~ window, data = df_sub, exact = FALSE)

  # group summaries
  summary <- df_sub %>%
    dplyr::group_by(window) %>%
    dplyr::summarise(
      n = dplyr::n(),
      median = median(logConc, na.rm = TRUE),
      mean = mean(logConc, na.rm = TRUE),
      .groups = "drop"
    )

  # ---------------------------
  # EFFECT SIZE (rank-biserial)
  # ---------------------------
  n_in  <- sum(df_sub$window == "in")
  n_out <- sum(df_sub$window == "out")

  w_stat <- as.numeric(w$statistic)

  rank_biserial <- (2 * w_stat) / (n_in * n_out) - 1

  prob_superiority <- (rank_biserial + 1) / 2

  list(
    summary = summary,
    test = data.frame(
      statistic = w_stat,
      p_value = w$p.value,
      rank_biserial = rank_biserial,
      prob_superiority = prob_superiority
    )
  )
}



all_species <- sort(unique(dfRawClean$species))
all_regions <- c("MAG", "PEI", "HAL", "BOF", "GOM")


collect_window_wilcox_results <- function(df_monthly, df_raw, threshold = 0.75) {

  all_species <- sort(unique(df_monthly$species))
  all_regions <- c("MAG", "PEI", "HAL", "BOF", "GOM")

  combos <- expand.grid(
    species = all_species,
    region = all_regions,
    stringsAsFactors = FALSE
  )

  results <- purrr::pmap_dfr(combos, function(species, region) {

    # -------------------------
    # monthly data for window calculation
    # -------------------------
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
    # window_res <- calc_window_optimized(
    #   plot_df,
    #   optimize = "mean_diff",
    #   min_width = 2,
    #   max_width = 6,
    #   must_include_peak = TRUE,
    #   tie_break = "wider"
    # )
    # window_res <- calc_window_plateaus(
    #   plot_df,
    #   min_width = 1,
    #   max_width = 11,
    #   width_weight = "none",
    #   edge_scale = "inside_high"
    # )
    if (is.null(window_res)) {
      print("window is null")
      return(NULL)
    }

    # -------------------------
    # raw data for Wilcoxon test
    # -------------------------
    wilcox_result <- test_window_wilcox(
      df_raw = df_raw,
      species = species,
      region = region,
      window_res = window_res
    )

    inside_months <- if (!window_res$wrap_around) {
      window_res$start_month:window_res$end_month
    } else {
      c(window_res$start_month:12, 1:window_res$end_month)
    }

    n_above_threshold_outside_window <- sum(
      plot_df$value >= threshold & !plot_df$month %in% inside_months,
      na.rm = TRUE
    )

    if (is.null(wilcox_result)) {
      print(paste0("wilcox is null at ", threshold, " for ", species, " in ", region))
      return(NULL)
    }

    summary_wide <- wilcox_result$summary %>%
      tidyr::pivot_wider(
        names_from = window,
        values_from = c(n, median, mean),
        names_glue = "{.value}_{window}"
      )

    dplyr::bind_cols(
      tibble::tibble(
        species = species,
        region = region,
        start_month = window_res$start_month,
        end_month = window_res$end_month,
        wrap_around = window_res$wrap_around,
        n_above_threshold_outside_window = n_above_threshold_outside_window
      ),
      summary_wide,
      tibble::as_tibble(wilcox_result$test)
    )
  })

  results
}




valid_combos_df <- dfRawClean %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    has_data = any(!is.na(logConc) & logConc > 0),
    .groups = "drop"
  ) %>%
  dplyr::filter(has_data)

total_valid_combos <- nrow(valid_combos_df)

threshold_avg_from_values_df <- threshold_values_df %>%
  dplyr::group_by(threshold) %>%
  dplyr::summarise(
    mean_prob_superiority = mean(prob_superiority, na.rm = TRUE),
    sd_prob_superiority   = sd(prob_superiority, na.rm = TRUE),
    n_valid = sum(!is.na(prob_superiority)),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    ymin = mean_prob_superiority - sd_prob_superiority,
    ymax = mean_prob_superiority + sd_prob_superiority,
    prop_valid = n_valid / total_valid_combos
  )







































threshold_values_df <- purrr::map_dfr(
  seq(0.50, 0.99, by = 0.01),
  function(thresh) {

    wilcox_results_df <- collect_window_wilcox_results(
      df_monthly = df,
      df_raw = dfRawClean,
      threshold = thresh
    )

    wilcox_results_df %>%
      dplyr::mutate(
        threshold = thresh,
        window_duration = dplyr::if_else(
          wrap_around,
          12 - start_month + end_month + 1,
          end_month - start_month + 1
        )
      ) %>%
      dplyr::select(
        threshold,
        species,
        region,
        prob_superiority,
        start_month,
        end_month,
        wrap_around,
        window_duration,
        n_above_threshold_outside_window
      )
  }
)

threshold_avg_from_values_df <- threshold_values_df %>%
  dplyr::group_by(threshold) %>%
  dplyr::summarise(
    mean_prob_superiority = mean(prob_superiority, na.rm = TRUE),
    sd_prob_superiority   = sd(prob_superiority, na.rm = TRUE),
    n_valid = sum(!is.na(prob_superiority)),

    mean_window_duration = mean(window_duration, na.rm = TRUE),
    sd_window_duration   = sd(window_duration, na.rm = TRUE),

    mean_above_outside = mean(n_above_threshold_outside_window, na.rm = TRUE),

    .groups = "drop"
  ) %>%
  dplyr::mutate(
    ymin = mean_prob_superiority - sd_prob_superiority,
    ymax = mean_prob_superiority + sd_prob_superiority,

    prop_valid = n_valid / total_valid_combos,

    dur_ymin = (mean_window_duration - sd_window_duration) / 12,
    dur_ymax = (mean_window_duration + sd_window_duration) / 12,

    above_outside_scaled = mean_above_outside / 12
  )

p <- ggplot(threshold_avg_from_values_df, aes(x = threshold)) +

  # --- prob superiority ribbon ---
  geom_ribbon(
    aes(
      ymin = ymin,
      ymax = ymax,
      fill = "Prob. superiority (± SD)"
    ),
    alpha = 0.25
  ) +

  geom_line(
    aes(
      y = mean_prob_superiority,
      color = "Prob. superiority (mean)"
    ),
    linewidth = 1
  ) +
  geom_point(
    aes(
      y = mean_prob_superiority,
      color = "Prob. superiority (mean)"
    ),
    size = 2
  ) +

  # --- above-threshold outside window ---
  geom_line(
    aes(
      y = above_outside_scaled,
      color = "Above-threshold outside window"
    ),
    linetype = "dotdash",
    linewidth = 1
  ) +
  geom_point(
    aes(
      y = above_outside_scaled,
      color = "Above-threshold outside window"
    ),
    shape = 4,
    size = 2
  ) +

  # --- prop valid ---
  geom_line(
    aes(
      y = prop_valid,
      color = "Proportion valid"
    ),
    linetype = "dashed",
    linewidth = 1
  ) +
  geom_point(
    aes(
      y = prop_valid,
      color = "Proportion valid"
    ),
    shape = 17,
    size = 2
  ) +

  # --- duration ribbon ---
  geom_ribbon(
    aes(
      ymin = dur_ymin,
      ymax = dur_ymax,
      fill = "Window duration (± SD)"
    ),
    alpha = 0.15
  ) +

  geom_line(
    aes(
      y = mean_window_duration / 12,
      color = "Window duration (mean)"
    ),
    linetype = "dotted",
    linewidth = 1
  ) +
  geom_point(
    aes(
      y = mean_window_duration / 12,
      color = "Window duration (mean)"
    ),
    shape = 15,
    size = 2
  ) +

  scale_color_manual(
    values = c(
      "Prob. superiority (mean)" = "black",
      "Proportion valid" = "red",
      "Window duration (mean)" = "blue",
      "Above-threshold outside window" = "purple"
    ),
    name = "Lines"
  ) +

  scale_fill_manual(
    values = c(
      "Prob. superiority (± SD)" = "grey70",
      "Window duration (± SD)" = "blue"
    ),
    name = "Ribbons"
  ) +

  scale_x_continuous(
    breaks = seq(0.50, 0.99, by = 0.05),
    labels = scales::percent_format(accuracy = 1)
  ) +

  scale_y_continuous(
    name = "Probability / proportion",
    sec.axis = sec_axis(~ . * 12, name = "Months")
  ) +

  theme_classic() +
  labs(
    title = "Threshold performance, coverage, and window duration",
    x = "Window threshold"
  )
p
ggsave("manuscript_figures/figure_B1(Appendix B).png", p, width = 9, height = 9, dpi = 300)







































