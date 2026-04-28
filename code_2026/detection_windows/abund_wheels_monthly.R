source('./AIS_eDNA_data_prep.R')


df <- dfMonths %>% select(region, species, month_original, normMeanLogConc)


prep_monthly_signal <- function(df) {

  df %>%
    dplyr::group_by(month_original) %>%
    dplyr::summarise(
      value = mean(normMeanLogConc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::rename(month = month_original) %>%
    dplyr::arrange(month)
}


interp_monthly_circular <- function(df) {

  full <- data.frame(month = 1:12)

  x <- full %>%
    dplyr::left_join(df, by = "month") %>%
    dplyr::arrange(month)

  # track missingness BEFORE interpolation
  x$is_observed <- !is.na(x$value)

  # no data case
  if (all(is.na(x$value))) {
    x$value <- NA
    x$is_imputed <- FALSE
    return(x)
  }

  # single observed point → flat line, but all imputed except original
  if (sum(x$is_observed) == 1) {
    val <- x$value[x$is_observed]
    x$value <- rep(val, 12)
    x$is_imputed <- !x$is_observed
    return(x)
  }

  idx <- which(x$is_observed)
  y <- x$value[idx]

  # circular wrap
  idx_ext <- c(idx - 12, idx, idx + 12)
  y_ext <- c(y, y, y)

  interp <- approx(
    x = idx_ext,
    y = y_ext,
    xout = 1:12,
    method = "linear",
    rule = 2
  )$y

  x$value <- interp

  # imputed = originally missing
  x$is_imputed <- !x$is_observed

  x
}

classify_months <- function(df, threshold = 0.75) {

  df %>%
    dplyr::mutate(
      status = dplyr::case_when(
        value >= threshold ~ "above",
        TRUE ~ "below"
      )
    )
}

plot_monthly_wheel <- function(df, region, title = NULL) {

  region_colors <- c(
    MAG = "#00A08A",
    PEI = "#446455",
    HAL = "#CCAA4F",
    BOF = "#5BBCD6",
    GOM = "#fb8072"
  )

  if (!region %in% names(region_colors)) {
    stop("Unknown region: ", region)
  }

  above_color <- region_colors[[region]]

  ggplot2::ggplot(df) +

    # -------------------------
  # RADIAL LINES (months)
  # -------------------------
  ggplot2::geom_vline(
    xintercept = 1:12,
    color = "black",
    linewidth = 0.3,
    alpha = 0.5
  ) +

    # -------------------------
  # CIRCULAR GRID LINES
  # -------------------------
  ggplot2::geom_hline(
    yintercept = c(0.25, 0.5, 0.75, 1),
    color = "black",
    linewidth = 0.3,
    alpha = 0.5
  ) +

    # -------------------------
  # BASE LAYER (all bars)
  # -------------------------
  ggplot2::geom_col(
    ggplot2::aes(
      x = month,
      y = value,
      fill = status
    ),
    width = 0.9
  ) +

    ggplot2::scale_fill_manual(
      values = c(
        above = above_color,
        below = "grey70"
      ),
      labels = c(
        above = "Above threshold",
        below = "Below threshold"
      )
    ) +

    # -------------------------
  # STRIPED OVERLAY (IMPUTED ONLY)
  # -------------------------
  ggpattern::geom_col_pattern(
    data = dplyr::filter(df, is_imputed == TRUE),
    ggplot2::aes(x = month, y = value),
    width = 0.9,
    fill = NA,
    pattern = "stripe",
    pattern_fill = "white",
    pattern_color = "white",
    pattern_density = 0.1,
    pattern_spacing = 0.02
  ) +

    ggplot2::coord_polar(start = 0) +

    ggplot2::scale_x_continuous(
      limits = c(0.5, 12.5),
      breaks = 1:12,
      labels = month.abb
    ) +

    ggplot2::ylim(0, 1) +

    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 32),
      legend.title = ggplot2::element_blank(),
      legend.position = "none"
    ) +

    ggplot2::labs(title = title)
}


plot_monthly_wheel_window <- function(df, region, title = NULL) {

  region_colors <- c(
    MAG = "#00A08A",
    PEI = "#446455",
    HAL = "#CCAA4F",
    BOF = "#5BBCD6",
    GOM = "#fb8072"
  )

  if (!region %in% names(region_colors)) {
    stop("Unknown region: ", region)
  }

  above_color <- region_colors[[region]]

  df <- df %>%
    dplyr::arrange(month)
  win <- calc_window_simple(df, threshold = 0.9)
  # calculate optimized window
  # win <- calc_window_optimized(
  #   df,
  #   value_col = "value",
  #   optimize = "mean_diff",
  #   min_width = 2,
  #   max_width = 6,
  #   must_include_peak = TRUE,
  #   tie_break = "wider"
  # )
  # win <- calc_window_plateaus(
  #   plot_df,
  #   min_width = 1,
  #   max_width = 11,
  #   width_weight = "none",
  #   edge_scale = "inside_high"
  # )

  # assign wedge status based on selected window
  df <- df %>%
    dplyr::mutate(
      status = dplyr::case_when(
        # inside window (same as before)
        !win$wrap_around &
          month >= win$start_month &
          month <= win$end_month ~ "inside",

        win$wrap_around &
          (month >= win$start_month | month <= win$end_month) ~ "inside",

        # NEW: outside but above threshold
        value >= 0.9 ~ "above_threshold",

        # everything else
        TRUE ~ "outside"
      )
    )

  ggplot2::ggplot(df) +

    # -------------------------
  # RADIAL LINES (months)
  # -------------------------
  ggplot2::geom_vline(
    xintercept = 1:12,
    color = "black",
    linewidth = 0.3,
    alpha = 0.5
  ) +

    # -------------------------
  # CIRCULAR GRID LINES
  # -------------------------
  ggplot2::geom_hline(
    yintercept = c(0.25, 0.5, 0.75, 1),
    color = "black",
    linewidth = 0.3,
    alpha = 0.5
  ) +

    # -------------------------
  # BASE LAYER (all bars)
  # -------------------------
  ggplot2::geom_col(
    ggplot2::aes(
      x = month,
      y = value,
      fill = status
    ),
    width = 0.9
  ) +

    ggplot2::scale_fill_manual(
      values = c(
        inside = above_color,
        above_threshold = "grey40",   # darker grey
        outside = "grey70"            # lighter grey
      ),
      labels = c(
        inside = "Optimal window",
        above_threshold = "Above threshold (outside window)",
        outside = "Below threshold"
      )
    ) +

    # -------------------------
  # STRIPED OVERLAY (IMPUTED ONLY)
  # -------------------------
  ggpattern::geom_col_pattern(
    data = dplyr::filter(df, is_imputed == TRUE),
    ggplot2::aes(x = month, y = value),
    width = 0.9,
    fill = NA,
    pattern = "stripe",
    pattern_fill = "white",
    pattern_color = "white",
    pattern_density = 0.1,
    pattern_spacing = 0.02
  ) +

    ggplot2::coord_polar(start = 0) +

    ggplot2::scale_x_continuous(
      limits = c(0.5, 12.5),
      breaks = 1:12,
      labels = month.abb
    ) +

    ggplot2::ylim(0, 1) +

    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(size = 32),
      legend.title = ggplot2::element_blank(),
      legend.position = "none"
    ) +

    ggplot2::labs(title = title)
}


# df_sub <- df %>%
#   dplyr::filter(
#     species == "Membranipora membranacea",
#     region == "MAG"
#   )
#
# plot_df <- df_sub %>%
#   prep_monthly_signal() %>%
#   interp_monthly_circular() %>%
#   classify_months(threshold = 0.75)
#
# plot_monthly_wheel(plot_df, "Membranipora membranacea – MAG")


if (!dir.exists("abundance_wheels_90_dark")) {
  dir.create("abundance_wheels_90_dark")
}


all_species <- unique(df$species)
all_regions <- unique(df$region)

for (sp in all_species) {
  for (reg in all_regions) {

    df_sub <- df %>%
      dplyr::filter(
        species == sp,
        region == reg
      )

    # skip empty combos
    if (nrow(df_sub) == 0) next

    # build plot data
    plot_df <- df_sub %>%
      prep_monthly_signal() %>%
      interp_monthly_circular() %>%
      classify_months(threshold = 0.75)

    # generate plot
    p <- plot_monthly_wheel(
      df = plot_df,
      region = reg
      # title = paste(sp, "-", reg)
    )

    # safe filename
    file_name <- paste0(
      "abundance_wheels/",
      gsub(" ", "_", sp),
      "__",
      gsub(" ", "_", reg),
      ".png"
    )

    # save
    ggsave(
      filename = file_name,
      plot = p,
      width = 6,
      height = 6,
      dpi = 300
    )
  }
}

##################################
#save window wheels instead of threshold wheels
##################################


if (!dir.exists("abundance_plateau_wheels")) {
  dir.create("abundance_plateau_wheels")
}


all_species <- unique(df$species)
all_regions <- unique(df$region)

for (sp in all_species) {
  for (reg in all_regions) {

    df_sub <- df %>%
      dplyr::filter(
        species == sp,
        region == reg
      )

    # skip empty combos
    if (nrow(df_sub) == 0) next

    # build plot data
    plot_df <- df_sub %>%
      prep_monthly_signal() %>%
      interp_monthly_circular() %>%
      classify_months(threshold = 0.90)

    # generate plot
    p <- plot_monthly_wheel_window(
      df = plot_df,
      region = reg
      # title = paste(sp, "-", reg)
    )

    # safe filename
    file_name <- paste0(
      "abundance_wheels_90_dark/",
      gsub(" ", "_", sp),
      "__",
      gsub(" ", "_", reg),
      ".png"
    )

    # save
    ggsave(
      filename = file_name,
      plot = p,
      width = 6,
      height = 6,
      dpi = 300
    )
  }
}


###################################
#CALC WINDOW
###################################


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


calc_window_optimized <- function(
    plot_df,
    value_col = "value",
    optimize = c("mean_diff", "prob_superiority"),
    min_width = 1,
    max_width = NULL,
    must_include_peak = TRUE,
    tie_break = c("wider", "narrower")
) {

  optimize <- match.arg(optimize)
  tie_break <- match.arg(tie_break)

  df <- plot_df %>%
    dplyr::arrange(month)

  x <- df[[value_col]]
  n <- length(x)

  if (is.null(max_width)) {
    max_width <- n - 1
  }

  if (n < 2) {
    stop("plot_df must contain at least 2 months.")
  }

  if (min_width < 1 || max_width >= n || min_width > max_width) {
    stop("Need 1 <= min_width <= max_width < number of months.")
  }

  peak_month <- which.max(x)

  # helper: circular sequence of indices
  get_window_idx <- function(start, width, n) {
    ((start - 1 + 0:(width - 1)) %% n) + 1
  }

  # helper: probability of superiority
  calc_prob_superiority <- function(in_vals, out_vals) {
    comps <- outer(in_vals, out_vals, FUN = "-")
    mean((comps > 0) + 0.5 * (comps == 0))
  }

  candidates <- list()
  k <- 1

  for (width in min_width:max_width) {
    for (start in 1:n) {

      idx_in <- get_window_idx(start, width, n)
      idx_out <- setdiff(seq_len(n), idx_in)

      if (must_include_peak && !(peak_month %in% idx_in)) {
        next
      }

      in_vals <- x[idx_in]
      out_vals <- x[idx_out]

      mean_in <- mean(in_vals, na.rm = TRUE)
      mean_out <- mean(out_vals, na.rm = TRUE)
      mean_diff <- mean_in - mean_out

      prob_superiority <- calc_prob_superiority(in_vals, out_vals)
      rank_biserial <- 2 * prob_superiority - 1

      score <- switch(
        optimize,
        mean_diff = mean_diff,
        prob_superiority = prob_superiority
      )

      end_month <- idx_in[length(idx_in)]

      candidates[[k]] <- data.frame(
        start_month = start,
        end_month = end_month,
        width = width,
        wrap_around = start > end_month,
        mean_in = mean_in,
        mean_out = mean_out,
        mean_diff = mean_diff,
        prob_superiority = prob_superiority,
        rank_biserial = rank_biserial,
        score = score
      )

      k <- k + 1
    }
  }

  candidate_table <- dplyr::bind_rows(candidates)

  # sort candidates
  candidate_table <- candidate_table %>%
    dplyr::arrange(
      dplyr::desc(score),
      if (tie_break == "wider") dplyr::desc(width) else width,
      start_month
    )

  best <- candidate_table[1, ]

  list(
    start_month = best$start_month,
    end_month = best$end_month,
    wrap_around = best$wrap_around,
    peak_month = peak_month,
    width = best$width,
    mean_in = best$mean_in,
    mean_out = best$mean_out,
    mean_diff = best$mean_diff,
    prob_superiority = best$prob_superiority,
    rank_biserial = best$rank_biserial,
    score = best$score,
    optimize = optimize,
    candidates = candidate_table
  )
}

calc_window_plateaus <- function(
    plot_df,
    value_col = "value",
    min_width = 2,
    max_width = NULL,
    must_include_peak = TRUE,
    width_weight = c("none", "sqrt", "linear"),
    tie_break = c("wider", "narrower"),
    edge_scale = c("none", "inside_high", "inside_high_sq")
) {

  width_weight <- match.arg(width_weight)
  tie_break <- match.arg(tie_break)
  edge_scale <- match.arg(edge_scale)

  df <- plot_df %>%
    dplyr::arrange(month)

  x <- df[[value_col]]
  n <- length(x)

  if (n < 2) {
    stop("plot_df must contain at least 2 rows.")
  }

  if (is.null(max_width)) {
    max_width <- n - 1
  }

  if (min_width < 1 || max_width >= n || min_width > max_width) {
    stop("Need 1 <= min_width <= max_width < number of rows.")
  }

  if (any(!is.finite(x))) {
    stop("All values in value_col must be finite.")
  }

  peak_month <- which.max(x)

  get_window_idx <- function(start, width, n) {
    ((start - 1 + 0:(width - 1)) %% n) + 1
  }

  prev_idx <- function(i, n) {
    ((i - 2) %% n) + 1
  }

  next_idx <- function(i, n) {
    (i %% n) + 1
  }

  width_multiplier <- function(width, method) {
    switch(
      method,
      none = 1,
      sqrt = sqrt(width),
      linear = width
    )
  }

  scale_edge <- function(drop, inside_val, method) {
    raw_drop <- max(0, drop)

    switch(
      method,
      none = raw_drop,
      inside_high = raw_drop * inside_val,
      inside_high_sq = raw_drop * inside_val^2
    )
  }

  candidates <- list()
  k <- 1

  for (width in min_width:max_width) {
    for (start in 1:n) {

      idx_in <- get_window_idx(start, width, n)
      end_month <- idx_in[length(idx_in)]

      if (must_include_peak && !(peak_month %in% idx_in)) {
        next
      }

      left_outside  <- prev_idx(start, n)
      right_outside <- next_idx(end_month, n)

      left_drop_raw  <- x[start] - x[left_outside]
      right_drop_raw <- x[end_month] - x[right_outside]

      left_drop  <- scale_edge(left_drop_raw,  x[start],     edge_scale)
      right_drop <- scale_edge(right_drop_raw, x[end_month], edge_scale)

      drop_sum <- left_drop + right_drop
      score <- drop_sum * width_multiplier(width, width_weight)

      candidates[[k]] <- data.frame(
        start_month = start,
        end_month = end_month,
        wrap_around = start > end_month,
        width = width,
        peak_month = peak_month,
        left_inside_value = x[start],
        left_outside_month = left_outside,
        left_outside_value = x[left_outside],
        left_drop_raw = left_drop_raw,
        left_drop = left_drop,
        right_inside_value = x[end_month],
        right_outside_month = right_outside,
        right_outside_value = x[right_outside],
        right_drop_raw = right_drop_raw,
        right_drop = right_drop,
        drop_sum = drop_sum,
        score = score
      )

      k <- k + 1
    }
  }

  candidate_table <- dplyr::bind_rows(candidates)

  candidate_table <- candidate_table %>%
    dplyr::arrange(
      dplyr::desc(score),
      dplyr::desc(drop_sum),
      if (tie_break == "wider") dplyr::desc(width) else width,
      start_month
    )

  best <- candidate_table[1, ]

  list(
    start_month = best$start_month,
    end_month = best$end_month,
    wrap_around = best$wrap_around,
    width = best$width,
    peak_month = best$peak_month,
    left_drop = best$left_drop,
    right_drop = best$right_drop,
    drop_sum = best$drop_sum,
    score = best$score,
    candidates = candidate_table
  )
}


selected_species <- "Botrylloides violaceus"
selected_region <- "MAG"

df_sub <- df %>%
  dplyr::filter(
    species == selected_species,
    region == selected_region
  )

plot_df <- df_sub %>%
  prep_monthly_signal() %>%
  interp_monthly_circular() %>%
  classify_months(threshold = 0.75)

window_res <- calc_window_simple(plot_df, threshold = 0.85)

test_window_wilcox(df, selected_species, selected_region, window_res)










##############################
#WILCOXON TEST BETWEEN OPTIMAL WINDOW OF MONTHS AND SUBOPTIMAL MONTHS
##############################

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


wilcox_result <- test_window_wilcox(dfRawClean, selected_species, selected_region, window_res)
wilcox_result

library(dplyr)
library(purrr)
library(tidyr)

all_species <- sort(unique(dfRawClean$species))
all_regions <- c("MAG", "PEI", "HAL", "BOF", "GOM")
library(dplyr)
library(purrr)
library(tidyr)
library(tibble)

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

    if (is.null(wilcox_result)) {
      print("wilcox is null")
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
        wrap_around = window_res$wrap_around
      ),
      summary_wide,
      tibble::as_tibble(wilcox_result$test)
    )
  })

  results
}

wilcox_results_df <- collect_window_wilcox_results(
  df_monthly = df,
  df_raw = dfRawClean,
  threshold = 0.8
)


region_colors <- c(
  MAG = "#00A08A",
  PEI = "#446455",
  HAL = "#CCAA4F",
  BOF = "#5BBCD6",
  GOM = "#fb8072"
)

region_order <- c("MAG", "PEI", "HAL", "BOF", "GOM")

species_order <- c(
  "Membranipora membranacea",
  "Botrylloides violaceus",
  "Didemnum vexillum",
  "Ciona intestinalis",
  "Carcinus maenas"
)

plot_df <- wilcox_results_df %>%
  filter(
    !(
      (species == "Didemnum vexillum" & region == "MAG")
    )
  ) %>%
  mutate(
    region = factor(region, levels = rev(region_order)),
    species = factor(species, levels = species_order)
  ) %>%
  arrange(species, region)

ggplot(plot_df, aes(x = prob_superiority, y = region, fill = region)) +
  geom_col(width = 0.8) +
  geom_vline(xintercept = 0.5, linetype = "dashed") +
  geom_vline(xintercept = 0.25, linetype = "dashed") +
  geom_vline(xintercept = 0.75, linetype = "dashed") +
  geom_vline(xintercept = 1.0, linetype = "dashed") +
  scale_fill_manual(values = region_colors) +
  scale_x_continuous(limits = c(0, 1)) +
  facet_grid(species ~ ., scales = "free_y", space = "free_y", switch = "y") +
  labs(
    x = "Probability of superiority",
    y = NULL,
    title = "Probability of superiority (optimal period vs suboptimal period)"
  ) +
  theme_minimal() +
  theme(
    legend.title = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    strip.placement = "outside",
    strip.text.y.left = element_text(
      angle = 0,
      face = "italic",
      hjust = 1
    ),
    panel.spacing.y = unit(0.6, "lines")
  )
































##############################
#STITCH TOGETHER PNGS
##############################

library(magick)

regions <- c("MAG", "PEI", "HAL", "BOF", "GOM")

species <- c(
  "Membranipora membranacea",
  "Botrylloides violaceus",
  "Didemnum vexillum",
  "Ciona intestinalis",
  "Carcinus maenas"
)

species_file <- gsub(" ", "_", species)

# -------------------------
# panel/header dimensions
# -------------------------
panel_w <- 1800
panel_h <- 1800

row_header_w <- 500
col_header_h <- 500
header_text_size <- 140
pad <- 40  # <-- adjust this (20–80 usually looks good)
# -------------------------
# column labels (top)
# -------------------------
species_labels <- lapply(species, function(sp) {
  image_blank(width = panel_w, height = col_header_h, color = "white") |>
    image_annotate(
      text = sp,
      size = header_text_size,
      gravity = "center"
    ) |>
    image_border("white", paste0(pad, "x0"))
})

top_row <- image_append(image_join(species_labels))

# -------------------------
# row labels (left side)
# IMPORTANT:
# make canvas wide-and-short first,
# rotate, then force exact final size
# -------------------------
region_labels <- lapply(regions, function(rg) {
  image_blank(width = panel_h, height = row_header_w, color = "white") |>
    image_annotate(
      text = rg,
      size = header_text_size,
      gravity = "center",
      weight = 700
    ) |>
    image_rotate(-90) |>
    image_extent(
      geometry = paste0(row_header_w, "x", panel_h),
      gravity = "center",
      color = "white"
    )
})

left_col <- image_append(image_join(region_labels), stack = TRUE)

# -------------------------
# top-left corner
# -------------------------
corner <- image_blank(width = row_header_w, height = col_header_h, color = "white")
top_full <- image_append(c(corner, top_row))

# -------------------------
# expected file layout
# -------------------------
expected <- expand.grid(
  species = species_file,
  region = regions,
  stringsAsFactors = FALSE
)

expected$filename <- paste0(expected$species, "__", expected$region, ".png")

files <- list.files("abundance_wheels_90_dark", pattern = "\\.png$", full.names = TRUE)
file_map <- setNames(files, basename(files))

expected$path <- file_map[expected$filename]

# blank placeholder for missing panels
blank <- image_blank(width = panel_w, height = panel_h, color = "white")

# ensure correct ordering
expected$region <- factor(expected$region, levels = regions)
expected$species <- factor(expected$species, levels = species_file)

expected <- expected[order(expected$region, expected$species), ]

# -------------------------
# build rows
# -------------------------


rows <- lapply(split(expected, expected$region), function(df_row) {
  imgs <- lapply(df_row$path, function(p) {
    if (is.na(p)) {
      blank
    } else {
      image_read(p) |>
        image_resize(paste0(panel_w, "x", panel_h, "!"))
    } |>
      image_border("white", paste0(pad, "x0"))  # horizontal padding only
  })
  image_append(image_join(imgs))
})

final <- image_append(image_join(rows), stack = TRUE)

# -------------------------
# combine everything
# -------------------------
grid_with_labels <- image_append(c(left_col, final))
final_labeled <- image_append(c(top_full, grid_with_labels), stack = TRUE)

image_write(final_labeled, "saved_figures/combined_monthly_wheels_90_dark.png")



















