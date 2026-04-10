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

plot_monthly_wheel <- function(df, title = NULL) {

  ggplot2::ggplot(df) +

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
        above = "#5B2C83",
        below = "grey70"
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
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_blank()
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
#

if (!dir.exists("abundance_wheels")) {
  dir.create("abundance_wheels")
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
      plot_df,
      title = paste(sp, "-", reg)
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




###################################
#CALC WINDOW
###################################


calc_window_simple <- function(plot_df, threshold = 0.75) {

  df <- plot_df %>%
    dplyr::arrange(month) %>%
    dplyr::mutate(above = value >= threshold)

  x <- df$above

  peak_month <- which.max(df$value)

  n <- length(x)

  expand_one_side <- function(start, step) {

    i <- start
    repeat {

      next_i <- ((i - 1 + step) %% n) + 1

      if (!x[next_i]) break
      if (next_i == peak_month) break

      i <- next_i

      if (i == start) break
    }

    return(i)
  }

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

  start_month <- left
  end_month <- right

  list(
    start_month = start_month,
    end_month = end_month,
    wrap_around = start_month > end_month
  )
}


df_sub <- df %>%
  dplyr::filter(
    species == "Membranipora membranacea",
    region == "PEI"
  )

plot_df <- df_sub %>%
  prep_monthly_signal() %>%
  interp_monthly_circular() %>%
  classify_months(threshold = 0.75)




calc_window_simple(plot_df, threshold = 0.75)
