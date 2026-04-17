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
    # window_res <- calc_window_plateaus(
    #   plot_df,
    #   min_width = 1,
    #   max_width = 11,
    #   width_weight = "none",
    #   edge_scale = "inside_high"
    # )
    if (is.null(window_res)) {
      return(NULL)
    }

    peak_month <- plot_df$month[which.max(plot_df$value)]

    tibble::tibble(
      species = species,
      region = region,
      start_month = window_res$start_month,
      end_month = window_res$end_month,
      wrap_around = window_res$wrap_around,
      peak_month = peak_month
    )
  })
}

###################################
# FIXED MONTH SHIFT: MARCH -> NEXT MARCH
###################################

shift_month_march <- function(m) {
  ifelse(m < 3, m + 12, m)
}

###################################
# PREP WINDOW SEGMENTS FOR PLOTTING
###################################

prep_window_segments <- function(window_df) {

  window_df %>%
    dplyr::mutate(
      start_plot = shift_month_march(start_month),
      end_plot   = shift_month_march(end_month),
      peak_plot  = shift_month_march(peak_month)
    ) %>%
    dplyr::mutate(
      end_plot = ifelse(end_plot < start_plot, end_plot + 12, end_plot),
      peak_plot = ifelse(peak_plot < start_plot, peak_plot + 12, peak_plot)
    ) %>%
    dplyr::mutate(
      start_plot = ifelse(
        species == "Ciona intestinalis" & region == "BOF",
        3,
        start_plot
      ),
      end_plot = ifelse(
        species == "Ciona intestinalis" & region == "BOF",
        15,
        end_plot
      ),
      peak_plot = ifelse(
        species == "Ciona intestinalis" & region == "BOF",
        shift_month_march(peak_month),
        peak_plot
      )
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

window_plot_df <- collect_window_plot_data(
  df_monthly = df,
  threshold = 0.8
) %>%
  dplyr::filter(
    !(species == "Didemnum vexillum" & region == "MAG")
  ) %>%
  dplyr::mutate(
    species = factor(species, levels = species_order),
    region = factor(region, levels = rev(region_order))
  )

window_plot_df2 <- prep_window_segments(window_plot_df)

x_min <- 3
x_max <- 15

###################################
# PLOT
###################################

ggplot(window_plot_df2) +
  geom_segment(
    aes(
      x = start_plot,
      xend = end_plot,
      y = region,
      yend = region,
      color = region
    ),
    linewidth = 4,
    lineend = "round"
  ) +
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
    limits = c(x_min, x_max),
    breaks = seq(x_min, x_max, by = 1),
    labels = shifted_month_labels,
    expand = expansion(mult = c(0.02, 0.02))
  ) +
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
  )






























##################################
#FRIEDMAN
##################################


circular_midpoint <- function(start, end, n = 12) {
  d <- ((end - start + n) %% n)
  ((start + d / 2 - 1) %% n) + 1
}


window_centers_df <- window_plot_df %>%
  dplyr::mutate(
    center_month = circular_midpoint(start_month, end_month)
  )


center_matrix <- window_centers_df %>%
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



rank_df <- window_centers_df %>%
  dplyr::mutate(
    center_shifted = shift_month_march(center_month)
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


rank_df_complete <- rank_df %>%
  group_by(species) %>%
  filter(n() == 5) %>%   # or whatever full region count is
  ungroup()


region_complete_rank_summary <- rank_df_complete %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    mean_rank = mean(rank, na.rm = TRUE),
    sd_rank = sd(rank, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(mean_rank)


















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
