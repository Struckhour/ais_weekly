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





selected_species <- "Membranipora membranacea"
selected_region <- "PEI"

df_sub <- df %>%
  dplyr::filter(
    species == selected_species,
    region == selected_region
  )

plot_df <- df_sub %>%
  prep_monthly_signal() %>%
  interp_monthly_circular() %>%
  classify_months(threshold = 0.75)

window_res <- calc_window_simple(plot_df, threshold = 0.75)












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

    if (is.null(window_res)) {
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
  threshold = 0.75
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
      (species == "Didemnum vexillum" & region == "MAG") |
       (species == "Ciona intestinalis" & region == "BOF")
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

files <- list.files("abundance_wheels", pattern = "\\.png$", full.names = TRUE)
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

image_write(final_labeled, "saved_figures/combined_monthly_wheels.png")



