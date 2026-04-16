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
# CHOOSE A GLOBAL CUTOFF MONTH
###################################
# We want to "cut" the circular year at the month that splits
# the fewest windows, so most windows stay continuous.

choose_cut_month <- function(window_df) {

  month_is_inside_window <- function(cut, start, end, wrap) {
    if (!wrap) {
      cut >= start && cut <= end
    } else {
      cut >= start || cut <= end
    }
  }

  cut_scores <- sapply(1:12, function(cut) {
    sum(
      purrr::pmap_lgl(
        list(window_df$start_month,
             window_df$end_month,
             window_df$wrap_around),
        ~ month_is_inside_window(cut, ..1, ..2, ..3)
      ),
      na.rm = TRUE
    )
  })

  which.min(cut_scores)
}

###################################
# SHIFT MONTHS TO A NON-WRAPPING AXIS
###################################

shift_month <- function(m, cut_month) {
  ifelse(m <= cut_month, m + 12, m)
}

prep_window_segments <- function(window_df, cut_month) {

  out <- window_df %>%
    mutate(
      start_plot = shift_month(start_month, cut_month),
      end_plot   = shift_month(end_month, cut_month),
      peak_plot  = shift_month(peak_month, cut_month)
    ) %>%
    mutate(
      # if end is left of start after shifting, push it forward one cycle
      end_plot = ifelse(end_plot < start_plot, end_plot + 12, end_plot),
      peak_plot = ifelse(peak_plot < start_plot, peak_plot + 12, peak_plot)
    )

  out
}

###################################
# LABEL FUNCTION FOR SHIFTED AXIS
###################################

shifted_month_labels <- function(x) {
  # map to month index 1–12
  m <- ((round(x) - 1) %% 12) + 1

  # single-letter labels (no collapsing of positions)
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
  threshold = 0.75
) %>%
  filter(
      !(species == "Didemnum vexillum" & region == "MAG")
  ) %>%
  mutate(
    species = factor(species, levels = species_order),
    region = factor(region, levels = rev(region_order))
  )

cut_month <- choose_cut_month(window_plot_df)

window_plot_df2 <- prep_window_segments(window_plot_df, cut_month)

# dynamic axis range
x_min <- cut_month + 1
x_max <- max(window_plot_df2$end_plot, na.rm = TRUE)

###################################
# PLOT
###################################

ggplot(window_plot_df2) +
  geom_segment(
    data = window_plot_df2 %>%
      filter(!(species == "Ciona intestinalis" & region == "BOF")),
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
    breaks = seq(x_min, floor(x_max), by = 1),
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
  threshold = 0.75
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
