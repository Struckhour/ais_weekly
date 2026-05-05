View(plateAbundance)


library(tidyverse)
library(lubridate)





life_state_colors <- c("#1b9e77",
                       "#d95f02",
                       "#7570b3",
                       "#e7298a",
                       "#66a61e")

alt_colors <- c("#e41a1c",
                "#ff7f00",
                "#4daf4a",
                "#377eb8",
                "#984ea3"
                )

# 1. Define global state order explicitly
# 1. Explicit order (TOP → BOTTOM)
state_levels <- c("R", "G", "B", "S", "D")

state_labels <- c(
  "R" = "Recruitment",
  "G" = "Growth",
  "B" = "Breeding",
  "S" = "Senescence",
  "D" = "Dieoff"
)

# 2. Prepare data
plateStates <- plateAbundance %>%
  filter(date >= as.Date("2022-06-01") & date <= as.Date("2028-05-31")) %>%
  separate_rows(State, sep = ",\\s*") %>%
  mutate(
    State = trimws(State),
    State = factor(State, levels = state_levels),  # <-- KEY
    species = factor(species, levels = unique(species))
  ) %>%
  filter(!is.na(State)) %>%
  mutate(
    state_row = as.numeric(State)
  )

# 3. Plot
ggplot(plateStates, aes(x = date, y = state_row, fill = State)) +
  geom_tile(width = 8, height = 0.8) +
  scale_fill_manual(
    values = alt_colors,
    breaks = state_levels,      # <-- legend order fixed
    labels = state_labels
  ) +
  facet_grid(region ~ species) +
  scale_y_reverse(              # <-- flips so first level is TOP
    breaks = 1:length(state_levels),
    labels = state_labels[state_levels]
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  theme_classic() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.spacing = unit(0.8, "lines")
  ) +
  labs(
    title = "Timeline of State Changes Over Time by Species and Region",
    x = "Month",
    y = "State",
    fill = "State"
  )





























library(patchwork)

normalize_to_max <- function(x) {
  mx <- max(x, na.rm = TRUE)
  if (!is.finite(mx) || mx == 0) return(rep(NA_real_, length(x)))
  x / mx
}

# ----------------------------
# Settings
# ----------------------------
target_species_vec <- c("Botrylloides violaceus", "Ciona intestinalis")
target_region_vec  <- c("MAG", "PEI", "HAL", "BOF")

state_levels <- c("R", "G", "B", "S", "D")

state_labels <- c(
  "R" = "Recruitment",
  "G" = "Growth",
  "B" = "Breeding",
  "S" = "Senescence",
  "D" = "Dieoff"
)

alt_colors <- c(
  "#e41a1c",
  "#ff7f00",
  "#4daf4a",
  "#377eb8",
  "#984ea3"
)

# ----------------------------
# Global ranges across all included data
# ----------------------------
all_plate_dates <- plateAbundance %>%
  filter(
    species %in% target_species_vec,
    region %in% target_region_vec,
    date >= as.Date("2022-06-01"),
    date <= as.Date("2028-05-31")
  ) %>%
  pull(date)

all_qpcr <- dfRawClean %>%
  filter(
    species %in% target_species_vec,
    region %in% target_region_vec
  )

all_qpcr_dates <- all_qpcr$date

global_x_min <- min(c(all_plate_dates, all_qpcr_dates), na.rm = TRUE)
global_x_max <- max(c(all_plate_dates, all_qpcr_dates), na.rm = TRUE)

global_y_min <- max(1, min(all_qpcr$concentration + 1, na.rm = TRUE))
global_y_max <- max(all_qpcr$concentration + 1, na.rm = TRUE)

# ----------------------------
# Helper: top species strip
# ----------------------------
make_species_strip <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, size = 4) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
}

# ----------------------------
# Helper: right-side region strip
# ----------------------------
make_region_strip <- function(label) {
  ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = label, angle = -90, size = 4) +
    xlim(0, 1) +
    ylim(0, 1) +
    theme_void() +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
    )
}

# ----------------------------
# Function for one species-region panel
# ----------------------------
make_species_region_panel <- function(target_species, target_region,
                                      show_x_text = FALSE,
                                      show_left_y = TRUE,
                                      show_legend = FALSE) {

  # ---- Life-stage data ----
  plateStates_one <- plateAbundance %>%
    filter(
      species == target_species,
      region == target_region,
      date >= as.Date("2022-06-01"),
      date <= as.Date("2028-05-31")
    ) %>%
    separate_rows(State, sep = ",\\s*") %>%
    mutate(
      State = trimws(State),
      State = factor(State, levels = state_levels),
      state_row = as.numeric(State)
    ) %>%
    filter(!is.na(State))

  # ---- qPCR data ----

  qpcr_one <- dfRawClean %>%
    filter(
      species == target_species,
      region == target_region
    ) %>%
    group_by(date) %>%
    summarise(
      value = mean(logConc, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      value_norm = normalize_to_max(value),
      source = "qPCR logConc"
    )

  # ---- Plate abundance data ----

  plate_one <- plateAbundance %>%
    filter(
      species == target_species,
      region == target_region,
      date >= as.Date("2022-06-01"),
      date <= as.Date("2028-05-31")
    ) %>%
    group_by(date) %>%
    summarise(
      value = mean(Avg, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      value_norm = normalize_to_max(value),
      source = "Plate Avg"
    )

  # ---- Combined abundance data ----
  abundance_one <- bind_rows(qpcr_one, plate_one)

  # ---- If both are empty, return blank ----
  if (nrow(plateStates_one) == 0 && nrow(qpcr_one) == 0) {
    return(
      plot_spacer() / plot_spacer() +
        plot_layout(heights = c(3, 2))
    )
  }



  # ---- Top plot: qPCR + plate abundance ----
  p_qpcr <- ggplot(abundance_one, aes(x = date, y = value_norm, color = source)) +
    geom_point(size = 1.2, alpha = 0.6) +
    geom_line(alpha = 0.8) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_x_date(
      limits = c(global_x_min, global_x_max),
      date_breaks = "1 month",
      date_labels = "%b"
    ) +
    theme_classic() +
    labs(
      x = NULL,
      y = if (show_left_y) "Relative\nabundance" else NULL,
      color = "Source"
    ) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.y = if (show_left_y) element_text() else element_blank(),
      axis.text.y  = if (show_left_y) element_text() else element_blank(),
      axis.ticks.y = if (show_left_y) element_line() else element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      legend.position = if (show_legend) "right" else "none"
    )

  # ---- Bottom plot: life stages
  # PEI Ciona gets a blank life-stage panel but keeps qPCR
  if (target_region == "PEI" && target_species == "Ciona intestinalis") {
    p_stage <- ggplot() +
      scale_y_reverse(
        breaks = 1:length(state_levels),
        labels = state_labels[state_levels],
        limits = c(length(state_levels) + 0.5, 0.5)
      ) +
      scale_x_date(
        limits = c(global_x_min, global_x_max),
        date_breaks = "1 month",
        date_labels = "%b"
      ) +
      theme_classic() +
      labs(
        x = if (show_x_text) "Month" else NULL,
        y = if (show_left_y) "State" else NULL
      ) +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text.x = if (show_x_text) element_text() else element_blank(),
        axis.ticks.x = if (show_x_text) element_line() else element_blank(),
        axis.title.y = if (show_left_y) element_text() else element_blank(),
        axis.text.y  = if (show_left_y) element_text() else element_blank(),
        axis.ticks.y = if (show_left_y) element_line() else element_blank(),
        legend.position = "none"
      )
  } else {
    p_stage <- ggplot(plateStates_one, aes(x = date, y = state_row, fill = State)) +
      geom_tile(width = 8, height = 0.8) +
      scale_fill_manual(
        values = alt_colors,
        breaks = state_levels,
        labels = state_labels
      ) +
      scale_y_reverse(
        breaks = 1:length(state_levels),
        labels = state_labels[state_levels],
        limits = c(length(state_levels) + 0.5, 0.5)
      ) +
      scale_x_date(
        limits = c(global_x_min, global_x_max),
        date_breaks = "1 month",
        date_labels = "%b"
      ) +
      theme_classic() +
      labs(
        x = if (show_x_text) "Month" else NULL,
        y = if (show_left_y) "State" else NULL,
        fill = "State"
      ) +
      theme(
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        axis.text.x = if (show_x_text) element_text() else element_blank(),
        axis.ticks.x = if (show_x_text) element_line() else element_blank(),
        axis.title.y = if (show_left_y) element_text() else element_blank(),
        axis.text.y  = if (show_left_y) element_text() else element_blank(),
        axis.ticks.y = if (show_left_y) element_line() else element_blank(),
        legend.position = if (show_legend) "right" else "none"
      )
  }

  p_qpcr / p_stage + plot_layout(heights = c(3, 2))
}

# ----------------------------
# Build all panels
# ----------------------------
panel_list <- list()

for (reg in target_region_vec) {
  for (i in seq_along(target_species_vec)) {
    sp <- target_species_vec[i]

    panel_list[[paste(reg, sp, sep = "_")]] <- make_species_region_panel(
      target_species = sp,
      target_region = reg,
      show_x_text = (reg == tail(target_region_vec, 1)),
      show_left_y = (i == 1),
      show_legend = (reg == "BOF" && sp == "Ciona intestinalis")
    )
  }
}

# ----------------------------
# Top strip row
# ----------------------------
top_row <- wrap_plots(
  make_species_strip("Botrylloides violaceus"),
  make_species_strip("Ciona intestinalis"),
  plot_spacer(),
  ncol = 3,
  widths = c(1, 1, 0.08)
)

# ----------------------------
# Region rows
# ----------------------------
row_MAG <- wrap_plots(
  panel_list[["MAG_Botrylloides violaceus"]],
  panel_list[["MAG_Ciona intestinalis"]],
  make_region_strip("MAG"),
  ncol = 3,
  widths = c(1, 1, 0.08)
)

row_PEI <- wrap_plots(
  panel_list[["PEI_Botrylloides violaceus"]],
  panel_list[["PEI_Ciona intestinalis"]],
  make_region_strip("PEI"),
  ncol = 3,
  widths = c(1, 1, 0.08)
)

row_HAL <- wrap_plots(
  panel_list[["HAL_Botrylloides violaceus"]],
  panel_list[["HAL_Ciona intestinalis"]],
  make_region_strip("HAL"),
  ncol = 3,
  widths = c(1, 1, 0.08)
)

row_BOF <- wrap_plots(
  panel_list[["BOF_Botrylloides violaceus"]],
  panel_list[["BOF_Ciona intestinalis"]],
  make_region_strip("BOF"),
  ncol = 3,
  widths = c(1, 1, 0.08)
)

# ----------------------------
# Final plot
# ----------------------------
final_plot <- wrap_plots(
  top_row,
  row_MAG,
  row_PEI,
  row_HAL,
  row_BOF,
  ncol = 1,
  heights = c(0.08, 1, 1, 1, 1),
  guides = "collect"
) &
  theme(legend.position = "right")

final_plot
