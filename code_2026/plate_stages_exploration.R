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
