View(plateAbundance)


library(tidyverse)
library(lubridate)

# 1. Define global state order (consistent everywhere)
all_states <- plateAbundance %>%
  separate_rows(State, sep = ",\\s*") %>%
  mutate(State = trimws(State)) %>%
  filter(State != "") %>%
  distinct(State) %>%
  arrange(State) %>%
  pull(State)

# 2. Prepare data
plateStates <- plateAbundance %>%
  filter(date >= as.Date("2022-06-01") & date <= as.Date("2028-05-31")) %>%
  separate_rows(State, sep = ",\\s*") %>%
  mutate(
    State = trimws(State),
    species = factor(species, levels = unique(species))
  ) %>%
  filter(State != "") %>%
  mutate(
    state_row = as.numeric(factor(State, levels = all_states))
  )

# 3. Plot
ggplot(plateStates, aes(x = date, y = state_row, fill = State)) +
  geom_tile(width = 8, height = 0.8) +
  facet_grid(region ~ species) +
  scale_fill_brewer(palette = "Set2") +
  scale_y_continuous(
    breaks = 1:length(all_states),
    labels = all_states
  ) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b") +
  guides(fill = guide_legend(reverse = TRUE)) +
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
