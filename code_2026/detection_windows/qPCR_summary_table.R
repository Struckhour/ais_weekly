source('./AIS_eDNA_data_prep.R')


############################
#SPECIES SD SUMMARY
############################
# --- restrict to species actually tested in each region ---
region_species_tested <- dfRawClean %>%
  dplyr::distinct(region, species)

all_samples <- dfRawClean %>%
  dplyr::distinct(materialSampleID, region, date, station, month)

sample_species_grid <- all_samples %>%
  dplyr::left_join(region_species_tested, by = "region")

# --- collapse to sample × species ---
df_full <- sample_species_grid %>%
  dplyr::left_join(
    dfRawClean %>%
      dplyr::group_by(materialSampleID, species) %>%
      dplyr::summarise(
        logConc = mean(logConc, na.rm = TRUE),
        sample_detected = as.integer(any(detected == 1)),
        .groups = "drop"
      ),
    by = c("materialSampleID", "species")
  ) %>%
  dplyr::mutate(
    sample_detected = tidyr::replace_na(sample_detected, 0),
    logConc = tidyr::replace_na(logConc, 0)
  )

# --- positive-only data ---
df_pos <- df_full %>%
  dplyr::filter(sample_detected == 1)

# --- sample-level detection by species ---
sample_detection_species <- df_full %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    n_obs = dplyr::n(),
    n_samples = dplyr::n_distinct(materialSampleID),
    n_regions = dplyr::n_distinct(region),
    n_sample_detected = sum(sample_detected == 1),
    sample_detection_rate = mean(sample_detected == 1),
    .groups = "drop"
  )

# --- replicate-level detection by species ---
rep_detection_species <- dfRawClean %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    n_reps = dplyr::n(),
    n_pos_reps = sum(detected == 1, na.rm = TRUE),
    rep_detection_rate = mean(detected == 1, na.rm = TRUE),
    .groups = "drop"
  )

# --- monthly SD by species, positive only ---
monthly_sd_species <- df_pos %>%
  dplyr::group_by(species, month) %>%
  dplyr::summarise(
    monthly_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    mean_monthly_sd = mean(monthly_sd, na.rm = TRUE),
    .groups = "drop"
  )

# --- weekly SD by species, positive only ---
weekly_sd_species <- df_pos %>%
  dplyr::mutate(
    week = lubridate::floor_date(date, "week")
  ) %>%
  dplyr::group_by(species, week) %>%
  dplyr::summarise(
    weekly_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    mean_weekly_sd = mean(weekly_sd, na.rm = TRUE),
    .groups = "drop"
  )

# --- within-sample SD by species, positive only ---
sample_sd_species <- dfRawClean %>%
  dplyr::filter(detected == 1) %>%
  dplyr::group_by(species, materialSampleID) %>%
  dplyr::summarise(
    sample_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    mean_sd_within_samples = mean(sample_sd, na.rm = TRUE),
    .groups = "drop"
  )

# --- concentration stats by species, positive only ---
summary_species_sd <- df_pos %>%
  dplyr::group_by(species) %>%
  dplyr::summarise(
    mean_log_conc = mean(logConc, na.rm = TRUE),
    median_log_conc = median(logConc, na.rm = TRUE),
    max_log_conc  = max(logConc, na.rm = TRUE),
    sd_log_conc   = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  )

# --- combine ---
summary_species_final <- sample_detection_species %>%
  dplyr::left_join(rep_detection_species, by = "species") %>%
  dplyr::left_join(summary_species_sd, by = "species") %>%
  dplyr::left_join(monthly_sd_species, by = "species") %>%
  dplyr::left_join(weekly_sd_species, by = "species") %>%
  dplyr::left_join(sample_sd_species, by = "species") %>%
  dplyr::arrange(species)

summary_species_final %>%
  dplyr::select(-n_samples) %>%
  dplyr::rename(
    n_samples = n_obs
  ) %>%
  dplyr::select(
    species,
    n_regions,
    n_samples,
    n_sample_detected,
    sample_detection_rate,
    n_reps,
    n_pos_reps,
    rep_detection_rate,
    mean_log_conc,
    median_log_conc,
    max_log_conc,
    sd_log_conc,
    mean_monthly_sd,
    mean_weekly_sd,
    mean_sd_within_samples
  ) %>%
  gt::gt() %>%
  gt::fmt_integer(
    columns = c(n_regions, n_samples, n_sample_detected, n_reps, n_pos_reps)
  ) %>%
  gt::fmt_percent(
    columns = c(sample_detection_rate, rep_detection_rate),
    decimals = 1
  ) %>%
  gt::fmt_number(
    columns = c(
      mean_log_conc, median_log_conc, max_log_conc, sd_log_conc,
      mean_monthly_sd, mean_weekly_sd, mean_sd_within_samples
    ),
    decimals = 2
  ) %>%
  gt::cols_label(
    species = "Species",
    n_regions = "N regions",
    n_samples = "N samples",
    n_sample_detected = "N sample detections",
    sample_detection_rate = "Sample detection rate",
    n_reps = "N replicates",
    n_pos_reps = "N replicate detections",
    rep_detection_rate = "Replicate detection rate",
    mean_log_conc = "Mean log conc. (positive only)",
    median_log_conc = "Median log conc. (positive only)",
    max_log_conc = "Max log conc.",
    sd_log_conc = "Annual SD",
    mean_monthly_sd = "Monthly SD",
    mean_weekly_sd = "Weekly SD",
    mean_sd_within_samples = "Within-sample SD"
  ) %>%
  gt::tab_header(
    title = "Summary of eDNA detection and concentration variability by species"
  ) %>%
  gt::tab_options(
    table.font.size = gt::px(13),
    heading.title.font.size = gt::px(16),
    column_labels.font.size = gt::px(14),
    data_row.padding = gt::px(6),
    column_labels.padding = gt::px(2),
    heading.padding = gt::px(4),
    table.width = gt::pct(90)
  )


library(dplyr)
library(writexl)

summary_species_final %>%
  select(-n_samples) %>%
  rename(n_samples = n_obs) %>%
  select(
    species,
    n_regions,
    n_samples,
    n_sample_detected,
    sample_detection_rate,
    n_reps,
    n_pos_reps,
    rep_detection_rate,
    mean_log_conc,
    median_log_conc,
    max_log_conc,
    sd_log_conc,
    mean_monthly_sd,
    mean_weekly_sd,
    mean_sd_within_samples
  ) %>%
  rename(
    Species = species,
    `N regions` = n_regions,
    `N samples` = n_samples,
    `N sample detections` = n_sample_detected,
    `Sample detection rate` = sample_detection_rate,
    `N replicates` = n_reps,
    `N replicate detections` = n_pos_reps,
    `Replicate detection rate` = rep_detection_rate,
    `Mean log conc. (positive only)` = mean_log_conc,
    `Median log conc. (positive only)` = median_log_conc,
    `Max log conc.` = max_log_conc,
    `Annual SD` = sd_log_conc,
    `Monthly SD` = mean_monthly_sd,
    `Weekly SD` = mean_weekly_sd,
    `Within-sample SD` = mean_sd_within_samples
  ) %>%
  writexl::write_xlsx("summary_species_final.xlsx")

#####################
#new version of region with samples x species
#####################


# --- restrict to species actually tested in each region ---
region_species_tested <- dfRawClean %>%
  distinct(region, species)

all_samples <- dfRawClean %>%
  distinct(materialSampleID, region, date, month)

sample_species_grid <- all_samples %>%
  left_join(region_species_tested, by = "region")

# --- collapse to sample × species ---
df_full <- sample_species_grid %>%
  left_join(
    dfRawClean %>%
      group_by(materialSampleID, species) %>%
      summarise(
        logConc = mean(logConc, na.rm = TRUE),
        sample_detected = as.integer(any(detected == 1)),
        .groups = "drop"
      ),
    by = c("materialSampleID", "species")
  ) %>%
  mutate(
    sample_detected = replace_na(sample_detected, 0),
    logConc = replace_na(logConc, 0)
  )

# --- positive-only data ---
df_pos <- df_full %>%
  filter(sample_detected == 1)

# --- sample-level detection ---
sample_detection_region <- df_full %>%
  group_by(region) %>%
  summarise(
    n_obs = n(),
    n_samples = n_distinct(materialSampleID),
    n_species = n_distinct(species),
    n_sample_detected = sum(sample_detected == 1),
    sample_detection_rate = mean(sample_detected == 1),
    .groups = "drop"
  )

# --- replicate-level detection ---
rep_detection_region <- dfRawClean %>%
  group_by(region) %>%
  summarise(
    n_reps = n(),
    n_pos_reps = sum(detected == 1, na.rm = TRUE),
    rep_detection_rate = mean(detected == 1, na.rm = TRUE),
    .groups = "drop"
  )

# --- monthly SD (positive only) ---
monthly_sd_region <- df_pos %>%
  dplyr::group_by(region, month) %>%
  dplyr::summarise(
    monthly_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    mean_monthly_sd = mean(monthly_sd, na.rm = TRUE),
    .groups = "drop"
  )

# --- weekly SD (positive only) ---
weekly_sd_region <- df_pos %>%
  dplyr::mutate(
    week = lubridate::floor_date(date, "week")
  ) %>%
  dplyr::group_by(region, week) %>%
  dplyr::summarise(
    weekly_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    mean_weekly_sd = mean(weekly_sd, na.rm = TRUE),
    .groups = "drop"
  )

# --- within-sample SD (positive only) ---
sample_sd_region <- dfRawClean %>%
  dplyr::filter(detected == 1) %>%
  dplyr::group_by(region, materialSampleID, species) %>%
  dplyr::summarise(
    sample_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(region) %>%
  dplyr::summarise(
    mean_sd_within_samples = mean(sample_sd, na.rm = TRUE),
    .groups = "drop"
  )
# --- concentration stats (positive only) ---
summary_region_sd <- df_pos %>%
  group_by(region) %>%
  summarise(
    mean_log_conc = mean(logConc, na.rm = TRUE),
    median_log_conc = median(logConc, na.rm = TRUE),
    max_log_conc  = max(logConc, na.rm = TRUE),
    sd_log_conc   = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  )

# --- combine ---
summary_region_final <- sample_detection_region %>%
  dplyr::left_join(rep_detection_region, by = "region") %>%
  dplyr::left_join(summary_region_sd, by = "region") %>%
  dplyr::left_join(monthly_sd_region, by = "region") %>%
  dplyr::left_join(weekly_sd_region, by = "region") %>%
  dplyr::left_join(sample_sd_region, by = "region") %>%
  dplyr::arrange(region)


summary_region_final %>%
  dplyr::select(-n_samples) %>%
  dplyr::rename(
    n_samples = n_obs
  ) %>%
  dplyr::select(
    region,
    n_species,
    n_samples,
    n_sample_detected,
    sample_detection_rate,
    n_reps,
    n_pos_reps,
    rep_detection_rate,
    mean_log_conc,
    median_log_conc,
    max_log_conc,
    sd_log_conc,
    mean_monthly_sd,
    mean_weekly_sd,
    mean_sd_within_samples
  ) %>%
  gt::gt() %>%
  gt::fmt_integer(
    columns = c(n_species, n_samples, n_sample_detected, n_reps, n_pos_reps)
  ) %>%
  gt::fmt_percent(
    columns = c(sample_detection_rate, rep_detection_rate),
    decimals = 1
  ) %>%
  gt::fmt_number(
    columns = c(
      mean_log_conc, median_log_conc, max_log_conc, sd_log_conc,
      mean_monthly_sd, mean_weekly_sd, mean_sd_within_samples
    ),
    decimals = 2
  ) %>%
  gt::cols_label(
    region = "Region",
    n_species = "N species",
    n_samples = "N samples",
    n_sample_detected = "N sample detections",
    sample_detection_rate = "Sample detection rate",
    n_reps = "N replicates",
    n_pos_reps = "N replicate detections",
    rep_detection_rate = "Replicate detection rate",
    mean_log_conc = "Mean log conc. (positive only)",
    median_log_conc = "Median log conc. (positive only)",
    max_log_conc = "Max log conc.",
    sd_log_conc = "Annual SD",
    mean_monthly_sd = "Monthly SD",
    mean_weekly_sd = "Weekly SD",
    mean_sd_within_samples = "Within-sample SD"
  ) %>%
  gt::tab_header(
    title = "Summary of eDNA detection and concentration variability by region"
  ) %>%
  gt::tab_options(
    table.font.size = gt::px(13),
    heading.title.font.size = gt::px(16),
    column_labels.font.size = gt::px(14),
    data_row.padding = gt::px(6),
    column_labels.padding = gt::px(2),
    heading.padding = gt::px(4),
    table.width = gt::pct(90)
  )


library(dplyr)
library(writexl)

summary_region_final %>%
  select(-n_samples) %>%
  rename(n_samples = n_obs) %>%
  select(
    region,
    n_species,
    n_samples,
    n_sample_detected,
    sample_detection_rate,
    n_reps,
    n_pos_reps,
    rep_detection_rate,
    mean_log_conc,
    median_log_conc,
    max_log_conc,
    sd_log_conc,
    mean_monthly_sd,
    mean_weekly_sd,
    mean_sd_within_samples
  ) %>%
  rename(
    Region = region,
    `N species` = n_species,
    `N samples` = n_samples,
    `N sample detections` = n_sample_detected,
    `Sample detection rate` = sample_detection_rate,
    `N replicates` = n_reps,
    `N replicate detections` = n_pos_reps,
    `Replicate detection rate` = rep_detection_rate,
    `Mean log conc. (positive only)` = mean_log_conc,
    `Median log conc. (positive only)` = median_log_conc,
    `Max log conc.` = max_log_conc,
    `Annual SD` = sd_log_conc,
    `Monthly SD` = mean_monthly_sd,
    `Weekly SD` = mean_weekly_sd,
    `Within-sample SD` = mean_sd_within_samples
  ) %>%
  writexl::write_xlsx("summary_region_final.xlsx")


###########################
#all species x region
###########################

# --- sample-level detection (species × region) ---
sample_detection_sr <- df_full %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    n_obs = dplyr::n(),
    n_samples = dplyr::n_distinct(materialSampleID),
    n_sample_detected = sum(sample_detected == 1),
    sample_detection_rate = mean(sample_detected == 1),
    .groups = "drop"
  )

# --- replicate-level detection ---
rep_detection_sr <- dfRawClean %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    n_reps = dplyr::n(),
    n_pos_reps = sum(detected == 1, na.rm = TRUE),
    rep_detection_rate = mean(detected == 1, na.rm = TRUE),
    .groups = "drop"
  )

# --- monthly SD ---
monthly_sd_sr <- df_pos %>%
  dplyr::group_by(species, region, month) %>%
  dplyr::summarise(
    monthly_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    mean_monthly_sd = mean(monthly_sd, na.rm = TRUE),
    .groups = "drop"
  )

# --- weekly SD ---
weekly_sd_sr <- df_pos %>%
  dplyr::mutate(
    week = lubridate::floor_date(date, "week")
  ) %>%
  dplyr::group_by(species, region, week) %>%
  dplyr::summarise(
    weekly_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    mean_weekly_sd = mean(weekly_sd, na.rm = TRUE),
    .groups = "drop"
  )

# --- within-sample SD ---
sample_sd_sr <- dfRawClean %>%
  dplyr::filter(detected == 1) %>%
  dplyr::group_by(species, region, materialSampleID) %>%
  dplyr::summarise(
    sample_sd = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    mean_sd_within_samples = mean(sample_sd, na.rm = TRUE),
    .groups = "drop"
  )

# --- concentration stats ---
summary_sd_sr <- df_pos %>%
  dplyr::group_by(species, region) %>%
  dplyr::summarise(
    mean_log_conc = mean(logConc, na.rm = TRUE),
    max_log_conc  = max(logConc, na.rm = TRUE),
    sd_log_conc   = sd(logConc, na.rm = TRUE),
    .groups = "drop"
  )

# --- combine ---
summary_sr_final <- sample_detection_sr %>%
  dplyr::left_join(rep_detection_sr, by = c("species", "region")) %>%
  dplyr::left_join(summary_sd_sr, by = c("species", "region")) %>%
  dplyr::left_join(monthly_sd_sr, by = c("species", "region")) %>%
  dplyr::left_join(weekly_sd_sr, by = c("species", "region")) %>%
  dplyr::left_join(sample_sd_sr, by = c("species", "region")) %>%
  dplyr::arrange(species, region)

summary_sr_final_clean <- summary_sr_final %>%
  dplyr::arrange(species, region) %>%
  dplyr::mutate(
    species = as.character(species),
    species_display = dplyr::if_else(
      species == dplyr::lag(species, default = ""),
      "",
      species
    )
  )

summary_sr_final_clean %>%
  dplyr::select(-n_samples) %>%
  dplyr::rename(n_samples = n_obs) %>%
  dplyr::select(
    species_display, region,
    n_samples,
    n_sample_detected,
    sample_detection_rate,
    n_reps,
    n_pos_reps,
    rep_detection_rate,
    mean_log_conc,
    max_log_conc,
    sd_log_conc,
    mean_monthly_sd,
    mean_weekly_sd,
    mean_sd_within_samples
  ) %>%
  gt::gt() %>%
  gt::fmt_integer(
    columns = c(n_samples, n_sample_detected, n_reps, n_pos_reps)
  ) %>%
  gt::cols_label(
    species_display = "Species",
    region = "Region",
    n_samples = "N samples",
    n_sample_detected = "N sample detections",
    sample_detection_rate = "Sample detection rate",
    n_reps = "N replicates",
    n_pos_reps = "N replicate detections",
    rep_detection_rate = "Replicate detection rate",
    mean_log_conc = "Mean log conc. (positive only)",
    max_log_conc = "Max log conc.",
    sd_log_conc = "Annual SD",
    mean_monthly_sd = "Monthly SD",
    mean_weekly_sd = "Weekly SD",
    mean_sd_within_samples = "Within-sample SD"
  ) %>%
  gt::cols_align(
    align = "center",
    everything()
  ) %>%
  gt::cols_width(
    everything() ~ gt::px(90)   # key line
  ) %>%
  gt::fmt_percent(
    columns = c(sample_detection_rate, rep_detection_rate),
    decimals = 1
  ) %>%
  gt::fmt_number(
    columns = c(
      mean_log_conc, max_log_conc, sd_log_conc,
      mean_monthly_sd, mean_weekly_sd, mean_sd_within_samples
    ),
    decimals = 2
  ) %>%
  gt::tab_header(
    title = "Summary of eDNA detection and concentration variability by species and region"
  ) %>%
  gt::cols_width(
    species_display ~ gt::px(95),
    region ~ gt::px(45),
    n_samples ~ gt::px(55),
    n_sample_detected ~ gt::px(65),
    sample_detection_rate ~ gt::px(70),
    n_reps ~ gt::px(55),
    n_pos_reps ~ gt::px(65),
    rep_detection_rate ~ gt::px(70),
    mean_log_conc ~ gt::px(75),
    max_log_conc ~ gt::px(55),
    sd_log_conc ~ gt::px(55),
    mean_monthly_sd ~ gt::px(55),
    mean_weekly_sd ~ gt::px(55),
    mean_sd_within_samples ~ gt::px(65)
  ) %>%
  gt::tab_options(
    table.font.size = gt::px(11),
    heading.title.font.size = gt::px(13),
    column_labels.font.size = gt::px(11),
    data_row.padding = gt::px(3),
    column_labels.padding = gt::px(1),
    heading.padding = gt::px(2),
    table.width = gt::pct(100)
  )


library(dplyr)
library(writexl)

summary_sr_final_clean %>%
  select(-n_samples) %>%
  rename(n_samples = n_obs) %>%
  select(
    species_display,
    region,
    n_samples,
    n_sample_detected,
    sample_detection_rate,
    n_reps,
    n_pos_reps,
    rep_detection_rate,
    mean_log_conc,
    max_log_conc,
    sd_log_conc,
    mean_monthly_sd,
    mean_weekly_sd,
    mean_sd_within_samples
  ) %>%
  rename(
    Species = species_display,
    Region = region,
    `N samples` = n_samples,
    `N sample detections` = n_sample_detected,
    `Sample detection rate` = sample_detection_rate,
    `N replicates` = n_reps,
    `N replicate detections` = n_pos_reps,
    `Replicate detection rate` = rep_detection_rate,
    `Mean log conc. (positive only)` = mean_log_conc,
    `Max log conc.` = max_log_conc,
    `Annual SD` = sd_log_conc,
    `Monthly SD` = mean_monthly_sd,
    `Weekly SD` = mean_weekly_sd,
    `Within-sample SD` = mean_sd_within_samples
  ) %>%
  writexl::write_xlsx("summary_species_region_final.xlsx")

























###############
#COMPARISONS without zeroes
###############

ggplot(df_pos, aes(x = region, y = logConc, fill = region)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~ species, scales = "free_y") +
  scale_fill_manual(values = hybrid_color) +
  theme_classic()


library(rstatix)

df_pos %>%
  group_by(species) %>%
  kruskal_test(logConc ~ region)

pairwise_results <- df_pos %>%
  group_by(species) %>%
  pairwise_wilcox_test(logConc ~ region, p.adjust.method = "BH")


###################
#COMPARISON WITH ZEROES
###################

p <- ggplot(df_full, aes(x = region, y = logConc, fill = region)) +
  geom_boxplot(outlier.alpha = 0.2) +
  facet_wrap(~ species, scales = "free_y") +
  scale_fill_manual(values = hybrid_color) +
  theme_classic()

ggsave("manuscript_figures/figure_2.png", p, width = 8, height = 6, dpi = 300)

ggplot(df_full, aes(x = region, y = logConc, fill = region)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  facet_wrap(~ species, scales = "free_y") +
  scale_fill_manual(values = hybrid_color) +
  theme_classic()

ggplot(df_full, aes(x = region, y = logConc, fill = region)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  facet_wrap(~ species, scales = "free_y") +
  scale_fill_manual(values = hybrid_color) +
  theme_classic()

ggplot(df_full, aes(x = region, y = logConc, fill = region)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  facet_wrap(~ species, scales = "free_y") +
  scale_fill_manual(values = c("black", "black", "black", "black", "black")) +
  theme_classic()

df_full %>%
  group_by(species) %>%
  kruskal_test(logConc ~ region)

pairwise_results <- df_full %>%
  group_by(species) %>%
  pairwise_wilcox_test(logConc ~ region, p.adjust.method = "BH")

print(pairwise_results, n=Inf, width=Inf)

################
#STATION
################

ggplot(df_full, aes(x = interaction(region, station), y = logConc, fill = region)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white") +
  facet_wrap(~ species, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = hybrid_color) +
  theme_classic() +
  labs(x = "Region / Station") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.spacing = unit(0.5, "lines"))

ggplot(df_full, aes(x = interaction(region, station), y = logConc, fill = region)) +
  geom_boxplot(width = 0.5, outlier.alpha = 0.3) +
  facet_wrap(~ species, scales = "free_y", ncol = 1) +
  scale_fill_manual(values = hybrid_color) +
  theme_classic() +
  labs(x = "Region / Station") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(0.5, "lines")
  )





library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)
library(lmerTest)
library(emmeans)

# --- add week ---
df_full_week <- df_full %>%
  mutate(
    week = floor_date(date, "week")
  )

# --- fit mixed models and run all pairwise region comparisons ---
pairwise_emmeans_results <- df_full_week %>%
  group_by(species) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lmerTest::lmer(logConc ~ region + (1 | week), data = .x)),
    emm = map(model, ~ emmeans::emmeans(.x, ~ region)),
    pairwise = map(emm, ~ pairs(.x, adjust = "tukey")),
    pairwise_table = map(pairwise, ~ as.data.frame(.x))
  ) %>%
  select(species, pairwise_table) %>%
  unnest(pairwise_table)

# --- view results ---
print(pairwise_emmeans_results, n = Inf, width = Inf)
