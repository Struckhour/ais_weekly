################################################################################
# Script:    AIS_eDNA_data_prep.R
# Purpose:   Clean, merge, and enrich weekly AIS eDNA and plate–abundance data
# Author:    Melissa K. Morrison
# Created:   20‑Jun‑2025
# Modified:  25‑Jun‑2025  (extensive annotation for long‑term storage)
#
# Overview ---------------------------------------------------------------------
# 1) Load raw metadata and sample tables
# 2) Merge metadata onto sample rows (by materialSampleID)
# 3) Re‑order months so plots start at each region’s first month of sampling
# 4) Remove known outliers from Halifax Oct‑19‑2023
# 5) Summarise eDNA by week and month
# 6) Load and tidy plate‑abundance data (for tunicates only)
# 7) Build tidy “stage” table with life‑stage & colony‑state decomposition
# 8) Prepare correlation‑ready data frames (raw / cleaned)
#
# NOTES:   ‑ All file paths are relative to project root (`./data/pub2/…`)
#          ‑ “scale_propMinAIS()” rescales concentrations to 0‑1 locally
################################################################################

# ───────────────────────────────────────────────────────────────────────────────
# 1. Load libraries -------------------------------------------------------------
#    • Keep tidyverse first (includes ggplot2, dplyr, tidyr, stringr, readr, …)
#    • Avoid loading full tidyverse twice
#    • mgcv is for later GAM modelling
# ───────────────────────────────────────────────────────────────────────────────
library(tidyverse)      # readr, dplyr, tidyr, ggplot2, stringr, etc.
library(readxl)         # read Excel if needed elsewhere
library(lubridate)      # dates
library(patchwork)      # plot layout
library(ggh4x)          # faceting helpers
library(ggExtra)        # marginal distributions
library(viridis)        # colour scales
library(grid)           # low‑level plotting
library(cowplot)        # ggplot utils
library(gridExtra)      # arrangeGrob, etc.
library(DescTools)      # miscellaneous stats
library(mgcv)           # generalized additive models

# ───────────────────────────────────────────────────────────────────────────────
# 2. Read metadata (sample‑level environmental data + coordinates) -------------
# ───────────────────────────────────────────────────────────────────────────────
metadata <- read.csv(
  "./data/AIS_weekly_metadata_BOF_adjusted.csv",
  check.names = FALSE,
  na.strings  = ""
) %>%
  mutate(
    date  = as.Date(eventDate, origin = "1899-12-30"),
    year  = year(date),
    month = month(date),

    waterTemp_C = dplyr::if_else(
      region == "PEI" &
        date == as.Date("2023-11-30") &
        trimws(waterTemp_C) == "27.5",
      NA_character_,
      waterTemp_C
    )
  ) %>%
  filter(is.na(controlType))

# ───────────────────────────────────────────────────────────────────────────────
# 2b. Read SBE data and make daily means for Halifax only ----------------------
# ───────────────────────────────────────────────────────────────────────────────

new_data <- read_xlsx(
  "./data/BIO_Stn24_SBE_CT_2023_2024_dateorderd.xlsx",
  na = ""
)

new_data <- new_data %>% select(date, temp, PSU, ts) %>%
  filter(ts >= as.Date("2023-04-27") & ts <= as.Date("2024-06-3"))

new_data_clean <- new_data %>%
  mutate(
    date = as.Date(date)
  )

sbe_daily <- new_data_clean %>%
  group_by(date) %>%
  summarise(
    temp = mean(temp, na.rm = TRUE),
    PSU = mean(PSU, na.rm = TRUE),
    .groups = "drop"
  )

sbe_daily <- new_data %>%
  mutate(date = as.Date(date)) %>%
  group_by(date) %>%
  summarise(
    sbe_temp = mean(temp, na.rm = TRUE),
    sbe_sal  = mean(PSU, na.rm = TRUE),
    .groups = "drop"
  )

# ───────────────────────────────────────────────────────────────────────────────
# 2c. Replace HAL metadata temp/salinity with SBE daily values -----------------
# ───────────────────────────────────────────────────────────────────────────────
metadata <- metadata %>%
  mutate(
    waterTemp_C  = as.numeric(waterTemp_C),
    salinity_ppt = as.numeric(salinity_ppt)
  ) %>%
  left_join(sbe_daily, by = "date") %>%
  mutate(
    waterTemp_C = dplyr::if_else(
      region == "HAL" & !is.na(sbe_temp),
      sbe_temp,
      waterTemp_C
    ),
    salinity_ppt = dplyr::if_else(
      region == "HAL" & !is.na(sbe_sal),
      sbe_sal,
      salinity_ppt
    )
  ) %>%
  select(-sbe_temp, -sbe_sal)

#CHANGE RIDICULOUS PH to NA
metadata <- metadata %>%
  dplyr::mutate(
    pH = as.numeric(pH),
    pH = dplyr::if_else(pH > 0 & pH < 10, pH, NA_real_)
  )

# ── OPTIONAL: save coordinates for static map (run once then comment out) ──
# metadata %>%
#   distinct(samplingStation, decimalLatitude, decimalLongitude) %>%
#   write.csv("./data/pub2/mapcoords.csv", row.names = FALSE)

# ───────────────────────────────────────────────────────────────────────────────
# 3. Read eDNA qPCR concentration table ----------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
samples <- read.csv(
  "./data/AIS_weekly_data.csv",
  check.names = FALSE,
  na.strings  = ""
)

# ───────────────────────────────────────────────────────────────────────────────
# 4. Helper: min‑max scaling to [0,1] ----------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
scale_propMinAIS <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (diff(rng) == 0) return(rep(0, length(x)))
  (x - rng[1]) / diff(rng)
}

# ───────────────────────────────────────────────────────────────────────────────
# 5. Merge metadata onto sample rows -------------------------------------------
#    • We match by materialSampleID
#    • Each column is filled from metadata where available
# ───────────────────────────────────────────────────────────────────────────────
dfRaw <- samples %>%
  mutate(
    # Vectorised look‑ups (match() returns row index in metadata)
    idx  = match(materialSampleID, metadata$materialSampleID),

    date            = metadata$date[idx],
    region          = metadata$region[idx],
    decimalLatitude = metadata$decimalLatitude[idx],
    decimalLongitude= metadata$decimalLongitude[idx],
    station         = metadata$samplingStation[idx],

    # Environmental variables
    temp = metadata$waterTemp_C[idx],
    sal  = metadata$salinity_ppt[idx],
    tds  = metadata$`TDS_mg/L`[idx],
    pH   = metadata$pH[idx],
    turb = metadata$turbidity_NTU[idx],
    chl  = metadata$chlorophyll[idx],
    tss  = metadata$`TSS_mg/L`[idx],

    # Time helpers
    year  = year(date),
    month = month(date),

    # eDNA columns
    concentration = as.numeric(concentration),
    species       = scientificName,
    primer        = "COI",
    domain        = "Eukaryota"
  ) %>%
  # Safeguards -----------------------------------------------------------------
  drop_na(date, decimalLatitude) %>%        # remove samples lacking essential info
  filter(!is.na(kingdom)) %>%             # discard taxonomy blanks
  # Remove metadata columns we do not need downstream
  select(-idx, -basisOfRecord, -recordedBy,
         -occurrenceStatus, -quantificationCycle, -concentrationUnit) %>%
  mutate(concentration = coalesce(concentration, 0),
         detected      = if_else(concentration > 0, 1, 0)) # treat any non‑zero as “detected”

# ───────────────────────────────────────────────────────────────────────────────
# 5. Re‑order months so each region’s “Month 1” = first sampling month ----------
# ───────────────────────────────────────────────────────────────────────────────
start_months <- c(MAG = 5, BOF = 5, PEI = 4, HAL = 4, GOM = 6)  # May, May, …

dfRaw <- dfRaw %>%
  mutate(
    start_month     = start_months[region],
    month_reordered = (month - start_month + 12) %% 12 + 1,
    calendar_month = month
  )

# ───────────────────────────────────────────────────────────────────────────────
# 6. Create rotated month‑label lookup (single‑letter abbreviations) -----------
# ───────────────────────────────────────────────────────────────────────────────
month_labels_df <- map_dfr(names(start_months), function(reg) {
  start  <- start_months[[reg]]
  labels <- month.abb[(start - 1 + 0:11) %% 12 + 1]
  tibble(region = reg,
         month_reordered = 1:12,
         month_plot_label = str_sub(labels, 1, 1))
})

# For quick per‑region look‑ups later
month_labels_split <- split(month_labels_df, ~region)

# ───────────────────────────────────────────────────────────────────────────────
# 8. Remove known Halifax outliers (19‑Oct-2023 samples) -----------------------
# ───────────────────────────────────────────────────────────────────────────────
outlierHal <- dfRaw %>%
  dplyr::filter(
    region == "HAL",
    date == as.Date("2023-10-19") | species == "Didemnum vexillum"
  )

dfRawClean <- dfRaw %>%
  dplyr::anti_join(
    outlierHal,
    by = c("region", "date", "species")
  ) %>%
  dplyr::filter(
    !(species == "Didemnum vexillum" & !region %in% c("BOF", "GOM")),
    station != "Causeway"
  )

# ───────────────────────────────────────────────────────────────────────────────
# 7. Colour palette matched to static map (for consistent figures) -------------
# ───────────────────────────────────────────────────────────────────────────────
mapCols <- setNames(
  c("#E15759", "#F28E2B", "#4E79A7", "#76B7B2", "#59A14F"),
  levels(dfRawClean$region)
)
# ───────────────────────────────────────────────────────────────────────────────
# 9. Log10‑transform concentration to stabilise variance -----------------------
# ───────────────────────────────────────────────────────────────────────────────
dfRawClean <- dfRawClean %>%
  mutate(logConc = log(concentration + 1))

# ───────────────────────────────────────────────────────────────────────────────
# 10. Force numeric type on all environmental variables ------------------------
# ───────────────────────────────────────────────────────────────────────────────
num_vars <- c("concentration", "temp", "sal", "turb", "pH", "tss", "tds", "chl")
dfRawClean[num_vars] <- lapply(dfRawClean[num_vars], as.numeric)

# ───────────────────────────────────────────────────────────────────────────────
# 11. Week‑scale summary -------------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
dfWeeks <- dfRawClean %>%
  select(region, species, date, concentration, logConc, month, temp, sal, pH, turb, tds, chl, tss, decimalLatitude, month_reordered) %>%
  group_by(region, species) %>%
  # Week index relative to first sample in that reg.spp group
  mutate(
    sampWeek = floor(as.period(date - min(date)) / weeks()) + 1,
    day_of_year = yday(date),
    week_of_year = floor((day_of_year - 1) / 7) + 1
  ) %>%
  ungroup() %>%
  group_by(region, species, week_of_year) %>%
  summarise(
    nReps      = n(),
    meanConc   = mean(concentration, na.rm = TRUE),
    meanLogConc = mean(logConc, na.rm = TRUE),
    sdConc     = sd(concentration,  na.rm = TRUE),
   # date       = first(date),                     # representative date
  #  month      = first(month),
    meanTemp   = mean(temp, na.rm = TRUE),
    meanSal    = mean(sal,  na.rm = TRUE),
    meanPH   = mean(pH, na.rm = TRUE),
    meanTurb    = mean(turb,  na.rm = TRUE),
    meanTds   = mean(tds, na.rm = TRUE),
    meanChl    = mean(chl,  na.rm = TRUE),
    meanTss    = mean(tss,  na.rm = TRUE),
    meanLat    = mean(decimalLatitude, na.rm = TRUE),
   # month_reordered = first(month_reordered),
    .groups = "drop"
  ) %>%
  # Scale concentration 0‑1 within each region × species
  group_by(region, species) %>%
  mutate(
    scaleConc = scale_propMinAIS(meanConc),
    scaleLogConc = scale_propMinAIS(meanLogConc),
    scaleTemp = scale_propMinAIS(meanTemp),
    scaleSal = scale_propMinAIS(meanSal),
    scalePH = scale_propMinAIS(meanPH)
    ) %>%
  ungroup()

# ───────────────────────────────────────────────────────────────────────────────
# 12. Month‑scale summary ------------------------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
dfMonths <- dfRawClean %>%
  group_by(region, species, month_reordered) %>%
  summarise(
    nReps    = n(),
    nDet     = sum(detected, na.rm = TRUE),
    det_rate = if_else(nReps > 0, nDet / nReps, NA_real_),
    meanConc = mean(concentration, na.rm = TRUE),
    meanPosConc = mean(concentration[concentration > 0], na.rm = TRUE),
    meanLogConc = mean(logConc, na.rm = TRUE),
    meanPosLogConc = mean(logConc[logConc > 0], na.rm = TRUE),
    month_original = first(month),
    sdConc   = sd(concentration,   na.rm = TRUE),
    .groups  = "drop"
  ) %>%
  # Normalize within region × species × year(?)
  group_by(region, species) %>%
  mutate(
    normMeanLogConc = if (all(is.na(meanLogConc))) {
      NA_real_
    } else {
      meanLogConc / max(meanLogConc, na.rm = TRUE)
    },
    normPosMeanLogConc = if (all(is.na(meanPosLogConc))) {
      NA_real_
    } else {
      meanPosLogConc / max(meanPosLogConc, na.rm = TRUE)
    },
    normMeanLogConc = replace_na(normMeanLogConc, 0),
    normPosMeanLogConc = replace_na(normPosMeanLogConc, 0),
    normMeanConc = if (all(is.na(meanConc))) {
      NA_real_
    } else {
      meanConc / max(meanConc, na.rm = TRUE)
    },
    normPosMeanConc = if (all(is.na(meanPosConc))) {
      NA_real_
    } else {
      meanPosConc / max(meanPosConc, na.rm = TRUE)
    },
    normMeanConc = replace_na(normMeanConc, 0),
    normPosMeanConc = replace_na(normPosMeanConc, 0)
  ) %>%
  ungroup() %>%
  # Join rotated month order for plotting
  left_join(month_labels_df, by = c("region", "month_reordered" = "month_reordered")) %>%
  group_by(region, species) %>%
  mutate(scaleConc = scale_propMinAIS(meanConc)) %>%
  ungroup() %>%
  mutate(
    # Handle edge cases where scaling returns NaN
    scaleConc = replace_na(scaleConc, 0),
    # Rescale replicate / detection counts so bubble sizes plot nicely
    meanReps  = if_else(nReps > 9, nReps / 2, nReps),
    meanDets  = if_else(nReps > 9, nDet  / 2, nDet)
  )

# ───────────────────────────────────────────────────────────────────────────────
# 14. Factor levels for species and regions (consistent across plots) -----------
# ───────────────────────────────────────────────────────────────────────────────
dfRawClean$species <- factor(
  dfRawClean$species,
  levels = c("Membranipora membranacea",
             "Botrylloides violaceus",
             "Didemnum vexillum",
             "Ciona intestinalis",
             "Carcinus maenas")
)

dfRawClean$region <- factor(
  dfRawClean$region,
  levels = c("MAG","PEI","HAL","BOF","GOM")
)
# ───────────────────────────────────────────────────────────────────────────────
# 15. Plate‑abundance data (for tunicates) -------------------------------------
# ───────────────────────────────────────────────────────────────────────────────
plateAbundance <- read.csv("./data/AISPlateAnalysis.csv", check.names = FALSE) %>%
  mutate(
    date   = as.Date(Date, origin = "1899-12-30"),
    year   = year(date),
    month  = month(date),
    region = factor(Region, levels = c("MAG", "PEI", "HAL", "BOF")),
    species= Species
  ) %>%
  filter(species %in% c("Botrylloides violaceus", "Ciona intestinalis")) %>%
  group_by(region, species) %>%
  mutate(
    sampWeek = floor(as.period(date - min(dfRaw$date)) / weeks()) + 1
  ) %>%
  ungroup()

# ───────────────────────────────────────────────────────────────────────────────
# 16. Life‑stage / colony state decomposition ----------------------------------
# ───────────────────────────────────────────────────────────────────────────────
# Expand comma-separated life stages and states into long format,
# aligning corresponding stages and states row-by-row.
dfStages <- plateAbundance %>%
  select(region, species, sampWeek, date, State, `Life stage`) %>%
  separate_rows(State, `Life stage`, sep = ",") %>%
  mutate(
    state = str_trim(State),
    stage = str_trim(`Life stage`)
  ) %>%
  filter(state != "" | stage != "") %>%
  mutate(
    # Convert state codes to descriptive factor levels
    state = recode(state,
                   "R" = "Recruitment",
                   "B" = "Breeding",
                   "G" = "Growth",
                   "S" = "Senesence",
                   "D" = "Dieback"),
    state = factor(state, levels = c("Recruitment", "Breeding", "Growth", "Senesence", "Dieback"))
  ) %>%
  select(region, species, sampWeek, date, state, stage)

# ───────────────────────────────────────────────────────────────────────────────
# 17. Correlation matrices ------------------------------------------------------
#    • Aggregated to station‑date mean values
# ───────────────────────────────────────────────────────────────────────────────
core_cols <- c("concentration", "temp", "sal", "turb", "pH", "tss", "tds", "chl")

summarise_env <- function(data) {
  data %>%
    group_by(region, station, species, date) %>%
    summarise(across(all_of(core_cols), mean, na.rm = TRUE), .groups = "drop") %>%
    mutate(     # Remove GOM because environmental variables were not obtained
      region  = factor(region, levels = c("MAG", "PEI", "HAL", "BOF")),
      species = factor(species, levels = levels(dfRawClean$species))
    )
}

dfCor  <- summarise_env(dfRawClean)
dfCorRegion <- split(dfCor, ~region)

# ───────────────────────────────── END OF FILE ─────────────────────────────────
