source('./AIS_eDNA_data_prep.R')


View(dfRawClean)

propdf <- dfRawClean %>% select(logConc, region, species, month, year, date, temp, sal, tds, pH, turb, chl, tss)
View(propdf)



selected_species <- "Carcinus maenas"
selected_region <- "MAG"
df <- propdf %>% filter(species == selected_species, region == selected_region)

##########################
#CHECK FOR SPECIES VS PROPERTY
##########################
df <- propdf %>% filter(species == selected_species)

ggplot(df, aes(x = temp, y = logConc)) +
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_minimal() +
  labs(x = "Temperature", y = "log Concentration")

cor.test(df$temp, df$logConc, method = "pearson", use = "complete.obs")
cor.test(df$temp, df$logConc, method = "spearman", use = "complete.obs")
model <- lm(logConc ~ temp, data = df)
summary(model)



################
#HEATMAP
################

get_correlations <- function(df, property) {

  df %>%
    group_by(region, species) %>%
    summarise(
      n = sum(!is.na(.data[[property]]) & !is.na(logConc)),
      r = if (n >= 3) {
        cor(.data[[property]], logConc, use = "complete.obs")
      } else {
        NA_real_
      },
      .groups = "drop"
    )
}

# cor_df <- get_correlations(propdf, "temp")
# cor_df <- get_correlations(propdf, "pH")
# cor_df <- get_correlations(propdf, "turb")
# cor_df <- get_correlations(propdf, "sal")
# cor_df <- get_correlations(propdf, "tds")
# cor_df <- get_correlations(propdf, "chl")
# cor_df <- get_correlations(propdf, "tss")
cor_df <- get_correlations(propdf, "pH") %>%
  group_by(region) %>%
  filter(!all(is.na(r))) %>%   # remove regions with no usable data at all
  ungroup() %>%
  mutate(region = factor(region, levels = rev(region_order)))

cor_df <- cor_df %>%
  tidyr::complete(region, species)

cor_df <- cor_df %>%
  mutate(region = factor(region, levels = rev(region_order)))

ggplot(cor_df, aes(x = species, y = region, fill = r)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "orange",
    midpoint = 0,
    limits = c(-1, 1),
    na.value = "grey80",
    name = "Pearson r"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "Species", y = "Region")
