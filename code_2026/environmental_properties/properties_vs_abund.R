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

      data_complete = list(
        na.omit(data.frame(
          x = .data[[property]],
          logConc = logConc
        ))
      ),

      n = nrow(data_complete[[1]]),

      r = if (n >= 3) {
        cor(
          data_complete[[1]]$x,
          data_complete[[1]]$logConc
        )
      } else {
        NA_real_
      },

      .groups = "drop"
    )
}

cor_df <- get_correlations(propdf, "temp")
cor_df <- get_correlations(propdf, "pH")
cor_df <- get_correlations(propdf, "turb")
cor_df <- get_correlations(propdf, "sal")
cor_df <- get_correlations(propdf, "tds")
cor_df <- get_correlations(propdf, "chl")
cor_df <- get_correlations(propdf, "tss")

ggplot(cor_df, aes(x = species, y = region, fill = r)) +
  geom_tile(color = "grey80") +
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "orange",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Pearson r"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  ) +
  labs(x = "Species", y = "Region")
