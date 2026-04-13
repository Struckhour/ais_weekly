

library(rpart)

source('./AIS_eDNA_data_prep.R')

sp <- "Membranipora membranacea"
sp <- "Botrylloides violaceus"
sp <- "Ciona intestinalis"

df_sp <- dfWeeks %>%
  filter(species == sp)


df_sp_clean <- df_sp %>%
  as.data.frame() %>%
  dplyr::filter(
    is.finite(scaleLogConc),
    is.finite(week_of_year),
    is.finite(meanTemp),
    is.finite(meanSal),
    is.finite(meanPH),
    is.finite(meanLat)
  )



tree <- rpart(
  scaleLogConc ~ meanTemp + meanSal + meanPH + meanLat,
  data = df_sp_clean
)

plot(tree)
text(tree)

df_sp_clean$pred_tree <- predict(tree, newdata = df_sp_clean)
ggplot(df_sp_clean, aes(x = pred_tree, y = scaleLogConc)) +
  geom_point(alpha = 0.4) +
  geom_abline(slope = 1, intercept = 0, color = "red") +
  labs(
    x = "Predicted (tree)",
    y = "Observed",
    title = paste("Regression Tree —", sp)
  ) +
  theme_minimal()

r2_tree <- cor(df_sp_clean$pred_tree, df_sp_clean$scaleLogConc)^2
rmse_tree <- sqrt(mean((df_sp_clean$scaleLogConc - df_sp_clean$pred_tree)^2))

r2_tree
rmse_tree
printcp(tree)
best_cp <- tree$cptable[which.min(tree$cptable[,"xerror"]), "CP"]

tree_pruned <- prune(tree, cp = best_cp)
plot(tree_pruned)
text(tree_pruned)
