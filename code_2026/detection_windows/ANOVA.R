source('./AIS_eDNA_data_prep.R')


library(dplyr)
library(broom)

anova_results <- dfRawClean %>%
  mutate(month = factor(month)) %>%
  group_by(species, region) %>%
  group_modify(~{

    if (n_distinct(.x$month) < 2) {
      return(tibble(
        F = NA_real_,
        p = NA_real_,
        df1 = NA_real_,
        df2 = NA_real_
      ))
    }

    fit <- aov(logConc ~ month, data = .x)
    out <- tidy(fit)

    tibble(
      F = out$statistic[1],
      p = out$p.value[1],
      df1 = out$df[1],
      df2 = out$df[2]
    )
  }) %>%
  ungroup()

anova_results
