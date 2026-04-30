deer_sst <- read.csv(
  "./data/Deer_Isle_SST.csv",
  row.names = 1,
  check.names = FALSE,
  na.strings  = ""
) %>%
  dplyr::mutate(
    station = dplyr::recode(
      site_name,
      "Deer Isle 1" = "Sand Beach",
      "Deer Isle2" = "Causeway",
      "Deer Isle3" = "Isle Haut ferry"
    )
  )

deer_sss_site1 <- read.csv(
  "./data/SSS_Deer_Island_site1.csv",
  check.names = FALSE,
  na.strings  = ""
) %>%
  mutate(station = "Sand Beach")

deer_sss_site2 <- read.csv(
  "./data/SSS_Deer_Island_site2.csv",
  check.names = FALSE,
  na.strings  = ""
) %>%
  mutate(station = "Causeway")

deer_sss_site3 <- read.csv(
  "./data/SSS_Deer_Island_site3.csv",
  check.names = FALSE,
  na.strings  = ""
) %>%
  mutate(station = "Isle Haut ferry")

deer_sss <- bind_rows(deer_sss_site1, deer_sss_site2, deer_sss_site3)
