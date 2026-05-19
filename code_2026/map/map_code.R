library(sf)
library(dplyr)
library(ggplot2)

# bioregions <- st_read("code_2026/map/FederalMarineBioregions_SHP/FederalMarineBioregions.shp")

bioregions <- st_read("code_2026/map/meow/meow_ecos.shp")

st_crs(bioregions)
names(bioregions)
plot(st_geometry(bioregions))


library(sf)
library(dplyr)

library(dplyr)

region_centroids <- dfRawClean %>%
  group_by(region) %>%
  summarise(
    mean_lat = mean(decimalLatitude, na.rm = TRUE),
    mean_lon = mean(decimalLongitude, na.rm = TRUE),
    .groups = "drop"
  )

points_sf <- st_as_sf(
  region_centroids,
  coords = c("mean_lon", "mean_lat"),
  crs = 4326
) %>%
  mutate(region = factor(region, levels = names(region_colors)))

points_sf <- points_sf %>%
  mutate(
    nudge_y = ifelse(region %in% c("MAG", "GOM", "PEI"), -0.2, 0.2)
  )

bioregions_4326 <- st_transform(bioregions, 4326)


library(ggrepel)

names(bioregions_4326)

bioregions_subset <- bioregions_4326 %>%
  dplyr::filter(ECOREGION %in% c("Scotian Shelf", "Gulf of St. Lawrence - Eastern Scotian Shelf", "Gulf of Maine/Bay of Fundy"))


library(rnaturalearth)
library(rnaturalearthdata)
library(sf)

land <- ne_countries(
  scale = "medium",
  returnclass = "sf"
) %>%
  st_transform(4326)

land_labels <- tibble::tribble(
  ~label,              ~lon,   ~lat,
  "Maine",             -68.7,  45.0,
  "New Brunswick",     -66.2,  46.2,
  "Prince Edward\nIsland", -63.3, 46.8,
  "Magdalen\nIslands", -62.4, 47.6,
  "Nova Scotia",       -65,  44.5
)

water_labels <- tibble::tribble(
  ~label,                  ~lon,   ~lat,
  "Gulf of Maine",         -68.2,  43.2,
  "Bay of Fundy",          -66.1,  44.9,
  "Scotian Shelf",         -63.5,  43.7,
  "Gulf of St. Lawrence",  -63.5,  48.1
)

bioregion_colors <- c(
  "Scotian Shelf" = "#b3e2cd",          # soft purple
  "Gulf of St. Lawrence - Eastern Scotian Shelf" = "#fdcdac",    # soft green
  "Gulf of Maine/Bay of Fundy" = "#cbd5e8"
)

bioregions_subset <- bioregions_subset %>%
  mutate(
    ECOREGION = factor(
      ECOREGION,
      levels = c(
        "Gulf of St. Lawrence - Eastern Scotian Shelf",
        "Scotian Shelf",
        "Gulf of Maine/Bay of Fundy"
      )
    )
  )

p <- ggplot() +
  geom_sf(
    data = land,
    fill = "grey90",
    color = "grey50",
    linewidth = 0.3
  ) +
  geom_sf(
    data = bioregions_subset,
    aes(fill = ECOREGION),
    color = "black",
    linetype = "dotted",
    linewidth = 0.6,   # optional: make dots more visible
    alpha = 0.2
  ) +
  scale_fill_manual(
    values = bioregion_colors
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 0.6)  # legend more opaque
    )
  ) +
  geom_text(
    data = land_labels,
    aes(x = lon, y = lat, label = label),
    fontface = "bold",
    size = 3.2,
    color = "grey20"
  ) +
  geom_text(
    data = water_labels,
    aes(x = lon, y = lat, label = label),
    fontface = "italic",
    size = 3.1,
    color = "grey25"
  ) +
  geom_sf(
    data = points_sf,
    aes(color = region),
    size = 4
  ) +
  geom_sf_text(
    data = points_sf,
    aes(label = region),
    color = "black",
    size = 5,
    fontface = "bold",
    nudge_y = points_sf$nudge_y,
    nudge_x = 0.1
  ) +
  scale_color_manual(values = region_colors) +
  coord_sf(
    xlim = c(-70, -61),
    ylim = c(43, 48.3),
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = NULL,
    fill = "Marine Ecoregions of the World Bioregion",
    color = "Region"
  )

p
ggsave("code_2026/map/figure_1_before_edit.png", p, width = 9, height = 9, dpi = 300)



#SMALL GOM PLOT
# SMALL GOM PLOT

gom_station_points <- dfRawClean %>%
  filter(region == "GOM") %>%
  group_by(station) %>%
  summarise(
    mean_lat = mean(decimalLatitude, na.rm = TRUE),
    mean_lon = mean(decimalLongitude, na.rm = TRUE),
    .groups = "drop"
  )

gom_points_sf <- st_as_sf(
  gom_station_points,
  coords = c("mean_lon", "mean_lat"),
  crs = 4326
)

gom_p <- ggplot() +
  geom_sf(
    data = land,
    fill = "grey90",
    color = "grey50",
    linewidth = 0.3
  ) +
  geom_sf(
    data = bioregions_subset,
    aes(fill = ECOREGION),
    color = "black",
    linetype = "dotted",
    linewidth = 0.6,
    alpha = 0.2
  ) +
  scale_fill_manual(values = bioregion_colors) +
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 0.6)
    )
  ) +
  geom_sf(
    data = gom_points_sf,
    color = region_colors["GOM"],
    size = 4
  ) +
  geom_sf_text(
    data = gom_points_sf,
    aes(label = station),
    color = "black",
    size = 4,
    fontface = "bold",
    nudge_y = 0.05
  ) +
  coord_sf(
    xlim = c(-68.9, -68.5),
    ylim = c(44.0, 44.7),
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    x = NULL,
    y = NULL,
    fill = "Marine Ecoregions of the World Bioregion"
  )

gom_p

gom_points_sf <- st_as_sf(
  gom_station_points,
  coords = c("mean_lon", "mean_lat"),
  crs = 4326
)

st_distance(gom_points_sf[1, ], gom_points_sf[2, ])

# ggsave("code_2026/map/gom_p.png", p, width = 9, height = 9, dpi = 300)
