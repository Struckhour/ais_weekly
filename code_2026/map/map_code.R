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







































sampling_locations <- dfRawClean %>%
  filter(
    !is.na(region),
    !is.na(station),
    !is.na(decimalLongitude),
    !is.na(decimalLatitude)
  ) %>%
  distinct(
    region,
    station,
    decimalLongitude,
    decimalLatitude
  ) %>%
  arrange(region, station)

sampling_locations

sampling_locations_sf <- sampling_locations %>%
  mutate(
    region = factor(region, levels = names(region_colors))
  ) %>%
  st_as_sf(
    coords = c("decimalLongitude", "decimalLatitude"),
    crs = 4326,
    remove = FALSE
  )

make_region_map <- function(
    region_code,
    sites_sf = sampling_locations_sf,
    land_sf = land,
    region_palette = region_colors,
    padding_lon = 0.15,
    padding_lat = 0.12
) {

  region_sites <- sites_sf %>%
    filter(region == region_code)

  if (nrow(region_sites) == 0) {
    stop(paste("No sampling locations found for", region_code))
  }

  site_coordinates <- st_coordinates(region_sites)

  x_range <- range(site_coordinates[, "X"], na.rm = TRUE)
  y_range <- range(site_coordinates[, "Y"], na.rm = TRUE)

  # Give maps a reasonable minimum width and height when sites
  # are very close together.
  # if (diff(x_range) < 0.05) {
  #   x_range <- mean(x_range) + c(-0.025, 0.025)
  # }
  #
  # if (diff(y_range) < 0.05) {
  #   y_range <- mean(y_range) + c(-0.025, 0.025)
  # }

  x_limits <- x_range + c(-padding_lon, padding_lon)
  y_limits <- y_range + c(-padding_lat, padding_lat)

  ggplot() +
    geom_sf(
      data = land_sf,
      fill = "grey90",
      color = "grey50",
      linewidth = 0.3
    ) +
    geom_sf(
      data = region_sites,
      aes(color = region),
      size = 3.5
    ) +
    geom_text_repel(
      data = region_sites,
      aes(
        x = decimalLongitude,
        y = decimalLatitude,
        label = station
      ),
      size = 3.2,
      fontface = "bold",
      min.segment.length = 0,
      box.padding = 0.4,
      point.padding = 0.3,
      seed = 1
    ) +
    scale_color_manual(
      values = region_palette,
      limits = names(region_palette)
    ) +
    coord_sf(
      xlim = x_limits,
      ylim = y_limits,
      expand = FALSE
    ) +
    labs(
      title = region_code,
      x = NULL,
      y = NULL
    ) +
    theme_classic() +
    theme(
      legend.position = "none",
      plot.title = element_text(
        face = "bold",
        hjust = 0.5
      )
    )
}

p_MAG <- make_region_map("MAG", padding_lon = 1.85, padding_lat = 1.82)
p_MAG
p_PEI <- make_region_map("PEI")
p_HAL <- make_region_map("HAL")
p_BOF <- make_region_map("BOF")
p_GOM <- make_region_map("GOM")






























library(ggspatial)
make_region_map <- function(
    region_code,
    sites_sf = sampling_locations_sf,
    padding_lon = 0.05,
    padding_lat = 0.04,
    tile_zoom = 1,
    merge_distance_m = 100
) {

  region_sites <- sites_sf %>%
    filter(as.character(region) == region_code)

  if (nrow(region_sites) == 0) {
    stop(paste("No sampling locations found for", region_code))
  }

  # ------------------------------------------------------------
  # Merge sampling locations within merge_distance_m
  # ------------------------------------------------------------

  region_sites_m <- region_sites %>%
    st_transform(3347)

  if (nrow(region_sites_m) > 1) {

    distance_matrix <- st_distance(region_sites_m)

    site_clusters <- hclust(
      as.dist(distance_matrix),
      method = "single"
    )

    region_sites_m$cluster_id <- cutree(
      site_clusters,
      h = merge_distance_m
    )

  } else {

    region_sites_m$cluster_id <- 1L
  }

  region_sites_merged <- region_sites_m %>%
    group_by(cluster_id) %>%
    summarise(
      region = first(region),
      station = paste(
        sort(unique(station)),
        collapse = " / "
      ),
      n_locations = n(),
      geometry = st_centroid(st_union(geometry)),
      .groups = "drop"
    ) %>%
    st_transform(4326)

  # ------------------------------------------------------------
  # Calculate map limits and midpoint labels
  # ------------------------------------------------------------

  site_coordinates <- st_coordinates(region_sites_merged)

  x_range <- range(site_coordinates[, "X"], na.rm = TRUE)
  y_range <- range(site_coordinates[, "Y"], na.rm = TRUE)

  x_limits <- x_range + c(-padding_lon, padding_lon)
  y_limits <- y_range + c(-padding_lat, padding_lat)

  x_mid <- mean(x_limits)
  y_mid <- mean(y_limits)

  # ------------------------------------------------------------
  # Draw map
  # ------------------------------------------------------------

  ggplot() +

    annotation_map_tile(
      type = "cartolight",
      zoomin = tile_zoom,
      progress = "none"
    ) +

    geom_sf(
      data = region_sites_merged,
      shape = 21,
      fill = "grey10",
      color = "black",
      stroke = 0.4,
      size = 2.8
    ) +

    scale_x_continuous(
      breaks = x_mid,
      labels = function(x) {
        paste0(
          formatC(abs(x), format = "f", digits = 3),
          "°W"
        )
      }
    ) +

    scale_y_continuous(
      breaks = y_mid,
      labels = function(y) {
        paste0(
          formatC(abs(y), format = "f", digits = 3),
          "°N"
        )
      }
    ) +

    coord_sf(
      xlim = x_limits,
      ylim = y_limits,
      expand = FALSE,
      crs = st_crs(4326),
      default_crs = st_crs(4326)
    ) +

    labs(
      title = region_code,
      x = NULL,
      y = NULL,
      caption = "© OpenStreetMap contributors"
    ) +

    theme_classic() +

    theme(
      legend.position = "none",

      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        size = 14
      ),

      plot.caption = element_text(
        size = 6,
        colour = "grey40"
      ),

      axis.text.x = element_text(
        colour = "grey20",
        size = 18,
        margin = margin(t = 5)
      ),

      axis.text.y = element_text(
        colour = "grey20",
        size = 18,
        angle = 90,
        vjust = 0.5,
        margin = margin(r = 5)
      ),

      axis.ticks = element_blank(),

      panel.border = element_rect(
        colour = "grey40",
        fill = NA
      )
    )
}

p_MAG <- make_region_map("MAG", padding_lon = .035, padding_lat = .02, tile_zoom = 0)
p_MAG
p_PEI <- make_region_map("PEI", padding_lon = .035, padding_lat = .02, tile_zoom = 0)
p_PEI
p_HAL <- make_region_map("HAL", padding_lon = .035, padding_lat = .02, tile_zoom = 0)
p_HAL
p_BOF <- make_region_map("BOF", padding_lon = .035, padding_lat = .02, tile_zoom = 0)
p_BOF
p_GOM <- make_region_map("GOM", padding_lon = .035, padding_lat = .02, tile_zoom = 0)
p_GOM



# Create an output folder if it does not already exist
dir.create(
  "code_2026/map/regional_maps",
  recursive = TRUE,
  showWarnings = FALSE
)

p_MAG <- make_region_map(
  "MAG",
  padding_lon = 0.035,
  padding_lat = 0.02,
  tile_zoom = 0
)

p_PEI <- make_region_map(
  "PEI",
  padding_lon = 0.035,
  padding_lat = 0.02,
  tile_zoom = 0
)

p_HAL <- make_region_map(
  "HAL",
  padding_lon = 0.035,
  padding_lat = 0.02,
  tile_zoom = 0
)

p_BOF <- make_region_map(
  "BOF",
  padding_lon = 0.035,
  padding_lat = 0.02,
  tile_zoom = 0
)

p_GOM <- make_region_map(
  "GOM",
  padding_lon = 0.035,
  padding_lat = 0.02,
  tile_zoom = 0
)

ggsave(
  filename = "code_2026/map/regional_maps/MAG_map.png",
  plot = p_MAG,
  width = 8,
  height = 5,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = "code_2026/map/regional_maps/PEI_map.png",
  plot = p_PEI,
  width = 8,
  height = 5,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = "code_2026/map/regional_maps/HAL_map.png",
  plot = p_HAL,
  width = 8,
  height = 5,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = "code_2026/map/regional_maps/BOF_map.png",
  plot = p_BOF,
  width = 8,
  height = 5,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = "code_2026/map/regional_maps/GOM_map.png",
  plot = p_GOM,
  width = 8,
  height = 5,
  units = "in",
  dpi = 400,
  bg = "white"
)
