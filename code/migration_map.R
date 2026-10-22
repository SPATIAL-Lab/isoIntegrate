library(dplyr)
library(ggplot2)
library(png)
library(sf)

# Read in Data
## Cleaned US Census Data 
## From: Table 1. State-to-State Migration Flows: 
### 2023 American Community Survey 1-Year Estimates
dest_orig <-read.csv('data/dest_orig.csv')

## US States Shapefile
states_shape <- st_read('data/states.shp', quiet = TRUE) %>%
  ### Isolate polygon and state abbreviation
  select(STUSPS, geometry)

## Enamel linear regression d18O isoscape
iso <- readPNG('data/enamLP.png')

# Create Bounding Box for States Shapefile
states_box <- st_bbox(states_shape)

# Crop Isoscape to Bounding Box
iso_mask <- apply(iso, c(1, 2), function(px) any(px < 0.99))
iso_rows <- which(apply(iso_mask, 1, any))
iso_cols <- which(apply(iso_mask, 2, any))
iso <- iso[min(iso_rows):max(iso_rows), min(iso_cols):max(iso_cols), , drop = FALSE]

# Isolate the 15 Largest Migration Flows
migration.15 <- dest_orig %>%
  arrange(desc(movers)) %>%
  slice_head(n=15)

# Define Centroid Location to Anchor Flow Arrows
state_centroids <- states_shape %>%
  mutate(
    abbr = STUSPS,
    centroid = st_point_on_surface(geometry),
    lon = st_coordinates(centroid)[, 1],
    lat = st_coordinates(centroid)[, 2]
  ) %>%
  st_drop_geometry()

# Join Centroid Coordinates to Migration Flows
flows <- migration.15 %>%
  left_join(state_centroids, by = c("origin" = "abbr")) %>%
  rename(lon_orig = lon, lat_orig = lat) %>%
  left_join(state_centroids, by = c("destination" = "abbr")) %>%
  rename(lon_dest = lon, lat_dest = lat)

# Map Creation
map <- ggplot() +
  ## PNG isoscape
  annotation_raster(
    iso,
    xmin = states_box["xmin"],
    xmax = states_box["xmax"],
    ymin = states_box["ymin"],
    ymax = states_box["ymax"]
  ) +
  ## State boundaries
  geom_sf(data = states_shape, fill = NA, color = "black", size = 0.4) +
  ## Migration flows as arrows
  geom_curve(
    data = flows,
    aes(x = lon_orig, y = lat_orig, xend = lon_dest, 
        yend = lat_dest, linewidth = movers),
    curvature = 0.33,
    color = "darkred",
    arrow = arrow(length = unit(0.08, "inches"), type = "closed")
  ) +
  scale_linewidth_continuous(range = c(0.1, 0.8), guide = "none") +
  coord_sf(
    xlim = c(states_box[["xmin"]], states_box[["xmax"]]),
    ylim = c(states_box[["ymin"]], states_box[["ymax"]]),
    expand = FALSE
  ) +
  ## Labels and Legend
  labs(
    x = NULL,
    y = NULL,
    linewidth = "Number of Movers"
  ) + 
  ## Styling
  theme_bw(base_size = 7) +
  theme(
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.key.width = unit(0.8, "lines"),
    legend.key.height = unit(0.8, "lines"),
    legend.position = "bottom",
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  )
