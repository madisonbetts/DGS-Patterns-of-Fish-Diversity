#-----------------
# COGU-1
# Cottus gulosus
# inland riffle sculpin
#-----------------

library(ggplot2)
library(dplyr)
library(geosphere)
library(maps)

# -----------------------------
# site IDs
# 1 Merced River (MER)
# 2 Mokelumne River (MOK)
# 3 Hot Springs Creek, Lassen (HSS)
# 4 Sacramento River at Cantara (SAC)
# -----------------------------

# -----------------------------
# coordinates
# exact sampling coordinates from Table 2
# -----------------------------
COGU_1_coords <- data.frame(
  site_id = 1:4,
  lat = c(
    37.6532,
    38.3159,
    40.435258,
    41.2661
  ),
  lon = c(
    -119.7825,
    -120.7104,
    -121.36544,
    -122.3078
  )
)

# -----------------------------
# pairwise FST matrix
# uncorrected pairwise FST between sampling locations
# extracted from authors' fst matrix / Supplementary Fig. S2
# -----------------------------
COGU_1_fst <- matrix(
  c(
    0.00, 0.45, 0.63, 0.60,
    0.45, 0.00, 0.80, 0.73,
    0.63, 0.80, 0.00, 0.71,
    0.60, 0.73, 0.71, 0.00
  ),
  nrow = 4,
  byrow = TRUE
)

rownames(COGU_1_fst) <- colnames(COGU_1_fst) <- as.character(COGU_1_coords$site_id)
COGU_1_fst[COGU_1_fst < 0] <- 0
diag(COGU_1_fst) <- 0

# -----------------------------
# plotting labels only
# -----------------------------
plot_df <- COGU_1_coords
plot_df$label <- c(
  "1. Merced River",
  "2. Mokelumne River",
  "3. Hot Springs Creek",
  "4. Sacramento River"
)

# -----------------------------
# site map
# -----------------------------
usa <- ggplot2::map_data("world", region = "USA")
can <- ggplot2::map_data("world", region = "Canada")
world_map <- bind_rows(usa, can)

site_map <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "grey90",
    color = "grey50",
    linewidth = 0.2
  ) +
  geom_point(
    data = plot_df,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = plot_df,
    aes(x = lon, y = lat, label = label),
    hjust = 0,
    nudge_x = 0.10,
    size = 3
  ) +
  coord_fixed(1.25, xlim = c(-124.5, -118.0), ylim = c(36.5, 42.3)) +
  theme_bw() +
  labs(
    title = "COGU-1 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

print(site_map)

# -----------------------------
# pairwise geographic distance
# -----------------------------
coords_mat <- as.matrix(COGU_1_coords[, c("lon", "lat")])
geo_dist_km <- geosphere::distm(coords_mat, fun = geosphere::distGeo) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- as.character(COGU_1_coords$site_id)

# -----------------------------
# IBD dataframe + plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(COGU_1_fst)[row(COGU_1_fst)[upper.tri(COGU_1_fst)]],
  site2   = colnames(COGU_1_fst)[col(COGU_1_fst)[upper.tri(COGU_1_fst)]],
  fst     = COGU_1_fst[upper.tri(COGU_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

ibd_plot <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    title = "COGU-1 isolation by distance",
    x = "Geographic distance (km)",
    y = expression(F[ST])
  )

print(ibd_plot)

# -----------------------------
# save RData
# -----------------------------
save(
  COGU_1_fst,
  COGU_1_coords,
  file = file.path("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/COGU-1/data/COGU-1.RData")
)
