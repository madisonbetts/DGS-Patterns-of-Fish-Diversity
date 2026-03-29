#-----------------
# COAS-1
# Cottus asper
# prickly sculpin
#-----------------

library(ggplot2)
library(dplyr)
library(geosphere)
library(maps)

# -----------------------------
# site IDs
# 1 McCloud Arm, Shasta Reservoir (McC M)
# 2 Napa River (NAP)
# 3 Waddell Creek (WAD)
# 4 Suisun Marsh (SUS)
# 5 Putah Creek, Russell Ranch (PUT RR)
# 6 Big Canyon Creek, Putah (PUT BC)
# -----------------------------

# -----------------------------
# coordinates
# exact sampling coordinates from Table 2
# -----------------------------
COAS_1_coords <- data.frame(
  site_id = 1:6,
  lat = c(
    40.932808,
    38.268000,
    37.099517,
    38.070000,
    38.533401,
    38.833417
  ),
  lon = c(
    -122.24928,
    -122.28400,
    -122.27549,
    -122.07000,
    -121.85026,
    -122.65006
  )
)

# -----------------------------
# pairwise FST matrix
# uncorrected pairwise FST between sampling locations
# extracted from authors' fst matrix / Supplementary Fig. S2
# -----------------------------
COAS_1_fst <- matrix(
  c(
    0.00, 0.81, 0.85, 0.84, 0.84, 0.91,
    0.81, 0.00, 0.22, 0.12, 0.18, 0.54,
    0.85, 0.22, 0.00, 0.28, 0.36, 0.65,
    0.84, 0.12, 0.28, 0.00, 0.21, 0.59,
    0.84, 0.18, 0.36, 0.21, 0.00, 0.46,
    0.91, 0.54, 0.65, 0.59, 0.46, 0.00
  ),
  nrow = 6,
  byrow = TRUE
)

rownames(COAS_1_fst) <- colnames(COAS_1_fst) <- as.character(COAS_1_coords$site_id)
COAS_1_fst[COAS_1_fst < 0] <- 0
diag(COAS_1_fst) <- 0

# -----------------------------
# plotting labels only
# -----------------------------
plot_df <- COAS_1_coords
plot_df$label <- c(
  "1. McCloud Arm",
  "2. Napa River",
  "3. Waddell Creek",
  "4. Suisun Marsh",
  "5. Putah RR",
  "6. Big Canyon Creek"
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
    nudge_x = 0.08,
    size = 3
  ) +
  coord_fixed(1.25, xlim = c(-123.3, -121.4), ylim = c(36.7, 41.3)) +
  theme_bw() +
  labs(
    title = "COAS-1 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

print(site_map)

# -----------------------------
# pairwise geographic distance
# -----------------------------
coords_mat <- as.matrix(COAS_1_coords[, c("lon", "lat")])
geo_dist_km <- geosphere::distm(coords_mat, fun = geosphere::distGeo) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- as.character(COAS_1_coords$site_id)

# -----------------------------
# IBD dataframe + plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(COAS_1_fst)[row(COAS_1_fst)[upper.tri(COAS_1_fst)]],
  site2   = colnames(COAS_1_fst)[col(COAS_1_fst)[upper.tri(COAS_1_fst)]],
  fst     = COAS_1_fst[upper.tri(COAS_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

ibd_plot <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    title = "COAS-1 isolation by distance",
    x = "Geographic distance (km)",
    y = expression(F[ST])
  )

print(ibd_plot)

# -----------------------------
# save RData
# -----------------------------
save(
  COAS_1_fst,
  COAS_1_coords,
  file = file.path("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/COAS-1/data/COAS-1.RData")
)
