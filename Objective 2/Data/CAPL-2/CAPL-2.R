# ============================================================
# CAPL-2 | Catostomus platyrhynchus / mountain sucker dataset
# Black Hills, South Dakota
# coordinates tightened to best-available stream-reach estimates
# ============================================================

# -----------------------------
# 0) libraries
# -----------------------------
library(ggplot2)
library(geosphere)
library(dplyr)
library(tidyr)
library(maps)
library(mapdata)

# -----------------------------
# 1) set working directory
# -----------------------------
setwd("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CAPL-2")

# create data subfolder if needed
if (!dir.exists("data")) dir.create("data")

# -----------------------------
# 2) site coordinates
# Annie Creek excluded because it was not used in the FST analyses
# site order matches Table 2 in the paper
# -----------------------------
CAPL_2_coords <- data.frame(
  site_id = 1:5,
  site_name = c(
    "Boxelder Creek",
    "Bear Butte Creek",
    "Elk Creek",
    "North Fork Rapid Creek",
    "Whitewood Creek"
  ),
  lat = c(
    44.13165,
    44.40650,
    44.29220,
    44.04890,
    44.45870
  ),
  lon = c(
    -103.29879,
    -103.37640,
    -103.47030,
    -103.33920,
    -103.71080
  )
)

# -----------------------------
# 3) pairwise FST matrix
# values transcribed from Table 2
# order:
# 1 Boxelder Creek
# 2 Bear Butte Creek
# 3 Elk Creek
# 4 North Fork Rapid Creek
# 5 Whitewood Creek
# -----------------------------
CAPL_2_fst <- matrix(c(
  0,      0.057, 0.039, 0.027, 0.049,
  0.057,  0,     0.048, 0.074, 0.037,
  0.039,  0.048, 0,     0.026, 0.022,
  0.027,  0.074, 0.026, 0,     0.040,
  0.049,  0.037, 0.022, 0.040, 0
), nrow = 5, byrow = TRUE)

rownames(CAPL_2_fst) <- CAPL_2_coords$site_id
colnames(CAPL_2_fst) <- CAPL_2_coords$site_id

# enforce workflow rule
CAPL_2_fst[CAPL_2_fst < 0] <- 0

# -----------------------------
# 4) geographic distance matrix
# great-circle distance in km
# -----------------------------
coord_mat <- as.matrix(CAPL_2_coords[, c("lon", "lat")])
CAPL_2_dist_km <- geosphere::distm(coord_mat, fun = geosphere::distHaversine) / 1000
rownames(CAPL_2_dist_km) <- CAPL_2_coords$site_id
colnames(CAPL_2_dist_km) <- CAPL_2_coords$site_id

# -----------------------------
# 5) long-format dataframe for IBD plot
# -----------------------------
fst_long <- as.data.frame(as.table(CAPL_2_fst), stringsAsFactors = FALSE)
names(fst_long) <- c("site1", "site2", "fst")

dist_long <- as.data.frame(as.table(CAPL_2_dist_km), stringsAsFactors = FALSE)
names(dist_long) <- c("site1", "site2", "dist_km")

ibd_df <- fst_long %>%
  left_join(dist_long, by = c("site1", "site2")) %>%
  mutate(
    site1 = as.integer(site1),
    site2 = as.integer(site2)
  ) %>%
  filter(site1 < site2)

# -----------------------------
# 6) map of sites
# US + Canada background
# zoom to extent of points
# do not save plot
# -----------------------------
us_map <- map_data("state")
world_map <- map_data("world")

x_pad <- 0.75
y_pad <- 0.50

x_rng <- range(CAPL_2_coords$lon)
y_rng <- range(CAPL_2_coords$lat)

ggplot() +
  geom_polygon(
    data = subset(world_map, region %in% c("USA", "Canada")),
    aes(x = long, y = lat, group = group),
    fill = "gray95",
    color = "gray70",
    linewidth = 0.3
  ) +
  geom_polygon(
    data = us_map,
    aes(x = long, y = lat, group = group),
    fill = "gray90",
    color = "gray65",
    linewidth = 0.3
  ) +
  geom_point(
    data = CAPL_2_coords,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = CAPL_2_coords,
    aes(x = lon, y = lat, label = site_name),
    nudge_y = 0.03,
    size = 3
  ) +
  coord_quickmap(
    xlim = c(x_rng[1] - x_pad, x_rng[2] + x_pad),
    ylim = c(y_rng[1] - y_pad, y_rng[2] + y_pad)
  ) +
  labs(
    title = "CAPL-2 sampling sites",
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_bw()

# -----------------------------
# 7) IBD plot
# do not save plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    title = "CAPL-2 IBD plot",
    x = "Geographic distance (km)",
    y = expression(F[ST])
  ) +
  theme_bw()

# -----------------------------
# 8) save objects
# -----------------------------
save(
  CAPL_2_fst,
  CAPL_2_coords,
  file = "data/CAPL-2.RData"
)