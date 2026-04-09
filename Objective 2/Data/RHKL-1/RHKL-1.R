#-----------------
# RHKL-1
# Rhinichthys klamathensis
# western speckled dace / Foskett dace system
#-----------------

# -----------------------------
# site IDs from the study
# -----------------------------
# 1 Foskett Spring
# 2 Dace Spring
# 3 Coleman Creek
# 4 Deep Creek
# 5 Twentymile Creek

# -----------------------------
# coordinates
# Foskett Spring anchored to DOGAMI spring locality
# the remaining four sites were digitized from the georeferenced
# Figure 3 map in Sidlauskas et al. (2024) and aligned to that
# exact Foskett Spring point
# -----------------------------
RHKL_1_coords <- data.frame(
  site_id = 1:5,
  lat = c(
    42.0697500,  # 1 Foskett Spring
    42.0558161,  # 2 Dace Spring
    41.9704686,  # 3 Coleman Creek
    42.1734399,  # 4 Deep Creek
    42.0631677   # 5 Twentymile Creek
  ),
  lon = c(
    -119.8396660, # 1 Foskett Spring
    -119.8393133, # 2 Dace Spring
    -119.7841430, # 3 Coleman Creek
    -119.9362790, # 4 Deep Creek
    -119.9552806  # 5 Twentymile Creek
  )
)

# -----------------------------
# pairwise FST matrix
# transcribed from Table 3
# table is presented as a lower triangle and rebuilt here
# as a full symmetric matrix
# -----------------------------
n_sites <- 5
RHKL_1_fst <- matrix(0, nrow = n_sites, ncol = n_sites)
rownames(RHKL_1_fst) <- colnames(RHKL_1_fst) <- as.character(1:5)

# site order in the matrix:
# 1 Foskett Spring
# 2 Dace Spring
# 3 Coleman Creek
# 4 Deep Creek
# 5 Twentymile Creek

lt <- list(
  `2` = c(0.0055),
  `3` = c(0.0671, 0.0704),
  `4` = c(0.1605, 0.1644, 0.1275),
  `5` = c(0.1467, 0.1499, 0.1128, 0.0604)
)

for (i in 2:n_sites) {
  RHKL_1_fst[i, 1:(i - 1)] <- lt[[as.character(i)]]
}

RHKL_1_fst <- RHKL_1_fst + t(RHKL_1_fst)
RHKL_1_fst[RHKL_1_fst < 0] <- 0
diag(RHKL_1_fst) <- 0

# -----------------------------
# plotting packages
# -----------------------------
library(ggplot2)
library(dplyr)
library(geosphere)

# -----------------------------
# site labels for plotting only
# not written into RHKL_1_coords
# -----------------------------
plot_labs <- c(
  "1. Foskett Spring",
  "2. Dace Spring",
  "3. Coleman Creek",
  "4. Deep Creek",
  "5. Twentymile Creek"
)

plot_df <- RHKL_1_coords
plot_df$label <- plot_labs
plot_df$nudge_x <- c(0.008, 0.008, 0.008, 0.008, -0.060)
plot_df$nudge_y <- c(0.008, -0.006, -0.012, 0.010, 0.012)
plot_df$hjust   <- c(0, 0, 0, 0, 0)

# -----------------------------
# base map: US + Canada
# -----------------------------
usa <- ggplot2::map_data("world", region = "USA")
can <- ggplot2::map_data("world", region = "Canada")
world_map <- bind_rows(usa, can)

# -----------------------------
# site map
# -----------------------------
ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "grey90", color = "grey50", linewidth = 0.2
  ) +
  geom_point(
    data = plot_df,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = plot_df,
    aes(
      x = lon + nudge_x,
      y = lat + nudge_y,
      label = label,
      hjust = hjust
    ),
    size = 3
  ) +
  coord_fixed(
    1.25,
    xlim = c(-120.05, -119.72),
    ylim = c(41.94, 42.20)
  ) +
  theme_bw() +
  labs(
    title = "RHKL-1 sampling sites",
    x = "Longitude",
    y = "Latitude"
  )

# -----------------------------
# pairwise geographic distance
# Euclidean / great-circle distance in km
# -----------------------------
coords_mat <- as.matrix(RHKL_1_coords[, c("lon", "lat")])
geo_dist_km <- geosphere::distm(coords_mat, fun = geosphere::distGeo) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- as.character(RHKL_1_coords$site_id)

# -----------------------------
# pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(RHKL_1_fst)[row(RHKL_1_fst)[upper.tri(RHKL_1_fst)]],
  site2   = colnames(RHKL_1_fst)[col(RHKL_1_fst)[upper.tri(RHKL_1_fst)]],
  fst     = RHKL_1_fst[upper.tri(RHKL_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_bw() +
  labs(
    title = "RHKL-1 isolation by distance",
    x = "Geographic distance (km)",
    y = expression(F[ST])
  )

# -----------------------------
# save RData
# -----------------------------
save(
  RHKL_1_fst,
  RHKL_1_coords,
  file = file.path("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/RHKL-1/data/RHKL-1.RData")
)
