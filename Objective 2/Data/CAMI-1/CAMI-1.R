# -----------------------------
# CAMI-1 Modoc sucker
# Catostomus microps
# site coordinates + FST matrix + IBD plot
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)

save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CAMI-1"

# -----------------------------
# 0) site coordinates
# best-available coordinates inferred from Figure 2 + named streams
# original sample order from paper:
# 1 Dutch Flat
# 2 Johnson
# 3 Thomas
# 4 Coffee Mill
# 5 Garden Gulch
# 6 Hulbert
# 7 Turner-mid_a
# 8 Turner-upp
# 9 Washington
# -----------------------------
CAMI_1_coords <- data.frame(
  site = as.character(1:9),
  #site_name = c("Dutch Flat", "Johnson", "Thomas", "Coffee Mill", "Garden Gulch", "Hulbert", "Turner-mid_a", "Turner-upp", "Washington"),
  lat = c(41.106, 41.115, 42.026, 41.217, 41.227, 41.246, 41.224, 41.241, 41.235),
  lon = c(-121.187, -121.164, -120.885, -121.139, -121.151, -121.177, -121.152, -121.142, -121.133),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
states_map <- map_data("state")
ca_or_map <- subset(states_map, region %in% c("california", "oregon"))

map_xlim <- c(-121.35, -120.75)
map_ylim <- c(41.00, 42.10)

ggplot() +
  geom_polygon(
    data = ca_or_map,
    aes(x = long, y = lat, group = group),
    fill = "grey92", color = "grey40", linewidth = 0.3
  ) +
  geom_point(
    data = CAMI_1_coords,
    aes(x = lon, y = lat),
    size = 2.5
  ) +
  geom_text(
    data = CAMI_1_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.015,
    size = 3.2
  ) +
  coord_fixed(xlim = map_xlim, ylim = map_ylim) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "CAMI-1 sampling sites"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(CAMI_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- CAMI_1_coords$site
colnames(geo_dist_km) <- CAMI_1_coords$site

# -----------------------------
# 3) pairwise FST matrix
# Appendix S1; above-diagonal FST values
# sample 10 was reassigned genetically to Sacramento sucker,
# so CAMI-1 uses samples 1-9 only
# -----------------------------
CAMI_1_fst <- matrix(c(
  0.000, 0.073, 0.226, 0.303, 0.386, 0.257, 0.275, 0.279, 0.276,
  0.073, 0.000, 0.198, 0.280, 0.376, 0.230, 0.242, 0.238, 0.251,
  0.226, 0.198, 0.000, 0.192, 0.269, 0.116, 0.185, 0.175, 0.162,
  0.303, 0.280, 0.192, 0.000, 0.324, 0.112, 0.069, 0.169, 0.028,
  0.386, 0.376, 0.269, 0.324, 0.000, 0.268, 0.311, 0.453, 0.274,
  0.257, 0.230, 0.116, 0.112, 0.268, 0.000, 0.121, 0.112, 0.085,
  0.275, 0.242, 0.185, 0.069, 0.311, 0.121, 0.000, 0.125, 0.037,
  0.279, 0.238, 0.175, 0.169, 0.453, 0.112, 0.125, 0.000, 0.130,
  0.276, 0.251, 0.162, 0.028, 0.274, 0.085, 0.037, 0.130, 0.000
), nrow = 9, byrow = TRUE)

rownames(CAMI_1_fst) <- colnames(CAMI_1_fst) <- CAMI_1_coords$site
CAMI_1_fst[CAMI_1_fst < 0] <- 0
diag(CAMI_1_fst) <- 0

# -----------------------------
# 4) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(CAMI_1_fst)[row(CAMI_1_fst)[upper.tri(CAMI_1_fst)]],
  site2   = colnames(CAMI_1_fst)[col(CAMI_1_fst)[upper.tri(CAMI_1_fst)]],
  fst     = CAMI_1_fst[upper.tri(CAMI_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 5) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "CAMI-1 isolation by distance"
  )

# -----------------------------
# 6) quick checks
# -----------------------------
stopifnot(identical(rownames(CAMI_1_fst), CAMI_1_coords$site))
stopifnot(identical(colnames(CAMI_1_fst), CAMI_1_coords$site))
stopifnot(isTRUE(all.equal(CAMI_1_fst, t(CAMI_1_fst))))

# -----------------------------
# 7) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  CAMI_1_fst,
  CAMI_1_coords,
  file = file.path(out_dir, "CAMI-1.RData")
)
