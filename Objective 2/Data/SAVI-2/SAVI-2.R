# -----------------------------
# SAVI-2 walleye (Sander vitreus)
# Stepien et al. 2009
# Population-level coordinates + pairwise FST matrix + map + IBD plot
# Save only SAVI_2_fst and SAVI_2_coords to SAVI-2.RData
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)
library(sf)
library(tigris)

options(tigris_use_cache = TRUE)

# ------------------------------
# paths
# ------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAVI-2"
data_dir <- file.path(base_dir, "data")

if (!dir.exists(data_dir)) dir.create(data_dir, recursive = TRUE)

# -----------------------------
# 0) site coordinates
# Table 2 reports population-level FST among 16 groups:
# A, B, C, D, E, F, G, H, I, J, K, L, M, N, O-P, Q
# Table 1 gives exact coordinates for A-J, M, N, O, P, Q.
# K, L, and O-P are composite groups in Table 2, so their
# coordinates are represented here by centroids of the
# component localities listed in Table 1.
# Numeric site IDs are used for matrix row/column names.
# -----------------------------

SAVI_2_coords <- data.frame(
  site = as.character(1:16),
  pop_code = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O-P", "Q"),
  locality = c(
    "Cedar Lake, MB",
    "McKim Lake, ON",
    "Mille Lacs, MN",
    "St. Louis River, MN",
    "Muskegon River, MI",
    "Alpena, MI",
    "Flint River, MI",
    "Moon & Musquash Rivers, ON",
    "Thames River, ON",
    "Detroit River, MI",
    "Western & Central Lake Erie basins centroid",
    "Eastern Lake Erie basin centroid",
    "Bay of Quinte, ON",
    "Oneida Lake, NY",
    "Ohio + Upper New River centroid",
    "North River, AL"
  ),
  lat = c(
    53.33,
    50.82,
    46.37,
    46.73,
    43.48,
    45.02,
    43.33,
    44.84,
    42.32,
    43.24,
    mean(c(42.09, 41.56, 41.63, 41.81, 41.72, 41.46, 41.78)),
    mean(c(42.46, 42.57, 42.81, 42.86)),
    44.16,
    43.28,
    mean(c(39.67, 36.83)),
    33.26
  ),
  lon = c(
    -100.10,
    -92.50,
    -93.19,
    -92.13,
    -85.83,
    -83.43,
    -84.05,
    -79.80,
    -82.45,
    -82.79,
    mean(c(-83.29, -83.65, -83.02, -82.79, -82.61, -82.89, -81.25)),
    mean(c(-79.41, -79.13, -78.86, -79.58)),
    -77.37,
    -75.44,
    mean(c(-80.87, -80.67)),
    -87.51
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# Plot full U.S. state outlines for context and include
# surrounding Canada/USA background across the study extent.
# -----------------------------
world_df <- ggplot2::map_data("world")
world_df <- world_df[world_df$region %in% c("USA", "Canada"), ]

states_sf <- tigris::states(cb = TRUE, year = 2022)
states_sf <- sf::st_transform(states_sf, 4326)

x_pad <- 6
y_pad <- 3
xlim_use <- range(SAVI_2_coords$lon) + c(-x_pad, x_pad)
ylim_use <- range(SAVI_2_coords$lat) + c(-y_pad, y_pad)

site_map <- ggplot() +
  geom_polygon(
    data = world_df,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey75",
    linewidth = 0.2
  ) +
  geom_sf(
    data = states_sf,
    fill = NA,
    color = "grey45",
    linewidth = 0.25,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = SAVI_2_coords,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text(
    data = SAVI_2_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.45,
    size = 3.1
  ) +
  coord_sf(xlim = xlim_use, ylim = ylim_use, expand = FALSE) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "SAVI-2 sampling locations"
  )

print(site_map)

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(SAVI_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- SAVI_2_coords$site
colnames(geo_dist_km) <- SAVI_2_coords$site

# -----------------------------
# 3) helper function:
# takes a vector of lower-triangle values (row-wise order)
# and reconstructs a full symmetric matrix
# -----------------------------
fill_sym_from_lower <- function(pops, vals, diag_val = 0) {
  n <- length(pops)
  stopifnot(length(vals) == n * (n - 1) / 2)

  mat <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(pops, pops)
  )

  k <- 1
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      mat[i, j] <- vals[k]
      mat[j, i] <- vals[k]
      k <- k + 1
    }
  }

  diag(mat) <- diag_val
  mat
}

# -----------------------------
# 4) pairwise FST matrix
# Values transcribed from Table 2 lower diagonal (hST/FST analogue)
# in row-wise order for the 16 Table 2 population groups.
# Any negative values are forced to 0.
# -----------------------------
fst_vals <- c(
  0.131,
  0.126, 0.200,
  0.083, 0.160, 0.110,
  0.123, 0.217, 0.162, 0.050,
  0.119, 0.253, 0.169, 0.069, 0.028,
  0.132, 0.231, 0.168, 0.058, 0.010, 0.010,
  0.188, 0.273, 0.214, 0.087, 0.048, 0.067, 0.046,
  0.138, 0.242, 0.145, 0.051, 0.032, 0.025, 0.013, 0.039,
  0.134, 0.252, 0.161, 0.062, 0.033, 0.022, 0.017, 0.052, 0.013,
  0.157, 0.228, 0.166, 0.067, 0.055, 0.060, 0.039, 0.072, 0.039, 0.042,
  0.131, 0.193, 0.133, 0.044, 0.040, 0.048, 0.027, 0.057, 0.024, 0.034, 0.008,
  0.126, 0.232, 0.138, 0.061, 0.039, 0.039, 0.028, 0.069, 0.026, 0.038, 0.059, 0.043,
  0.137, 0.228, 0.163, 0.056, 0.037, 0.064, 0.029, 0.060, 0.035, 0.062, 0.060, 0.040, 0.049,
  0.122, 0.216, 0.138, 0.050, 0.027, 0.038, 0.019, 0.050, 0.010, 0.028, 0.043, 0.023, 0.032, 0.018,
  0.171, 0.295, 0.206, 0.102, 0.067, 0.081, 0.068, 0.146, 0.114, 0.096, 0.096, 0.069, 0.085, 0.122, 0.074
)

fst_vals[fst_vals < 0] <- 0

SAVI_2_fst <- fill_sym_from_lower(
  pops = SAVI_2_coords$site,
  vals = fst_vals,
  diag_val = 0
)

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(SAVI_2_fst)[row(SAVI_2_fst)[upper.tri(SAVI_2_fst)]],
  site2   = colnames(SAVI_2_fst)[col(SAVI_2_fst)[upper.tri(SAVI_2_fst)]],
  fst     = SAVI_2_fst[upper.tri(SAVI_2_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 6) IBD plot
# -----------------------------
ibd_plot <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "SAVI-2 isolation by distance"
  )

print(ibd_plot)

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(SAVI_2_fst), SAVI_2_coords$site))
stopifnot(identical(colnames(SAVI_2_fst), SAVI_2_coords$site))
stopifnot(all(SAVI_2_fst >= 0))

# -----------------------------
# 8) save only requested objects
# -----------------------------
save(
  SAVI_2_fst,
  SAVI_2_coords,
  file = file.path(data_dir, "SAVI-2.RData")
)
