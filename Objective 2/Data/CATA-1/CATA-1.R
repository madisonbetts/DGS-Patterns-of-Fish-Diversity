# -----------------------------
# CATA-1 Tahoe sucker
# Truckee River, 2008
# pairwise FST matrix + site coordinates + map + IBD plot
# -----------------------------

library(ggplot2)
library(ggrepel)
library(geosphere)
library(maps)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CATA-1"

# -----------------------------
# 0) site coordinates
# Table 1 in Peacock et al. 2016
# site numbering matches FST matrix order below
# -----------------------------
CATA_1_coords <- data.frame(
  site = as.character(1:14),
  site_name = c(
    "VerdiB", "RenoC", "RenoB", "RenoA",
    "McRC", "McRB", "McRA", "McRD",
    "WadsA", "WadsB", "WadsC",
    "NixonA", "NixonB", "NixonC"
  ),
  lat = c(
    39.507, 39.524, 39.530, 39.514,
    39.526, 39.545, 39.547, 39.565,
    39.585, 39.591, 39.613,
    39.727, 39.818, 39.854
  ),
  lon = c(
    -119.902, -119.819, -119.795, -119.736,
    -119.610, -119.587, -119.573, -119.487,
    -119.444, -119.368, -119.306,
    -119.319, -119.350, -119.394
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# Plot US + Canada, zoom to point extent
# Do not save plot
# -----------------------------
usa <- map_data("state")
can <- map_data("world", region = "Canada")

xpad <- 1.0
ypad <- 0.6

xlim_map <- c(
  min(CATA_1_coords$lon) - xpad,
  max(CATA_1_coords$lon) + xpad
)

ylim_map <- c(
  min(CATA_1_coords$lat) - ypad,
  max(CATA_1_coords$lat) + ypad
)

ggplot() +
  geom_polygon(
    data = usa,
    aes(x = long, y = lat, group = group),
    fill = "grey92", color = "black", linewidth = 0.25
  ) +
  geom_polygon(
    data = can,
    aes(x = long, y = lat, group = group),
    fill = "grey86", color = "black", linewidth = 0.25
  ) +
  geom_point(
    data = CATA_1_coords,
    aes(x = lon, y = lat),
    size = 2.7
  ) +
  ggrepel::geom_text_repel(
    data = CATA_1_coords,
    aes(x = lon, y = lat, label = site),
    size = 3.4,
    seed = 1,
    min.segment.length = 0
  ) +
  coord_fixed(
    xlim = xlim_map,
    ylim = ylim_map
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "CATA-1 sampling sites"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(CATA_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- CATA_1_coords$site
colnames(geo_dist_km) <- CATA_1_coords$site

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
# S5 Table (B): Tahoe sucker, 2008
# lower triangle, row-wise order
#
# site order:
# 1 VerdiB
# 2 RenoC
# 3 RenoB
# 4 RenoA
# 5 McRC
# 6 McRB
# 7 McRA
# 8 McRD
# 9 WadsA
# 10 WadsB
# 11 WadsC
# 12 NixonA
# 13 NixonB
# 14 NixonC
#
# NixonB vs NixonC cell was not visible in the parsed table;
# set to 0 per user instruction.
# -----------------------------
fst_vals <- c(

  # 2
  -0.036,

  # 3
  -0.031, -0.039,

  # 4
  -0.008, -0.003, -0.016,

  # 5
  -0.006, 0.026, -0.013, 0.023,

  # 6
  0.017, 0.029, 0.003, 0.040, 0.008,

  # 7
  0.009, 0.020, 0.003, 0.029, 0.009, 0.004,

  # 8
  -0.003, 0.000, -0.009, 0.003, 0.015, 0.038, 0.027,

  # 9
  0.020, 0.028, 0.011, 0.040, 0.006, 0.003, -0.001, 0.039,

  # 10
  0.028, 0.036, 0.020, 0.038, 0.011, 0.003, 0.006, 0.040, 0.001,

  # 11
  0.020, 0.030, 0.013, 0.041, 0.007, -0.003, 0.004, 0.037, 0.003, 0.002,

  # 12
  -0.002, 0.022, 0.000, 0.031, 0.012, 0.002, 0.001, 0.035, 0.000, 0.007, 0.006,

  # 13
  0.014, 0.019, -0.004, 0.032, 0.009, -0.002, 0.001, 0.029, 0.006, -0.001, 0.001, 0.004,

  # 14
  -0.015, 0.001, 0.001, 0.019, 0.010, 0.002, -0.001, 0.020, 0.002, 0.009, 0.000, -0.001, 0
)

CATA_1_fst <- fill_sym_from_lower(
  pops = CATA_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

# always make fst < 0 equal to 0
CATA_1_fst[CATA_1_fst < 0] <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(CATA_1_fst)[row(CATA_1_fst)[upper.tri(CATA_1_fst)]],
  site2   = colnames(CATA_1_fst)[col(CATA_1_fst)[upper.tri(CATA_1_fst)]],
  fst     = CATA_1_fst[upper.tri(CATA_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
)

# -----------------------------
# 6) IBD plot
# Do not save plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "CATA-1 isolation by distance"
  )

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(CATA_1_fst), CATA_1_coords$site))
stopifnot(identical(colnames(CATA_1_fst), CATA_1_coords$site))
stopifnot(isTRUE(all.equal(CATA_1_fst, t(CATA_1_fst))))
stopifnot(length(fst_vals) == nrow(CATA_1_coords) * (nrow(CATA_1_coords) - 1) / 2)

# -----------------------------
# 8) save RData
# final saved coords df only keeps site + lat/lon
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

CATA_1_coords <- CATA_1_coords[, c("site", "lat", "lon")]

save(
  CATA_1_fst,
  CATA_1_coords,
  file = file.path(out_dir, "CATA-1.RData")
)
