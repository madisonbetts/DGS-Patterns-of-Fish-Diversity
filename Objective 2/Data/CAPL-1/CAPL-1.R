# -----------------------------
# CAPL-1 Mountain sucker
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
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CAPL-1"

# -----------------------------
# 0) site coordinates
# Table 1 in Peacock et al. 2016
# site numbering matches FST matrix order below
# -----------------------------
CAPL_1_coords <- data.frame(
  site = as.character(1:13),
  site_name = c(
    "VerdiC", "VerdiB", "RenoC", "RenoB", "RenoA",
    "McRC", "McRB", "McRA",
    "WadsA", "WadsB", "WadsC",
    "NixonB", "NixonC"
  ),
  lat = c(
    39.509, 39.507, 39.524, 39.530, 39.514,
    39.526, 39.545, 39.547,
    39.585, 39.591, 39.613,
    39.818, 39.854
  ),
  lon = c(
    -119.996, -119.902, -119.819, -119.795, -119.736,
    -119.610, -119.587, -119.573,
    -119.444, -119.368, -119.306,
    -119.350, -119.394
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

xlim_map <- c(min(CAPL_1_coords$lon) - xpad, max(CAPL_1_coords$lon) + xpad)
ylim_map <- c(min(CAPL_1_coords$lat) - ypad, max(CAPL_1_coords$lat) + ypad)

ggplot() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group),
               fill = "grey92", color = "black", linewidth = 0.25) +
  geom_polygon(data = can, aes(x = long, y = lat, group = group),
               fill = "grey86", color = "black", linewidth = 0.25) +
  geom_point(data = CAPL_1_coords, aes(x = lon, y = lat), size = 2.7) +
  ggrepel::geom_text_repel(
    data = CAPL_1_coords,
    aes(x = lon, y = lat, label = site),
    size = 3.4,
    seed = 1,
    min.segment.length = 0
  ) +
  coord_fixed(xlim = xlim_map, ylim = ylim_map) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", title = "CAPL-1 sampling sites")

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(CAPL_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- CAPL_1_coords$site
colnames(geo_dist_km) <- CAPL_1_coords$site

# -----------------------------
# 3) helper function
# -----------------------------
fill_sym_from_lower <- function(pops, vals, diag_val = 0) {
  n <- length(pops)
  stopifnot(length(vals) == n * (n - 1) / 2)

  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(pops, pops))

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
# S5 Table (A): Mountain sucker, 2008
# lower triangle, row-wise order
# -----------------------------
fst_vals <- c(
  # 2
  0.012,

  # 3
  0.016, -0.002,

  # 4
  0.013, 0.004, 0.003,

  # 5
  0.028, 0.023, 0.022, 0.022,

  # 6
  0.040, 0.028, 0.027, 0.027, -0.003,

  # 7
  0.038, 0.032, 0.027, 0.031, 0.004, 0.001,

  # 8
  0.045, 0.039, 0.034, 0.040, 0.005, -0.004, -0.001,

  # 9
  0.036, 0.031, 0.029, 0.030, -0.003, -0.001, 0.002, -0.003,

  # 10
  0.000, 0.020, 0.020, 0.022, -0.005, -0.007, 0.001, -0.003, -0.001,

  # 11
  0.038, 0.031, 0.030, 0.034, 0.001, -0.003, 0.003, -0.003, 0.001, -0.006,

  # 12
  0.029, 0.031, 0.030, 0.029, -0.001, 0.005, 0.014, 0.011, 0.008, -0.001, 0.003,

  # 13
  0.048, 0.040, 0.035, 0.042, -0.001, -0.004, 0.003, -0.006, 0.002, -0.020, -0.006, -0.035
)

CAPL_1_fst <- fill_sym_from_lower(
  pops = CAPL_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

CAPL_1_fst[CAPL_1_fst < 0] <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(CAPL_1_fst)[row(CAPL_1_fst)[upper.tri(CAPL_1_fst)]],
  site2   = colnames(CAPL_1_fst)[col(CAPL_1_fst)[upper.tri(CAPL_1_fst)]],
  fst     = CAPL_1_fst[upper.tri(CAPL_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
)

# -----------------------------
# 6) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "CAPL-1 isolation by distance"
  )

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(CAPL_1_fst), CAPL_1_coords$site))
stopifnot(identical(colnames(CAPL_1_fst), CAPL_1_coords$site))
stopifnot(isTRUE(all.equal(CAPL_1_fst, t(CAPL_1_fst))))
stopifnot(length(fst_vals) == nrow(CAPL_1_coords) * (nrow(CAPL_1_coords) - 1) / 2)

# -----------------------------
# 8) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

CAPL_1_coords <- CAPL_1_coords[, c("site", "lat", "lon")]

save(
  CAPL_1_fst,
  CAPL_1_coords,
  file = file.path(out_dir, "CAPL-1.RData")
)
