# -----------------------------
# RHOS-5 Lahontan speckled dace
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
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/RHOS-5"

# -----------------------------
# 0) site coordinates
# Table 1 in Peacock et al. 2016
# site numbering matches FST matrix order below
# -----------------------------
RHOS_5_coords <- data.frame(
  site = as.character(1:17),
  site_name = c(
    "TahoeB", "TahoeC", "TahoeA", "TrkC", "FaradC",
    "VerdiC", "VerdiB",
    "RenoC", "RenoB", "RenoA",
    "McRC", "McRB", "McRA",
    "WadsC", "NixonA", "NixonB", "NixonC"
  ),
  lat = c(
    39.162, 39.208, 39.274, 39.334, 39.450,
    39.509, 39.507,
    39.524, 39.530, 39.514,
    39.526, 39.545, 39.547,
    39.613, 39.727, 39.818, 39.854
  ),
  lon = c(
    -120.159, -120.198, -120.206, -120.164, -120.005,
    -119.996, -119.902,
    -119.819, -119.795, -119.736,
    -119.610, -119.587, -119.573,
    -119.306, -119.319, -119.350, -119.394
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map
# -----------------------------
usa <- map_data("state")
can <- map_data("world", region = "Canada")

xpad <- 1.0
ypad <- 0.8

xlim_map <- c(min(RHOS_5_coords$lon) - xpad, max(RHOS_5_coords$lon) + xpad)
ylim_map <- c(min(RHOS_5_coords$lat) - ypad, max(RHOS_5_coords$lat) + ypad)

ggplot() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group),
               fill = "grey92", color = "black", linewidth = 0.25) +
  geom_polygon(data = can, aes(x = long, y = lat, group = group),
               fill = "grey86", color = "black", linewidth = 0.25) +
  geom_point(data = RHOS_5_coords, aes(x = lon, y = lat), size = 2.7) +
  ggrepel::geom_text_repel(
    data = RHOS_5_coords,
    aes(x = lon, y = lat, label = site),
    size = 3.4,
    seed = 1,
    min.segment.length = 0
  ) +
  coord_fixed(xlim = xlim_map, ylim = ylim_map) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", title = "RHOS-5 sampling sites")

# -----------------------------
# 2) geographic distance matrix
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(RHOS_5_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- RHOS_5_coords$site
colnames(geo_dist_km) <- RHOS_5_coords$site

# -----------------------------
# 3) helper
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
# S5 Table (D): Lahontan speckled dace, 2008
# lower triangle, row-wise order
# note: parsed header omitted McRC, but row structure clearly includes it
# -----------------------------
fst_vals <- c(
  # 2
  0.007,

  # 3
  -0.007, 0.005,

  # 4
  0.001, 0.003, 0.004,

  # 5
  0.081, 0.026, 0.074, 0.005,

  # 6
  0.043, 0.033, 0.060, -0.003, -0.016,

  # 7
  0.034, 0.019, 0.042, -0.021, -0.071, -0.038,

  # 8
  0.036, 0.039, 0.044, -0.003, 0.025, 0.025, 0.008,

  # 9
  0.000, -0.040, -0.002, -0.022, -0.020, 0.015, -0.024, -0.008,

  # 10
  0.034, 0.022, 0.043, -0.002, 0.020, -0.012, -0.012, 0.008, -0.009,

  # 11
  0.028, 0.024, 0.044, -0.002, 0.028, -0.001, -0.016, -0.002, -0.003, -0.005,

  # 12
  0.069, 0.034, 0.081, 0.021, -0.015, -0.003, -0.037, 0.015, -0.009, 0.002, 0.017,

  # 13
  0.037, 0.035, 0.036, 0.006, 0.026, -0.003, -0.034, 0.002, -0.022, -0.005, 0.000, 0.004,

  # 14
  0.030, 0.019, 0.049, -0.003, 0.013, -0.013, -0.009, -0.010, -0.014, 0.001, -0.004, 0.000, 0.004,

  # 15
  0.023, 0.029, 0.001, -0.015, -0.018, 0.048, -0.052, -0.007, -0.018, -0.019, -0.015, 0.008, -0.042, 0.000,

  # 16
  0.053, 0.033, 0.079, 0.013, 0.026, -0.007, 0.002, 0.018, 0.010, 0.000, 0.005, 0.014, 0.015, 0.005, 0.000,

  # 17
  0.025, 0.018, 0.016, -0.005, 0.002, -0.006, -0.030, -0.006, -0.013, -0.005, -0.004, 0.007, -0.006, -0.043, 0.005, 0.005
)

RHOS_5_fst <- fill_sym_from_lower(
  pops = RHOS_5_coords$site,
  vals = fst_vals,
  diag_val = 0
)

RHOS_5_fst[RHOS_5_fst < 0] <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(RHOS_5_fst)[row(RHOS_5_fst)[upper.tri(RHOS_5_fst)]],
  site2   = colnames(RHOS_5_fst)[col(RHOS_5_fst)[upper.tri(RHOS_5_fst)]],
  fst     = RHOS_5_fst[upper.tri(RHOS_5_fst)],
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
    title = "RHOS-5 isolation by distance"
  )

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(RHOS_5_fst), RHOS_5_coords$site))
stopifnot(identical(colnames(RHOS_5_fst), RHOS_5_coords$site))
stopifnot(isTRUE(all.equal(RHOS_5_fst, t(RHOS_5_fst))))
stopifnot(length(fst_vals) == nrow(RHOS_5_coords) * (nrow(RHOS_5_coords) - 1) / 2)

# -----------------------------
# 8) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

RHOS_5_coords <- RHOS_5_coords[, c("site", "lat", "lon")]

save(
  RHOS_5_fst,
  RHOS_5_coords,
  file = file.path(out_dir, "RHOS-5.RData")
)
