# -----------------------------
# RHEG-1 Lahontan redside
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
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/RHEG-1"

# -----------------------------
# 0) site coordinates
# Table 1 in Peacock et al. 2016
# site numbering matches FST matrix order below
# -----------------------------
RHEG_1_coords <- data.frame(
  site = as.character(1:12),
  site_name = c(
    "TahoeB", "RenoA", "McRC", "McRB", "McRA", "McRD",
    "WadsA", "WadsB", "WadsC", "NixonA", "NixonB", "NixonC"
  ),
  lat = c(
    39.162, 39.514, 39.526, 39.545, 39.547, 39.565,
    39.585, 39.591, 39.613, 39.727, 39.818, 39.854
  ),
  lon = c(
    -120.159, -119.736, -119.610, -119.587, -119.573, -119.487,
    -119.444, -119.368, -119.306, -119.319, -119.350, -119.394
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

xlim_map <- c(min(RHEG_1_coords$lon) - xpad, max(RHEG_1_coords$lon) + xpad)
ylim_map <- c(min(RHEG_1_coords$lat) - ypad, max(RHEG_1_coords$lat) + ypad)

ggplot() +
  geom_polygon(data = usa, aes(x = long, y = lat, group = group),
               fill = "grey92", color = "black", linewidth = 0.25) +
  geom_polygon(data = can, aes(x = long, y = lat, group = group),
               fill = "grey86", color = "black", linewidth = 0.25) +
  geom_point(data = RHEG_1_coords, aes(x = lon, y = lat), size = 2.7) +
  ggrepel::geom_text_repel(
    data = RHEG_1_coords,
    aes(x = lon, y = lat, label = site),
    size = 3.4,
    seed = 1,
    min.segment.length = 0
  ) +
  coord_fixed(xlim = xlim_map, ylim = ylim_map) +
  theme_classic() +
  labs(x = "Longitude", y = "Latitude", title = "RHEG-1 sampling sites")

# -----------------------------
# 2) geographic distance matrix
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(RHEG_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- RHEG_1_coords$site
colnames(geo_dist_km) <- RHEG_1_coords$site

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
# S5 Table (C): Lahontan redside, 2008
# lower triangle, row-wise order
# -----------------------------
fst_vals <- c(
  # 2
  -0.0053,

  # 3
  0.0034, -0.004,

  # 4
  0.0354, 0.0224, 0.0174,

  # 5
  0.0606, 0.0380, 0.0316, 0.0036,

  # 6
  -0.0071, -0.0025, -0.0039, 0.0204, 0.0371,

  # 7
  0.0532, 0.0384, 0.0299, 0.0044, 0.0011, 0.0358,

  # 8
  0.0660, 0.0438, 0.0338, 0.0039, -0.0010, 0.0420, -0.0024,

  # 9
  0.0541, 0.0456, 0.0391, 0.0056, 0.0058, 0.0436, -0.0007, 0.0011,

  # 10
  0.0626, 0.0350, 0.0272, 0.0096, 0.0037, 0.0352, -0.0031, 0.0034, 0.0072,

  # 11
  0.0460, 0.0309, 0.0303, 0.0122, 0.0116, 0.0365, 0.0086, 0.0100, 0.0118, 0.0044,

  # 12
  0.0436, 0.0260, 0.0204, 0.0058, 0.0022, 0.0244, 0.0060, 0.0070, 0.0141, -0.0063, 0.0008
)

RHEG_1_fst <- fill_sym_from_lower(
  pops = RHEG_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

RHEG_1_fst[RHEG_1_fst < 0] <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(RHEG_1_fst)[row(RHEG_1_fst)[upper.tri(RHEG_1_fst)]],
  site2   = colnames(RHEG_1_fst)[col(RHEG_1_fst)[upper.tri(RHEG_1_fst)]],
  fst     = RHEG_1_fst[upper.tri(RHEG_1_fst)],
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
    title = "RHEG-1 isolation by distance"
  )

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(RHEG_1_fst), RHEG_1_coords$site))
stopifnot(identical(colnames(RHEG_1_fst), RHEG_1_coords$site))
stopifnot(isTRUE(all.equal(RHEG_1_fst, t(RHEG_1_fst))))
stopifnot(length(fst_vals) == nrow(RHEG_1_coords) * (nrow(RHEG_1_coords) - 1) / 2)

# -----------------------------
# 8) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

RHEG_1_coords <- RHEG_1_coords[, c("site", "lat", "lon")]

save(
  RHEG_1_fst,
  RHEG_1_coords,
  file = file.path(out_dir, "RHEG-1.RData")
)
