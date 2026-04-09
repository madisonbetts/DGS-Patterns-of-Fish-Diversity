# -----------------------------
# PEFL-4: yellow perch
# Perca flavescens
# -----------------------------

library(dplyr)
library(ggplot2)
library(geosphere)
library(maps)

# -----------------------------
# 0) paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/PEFL-4"
data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) site coordinates
# from Table 1 / Figure 1
# final coords df only has: site, lat, lon
# -----------------------------
PEFL_4_coords <- data.frame(
  site = as.character(1:17),
  lat = c(
    48.0,   # 1 Devils Lake, ND (approx from Fig. 1A)
    46.95,  # 2 Bad River, WI
    46.01,  # 3 Lac du Flambeau, WI
    45.70,  # 4 Little Tail Point, WI
    44.12,  # 5 Lake Winnebago, WI
    43.00,  # 6 Lake Michigan 1998
    43.00,  # 7 Lake Michigan 2002
    43.25,  # 8 Lake Ontario
    39.00,  # 9 Severn River, MD
    39.12,  # 10 Bush River, MD
    38.63,  # 11 Choptank River, MD
    38.38,  # 12 Nanticoke River, MD
    36.17,  # 13 Perquimans River, NC
    36.20,  # 14 Little River, NC
    36.12,  # 15 Pasquotank River, NC
    36.20,  # 16 North River, NC
    35.90   # 17 Scuppernong River, NC
  ),
  lon = c(
    -99.0,   # Devils Lake
    -90.85,  # Bad River
    -89.89,  # Lac du Flambeau
    -87.90,  # Little Tail Point / Green Bay
    -88.50,  # Lake Winnebago
    -87.00,  # Lake Michigan
    -87.00,  # Lake Michigan
    -77.90,  # Lake Ontario west basin sample shown in Fig. 1
    -76.60,  # Severn River
    -76.20,  # Bush River
    -76.00,  # Choptank River
    -75.90,  # Nanticoke River
    -76.35,  # Perquimans River
    -76.15,  # Little River
    -76.05,  # Pasquotank River
    -75.95,  # North River
    -76.15   # Scuppernong River
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) helper: fill symmetric matrix from lower triangle
# row-wise order:
# row2 has 1 value, row3 has 2 values, etc.
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
# 3) pairwise FST matrix
# Table 3, hST values below diagonal
# row-wise lower triangle
# -----------------------------
fst_vals <- c(
  # row 2
  0.172,
  
  # row 3
  0.111, 0.152,
  
  # row 4
  0.242, 0.137, 0.290,
  
  # row 5
  0.198, 0.131, 0.250, 0.072,
  
  # row 6
  0.342, 0.165, 0.386, 0.096, 0.276,
  
  # row 7
  0.337, 0.144, 0.384, 0.065, 0.249, 0.000,
  
  # row 8
  0.460, 0.113, 0.442, 0.338, 0.399, 0.298, 0.274,
  
  # row 9
  0.841, 0.633, 0.835, 0.810, 0.835, 0.786, 0.775, 0.533,
  
  # row 10
  0.817, 0.613, 0.809, 0.778, 0.806, 0.750, 0.736, 0.501, 0.051,
  
  # row 11
  0.852, 0.653, 0.848, 0.821, 0.845, 0.800, 0.789, 0.564, 0.043, 0.112,
  
  # row 12
  0.819, 0.615, 0.817, 0.777, 0.804, 0.752, 0.737, 0.503, 0.226, 0.220, 0.166,
  
  # row 13
  0.870, 0.685, 0.867, 0.849, 0.867, 0.831, 0.823, 0.615, 0.113, 0.234, 0.159, 0.304,
  
  # row 14
  0.889, 0.668, 0.888, 0.863, 0.885, 0.846, 0.836, 0.588, 0.053, 0.166, 0.102, 0.284, 0.024,
  
  # row 15
  0.867, 0.657, 0.866, 0.842, 0.864, 0.823, 0.811, 0.575, 0.088, 0.193, 0.152, 0.267, 0.032, 0.005,
  
  # row 16
  0.882, 0.693, 0.880, 0.861, 0.880, 0.844, 0.836, 0.624, 0.111, 0.224, 0.189, 0.366, 0.034, 0.019, 0.015,
  
  # row 17
  0.902, 0.718, 0.899, 0.881, 0.899, 0.867, 0.862, 0.655, 0.121, 0.251, 0.173, 0.368, 0.051, 0.000, 0.077, 0.070
)

PEFL_4_fst <- fill_sym_from_lower(
  pops = PEFL_4_coords$site,
  vals = fst_vals,
  diag_val = 0
)

PEFL_4_fst[PEFL_4_fst < 0] <- 0

# -----------------------------
# 4) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(PEFL_4_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- PEFL_4_coords$site
colnames(geo_dist_km) <- PEFL_4_coords$site

# -----------------------------
# 5) all-sites IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(PEFL_4_fst)[row(PEFL_4_fst)[upper.tri(PEFL_4_fst)]],
  site2   = colnames(PEFL_4_fst)[col(PEFL_4_fst)[upper.tri(PEFL_4_fst)]],
  fst     = PEFL_4_fst[upper.tri(PEFL_4_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 6) map
# include US + Canada, zoom to point extent
# plot only in RStudio; do not save
# -----------------------------
usa <- map_data("state")
can <- map_data("world", region = "Canada")

xpad <- 2.0
ypad <- 1.5

xlim_map <- c(
  min(PEFL_4_coords$lon) - xpad,
  max(PEFL_4_coords$lon) + xpad
)

ylim_map <- c(
  min(PEFL_4_coords$lat) - ypad,
  max(PEFL_4_coords$lat) + ypad
)

map_plot <- ggplot() +
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
    data = PEFL_4_coords,
    aes(x = lon, y = lat),
    size = 2.7
  ) +
  geom_text(
    data = PEFL_4_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.20,
    size = 3.2
  ) +
  coord_fixed(
    xlim = xlim_map,
    ylim = ylim_map
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "PEFL-4 sampling sites"
  )

print(map_plot)

# -----------------------------
# 7) all-sites IBD plot
# -----------------------------
ibd_plot <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(h[ST]),
    title = "PEFL-4 isolation by distance"
  )

print(ibd_plot)

# -----------------------------
# 8) Chesapeake-only IBD
# sites 9:12
# -----------------------------
ches_idx <- 9:12

ches_coords <- PEFL_4_coords[ches_idx, ]
ches_fst <- PEFL_4_fst[ches_idx, ches_idx]

ches_geo_dist_km <- geosphere::distm(
  as.matrix(ches_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(ches_geo_dist_km) <- ches_coords$site
colnames(ches_geo_dist_km) <- ches_coords$site

ches_ibd_df <- data.frame(
  site1   = rownames(ches_fst)[row(ches_fst)[upper.tri(ches_fst)]],
  site2   = colnames(ches_fst)[col(ches_fst)[upper.tri(ches_fst)]],
  fst     = ches_fst[upper.tri(ches_fst)],
  dist_km = ches_geo_dist_km[upper.tri(ches_geo_dist_km)]
)

ches_ibd_plot <- ggplot(ches_ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dotted") +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(h[ST]),
    title = "PEFL-4 Chesapeake Bay IBD"
  )

print(ches_ibd_plot)

# -----------------------------
# 9) Albemarle-only IBD
# sites 13:17
# -----------------------------
alb_idx <- 13:17

alb_coords <- PEFL_4_coords[alb_idx, ]
alb_fst <- PEFL_4_fst[alb_idx, alb_idx]

alb_geo_dist_km <- geosphere::distm(
  as.matrix(alb_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(alb_geo_dist_km) <- alb_coords$site
colnames(alb_geo_dist_km) <- alb_coords$site

alb_ibd_df <- data.frame(
  site1   = rownames(alb_fst)[row(alb_fst)[upper.tri(alb_fst)]],
  site2   = colnames(alb_fst)[col(alb_fst)[upper.tri(alb_fst)]],
  fst     = alb_fst[upper.tri(alb_fst)],
  dist_km = alb_geo_dist_km[upper.tri(alb_geo_dist_km)]
)

alb_ibd_plot <- ggplot(alb_ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dotted") +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(h[ST]),
    title = "PEFL-4 Albemarle Sound IBD"
  )

print(alb_ibd_plot)

# -----------------------------
# 10) quick checks
# -----------------------------
stopifnot(nrow(PEFL_4_coords) == 17)
stopifnot(identical(rownames(PEFL_4_fst), PEFL_4_coords$site))
stopifnot(identical(colnames(PEFL_4_fst), PEFL_4_coords$site))
stopifnot(isTRUE(all.equal(PEFL_4_fst, t(PEFL_4_fst))))
stopifnot(length(fst_vals) == nrow(PEFL_4_coords) * (nrow(PEFL_4_coords) - 1) / 2)

# -----------------------------
# 11) save RData
# -----------------------------
save(
  PEFL_4_fst,
  PEFL_4_coords,
  file = file.path(data_dir, "PEFL-4.RData")
)

# -----------------------------
# 12) export coords csv
# -----------------------------
write.csv(
  PEFL_4_coords,
  file = file.path(data_dir, "PEFL-4_coords.csv"),
  row.names = FALSE
)