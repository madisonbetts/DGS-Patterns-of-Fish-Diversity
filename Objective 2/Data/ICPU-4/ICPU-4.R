# -----------------------------
# ICPU-4: Channel Catfish
# Ictalurus punctatus
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)

# -----------------------------
# 0) paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ICPU-4"
data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) site coordinates
# order follows Table 2 FST order
# -----------------------------
ICPU_4_coords <- data.frame(
  site = as.character(1:10),
  site_name = c(
    "Mississippi Pool 2, MN",
    "St. Croix River, MN",
    "Lake Oahe, ND",
    "Lake Sakakawea, ND",
    "Red River - Wahpeton",
    "Red River - Fargo",
    "Red River - Grand Forks",
    "Red River - Drayton",
    "Red River - St. Andrews, MB",
    "James River, ND"
  ),
  lat = c(
    44.9050,
    44.9470,
    46.8080,
    47.4910,
    46.2652,  # Wahpeton
    46.8772,  # Fargo
    47.9253,  # Grand Forks
    48.5714,  # Drayton
    50.0270,  # St. Andrews
    46.9105
  ),
  lon = c(
    -93.0780,
    -92.7830,
    -100.7830,
    -101.4100,
    -96.6056,
    -96.7898,
    -97.0329,
    -97.1779,
    -96.9860,
    -98.7081
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) helper: fill symmetric matrix from lower triangle
# vals entered row-by-row:
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
# Table 2, values below diagonal, row-by-row
# -----------------------------
fst_vals <- c(
  
  # row 2: StC
  0.008,
  
  # row 3: Oahe
  0.009, 0.010,
  
  # row 4: Sak
  0.011, 0.012, 0.006,
  
  # row 5: 1Red
  0.034, 0.031, 0.033, 0.034,
  
  # row 6: 2Red
  0.014, 0.015, 0.011, 0.012, 0.026,
  
  # row 7: 3Red
  0.017, 0.023, 0.018, 0.017, 0.040, 0.017,
  
  # row 8: 4Red
  0.019, 0.023, 0.020, 0.019, 0.037, 0.018, 0.017,
  
  # row 9: StA
  0.026, 0.025, 0.026, 0.025, 0.039, 0.018, 0.021, 0.016,
  
  # row 10: James
  0.020, 0.023, 0.017, 0.021, 0.040, 0.027, 0.031, 0.032, 0.037
)

ICPU_4_fst <- fill_sym_from_lower(
  pops = ICPU_4_coords$site,
  vals = fst_vals,
  diag_val = 0
)

ICPU_4_fst[ICPU_4_fst < 0] <- 0

# -----------------------------
# 4) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(ICPU_4_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ICPU_4_coords$site
colnames(geo_dist_km) <- ICPU_4_coords$site

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(ICPU_4_fst)[row(ICPU_4_fst)[upper.tri(ICPU_4_fst)]],
  site2   = colnames(ICPU_4_fst)[col(ICPU_4_fst)[upper.tri(ICPU_4_fst)]],
  fst     = ICPU_4_fst[upper.tri(ICPU_4_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 6) map
# include US + Canada, zoom to point extent
# plot only in RStudio; do not save
# -----------------------------
usa <- map_data("state")
can <- map_data("world", region = "Canada")

xpad <- 1.5
ypad <- 1.0

xlim_map <- c(
  min(ICPU_4_coords$lon) - xpad,
  max(ICPU_4_coords$lon) + xpad
)

ylim_map <- c(
  min(ICPU_4_coords$lat) - ypad,
  max(ICPU_4_coords$lat) + ypad
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
    data = ICPU_4_coords,
    aes(x = lon, y = lat),
    size = 2.7
  ) +
  geom_text(
    data = ICPU_4_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.12,
    size = 3.3
  ) +
  coord_fixed(
    xlim = xlim_map,
    ylim = ylim_map
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ICPU-4 sampling sites"
  )

print(map_plot)

# -----------------------------
# 7) IBD plot
# plot only in RStudio; do not save
# -----------------------------
ibd_plot <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "ICPU-4 isolation by distance"
  )

print(ibd_plot)

# -----------------------------
# 7b) Red River-only IBD plot
# sites 5:9 = Red River system
# -----------------------------
red_idx <- 5:9

red_coords <- ICPU_4_coords[red_idx, ]
red_fst <- ICPU_4_fst[red_idx, red_idx]

# Euclidean distances (for consistency with main workflow)
red_geo_dist_km <- geosphere::distm(
  as.matrix(red_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(red_geo_dist_km) <- red_coords$site
colnames(red_geo_dist_km) <- red_coords$site

red_ibd_df <- data.frame(
  site1   = rownames(red_fst)[row(red_fst)[upper.tri(red_fst)]],
  site2   = colnames(red_fst)[col(red_fst)[upper.tri(red_fst)]],
  fst     = red_fst[upper.tri(red_fst)],
  dist_km = red_geo_dist_km[upper.tri(red_geo_dist_km)]
)

# plot
red_ibd_plot <- ggplot(red_ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.8, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dotted") +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "ICPU-4 Red River IBD"
  )

print(red_ibd_plot)

# -----------------------------
# 8) quick checks
# -----------------------------
stopifnot(identical(rownames(ICPU_4_fst), ICPU_4_coords$site))
stopifnot(identical(colnames(ICPU_4_fst), ICPU_4_coords$site))
stopifnot(isTRUE(all.equal(ICPU_4_fst, t(ICPU_4_fst))))
stopifnot(length(fst_vals) == nrow(ICPU_4_coords) * (nrow(ICPU_4_coords) - 1) / 2)

# -----------------------------
# 9) save RData
# -----------------------------
ICPU_4_coords <- ICPU_4_coords[, c("site", "lat", "lon")]
save(
  ICPU_4_fst,
  ICPU_4_coords,
  file = file.path(data_dir, "ICPU-4.RData")
)
