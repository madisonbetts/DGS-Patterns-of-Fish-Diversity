# -----------------------------
# ETPA-1: paleback darter
# Etheostoma pallididorsum
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)

# -----------------------------
# 0) paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETPA-1"
data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) site coordinates
# hard-coded from final raw coordinates
# order follows FST table order
# 1 Caddo River
# 2 Polk Creek
# 3 Lick Creek
# 4 Collier Creek
# 5 Big Fork Creek
# 6 Kate's Creek
# 7 Big Hill Creek
# 8 Mazarn Creek
# -----------------------------
ETPA_1_coords <- data.frame(
  site = as.character(1:8),
  lat = c(
    34.4864,
    34.4444,
    34.5037,
    34.4518,
    34.5012,
    34.5284,
    34.5260,
    34.4963
  ),
  lon = c(
    -93.8257,
    -93.7486,
    -93.7244,
    -93.6297,
    -93.9248,
    -93.8940,
    -93.8543,
    -93.4536
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
# Table 4, microsatellite FST below diagonal
# entered row-by-row from lower triangle
# -----------------------------
fst_vals <- c(
  
  # row 2: Polk Creek
  0.02,
  
  # row 3: Lick Creek
  0.02, 0.02,
  
  # row 4: Collier Creek
  0.02, 0.01, 0.02,
  
  # row 5: Big Fork Creek
  0.18, 0.18, 0.20, 0.18,
  
  # row 6: Kate's Creek
  0.24, 0.21, 0.24, 0.23, 0.33,
  
  # row 7: Big Hill Creek
  0.16, 0.13, 0.16, 0.14, 0.27, 0.13,
  
  # row 8: Mazarn Creek
  0.19, 0.22, 0.24, 0.23, 0.33, 0.36, 0.30
)

ETPA_1_fst <- fill_sym_from_lower(
  pops = ETPA_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

ETPA_1_fst[ETPA_1_fst < 0] <- 0

# -----------------------------
# 4) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(ETPA_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ETPA_1_coords$site
colnames(geo_dist_km) <- ETPA_1_coords$site

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(ETPA_1_fst)[row(ETPA_1_fst)[upper.tri(ETPA_1_fst)]],
  site2   = colnames(ETPA_1_fst)[col(ETPA_1_fst)[upper.tri(ETPA_1_fst)]],
  fst     = ETPA_1_fst[upper.tri(ETPA_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 6) map
# include US + Canada, but zoom to point extent
# plot only in RStudio; do not save
# -----------------------------
usa <- map_data("state")
can <- map_data("world", region = "Canada")

xpad <- 0.08
ypad <- 0.05

xlim_map <- c(
  min(ETPA_1_coords$lon) - xpad,
  max(ETPA_1_coords$lon) + xpad
)

ylim_map <- c(
  min(ETPA_1_coords$lat) - ypad,
  max(ETPA_1_coords$lat) + ypad
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
    data = ETPA_1_coords,
    aes(x = lon, y = lat),
    size = 2.7
  ) +
  geom_text(
    data = ETPA_1_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.008,
    size = 3.4
  ) +
  coord_fixed(
    xlim = xlim_map,
    ylim = ylim_map
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ETPA-1 sampling sites"
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
    title = "ETPA-1 isolation by distance"
  )

print(ibd_plot)

# -----------------------------
# 8) quick checks
# -----------------------------
stopifnot(identical(rownames(ETPA_1_fst), ETPA_1_coords$site))
stopifnot(identical(colnames(ETPA_1_fst), ETPA_1_coords$site))
stopifnot(isTRUE(all.equal(ETPA_1_fst, t(ETPA_1_fst))))
stopifnot(length(fst_vals) == nrow(ETPA_1_coords) * (nrow(ETPA_1_coords) - 1) / 2)

# -----------------------------
# 9) save RData
# -----------------------------
save(
  ETPA_1_fst,
  ETPA_1_coords,
  file = file.path(data_dir, "ETPA-1.RData")
)

# -----------------------------
# 10) export coords csv
# -----------------------------
write.csv(
  ETPA_1_coords,
  file = file.path(data_dir, "ETPA-1_coords.csv"),
  row.names = FALSE
)