# -----------------------------
# ETNI-2: Johnny darter (Etheostoma nigrum)
# 8 sites, FST from table row-by-row
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)

# -----------------------------
# 0) paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETNI-2"
data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) site coordinates
# best-available coordinates based on study figures / local crossings
# update these if you tighten them further
# -----------------------------
ETNI_2_coords <- data.frame(
  site_id = 1:8,
  site_name = c(
    "Beaver Creek",
    "Pipestem Creek",
    "Forest River",
    "Turtle River",
    "Lake Ida",
    "Shell River",
    "Fishhook River",
    "Coffeepot Landing"
  ),
  lat = c(
    46.8100,  # Beaver
    46.9300,  # Pipestem
    48.9120,  # Forest
    48.3920,  # Turtle
    46.7390,  # Lake Ida
    46.9140,  # Shell
    46.9000,  # Fishhook
    47.1100   # Coffeepot
  ),
  lon = c(
    -98.0600,  # Beaver
    -98.4200,  # Pipestem
    -97.5250,  # Forest
    -97.9500,  # Turtle
    -95.9350,  # Lake Ida
    -95.8660,  # Shell
    -95.1500,  # Fishhook
    -95.2000   # Coffeepot
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) helper: fill symmetric matrix
# vals must be lower triangle in row-wise order:
# row2(1), row3(2), row4(3), ...
# -----------------------------
fill_sym_from_lower <- function(n, vals, diag_val = 0) {
  stopifnot(length(vals) == n * (n - 1) / 2)
  
  mat <- matrix(0, nrow = n, ncol = n)
  
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
# 3) FST values from Table 1
# entered row by row from lower triangle
# -----------------------------
fst_vals <- c(
  # row 2
  0.02137,
  
  # row 3
  0.41000, 0.46890,
  
  # row 4
  0.60343, 0.63574, 0.21795,
  
  # row 5
  0.48222, 0.51667, 0.21016, 0.40701,
  
  # row 6
  0.56211, 0.60209, 0.28296, 0.32696, 0.23935,
  
  # row 7
  0.48374, 0.52020, 0.18324, 0.33357, 0.01903, 0.14189,
  
  # row 8
  0.50149, 0.55048, 0.27268, 0.36301, 0.16616, 0.06192, 0.10258
)

ETNI_2_fst <- fill_sym_from_lower(
  n = 8,
  vals = fst_vals,
  diag_val = 0
)

ETNI_2_fst[ETNI_2_fst < 0] <- 0

rownames(ETNI_2_fst) <- ETNI_2_coords$site_id
colnames(ETNI_2_fst) <- ETNI_2_coords$site_id

# -----------------------------
# 4) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(ETNI_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ETNI_2_coords$site_id
colnames(geo_dist_km) <- ETNI_2_coords$site_id

# -----------------------------
# 5) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(ETNI_2_fst)[row(ETNI_2_fst)[upper.tri(ETNI_2_fst)]],
  site2   = colnames(ETNI_2_fst)[col(ETNI_2_fst)[upper.tri(ETNI_2_fst)]],
  fst     = ETNI_2_fst[upper.tri(ETNI_2_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 6) map
# zoom to extent of points with a small buffer
# US + Canada included, but cropped to data extent
# -----------------------------
usa <- map_data("state")
can <- map_data("world", region = "Canada")

xpad <- 1.5
ypad <- 1.0

xlim_map <- c(min(ETNI_2_coords$lon) - xpad, max(ETNI_2_coords$lon) + xpad)
ylim_map <- c(min(ETNI_2_coords$lat) - ypad, max(ETNI_2_coords$lat) + ypad)

map_plot <- ggplot() +
  geom_polygon(
    data = usa,
    aes(x = long, y = lat, group = group),
    fill = "grey90", color = "black", linewidth = 0.3
  ) +
  geom_polygon(
    data = can,
    aes(x = long, y = lat, group = group),
    fill = "grey85", color = "black", linewidth = 0.3
  ) +
  geom_point(
    data = ETNI_2_coords,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text(
    data = ETNI_2_coords,
    aes(x = lon, y = lat, label = site_id),
    nudge_y = 0.10,
    size = 3.5
  ) +
  coord_fixed(
    xlim = xlim_map,
    ylim = ylim_map
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ETNI-2 Sampling Sites"
  )

print(map_plot)

# -----------------------------
# 7) IBD plot
# just print in RStudio, do not save
# -----------------------------
ibd_plot <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.8, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "ETNI-2 Isolation by Distance"
  )

print(ibd_plot)

# -----------------------------
# 8) quick checks
# -----------------------------
stopifnot(identical(rownames(ETNI_2_fst), as.character(ETNI_2_coords$site_id)))
stopifnot(identical(colnames(ETNI_2_fst), as.character(ETNI_2_coords$site_id)))
stopifnot(isTRUE(all.equal(ETNI_2_fst, t(ETNI_2_fst))))
stopifnot(length(fst_vals) == nrow(ETNI_2_coords) * (nrow(ETNI_2_coords) - 1) / 2)

# -----------------------------
# 9) save RData only
# -----------------------------
save(
  ETNI_2_fst,
  ETNI_2_coords,
  file = file.path(data_dir, "ETNI-2.RData")
)