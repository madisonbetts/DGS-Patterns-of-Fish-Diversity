# -----------------------------
# POSP-3: paddlefish
# Polyodon spathula
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)

# -----------------------------
# 0) paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/POSP-3"
data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) grouped population coordinates
# paper analyzed 4 grouped populations:
# 1 Ohio River
# 2 Red River
# 3 Yellowstone/Missouri River
# 4 Alabama River
#
# coordinates are best-available grouped reach centroids
# based on the site descriptions and Fig. 1
# -----------------------------
POSP_3_coords <- data.frame(
  site = as.character(1:4),
  lat = c(
    38.27,  # Ohio River, Jefferson Co., Kentucky
    33.28,  # Red River pooled OK + LA sites
    47.80,  # Yellowstone/Missouri River, McKenzie Co., ND
    32.63   # Alabama River pooled AL counties
  ),
  lon = c(
    -85.74,  # Louisville / Jefferson County reach
    -93.35,  # pooled Red River reach centroid
    -103.35, # lower Yellowstone / upper Missouri reach
    -86.85   # central Alabama River pooled reach
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
# Table 3, FST values above diagonal
# paper order:
# Ohio, Red, Yellowstone/Missouri, Alabama
# entered row-wise as lower triangle
# -----------------------------
fst_vals <- c(
  # row 2: Red
  0.0161,
  
  # row 3: Yellowstone/Missouri
  0.0170, 0.0094,
  
  # row 4: Alabama
  0.1381, 0.1194, 0.0952
)

POSP_3_fst <- fill_sym_from_lower(
  pops = POSP_3_coords$site,
  vals = fst_vals,
  diag_val = 0
)

POSP_3_fst[POSP_3_fst < 0] <- 0

# -----------------------------
# 4) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(POSP_3_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- POSP_3_coords$site
colnames(geo_dist_km) <- POSP_3_coords$site

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(POSP_3_fst)[row(POSP_3_fst)[upper.tri(POSP_3_fst)]],
  site2   = colnames(POSP_3_fst)[col(POSP_3_fst)[upper.tri(POSP_3_fst)]],
  fst     = POSP_3_fst[upper.tri(POSP_3_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 6) map
# include US + Canada, zoom to point extent
# plot only in RStudio; do not save
# -----------------------------
usa <- map_data("state")
can <- map_data("world", region = "Canada")

xpad <- 3
ypad <- 2

xlim_map <- c(
  min(POSP_3_coords$lon) - xpad,
  max(POSP_3_coords$lon) + xpad
)

ylim_map <- c(
  min(POSP_3_coords$lat) - ypad,
  max(POSP_3_coords$lat) + ypad
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
    data = POSP_3_coords,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text(
    data = POSP_3_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.35,
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
    title = "POSP-3 sampling sites"
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
    title = "POSP-3 isolation by distance"
  )

print(ibd_plot)

# -----------------------------
# 8) quick checks
# -----------------------------
stopifnot(nrow(POSP_3_coords) == 4)
stopifnot(identical(rownames(POSP_3_fst), POSP_3_coords$site))
stopifnot(identical(colnames(POSP_3_fst), POSP_3_coords$site))
stopifnot(isTRUE(all.equal(POSP_3_fst, t(POSP_3_fst))))
stopifnot(length(fst_vals) == nrow(POSP_3_coords) * (nrow(POSP_3_coords) - 1) / 2)

# -----------------------------
# 9) save RData
# -----------------------------
save(
  POSP_3_fst,
  POSP_3_coords,
  file = file.path(data_dir, "POSP-3.RData")
)

# -----------------------------
# 10) export coords csv
# -----------------------------
write.csv(
  POSP_3_coords,
  file = file.path(data_dir, "POSP-3_coords.csv"),
  row.names = FALSE
)