# -----------------------------
# MIDO-4: smallmouth bass
# Micropterus dolomieu
# -----------------------------

library(dplyr)
library(ggplot2)
library(geosphere)
library(maps)

# -----------------------------
# 0) paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/MIDO-4"
data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) posted Table I coordinates
# collapse locality coordinates to the 9 Table II regional groups
# final coords df only has: site, lat, lon
# site order matches Table II exactly
# -----------------------------
site_coords_raw <- data.frame(
  locality = c(
    "Cannon River, MN",
    "Apple River, WI",
    "St. Croix River, MN",
    "St Louis Bay, MN",
    "Whitefish Bay, MI",
    "Saginaw Bay, MI",
    "Georgian Bay, ON",
    "Anchor Bay, MI",
    "Sandusky Bay, OH",
    "Chagrin River, OH",
    "Grand River, OH",
    "Cattaraugus Creek, NY",
    "Long Point Bay, ON",
    "Fairbanks Point, NY",
    "Cranberry Lake, NY",
    "Little Moose Lake, NY",
    "Hudson River, NY",
    "Paint Creek, OH"
  ),
  region = c(
    "Upper Mississippi River",
    "Upper Mississippi River",
    "Upper Mississippi River",
    "Lake Superior",
    "Lake Superior",
    "Lake Huron",
    "Lake Huron",
    "Lake St Clair",
    "Lake Erie",
    "Lake Erie",
    "Lake Erie",
    "Lake Erie",
    "Lake Erie",
    "Lake Ontario",
    "St Lawrence River",
    "St Lawrence River",
    "Hudson River",
    "Ohio River"
  ),
  lat = c(
    44.52, 45.38, 45.13,
    46.72, 46.52,
    43.55, 44.59,
    42.63,
    41.45, 41.41, 41.75, 42.57, 42.66,
    43.30,
    44.21, 43.69,
    43.24,
    39.32
  ),
  lon = c(
    92.88, 92.44, 92.44,
    92.16, 85.01,
    83.90, 80.94,
    82.77,
    82.96, 81.41, 80.99, 79.11, 80.26,
    77.18,
    74.84, 74.92,
    73.83,
    83.08
  ),
  stringsAsFactors = FALSE
)

# convert to western hemisphere longitudes
site_coords_raw$lon <- -site_coords_raw$lon

region_order <- c(
  "Upper Mississippi River",
  "Lake Superior",
  "Lake Huron",
  "Lake Erie",
  "Lake St Clair",
  "Lake Ontario",
  "St Lawrence River",
  "Hudson River",
  "Ohio River"
)

MIDO_4_coords <- site_coords_raw %>%
  group_by(region) %>%
  summarise(
    lat = mean(lat),
    lon = mean(lon),
    .groups = "drop"
  ) %>%
  mutate(region = factor(region, levels = region_order)) %>%
  arrange(region) %>%
  mutate(site = as.character(1:n())) %>%
  select(site, lat, lon)

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
# 3) Table II pairwise theta_ST
# below diagonal = mtDNA
# above diagonal = microsatellites
#
# Objective 2 workflow uses one matrix.
# Here I use the microsatellite theta_ST values from ABOVE diagonal,
# entered row-by-row as lower triangle in Table II order:
# 1 Upper Mississippi River
# 2 Lake Superior
# 3 Lake Huron
# 4 Lake Erie
# 5 Lake St Clair
# 6 Lake Ontario
# 7 St Lawrence River
# 8 Hudson River
# 9 Ohio River
# -----------------------------
fst_vals <- c(
  # row 2
  0.11,
  
  # row 3
  0.22, 0.09,
  
  # row 4
  0.16, 0.23, 0.25,
  
  # row 5
  0.17, 0.10, 0.06, 0.21,
  
  # row 6
  0.17, 0.22, 0.27, 0.18, 0.15,
  
  # row 7
  0.24, 0.19, 0.19, 0.26, 0.13, 0.25,
  
  # row 8
  0.25, 0.26, 0.33, 0.34, 0.25, 0.38, 0.31,
  
  # row 9
  0.12, 0.15, 0.21, 0.19, 0.15, 0.19, 0.18, 0.20
)

MIDO_4_fst <- fill_sym_from_lower(
  pops = MIDO_4_coords$site,
  vals = fst_vals,
  diag_val = 0
)

MIDO_4_fst[MIDO_4_fst < 0] <- 0

# -----------------------------
# 4) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(MIDO_4_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- MIDO_4_coords$site
colnames(geo_dist_km) <- MIDO_4_coords$site

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(MIDO_4_fst)[row(MIDO_4_fst)[upper.tri(MIDO_4_fst)]],
  site2   = colnames(MIDO_4_fst)[col(MIDO_4_fst)[upper.tri(MIDO_4_fst)]],
  fst     = MIDO_4_fst[upper.tri(MIDO_4_fst)],
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
  min(MIDO_4_coords$lon) - xpad,
  max(MIDO_4_coords$lon) + xpad
)

ylim_map <- c(
  min(MIDO_4_coords$lat) - ypad,
  max(MIDO_4_coords$lat) + ypad
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
    data = MIDO_4_coords,
    aes(x = lon, y = lat),
    size = 2.7
  ) +
  geom_text(
    data = MIDO_4_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.18,
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
    title = "MIDO-4 sampling sites"
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
    y = expression(theta[ST]),
    title = "MIDO-4 isolation by distance"
  )

print(ibd_plot)

# -----------------------------
# 8) quick checks
# -----------------------------
stopifnot(nrow(MIDO_4_coords) == 9)
stopifnot(identical(rownames(MIDO_4_fst), MIDO_4_coords$site))
stopifnot(identical(colnames(MIDO_4_fst), MIDO_4_coords$site))
stopifnot(isTRUE(all.equal(MIDO_4_fst, t(MIDO_4_fst))))
stopifnot(length(fst_vals) == nrow(MIDO_4_coords) * (nrow(MIDO_4_coords) - 1) / 2)

# -----------------------------
# 9) save RData
# -----------------------------
save(
  MIDO_4_fst,
  MIDO_4_coords,
  file = file.path(data_dir, "MIDO-4.RData")
)

# -----------------------------
# 10) export coords csv
# -----------------------------
write.csv(
  MIDO_4_coords,
  file = file.path(data_dir, "MIDO-4_coords.csv"),
  row.names = FALSE
)