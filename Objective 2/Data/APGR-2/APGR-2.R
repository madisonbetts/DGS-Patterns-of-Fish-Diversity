# -----------------------------
# APGR-2: Freshwater Drum
# Aplodinotus grunniens
# Alabama River four-section dataset
# -----------------------------

library(ggplot2)
library(geosphere)
library(maps)

# -----------------------------
# 0) paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/APGR-2"
data_dir <- file.path(base_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) section coordinates
# site order:
# 1 = JBR
# 2 = MFR
# 3 = CLL
# 4 = LAR
#
# built from thesis Table 1:
# JBR = ALR270
# MFR = mean(ALR236, ALR221, ALR181)
# CLL = mean(ALR132, ALR109)
# LAR = mean(ALR72, ALR47, ALR3)
# -----------------------------
APGR_2_coords <- data.frame(
  site = as.character(1:4),
  lat = c(
    32.354761,
    32.341152,
    32.008633,
    31.398635
  ),
  lon = c(
    -86.493584,
    -86.929366,
    -87.417859,
    -87.704836
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) helper: fill symmetric matrix from lower triangle
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
# Rotar et al. 2025 Table 4
# Freshwater Drum values are ABOVE the diagonal
#
# original table order:
# LAR, CLL, MFR, JBR
#
# workflow order here:
# JBR, MFR, CLL, LAR
# -----------------------------
fst_vals <- c(
  # row 2: MFR vs JBR
  0.003,
  
  # row 3: CLL vs JBR, CLL vs MFR
  0.002, 0.004,
  
  # row 4: LAR vs JBR, LAR vs MFR, LAR vs CLL
  0.003, 0.004, 0.001
)

APGR_2_fst <- fill_sym_from_lower(
  pops = APGR_2_coords$site,
  vals = fst_vals,
  diag_val = 0
)

APGR_2_fst[APGR_2_fst < 0] <- 0

# -----------------------------
# 4) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(APGR_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- APGR_2_coords$site
colnames(geo_dist_km) <- APGR_2_coords$site

# -----------------------------
# 5) IBD dataframe
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(APGR_2_fst)[row(APGR_2_fst)[upper.tri(APGR_2_fst)]],
  site2   = colnames(APGR_2_fst)[col(APGR_2_fst)[upper.tri(APGR_2_fst)]],
  fst     = APGR_2_fst[upper.tri(APGR_2_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 6) map
# plot only in RStudio; do not save
# -----------------------------
usa <- map_data("state")
can <- map_data("world", region = "Canada")

xpad <- 1.2
ypad <- 0.9

xlim_map <- c(
  min(APGR_2_coords$lon) - xpad,
  max(APGR_2_coords$lon) + xpad
)

ylim_map <- c(
  min(APGR_2_coords$lat) - ypad,
  max(APGR_2_coords$lat) + ypad
)

map_plot <- ggplot() +
  geom_polygon(
    data = usa,
    aes(x = long, y = lat, group = group),
    fill = "grey92", color = "black", linewidth = 0.2
  ) +
  geom_polygon(
    data = can,
    aes(x = long, y = lat, group = group),
    fill = "grey86", color = "black", linewidth = 0.2
  ) +
  geom_point(
    data = APGR_2_coords,
    aes(x = lon, y = lat),
    size = 2.7
  ) +
  geom_text(
    data = APGR_2_coords,
    aes(x = lon, y = lat, label = site),
    nudge_y = 0.05,
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
    title = "APGR-2 sampling sites"
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
    title = "APGR-2 isolation by distance"
  )

print(ibd_plot)

# -----------------------------
# 8) quick checks
# -----------------------------
stopifnot(nrow(APGR_2_coords) == 4)
stopifnot(identical(rownames(APGR_2_fst), APGR_2_coords$site))
stopifnot(identical(colnames(APGR_2_fst), APGR_2_coords$site))
stopifnot(isTRUE(all.equal(APGR_2_fst, t(APGR_2_fst))))
stopifnot(length(fst_vals) == nrow(APGR_2_coords) * (nrow(APGR_2_coords) - 1) / 2)

# -----------------------------
# 9) save RData
# -----------------------------
save(
  APGR_2_fst,
  APGR_2_coords,
  file = file.path(data_dir, "APGR-2.RData")
)

# -----------------------------
# 10) export coords csv
# -----------------------------
write.csv(
  APGR_2_coords,
  file = file.path(data_dir, "APGR-2_coords.csv"),
  row.names = FALSE
)