# -----------------------------
# ETLE-1 Tuxedo Darter
# Washburn et al. 2020
# pairwise FST from Table 2 entered row-by-row
# -----------------------------

library(ggplot2)
library(ggrepel)
library(dplyr)
library(sf)
library(tigris)
library(geosphere)

options(tigris_use_cache = TRUE)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETLE-1"
data_dir <- file.path(save_dir, "data")

# -----------------------------
# 0) site coordinates
# numeric site IDs match the FST matrix order
# XYs are inferred reach centroids from the named Big South Fork localities
# -----------------------------
ETLE_1_coords <- data.frame(
  site = as.character(1:10),
  site_code = paste0("R", 1:10),
  site_name = c(
    "Station Camp Creek",
    "Big Island / Hurricane Branch",
    "Upstream Oil Well Branch",
    "Blue Heron / upstream Devil's Creek",
    "Roaring Paunch corridor",
    "Upstream Stover Branch",
    "Lower Stover-Rock corridor",
    "Rock Creek / Devil's Jump corridor",
    "Upper Yamacraw corridor",
    "Downstream KY-92 bridge / Yamacraw"
  ),
  lat = c(
    36.5462,
    36.5909,
    36.6242,
    36.6716,
    36.6798,
    36.7048,
    36.7110,
    36.7173,
    36.7205,
    36.7235
  ),
  lon = c(
    -84.6652,
    -84.6083,
    -84.5694,
    -84.5477,
    -84.5388,
    -84.5308,
    -84.5370,
    -84.5436,
    -84.5425,
    -84.5415
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) helper function:
# takes a vector of upper-triangle values entered row-wise
# and reconstructs a full symmetric matrix
# -----------------------------
fill_sym_from_upper <- function(pops, vals, diag_val = 0) {
  n <- length(pops)
  stopifnot(length(vals) == n * (n - 1) / 2)
  
  mat <- matrix(
    0,
    nrow = n,
    ncol = n,
    dimnames = list(pops, pops)
  )
  
  k <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      mat[i, j] <- vals[k]
      mat[j, i] <- vals[k]
      k <- k + 1
    }
  }
  
  diag(mat) <- diag_val
  mat
}

# -----------------------------
# 2) pairwise FST matrix
# Table 2 FST values are ABOVE the diagonal in the paper.
# Entered row-by-row:
# row R1 has 9 values, row R2 has 8 values, ..., row R9 has 1 value.
# -----------------------------
fst_vals <- c(
  # R1 vs R2:R10
  0.021, 0.021, 0.040, 0.025, 0.028, 0.021, 0.037, 0.021, 0.033,
  
  # R2 vs R3:R10
  0.002, 0.018, 0.014, 0.008, 0.005, 0.008, 0.006, 0.011,
  
  # R3 vs R4:R10
  0.010, 0.014, 0.015, 0.020, 0.012, 0.012, 0.023,
  
  # R4 vs R5:R10
  0.010, 0.016, 0.019, 0.012, 0.008, 0.023,
  
  # R5 vs R6:R10
  0.005, 0.012, 0.020, 0.007, 0.010,
  
  # R6 vs R7:R10
  0.015, 0.012, 0.002, 0.003,
  
  # R7 vs R8:R10
  0.004, 0.000, 0.014,
  
  # R8 vs R9:R10
  0.005, 0.014,
  
  # R9 vs R10
  0.001
)

ETLE_2_fst <- fill_sym_from_upper(
  pops = ETLE_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

# workflow rule
ETLE_2_fst[ETLE_2_fst < 0] <- 0
diag(ETLE_2_fst) <- 0

# -----------------------------
# 3) inspect matrix
# -----------------------------
print(round(ETLE_2_fst, 3))

# -----------------------------
# 4) map of sampling locations with US + Canada
# -----------------------------
states_sf <- tigris::states(cb = TRUE, year = 2022) |>
  sf::st_transform(4326)

canada_sf <- rnaturalearth::ne_countries(
  country = "Canada",
  scale = "medium",
  returnclass = "sf"
) |>
  sf::st_transform(4326)

x_pad <- 0.3
y_pad <- 0.2

xlim_use <- range(ETLE_1_coords$lon) + c(-x_pad, x_pad)
ylim_use <- range(ETLE_1_coords$lat) + c(-y_pad, y_pad)

p_map <- ggplot() +
  geom_sf(data = canada_sf, fill = "gray97", color = "gray75", linewidth = 0.2) +
  geom_sf(data = states_sf, fill = "gray95", color = "gray70", linewidth = 0.2) +
  geom_point(
    data = ETLE_1_coords,
    aes(x = lon, y = lat),
    size = 2.8
  ) +
  geom_text_repel(
    data = ETLE_1_coords,
    aes(x = lon, y = lat, label = paste0(site, ": ", site_code)),
    size = 3.2,
    max.overlaps = 100
  ) +
  coord_sf(
    xlim = xlim_use,
    ylim = ylim_use,
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    title = "ETLE-2 sampling reaches",
    x = "Longitude",
    y = "Latitude"
  )

print(p_map)

# -----------------------------
# 5) IBD plot
# use straight-line distance among reach centroids
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(ETLE_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ETLE_1_coords$site
colnames(geo_dist_km) <- ETLE_1_coords$site

ibd_df <- data.frame(
  site1   = rownames(ETLE_2_fst)[row(ETLE_2_fst)[upper.tri(ETLE_2_fst)]],
  site2   = colnames(ETLE_2_fst)[col(ETLE_2_fst)[upper.tri(ETLE_2_fst)]],
  fst     = ETLE_2_fst[upper.tri(ETLE_2_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    title = "ETLE-1 isolation by distance",
    x = "Euclidean distance among reach centroids (km)",
    y = expression(F[ST])
  )

print(p_ibd)

# -----------------------------
# 6) save
# -----------------------------
save(
  ETLE_2_fst,
  ETLE_1_coords,
  file = file.path(data_dir, "ETLE-1.RData")
)