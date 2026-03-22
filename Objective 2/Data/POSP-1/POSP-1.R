# -----------------------------
# POSP-1 paddlefish (Polyodon spathula)
# Heist & Mustapha 2008, TAFS 137:909-915
# Microsatellite pairwise FST + site coordinates + IBD plot
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
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/POSP-1"
data_dir <- file.path(save_dir, "data")

# -----------------------------
# 0) site metadata + coordinates
# numeric site IDs match the FST matrix order
# Exact sample XYs were not provided in the paper.
# These coordinates are representative points inferred from the named
# localities / dams / gages closest to the reported sampling areas.
# -----------------------------
POSP_1_coords <- data.frame(
  site = as.character(1:12),
  site_code = c("BN", "LA", "OK", "GL", "IL", "TB", "WV", "TL", "AR", "WR", "ND", "TN"),
  #locality = c(
  #  "Bayou Nezpique, Mermentau River, Louisiana",
  #  "Red River near Natchitoches, Louisiana",
  #  "Red River below Lake Texoma, Oklahoma",
  #  "Grand Lake / Neosho-Grand River, Oklahoma",
  #  "Big Muddy River near Mississippi confluence, Illinois",
  #  "Tombigbee River, Alabama",
  #  "Ohio River, R. C. Byrd tailwater, West Virginia",
  #  "Truman Lake, Missouri",
  #  "Arkansas River near Fort Smith, Arkansas",
  #  "White River, Arkansas",
  #  "Missouri River below Yellowstone confluence, North Dakota",
  #  "Tennessee River, Tennessee"
  #),
  #state = c("LA", "LA", "OK", "OK", "IL", "AL", "WV", "MO", "AR", "AR", "ND", "TN"),
  lat = c(
    30.2299751,
    31.8182,
    33.8182,
    36.4721250,
    37.7486100,
    32.8540,
    38.6818700,
    38.2664000,
    35.3917577,
    35.6053,
    47.9783300,
    35.3875
  ),
  lon = c(
    -92.6270630,
    -93.0849,
    -96.5725,
    -95.0392583,
    -89.3458300,
    -88.1580,
    -82.1854250,
    -93.3980650,
    -94.4324373,
    -91.2889,
    -103.9822200,
    -87.9861
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) helper function:
# takes a vector of lower-triangle values (row-wise order)
# and reconstructs a full symmetric matrix
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
# 2) pairwise FST matrix
# Table 2 values entered from the LOWER triangle row-wise.
# Original site order in the paper:
# BN, LA, OK, GL, IL, TB, WV, TL, AR, WR, ND, TN
# -----------------------------
fst_vals <- c(
  # LA
  0.0316,
  
  # OK
  0.0475, 0.0007,
  
  # GL
  0.1183, 0.0533, 0.0579,
  
  # IL
  0.0485, 0.0000, 0.0076, 0.0649,
  
  # TB
  0.1872, 0.0956, 0.1036, 0.1033, 0.1067,
  
  # WV
  0.0510, 0.0189, 0.0324, 0.0782, 0.0139, 0.1379,
  
  # TL (printed as TR in parsed table, locality is Truman Lake)
  0.0784, 0.0274, 0.0379, 0.0871, 0.0207, 0.1271, 0.0122,
  
  # AR
  0.0653, 0.0100, 0.0219, 0.0474, 0.0102, 0.0971, 0.0320, 0.0359,
  
  # WR
  0.0950, 0.0319, 0.0298, 0.1172, 0.0160, 0.1531, 0.0170, 0.0122, 0.0414,
  
  # ND
  0.0580, 0.0041, 0.0084, 0.0529, 0.0043, 0.0701, 0.0273, 0.0263, 0.0142, 0.0349,
  
  # TN
  0.0396, 0.0010, 0.0055, 0.0582, 0.0000, 0.1094, 0.0222, 0.0274, 0.0136, 0.0362, 0.0066
)

POSP_1_fst <- fill_sym_from_lower(
  pops = POSP_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

# enforce workflow rule
POSP_1_fst[POSP_1_fst < 0] <- 0

# -----------------------------
# 3) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(POSP_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- POSP_1_coords$site
colnames(geo_dist_km) <- POSP_1_coords$site

# -----------------------------
# 4) map of sampling locations using relevant states
# -----------------------------
relevant_states <- c(
  "Louisiana", "Oklahoma", "Illinois", "Alabama", "West Virginia",
  "Missouri", "Arkansas", "North Dakota", "Tennessee"
)

state_sf <- tigris::states(cb = TRUE, year = 2022, class = "sf") |>
  dplyr::filter(NAME %in% relevant_states)

site_sf <- st_as_sf(POSP_1_coords, coords = c("lon", "lat"), crs = 4326)

map_plot <- ggplot() +
  geom_sf(data = state_sf, fill = "grey95", color = "grey45", linewidth = 0.3) +
  geom_sf(data = site_sf, size = 2.4) +
  geom_text_repel(
    data = POSP_1_coords,
    aes(x = lon, y = lat, label = paste0(site, ": ", site_code)),
    size = 3.2,
    min.segment.length = 0,
    seed = 1
  ) +
  coord_sf(
    xlim = c(min(POSP_1_coords$lon) - 3, max(POSP_1_coords$lon) + 3),
    ylim = c(min(POSP_1_coords$lat) - 2, max(POSP_1_coords$lat) + 2),
    expand = FALSE
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "POSP-1 sampling locations",
    subtitle = "Labels show numeric site ID and paper abbreviation"
  )

print(map_plot)

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
POSP_1_ibd_df <- data.frame(
  site1 = rownames(POSP_1_fst)[row(POSP_1_fst)[upper.tri(POSP_1_fst)]],
  site2 = colnames(POSP_1_fst)[col(POSP_1_fst)[upper.tri(POSP_1_fst)]],
  fst = POSP_1_fst[upper.tri(POSP_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)],
  stringsAsFactors = FALSE
) |>
  left_join(POSP_1_coords[, c("site", "site_code")], by = c("site1" = "site")) |>
  rename(site1_code = site_code) |>
  left_join(POSP_1_coords[, c("site", "site_code")], by = c("site2" = "site")) |>
  rename(site2_code = site_code)

# -----------------------------
# 6) IBD plot
# -----------------------------
ibd_plot <- ggplot(POSP_1_ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.7, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(pairwise~F[ST]),
    title = "POSP-1 isolation by distance"
  )

print(ibd_plot)

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(POSP_1_fst), POSP_1_coords$site))
stopifnot(identical(colnames(POSP_1_fst), POSP_1_coords$site))
stopifnot(isTRUE(all.equal(POSP_1_fst, t(POSP_1_fst))))
stopifnot(all(diag(POSP_1_fst) == 0))
stopifnot(all(POSP_1_fst >= 0))

# -----------------------------
# 8) save .RData to data subfolder
# -----------------------------
save(
  POSP_1_fst,
  POSP_1_coords,
  file = file.path(data_dir, "POSP-1.RData")
)
