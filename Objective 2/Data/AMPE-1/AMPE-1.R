############################
## AMPE_1 workflow
## 1) site coordinates
## 2) pairwise FST matrix
## 3) long-format site-pair data frame
## 4) optional IBD plot using straight-line distance
############################

library(ggplot2)
library(ggrepel)
library(geosphere)

# -------------------------
# 1) Site coordinates
# -------------------------
# These are the sampling coordinates for the 16 rivers.
# Object name requested: AMPE_1_coords

AMPE_1_coords <- data.frame(
  site = c("ER","EF","BC","DC","LK","RED","LM","HR","SC","MA","Syd","TH","GR","RAS","RR","CC"),
  full_name = c(
    "Eel River",
    "East Fork White River",
    "Big Creek",
    "Deer Creek",
    "Licking River",
    "Red River",
    "Little Muskingum River",
    "Hocking River main",
    "Salt Creek",
    "Maumee River main",
    "Sydenham River",
    "Thames River upper",
    "Grand River upper",
    "Riviere au Saumon",
    "Richelieu River",
    "Champlain Canal"
  ),
  lat = c(
    40.828056,
    39.138611,
    38.809167,
    39.500556,
    38.208333,
    37.819722,
    39.411667,
    39.300833,
    39.433333,
    41.084167,
    42.646944,
    42.931944,
    43.127778,
    44.999167,
    45.635000,
    43.352500
  ),
  lon = c(
    -86.113889,
    -85.893889,
    -85.643889,
    -86.930278,
    -83.680278,
    -83.575833,
    -81.358611,
    -81.963889,
    -82.680000,
    -85.019722,
    -82.009722,
    -81.426389,
    -80.199167,
    -74.510556,
    -73.190556,
    -73.495556
  ),
  stringsAsFactors = FALSE
)

AMPE_1_coords

# -------------------------
# Quick map of sampling sites
# -------------------------
ggplot(AMPE_1_coords, aes(lon, lat)) +
  geom_point(size = 2.5) +
  geom_text_repel(aes(label = site), size = 3.5) +
  coord_fixed() +
  theme_bw() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "AMPE_1 sampling localities"
  )

# -------------------------
# 2) Pairwise FST matrix
# -------------------------
# Values below are entered ROW-WISE from the published table
# (row 2, then row 3, then row 4, etc.)

sites <- AMPE_1_coords$site

fst_vals <- c(
  0.075,
  0.085, 0.011,
  0.076, 0.024, 0.009,
  0.160, 0.103, 0.081, 0.078,
  0.144, 0.089, 0.069, 0.063, 0.032,
  0.103, 0.072, 0.063, 0.042, 0.075, 0.049,
  0.164, 0.119, 0.085, 0.073, 0.080, 0.046, 0.053,
  0.153, 0.139, 0.123, 0.112, 0.069, 0.060, 0.075, 0.081,
  0.081, 0.047, 0.058, 0.077, 0.148, 0.145, 0.120, 0.165, 0.162,
  0.062, 0.071, 0.084, 0.084, 0.172, 0.159, 0.121, 0.175, 0.154, 0.054,
  0.053, 0.047, 0.054, 0.053, 0.123, 0.110, 0.083, 0.126, 0.134, 0.050, 0.021,
  0.099, 0.077, 0.090, 0.088, 0.156, 0.149, 0.109, 0.168, 0.165, 0.090, 0.044, 0.055,
  0.114, 0.070, 0.056, 0.060, 0.159, 0.147, 0.115, 0.130, 0.171, 0.096, 0.116, 0.081, 0.105,
  0.148, 0.096, 0.098, 0.086, 0.184, 0.170, 0.118, 0.146, 0.190, 0.125, 0.143, 0.098, 0.093, 0.060,
  0.259, 0.170, 0.184, 0.193, 0.279, 0.267, 0.224, 0.237, 0.281, 0.243, 0.289, 0.204, 0.205, 0.155, 0.175
)

AMPE_1_fst <- matrix(
  0,
  nrow = length(sites),
  ncol = length(sites),
  dimnames = list(sites, sites)
)

# lower-triangle coordinates in ROW-WISE order
lt_idx <- which(lower.tri(AMPE_1_fst), arr.ind = TRUE)
lt_idx <- lt_idx[order(lt_idx[, 1], lt_idx[, 2]), ]

# fill lower triangle row-wise
AMPE_1_fst[lt_idx] <- fst_vals

# mirror to upper triangle
AMPE_1_fst[upper.tri(AMPE_1_fst)] <- t(AMPE_1_fst)[upper.tri(AMPE_1_fst)]

AMPE_1_fst  # QC
round(AMPE_1_fst, 3)

# -------------------------
# 3) Convert FST matrix to long format
# -------------------------
# Keep only one copy of each pair (lower triangle).

pair_idx <- lower.tri(AMPE_1_fst)

AMPE_1_ibd_df <- data.frame(
  site1 = rownames(AMPE_1_fst)[row(AMPE_1_fst)[pair_idx]],
  site2 = colnames(AMPE_1_fst)[col(AMPE_1_fst)[pair_idx]],
  fst = AMPE_1_fst[pair_idx],
  stringsAsFactors = FALSE
)

head(AMPE_1_ibd_df)

# -------------------------
# 4) Optional straight-line distance + IBD plot
# -------------------------
# No river-distance matrix is available for this dataset,
# so this uses great-circle geographic distance (km).

AMPE_1_geodists <- geosphere::distm(
  x = AMPE_1_coords[, c("lon", "lat")],
  fun = geosphere::distGeo
) / 1000

rownames(AMPE_1_geodists) <- AMPE_1_coords$site
colnames(AMPE_1_geodists) <- AMPE_1_coords$site

AMPE_1_ibd_df$geo_distance_km <- AMPE_1_geodists[pair_idx]

ggplot(AMPE_1_ibd_df, aes(geo_distance_km, fst)) +
  geom_point(size = 2.3) +
  geom_smooth(method = "lm", se = TRUE) +
  theme_bw() +
  labs(
    x = "Straight-line geographic distance (km)",
    y = expression(pairwise~F[ST]),
    title = "Isolation by distance for AMPE_1"
  )

# -------------------------
# Save outputs
# -------------------------

# path
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/AMPE-1/data"

# save fst, coords, geographic distances, and long-format df
save(
  AMPE_1_fst,
  AMPE_1_coords,
  #AMPE_1_geodists,
  #AMPE_1_ibd_df,
  file = file.path(save_dir, "AMPE_1.RData")
)