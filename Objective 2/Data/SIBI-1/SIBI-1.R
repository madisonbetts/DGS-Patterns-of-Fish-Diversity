# -----------------------------
# SIBI-1 Lahontan tui chub
# Siphateles bicolor ssp.
# sites + pairwise FST + IBD plot
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SIBI-1"

# -----------------------------
# 0) site coordinates
# estimated from named localities + paper map
# codes/order follow Table 1 / Table 2
# -----------------------------
SIBI_1_coords <- data.frame(
  site = as.character(1:11),
  code = c("TPZ", "LSL", "SPN", "TKS", "PYR", "TWN", "WLK", "STW", "DXV", "EWR", "SFH"),
  location = c(
    "Topaz Lake, NV",
    "Little Soda Lake, NV",
    "Spooner Lake, NV",
    "Tahoe Keys, CA",
    "Pyramid Lake, NV",
    "Twin Lakes, CA",
    "Walker Lake, NV",
    "Stillwater National Wildlife Refuge, NV",
    "Dixie Valley (Casey Pond), NV",
    "East Fork Walker River, CA",
    "South Fork Reservoir, NV"
  ),
  lat = c(
    38.6805,   # Topaz Lake
    39.5200,   # Little Soda Lake / Soda Lakes
    39.1060,   # Spooner Lake
    38.9257,   # Tahoe Keys
    39.9846,   # Pyramid Lake near Nixon
    38.1580,   # Twin Lakes (Bridgeport area, south of Bridgeport)
    38.6922,   # Walker Lake
    39.6560,   # Stillwater NWR / Fallon marsh area
    39.8200,   # Casey Pond / Dixie Valley (estimated from map)
    38.4140,   # East Fork Walker River near CA/NV state line
    40.6900    # South Fork Reservoir
  ),
  lon = c(
    -119.5340, # Topaz Lake
    -118.8800, # Little Soda Lake
    -119.9020, # Spooner Lake
    -120.0082, # Tahoe Keys
    -119.5010, # Pyramid Lake
    -119.3300, # Twin Lakes
    -118.7361, # Walker Lake
    -118.4010, # Stillwater NWR
    -117.9200, # Casey Pond / Dixie Valley (estimated)
    -119.1660, # East Fork Walker River
    -115.7800  # South Fork Reservoir
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
ggplot(SIBI_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.12, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "SIBI-1 sampling locations"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(SIBI_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- SIBI_1_coords$site
colnames(geo_dist_km) <- SIBI_1_coords$site
diag(geo_dist_km) <- 0

# -----------------------------
# 3) helper function:
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
# 4) pairwise FST matrix
# values transcribed from Table 2 LOWER triangle
# order: TPZ LSL SPN TKS PYR TWN WLK STW DXV EWR SFH
# -----------------------------
fst_vals <- c(
  
  # LSL
  0.080,
  
  # SPN
  0.082, 0.141,
  
  # TKS
  0.151, 0.193, 0.063,
  
  # PYR
  0.116, 0.177, 0.070, 0.024,
  
  # TWN
  0.070, 0.133, 0.079, 0.075, 0.062,
  
  # WLK
  0.031, 0.079, 0.054, 0.059, 0.044, 0.040,
  
  # STW
  0.065, 0.116, 0.067, 0.095, 0.081, 0.088, 0.038,
  
  # DXV
  0.098, 0.120, 0.149, 0.217, 0.202, 0.174, 0.103, 0.164,
  
  # EWR
  0.051, 0.099, 0.122, 0.144, 0.126, 0.038, 0.042, 0.093, 0.157,
  
  # SFH
  0.056, 0.138, 0.059, 0.083, 0.058, 0.071, 0.038, 0.048, 0.151, 0.095
)

SIBI_1_fst <- fill_sym_from_lower(
  pops = SIBI_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

# set negative FST to 0 and force diagonal = 0
SIBI_1_fst <- pmax(SIBI_1_fst, 0)
diag(SIBI_1_fst) <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(SIBI_1_fst)[row(SIBI_1_fst)[upper.tri(SIBI_1_fst)]],
  site2   = colnames(SIBI_1_fst)[col(SIBI_1_fst)[upper.tri(SIBI_1_fst)]],
  fst     = SIBI_1_fst[upper.tri(SIBI_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 6) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "SIBI-1 isolation by distance"
  )

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(SIBI_1_fst), SIBI_1_coords$site))
stopifnot(identical(colnames(SIBI_1_fst), SIBI_1_coords$site))
stopifnot(isTRUE(all.equal(SIBI_1_fst, t(SIBI_1_fst))))
stopifnot(all(diag(SIBI_1_fst) == 0))
stopifnot(all(diag(geo_dist_km) == 0))

# -----------------------------
# 8) save RData
# -----------------------------
out_dir <- file.path(base_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  SIBI_1_fst,
  SIBI_1_coords,
  file = file.path(out_dir, "SIBI-1.RData")
)