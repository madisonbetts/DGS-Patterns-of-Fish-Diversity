# -----------------------------
# PEFL-2 Yellow perch
# Perca flavescens
# sites + pairwise FST + IBD plot
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/PEFL-2"

# -----------------------------
# 0) site coordinates
# approximate from Fig. 1 + named localities
# order follows Table 5 exactly:
# EMM, PP, AB, SL, BDW, CHP, WP, CI, CP, MC
# -----------------------------
PEFL_2_coords <- data.frame(
  site = as.character(1:10),
  code = c("EMM", "PP", "AB", "SL", "BDW", "CHP", "WP", "CI", "CP", "MC"),
  location = c(
    "11 Mile Marsh, MI",
    "Palmer's Point, ON",
    "Ashmun Bay, MI",
    "Soldier Lake, MI",
    "Bai de Wasai, ON",
    "Churchville Point, ON",
    "Whipple Point, ON",
    "Cook Island, MI/ON",
    "Cedar Point, MI",
    "Mission Creek, MI"
  ),
  lat = c(
    46.27,  # 11 Mile Marsh
    46.50,  # Palmer's Point
    46.49,  # Ashmun Bay
    46.46,  # Soldier Lake
    46.43,  # Bai de Wasai
    46.50,  # Churchville Point
    46.37,  # Whipple Point
    46.52,  # Cook Island
    46.45,  # Cedar Point
    46.48   # Mission Creek
  ),
  lon = c(
    -84.10, # 11 Mile Marsh
    -83.94, # Palmer's Point
    -84.23, # Ashmun Bay
    -84.45, # Soldier Lake
    -84.00, # Bai de Wasai
    -83.85, # Churchville Point
    -83.96, # Whipple Point
    -84.08, # Cook Island
    -84.33, # Cedar Point
    -84.18  # Mission Creek / lower shipping channel
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
ggplot(PEFL_2_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.04, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "PEFL-2 sampling locations"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(PEFL_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- PEFL_2_coords$site
colnames(geo_dist_km) <- PEFL_2_coords$site
diag(geo_dist_km) <- 0

# -----------------------------
# 3) helper function
# lower-triangle values entered row-wise
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
# Table 5 lower triangle, transcribed row-wise
# order: EMM, PP, AB, SL, BDW, CHP, WP, CI, CP, MC
# negative FST values truncated to 0 afterward
# -----------------------------
fst_vals <- c(
  
  # PP
  0.0005,
  
  # AB
  0.0470, 0.0208,
  
  # SL
  0.1918, 0.2224, 0.2903,
  
  # BDW
  0.0107, 0.0140, 0.0258, 0.2492,
  
  # CHP
  0.1237, 0.1540, 0.2188, 0.2035, 0.1814,
  
  # WP
  0.0930, 0.1477, 0.2314, 0.2250, 0.1688, 0.0547,
  
  # CI
  -0.0248, 0.0033, 0.0265, 0.2207, -0.0272, 0.1575, 0.1189,
  
  # CP
  -0.0148, 0.0270, 0.0929, 0.2000, 0.0157, 0.0772, 0.0572, -0.0215,
  
  # MC
  0.0295, -0.0107, 0.0289, 0.2505, 0.0338, 0.1614, 0.1586, 0.0403, 0.0421
)

PEFL_2_fst <- fill_sym_from_lower(
  pops = PEFL_2_coords$site,
  vals = fst_vals,
  diag_val = 0
)

PEFL_2_fst <- pmax(PEFL_2_fst, 0)
diag(PEFL_2_fst) <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(PEFL_2_fst)[row(PEFL_2_fst)[upper.tri(PEFL_2_fst)]],
  site2   = colnames(PEFL_2_fst)[col(PEFL_2_fst)[upper.tri(PEFL_2_fst)]],
  fst     = PEFL_2_fst[upper.tri(PEFL_2_fst)],
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
    title = "PEFL-2 isolation by distance"
  )

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(PEFL_2_fst), PEFL_2_coords$site))
stopifnot(identical(colnames(PEFL_2_fst), PEFL_2_coords$site))
stopifnot(isTRUE(all.equal(PEFL_2_fst, t(PEFL_2_fst))))
stopifnot(all(diag(PEFL_2_fst) == 0))
stopifnot(all(diag(geo_dist_km) == 0))

# -----------------------------
# 8) save RData
# -----------------------------
out_dir <- file.path(base_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  PEFL_2_fst,
  PEFL_2_coords,
  file = file.path(out_dir, "PEFL-2.RData")
)