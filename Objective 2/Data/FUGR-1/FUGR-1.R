# -----------------------------
# FUGR-1 Gulf killifish
# site coordinates + FST matrix + IBD plot
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/FUGR-1"

# -----------------------------
# 0) site coordinates
# estimated from named localities in Fig. 1 / Table 1
# -----------------------------
FUGR_1_coords <- data.frame(
  site = c("AR", "OC", "LC", "LV", "KC", "DI", "MB", "TB", "SS", "BB"),
  lat = c(
    27.84,  # Port Aransas, TX
    28.45,  # Port O'Connor, TX
    30.21,  # Lake Charles, LA
    29.25,  # Leeville, LA
    30.41,  # Kiln, MS
    30.25,  # Dauphin Island, AL
    30.68,  # Mobile Bay, AL
    27.95,  # Tampa Bay, FL
    27.34,  # Sarasota, FL
    26.34   # Bonita Beach, FL
  ),
  lon = c(
    -97.07,  # Port Aransas, TX
    -96.41,  # Port O'Connor, TX
    -93.22,  # Lake Charles, LA
    -90.21,  # Leeville, LA
    -89.44,  # Kiln, MS
    -88.11,  # Dauphin Island, AL
    -88.04,  # Mobile Bay, AL
    -82.46,  # Tampa Bay, FL
    -82.54,  # Sarasota, FL
    -81.84   # Bonita Beach, FL
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
ggplot(FUGR_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.18, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "FUGR-1 sampling locations"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(FUGR_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- FUGR_1_coords$site
colnames(geo_dist_km) <- FUGR_1_coords$site

# -----------------------------
# 3) helper function:
# takes a vector of upper-triangle values (row-wise order)
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
# 4) pairwise FST matrix
# values transcribed from the UPPER triangle of Table 2
# row-wise: row 1 has 9 values, row 2 has 8 values, etc.
# -----------------------------
fst_vals <- c(
  
  # AR vs others
  0.003, 0.014, 0.022, 0.025, 0.055, 0.044, 0.075, 0.061, 0.049,
  
  # OC vs others
  0.013, 0.014, 0.019, 0.046, 0.042, 0.079, 0.065, 0.062,
  
  # LC vs others
  0.004, 0.014, 0.051, 0.045, 0.106, 0.087, 0.078,
  
  # LV vs others
  0.010, 0.036, 0.039, 0.103, 0.086, 0.078,
  
  # KC vs others
  0.035, 0.025, 0.092, 0.084, 0.077,
  
  # DI vs others
  0.013, 0.086, 0.075, 0.073,
  
  # MB vs others
  0.069, 0.072, 0.054,
  
  # TB vs others
  0.012, 0.035,
  
  # SS vs BB
  0.028
)

FUGR_1_fst <- fill_sym_from_upper(
  pops = FUGR_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(FUGR_1_fst)[row(FUGR_1_fst)[upper.tri(FUGR_1_fst)]],
  site2   = colnames(FUGR_1_fst)[col(FUGR_1_fst)[upper.tri(FUGR_1_fst)]],
  fst     = FUGR_1_fst[upper.tri(FUGR_1_fst)],
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
    title = "FUGR-1 isolation by distance"
  )

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(FUGR_1_fst), FUGR_1_coords$site))
stopifnot(identical(colnames(FUGR_1_fst), FUGR_1_coords$site))
stopifnot(isTRUE(all.equal(FUGR_1_fst, t(FUGR_1_fst))))
stopifnot(length(fst_vals) == nrow(FUGR_1_coords) * (nrow(FUGR_1_coords) - 1) / 2)

# -----------------------------
# 8) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  FUGR_1_fst,
  FUGR_1_coords,
  file = file.path(out_dir, "FUGR-1.RData")
)