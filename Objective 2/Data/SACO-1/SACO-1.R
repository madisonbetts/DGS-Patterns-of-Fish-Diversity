# -----------------------------
# SACL-1 bull trout
# site coordinates + FST matrix + IBD plot
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SACL-1/"

# -----------------------------
# 0) site coordinates
# -----------------------------
SACL_1_coords <- data.frame(
  site = c("UK","KI","AK","BO","CE","QU","MQ","LQ",
           "LO","AR","TR","MC","LI","HA","IS","UI"),
  lat = c(48.975420,48.963626,48.880076,48.873800,
          48.872593,48.829246,48.822867,48.807212,
          48.755303,48.706587,48.680406,48.583148,
          48.590829,48.517693,48.422121,48.421006),
  lon = c(-114.177208,-114.302978,-114.198573,-114.154136,
          -114.056863,-114.095342,-114.141250,-114.172885,
          -114.077180,-113.885388,-113.909754,-113.917419,
          -113.770487,-113.770215,-113.493614,-113.507906),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
ggplot(SACL_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.01, size = 4) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "SACL-1 sampling locations"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(SACL_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- SACL_1_coords$site
colnames(geo_dist_km) <- SACL_1_coords$site

# -----------------------------
# 3) pairwise FST matrix
# -----------------------------
SACL_1_fst <- matrix(
  0,
  nrow = nrow(SACL_1_coords),
  ncol = nrow(SACL_1_coords),
  dimnames = list(SACL_1_coords$site, SACL_1_coords$site)
)

SACL_1_fst["KI", "UK"] <- 0.241

SACL_1_fst["AK", c("UK","KI")] <- c(0.385, 0.081)

SACL_1_fst["BO", c("UK","KI","AK")] <- c(0.377, 0.068, 0.138)

SACL_1_fst["CE", c("UK","KI","AK","BO")] <- c(0.485, 0.175, 0.244, 0.212)

SACL_1_fst["QU", c("UK","KI","AK","BO","CE")] <- c(0.450, 0.145, 0.202, 0.163, 0.005)

SACL_1_fst["MQ", c("UK","KI","AK","BO","CE","QU")] <- c(0.514, 0.175, 0.243, 0.210, 0.001, 0.012)

SACL_1_fst["LQ", c("UK","KI","AK","BO","CE","QU","MQ")] <- c(0.411, 0.104, 0.167, 0.126, 0.058, 0.048, 0.063)

SACL_1_fst["LO", c("UK","KI","AK","BO","CE","QU","MQ","LQ")] <- c(0.396, 0.088, 0.165, 0.148, 0.154, 0.121, 0.168, 0.124)

SACL_1_fst["AR", c("UK","KI","AK","BO","CE","QU","MQ","LQ","LO")] <- c(0.615, 0.357, 0.368, 0.335, 0.524, 0.482, 0.568, 0.413, 0.454)

SACL_1_fst["TR", c("UK","KI","AK","BO","CE","QU","MQ","LQ","LO","AR")] <- c(0.561, 0.333, 0.352, 0.313, 0.503, 0.458, 0.541, 0.392, 0.430, 0.015)

SACL_1_fst["MC", c("UK","KI","AK","BO","CE","QU","MQ","LQ","LO","AR","TR")] <- c(0.297, 0.006, 0.106, 0.073, 0.192, 0.159, 0.202, 0.118, 0.116, 0.385, 0.364)

SACL_1_fst["LI", c("UK","KI","AK","BO","CE","QU","MQ","LQ","LO","AR","TR","MC")] <- c(0.349, 0.059, 0.132, 0.162, 0.227, 0.198, 0.245, 0.172, 0.159, 0.454, 0.429, 0.078)

SACL_1_fst["HA", c("UK","KI","AK","BO","CE","QU","MQ","LQ","LO","AR","TR","MC","LI")] <- c(0.550, 0.301, 0.345, 0.322, 0.396, 0.351, 0.430, 0.340, 0.326, 0.641, 0.621, 0.318, 0.311)

SACL_1_fst["IS", c("UK","KI","AK","BO","CE","QU","MQ","LQ","LO","AR","TR","MC","LI","HA")] <- c(0.421, 0.205, 0.265, 0.319, 0.343, 0.307, 0.371, 0.284, 0.283, 0.562, 0.542, 0.219, 0.226, 0.431)

SACL_1_fst["UI", c("UK","KI","AK","BO","CE","QU","MQ","LQ","LO","AR","TR","MC","LI","HA","IS")] <- c(0.561, 0.275, 0.333, 0.363, 0.375, 0.338, 0.412, 0.322, 0.333, 0.658, 0.625, 0.297, 0.306, 0.520, 0.216)

SACL_1_fst[upper.tri(SACL_1_fst)] <- t(SACL_1_fst)[upper.tri(SACL_1_fst)]
diag(SACL_1_fst) <- 0
SACL_1_fst[SACL_1_fst < 0] <- 0

# -----------------------------
# river distance matrix (km)
# same population order as SACL_1_fst
# values entered row-wise from lower triangle
# -----------------------------

# helper function:
# takes a vector of lower-triangle values (row-wise order)
# and reconstructs a full symmetric matrix
fill_sym_from_lower <- function(pops, vals, diag_val = 0) {
  
  n <- length(pops)
  
  # sanity check:
  # number of values must equal n(n-1)/2 for a full lower triangle
  stopifnot(length(vals) == n * (n - 1) / 2)
  
  # initialize empty matrix
  mat <- matrix(0,
                nrow = n,
                ncol = n,
                dimnames = list(pops, pops))
  
  # index through the vector of values
  k <- 1
  
  # fill LOWER triangle row-wise:
  # row 2 gets 1 value, row 3 gets 2 values, etc.
  for (i in 2:n) {
    for (j in 1:(i - 1)) {
      
      mat[i, j] <- vals[k]   # lower triangle
      mat[j, i] <- vals[k]   # mirror to upper triangle
      
      k <- k + 1             # move to next value
    }
  }
  
  # set diagonal (distance of site to itself)
  diag(mat) <- diag_val
  
  return(mat)
}


# -----------------------------
# vector of river distances (km)
# -----------------------------
rivdist_vals <- c(
  
  # KI
  3.8,
  
  # AK
  50.4, 46.6,
  
  # BO
  42.7, 38.9, 27.9,
  
  # CE
  63.9, 60.1, 49.1, 40.5,
  
  # QU
  60.8, 57.1, 46.1, 37.5, 3.0,
  
  # MQ
  60.4, 56.6, 45.6, 37.0, 3.5, 0.4,
  
  # LQ
  58.4, 54.7, 43.7, 35.1, 5.4, 2.4, 1.9,
  
  # LO
  55.9, 52.2, 41.2, 32.6, 31.6, 28.5, 28.1, 26.2,
  
  # AR
  82.6, 78.8, 67.9, 59.2, 58.3, 55.2, 54.8, 52.8, 42.9,
  
  # TR
  80.2, 76.5, 65.5, 56.8, 55.9, 52.8, 52.4, 50.4, 40.5, 2.4,
  
  # MC
  98.2, 94.5, 83.5, 74.9, 73.9, 70.8, 70.4, 68.5, 58.5, 66.7, 64.3,
  
  # LI
  121.9,118.2,107.2, 98.6, 97.6, 94.5, 94.1, 92.1, 82.2, 90.4, 88.0, 31.2,
  
  # HA
  115.3,111.6,100.6, 91.9, 91.0, 87.9, 87.5, 85.5, 75.6, 83.8, 81.4, 24.6, 25.3,
  
  # IS
  169.1,165.3,154.3,145.7,144.7,141.7,141.3,139.3,129.4,137.5,135.1, 78.4, 79.1, 66.7,
  
  # UI
  169.8,166.1,155.1,146.5,145.5,142.4,142.0,140.1,130.1,138.3,135.9, 79.1, 79.8, 67.4, 0.8
)


# -----------------------------
# build symmetric river distance matrix
# -----------------------------
SACO_1_rivdists <- fill_sym_from_lower(
  pops = rownames(SACL_1_fst),  # CRITICAL: ensures identical ordering
  vals = rivdist_vals,
  diag_val = 0                  # distance to self = 0
)

# -----------------------------
# integrity checks (do not skip)
# -----------------------------

# same ordering as FST matrix
stopifnot(identical(rownames(SACO_1_rivdists), rownames(SACL_1_fst)))
stopifnot(identical(colnames(SACO_1_rivdists), colnames(SACL_1_fst)))

# symmetry check (distance matrix must be symmetric)
stopifnot(isTRUE(all.equal(SACO_1_rivdists, t(SACO_1_rivdists))))


# -----------------------------
# 4) pairwise dataframe for IBD plot
# use in-river distances
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(SACL_1_fst)[row(SACL_1_fst)[upper.tri(SACL_1_fst)]],
  site2   = colnames(SACL_1_fst)[col(SACL_1_fst)[upper.tri(SACL_1_fst)]],
  fst     = SACL_1_fst[upper.tri(SACL_1_fst)],
  dist_km = SACO_1_rivdists[upper.tri(SACO_1_rivdists)]
)

# -----------------------------
# 5) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "In-river distance (km)",
    y = expression(F[ST]),
    title = "SACL-1 isolation by distance"
  )

# -----------------------------
# 6) quick checks
# -----------------------------
stopifnot(identical(rownames(SACL_1_fst), SACL_1_coords$site))
stopifnot(identical(colnames(SACL_1_fst), SACL_1_coords$site))
stopifnot(isTRUE(all.equal(SACL_1_fst, t(SACL_1_fst))))

# -----------------------------
# 7) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  SACL_1_fst,
  SACL_1_coords,
  file = file.path(out_dir, "SACL-1.RData")
)