
# -----------------------------
# RHOS-4 Speckled Dace
# basin coordinates + FST matrix + IBD plot
# Su et al. 2022, Transactions of the American Fisheries Society
# Pairwise FST from Supplemental Table S.2
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
save_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/RHOS-4"

# -----------------------------
# 0) site coordinates
# basin-level representative coordinates
# order follows Supplemental Table S.2
# -----------------------------
RHOS_4_coords <- data.frame(
  site = as.character(1:10),
  #site_name = c(
  #  "Santa Ana",
  #  "Amargosa",
  #  "Owens",
  #  "Long Valley",
  #  "Lahontan",
  #  "Sacramento",
  #  "Central California",
  #  "Klamath",
  #  "Butte Lake",
  #  "Warner"
  #),
  lat = c(
    33.98, 36.30, 37.30, 37.70, 39.50,
    38.60, 36.80, 42.20, 40.90, 41.80
  ),
  lon = c(
    -117.33, -116.30, -118.40, -118.80, -119.80,
    -121.50, -120.30, -123.80, -121.30, -120.20
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling basins
# -----------------------------
ggplot(RHOS_4_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.20, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "RHOS-4 sampling basins"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(RHOS_4_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- RHOS_4_coords$site
colnames(geo_dist_km) <- RHOS_4_coords$site

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
# Supplemental Table S.2
# order:
# 1 Santa Ana
# 2 Amargosa
# 3 Owens
# 4 Long Valley
# 5 Lahontan
# 6 Sacramento
# 7 Central California
# 8 Klamath
# 9 Butte Lake
# 10 Warner
# -----------------------------
fst_vals <- c(

  # 2 Amargosa
  0.683038,

  # 3 Owens
  0.655648, 0.163662,

  # 4 Long Valley
  0.495095, 0.382165, 0.300142,

  # 5 Lahontan
  0.549713, 0.328260, 0.246113, 0.259542,

  # 6 Sacramento
  0.645057, 0.555592, 0.511999, 0.443290, 0.412116,

  # 7 Central California
  0.753472, 0.620699, 0.586420, 0.488401, 0.463870, 0.194572,

  # 8 Klamath
  0.549713, 0.571066, 0.535016, 0.456907, 0.428349, 0.165698, 0.326771,

  # 9 Butte Lake
  0.679071, 0.576948, 0.533901, 0.454898, 0.424014, 0.083905, 0.228633, 0.217210,

  # 10 Warner
  0.679669, 0.559563, 0.507790, 0.418762, 0.388788, 0.288543, 0.429950, 0.319681, 0.329080
)

RHOS_4_fst <- fill_sym_from_lower(
  pops = RHOS_4_coords$site,
  vals = fst_vals,
  diag_val = 0
)

RHOS_4_fst[RHOS_4_fst < 0] <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(RHOS_4_fst)[row(RHOS_4_fst)[upper.tri(RHOS_4_fst)]],
  site2   = colnames(RHOS_4_fst)[col(RHOS_4_fst)[upper.tri(RHOS_4_fst)]],
  fst     = RHOS_4_fst[upper.tri(RHOS_4_fst)],
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
    title = "RHOS-4 isolation by distance"
  )

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(RHOS_4_fst), RHOS_4_coords$site))
stopifnot(identical(colnames(RHOS_4_fst), RHOS_4_coords$site))
stopifnot(isTRUE(all.equal(RHOS_4_fst, t(RHOS_4_fst))))
stopifnot(length(fst_vals) == nrow(RHOS_4_coords) * (nrow(RHOS_4_coords) - 1) / 2)

# -----------------------------
# 8) save RData
# -----------------------------
out_dir <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  RHOS_4_fst,
  RHOS_4_coords,
  file = file.path(out_dir, "RHOS-4.RData")
)
