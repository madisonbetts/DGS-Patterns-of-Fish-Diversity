# -----------------------------
# GINI-1 Gila nigra (headwater chub)
# sites + pairwise FST + IBD plot
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/GINI-1"

# -----------------------------
# 0) site coordinates
# approximate, anchored to named Gila locations
# -----------------------------
GINI_1_coords <- data.frame(
  site = c("1", "2", "3", "4"),
  code = c("WF1", "MF", "EF", "GM1"),
  lat = c(
    33.2294640,   # WF1 West Fork Gila at Cliff Dwellings
    33.2378021,   # MF Middle Fork / Woody's Corral
    33.1769800,   # EF East Fork above Gila River
    33.1794444    # GM1 Gila River below East Fork
  ),
  lon = c(
    -108.2655590,
    -108.1690382,
    -108.2010000,
    -108.2061111
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
ggplot(GINI_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.01, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "GINI-1 sampling locations"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(GINI_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- GINI_1_coords$site
colnames(geo_dist_km) <- GINI_1_coords$site
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
# 2010 table only
# order: WF1, MF, EF, GM1
# -----------------------------
fst_vals <- c(
  
  # MF
  0.028,
  
  # EF
  0.018, 0.023,
  
  # GM1
  0.017, 0.040, 0.037
)

GINI_1_fst <- fill_sym_from_lower(
  pops = GINI_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

GINI_1_fst <- pmax(GINI_1_fst, 0)
diag(GINI_1_fst) <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(GINI_1_fst)[row(GINI_1_fst)[upper.tri(GINI_1_fst)]],
  site2   = colnames(GINI_1_fst)[col(GINI_1_fst)[upper.tri(GINI_1_fst)]],
  fst     = GINI_1_fst[upper.tri(GINI_1_fst)],
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
    title = "GINI-1 isolation by distance"
  )

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(GINI_1_fst), GINI_1_coords$site))
stopifnot(identical(colnames(GINI_1_fst), GINI_1_coords$site))
stopifnot(isTRUE(all.equal(GINI_1_fst, t(GINI_1_fst))))
stopifnot(all(diag(GINI_1_fst) == 0))
stopifnot(all(diag(geo_dist_km) == 0))

# -----------------------------
# 8) save RData
# -----------------------------
out_dir <- file.path(base_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  GINI_1_fst,
  GINI_1_coords,
  file = file.path(out_dir, "GINI-1.RData")
)