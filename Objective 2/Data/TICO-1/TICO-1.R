# -----------------------------
# TICO-1 Tiaroga cobitis (loach minnow)
# sites + pairwise FST + IBD plot
# mutual 2009/2010 site-pairs averaged
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/TICO-1"

# -----------------------------
# 0) site coordinates
# approximate, anchored to named Gila locations
# -----------------------------
TICO_1_coords <- data.frame(
  site = c("1", "2", "3", "4"),
  code = c("GM1", "GM2", "GM3", "GM4"),
  lat = c(
    33.1794444,   # GM1 below East Fork
    33.0761800,   # GM2 above Turkey Creek
    33.0615000,   # GM3 above Mogollon Creek / near Gila
    33.0429800    # GM4 below Mogollon Creek
  ),
  lon = c(
    -108.2061111,
    -108.4880000,
    -108.5370000,
    -108.5282500
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
ggplot(TICO_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.01, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "TICO-1 sampling locations"
  )

# -----------------------------
# 2) geographic distance matrix
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(TICO_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- TICO_1_coords$site
colnames(geo_dist_km) <- TICO_1_coords$site
diag(geo_dist_km) <- 0

# -----------------------------
# 3) helper function
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
# order: GM1, GM2, GM3, GM4
# overlapping 2009/2010 pairs averaged:
# GM1-GM2 = mean(0.024, 0.011) = 0.0175
# GM1-GM3 = mean(0.021, 0.014) = 0.0175
# GM2-GM3 = mean(-0.002, 0.004) = 0.0010
# -----------------------------
fst_vals <- c(
  
  # GM2
  0.0175,
  
  # GM3
  0.0175, 0.0010,
  
  # GM4
  0.0150, 0.0050, 0.0050
)

TICO_1_fst <- fill_sym_from_lower(
  pops = TICO_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

TICO_1_fst <- pmax(TICO_1_fst, 0)
diag(TICO_1_fst) <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(TICO_1_fst)[row(TICO_1_fst)[upper.tri(TICO_1_fst)]],
  site2   = colnames(TICO_1_fst)[col(TICO_1_fst)[upper.tri(TICO_1_fst)]],
  fst     = TICO_1_fst[upper.tri(TICO_1_fst)],
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
    title = "TICO-1 isolation by distance"
  )

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(TICO_1_fst), TICO_1_coords$site))
stopifnot(identical(colnames(TICO_1_fst), TICO_1_coords$site))
stopifnot(isTRUE(all.equal(TICO_1_fst, t(TICO_1_fst))))
stopifnot(all(diag(TICO_1_fst) == 0))
stopifnot(all(diag(geo_dist_km) == 0))

# -----------------------------
# 8) save RData
# -----------------------------
out_dir <- file.path(base_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  TICO_1_fst,
  TICO_1_coords,
  file = file.path(out_dir, "TICO-1.RData")
)