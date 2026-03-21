# -----------------------------
# MEFU-1 Meda fulgida (spikedace)
# sites + pairwise FST + IBD plot
# mutual 2009/2010 site-pairs averaged
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/MEFU-1"

# -----------------------------
# 0) site coordinates
# approximate, anchored to named Gila locations
# -----------------------------
MEFU_1_coords <- data.frame(
  site = c("1", "2", "3", "4", "5", "6"),
  code = c("WF1", "MF", "WF2", "GM2", "GM3", "GM4"),
  lat = c(
    33.2294640,   # WF1
    33.2378021,   # MF
    33.2218611,   # WF2
    33.0761800,   # GM2 above Turkey Creek
    33.0615000,   # GM3 above Mogollon Creek / near Gila
    33.0429800    # GM4 below Mogollon Creek
  ),
  lon = c(
    -108.2655590,
    -108.1690382,
    -108.2420000,
    -108.4880000,
    -108.5370000,
    -108.5282500
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
ggplot(MEFU_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.01, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "MEFU-1 sampling locations"
  )

# -----------------------------
# 2) geographic distance matrix
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(MEFU_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- MEFU_1_coords$site
colnames(geo_dist_km) <- MEFU_1_coords$site
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
# order: WF1, MF, WF2, GM2, GM3, GM4
# overlapping 2009/2010 pairs averaged:
# WF2-GM2 = mean(0.031, 0.036) = 0.0335
# WF2-GM3 = mean(0.037, 0.027) = 0.0320
# GM2-GM3 = mean(-0.002, 0.007) = 0.0025
# -----------------------------
fst_vals <- c(
  
  # MF
  0.006,
  
  # WF2
  0.013, 0.000,
  
  # GM2
  0.042, 0.033, 0.0335,
  
  # GM3
  0.032, 0.027, 0.0320, 0.0025,
  
  # GM4
  0.036, 0.032, 0.0320, 0.0050, 0.000
)

MEFU_1_fst <- fill_sym_from_lower(
  pops = MEFU_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

MEFU_1_fst <- pmax(MEFU_1_fst, 0)
diag(MEFU_1_fst) <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(MEFU_1_fst)[row(MEFU_1_fst)[upper.tri(MEFU_1_fst)]],
  site2   = colnames(MEFU_1_fst)[col(MEFU_1_fst)[upper.tri(MEFU_1_fst)]],
  fst     = MEFU_1_fst[upper.tri(MEFU_1_fst)],
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
    title = "MEFU-1 isolation by distance"
  )

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(MEFU_1_fst), MEFU_1_coords$site))
stopifnot(identical(colnames(MEFU_1_fst), MEFU_1_coords$site))
stopifnot(isTRUE(all.equal(MEFU_1_fst, t(MEFU_1_fst))))
stopifnot(all(diag(MEFU_1_fst) == 0))
stopifnot(all(diag(geo_dist_km) == 0))

# -----------------------------
# 8) save RData
# -----------------------------
out_dir <- file.path(base_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  MEFU_1_fst,
  MEFU_1_coords,
  file = file.path(out_dir, "MEFU-1.RData")
)