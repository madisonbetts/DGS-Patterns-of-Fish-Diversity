# -----------------------------
# SACO-2 bull trout
# Salvelinus confluentus
# Warnock et al. 2010
# sites + pairwise FST + IBD plot
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SACO-2"

# -----------------------------
# 0) site coordinates
# order follows Table 1 exactly:
# Cb-CR, Ga-CR, Lo-CR, Wca-CR, Sca-CR, Mi-CR,
# Ra-OR, Sra-OR, Nra-OR, Du-OR, Hi-OR, Lli-OR, Uli-OR
#
# The paper does not provide exact sampling coordinates.
# These are approximate stream reference points based on
# named waterbodies / critical habitat points in the same
# drainages, plus the Livingstone Falls split for lower vs upper
# Livingstone.
# -----------------------------
SACO_2_coords <- data.frame(
  site = as.character(1:13),
  code = c(
    "Cb-CR", "Ga-CR", "Lo-CR", "Wca-CR", "Sca-CR", "Mi-CR",
    "Ra-OR", "Sra-OR", "Nra-OR", "Du-OR", "Hi-OR", "Lli-OR", "Uli-OR"
  ),
  location = c(
    "Carbondale River",
    "Gardiner Creek",
    "Lost Creek",
    "West Castle River",
    "South Castle River",
    "Mill Creek",
    "Racehorse Creek",
    "South Racehorse Creek",
    "North Racehorse Creek",
    "Dutch Creek",
    "Hidden Creek",
    "Lower Livingstone River",
    "Upper Livingstone River"
  ),
  lat = c(
    49.38654309020,
    49.37144943270,
    49.39888385160,
    49.20009722000,
    49.22233722222,
    49.24034401990,
    49.83089658140,
    49.75547785090,
    49.83781484500,
    49.96709838730,
    49.97132085210,
    50.10141111111,
    50.18327611111
  ),
  lon = c(
    -114.57512465300,
    -114.52875620700,
    -114.58250840600,
    -114.34369711500,
    -114.22821111111,
    -114.16792365800,
    -114.47506394400,
    -114.62825033600,
    -114.64854979400,
    -114.67267947700,
    -114.64932796100,
    -114.44437222222,
    -114.47625972222
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
ggplot(SACO_2_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.03, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "SACO-2 sampling locations"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(SACO_2_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- SACO_2_coords$site
colnames(geo_dist_km) <- SACO_2_coords$site
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
# Table 1 lower triangle, transcribed row-wise
# negative FST values truncated to 0 afterward
# -----------------------------
fst_vals <- c(

  # Ga-CR
  0.0505,

  # Lo-CR
  0.0497, 0.0742,

  # Wca-CR
  0.1140, 0.0787, 0.1051,

  # Sca-CR
  0.0613, 0.0484, 0.0452, 0.0672,

  # Mi-CR
  0.1455, 0.1327, 0.0688, 0.1686, 0.0622,

  # Ra-OR
  0.1962, 0.1873, 0.1627, 0.2878, 0.1646, 0.1907,

  # Sra-OR
  0.1869, 0.1891, 0.1645, 0.2827, 0.1612, 0.1930, -0.0080,

  # Nra-OR
  0.2398, 0.2359, 0.1883, 0.3057, 0.1907, 0.2109, 0.0559, 0.0601,

  # Du-OR
  0.0998, 0.1357, 0.0561, 0.1729, 0.1003, 0.1349, 0.1159, 0.1146, 0.1502,

  # Hi-OR
  0.1232, 0.1576, 0.0842, 0.2181, 0.1304, 0.1612, 0.0909, 0.1018, 0.1527, 0.0344,

  # Lli-OR
  0.1674, 0.1698, 0.1145, 0.2053, 0.1685, 0.2284, 0.2625, 0.2750, 0.3174, 0.1590, 0.1547,

  # Uli-OR
  0.1775, 0.1992, 0.1142, 0.2320, 0.1898, 0.2396, 0.2863, 0.2968, 0.3129, 0.1577, 0.1483, 0.0419
)

SACO_2_fst <- fill_sym_from_lower(
  pops = SACO_2_coords$site,
  vals = fst_vals,
  diag_val = 0
)

SACO_2_fst <- pmax(SACO_2_fst, 0)
diag(SACO_2_fst) <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(SACO_2_fst)[row(SACO_2_fst)[upper.tri(SACO_2_fst)]],
  site2   = colnames(SACO_2_fst)[col(SACO_2_fst)[upper.tri(SACO_2_fst)]],
  fst     = SACO_2_fst[upper.tri(SACO_2_fst)],
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
    title = "SACO-2 isolation by distance"
  )

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(SACO_2_fst), SACO_2_coords$site))
stopifnot(identical(colnames(SACO_2_fst), SACO_2_coords$site))
stopifnot(isTRUE(all.equal(SACO_2_fst, t(SACO_2_fst))))
stopifnot(all(diag(SACO_2_fst) == 0))
stopifnot(all(diag(geo_dist_km) == 0))

# -----------------------------
# 8) save RData
# -----------------------------
out_dir <- file.path(base_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  SACO_2_fst,
  SACO_2_coords,
  file = file.path(out_dir, "SACO-2.RData")
)
