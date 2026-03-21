# -----------------------------
# SAVI-1 Walleye
# Sander vitreus
# sites + pairwise FST + IBD plot
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAVI-1"

# -----------------------------
# 0) site coordinates
# approximate from Fig. 1 + named localities
# order follows Table 3 exactly:
# BDN, BR, CR, ML, MB, PR, WR
# -----------------------------
SAVI_1_coords <- data.frame(
  site = as.character(1:7),
  code = c("BDN", "BR", "CR", "ML", "MB", "PR", "WR"),
  location = c(
    "Bay de Noc, MI",
    "Bar River, ON",
    "Charlotte River, MI",
    "Manistique Lake, MI",
    "Munuscong Bay, MI",
    "Potogannissing River, MI",
    "Waiska River, MI"
  ),
  lat = c(
    45.72,   # Bay de Noc
    46.32,   # Bar River
    46.26,   # Charlotte River
    46.10,   # Manistique Lake
    46.18,   # Munuscong Bay
    45.92,   # Potogannissing River
    46.44    # Waiska River
  ),
  lon = c(
    -86.90,  # Bay de Noc
    -83.95,  # Bar River / Lake George area
    -84.16,  # Charlotte River
    -85.95,  # Manistique Lake
    -84.02,  # Munuscong Bay
    -83.85,  # Potogannissing River / Drummond Island area
    -84.37   # Waiska River
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
ggplot(SAVI_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.08, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "SAVI-1 sampling locations"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(SAVI_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- SAVI_1_coords$site
colnames(geo_dist_km) <- SAVI_1_coords$site
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
# Table 3 lower triangle, transcribed row-wise
# order: BDN, BR, CR, ML, MB, PR, WR
# -----------------------------
fst_vals <- c(
  
  # BR
  0.1009,
  
  # CR
  0.1052, 0.0629,
  
  # ML
  0.0311, 0.0886, 0.1018,
  
  # MB
  0.0410, 0.0399, 0.0469, 0.0327,
  
  # PR
  0.0167, 0.0868, 0.0728, 0.0245, 0.0157,
  
  # WR
  0.0273, 0.0676, 0.0661, 0.0354, 0.0051, 0.0014
)

SAVI_1_fst <- fill_sym_from_lower(
  pops = SAVI_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

SAVI_1_fst <- pmax(SAVI_1_fst, 0)
diag(SAVI_1_fst) <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(SAVI_1_fst)[row(SAVI_1_fst)[upper.tri(SAVI_1_fst)]],
  site2   = colnames(SAVI_1_fst)[col(SAVI_1_fst)[upper.tri(SAVI_1_fst)]],
  fst     = SAVI_1_fst[upper.tri(SAVI_1_fst)],
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
    title = "SAVI-1 isolation by distance"
  )

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(SAVI_1_fst), SAVI_1_coords$site))
stopifnot(identical(colnames(SAVI_1_fst), SAVI_1_coords$site))
stopifnot(isTRUE(all.equal(SAVI_1_fst, t(SAVI_1_fst))))
stopifnot(all(diag(SAVI_1_fst) == 0))
stopifnot(all(diag(geo_dist_km) == 0))

# -----------------------------
# 8) save RData
# -----------------------------
out_dir <- file.path(base_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  SAVI_1_fst,
  SAVI_1_coords,
  file = file.path(out_dir, "SAVI-1.RData")
)