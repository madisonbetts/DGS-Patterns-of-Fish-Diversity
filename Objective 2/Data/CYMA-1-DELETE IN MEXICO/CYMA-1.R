# -----------------------------
# CYMA-1 Desert pupfish
# Cyprinodon macularius
# site coordinates + FST matrix + IBD plot
# Source: Loftis et al. 2009 Conserv Genet 10:453-463
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CYMA-1"

# -----------------------------
# 0) site coordinates
# order follows Table 3 exactly for C. macularius:
# CLD, SFC, SPSS, CP, LS, CST, FDD, ED1, ED2
# coordinates are given in Methods / Fig. 1
# note: locality 1 and 3 longitudes are printed as 106 in the pdf,
# but given the Imperial County / Salton Sea locality this is clearly 116 W
# -----------------------------
CYMA_1_coords <- data.frame(
  site = as.character(1:9),
  code = c("CLD", "SFC", "SPSS", "CP", "LS", "CST", "FDD", "ED1", "ED2"),
  location = c(
    "County-line drain, Salton Sea, CA",
    "San Felipe Creek, Imperial County, CA",
    "Shoreline pool of Salton Sea near Trifolium 20A drain, CA",
    "Cerro Prieto slough, Baja California, MX",
    "Pozo del Tules, Laguna Salada, Baja California, MX",
    "Santa Clara Slough / Canal Sanchez Taboada terminus, Sonora, MX",
    "Flor del Desierto canal, Sonora, MX",
    "El Doctor spring 1, Sonora, MX",
    "El Doctor spring 2, Sonora, MX"
  ),
  lat = c(
    33 + 29/60 + 56/3600,
    33 +  6/60 +  1/3600,
    33 + 29/60 + 56/3600,
    32 + 25/60 + 26/3600,
    32 + 26/60 + 27/3600,
    32 +  3/60 + 28/3600,
    32 +  1/60 + 56/3600,
    31 + 56/60 + 12/3600,
    31 + 56/60 + 47/3600
  ),
  lon = c(
    -(116 +  4/60 + 45/3600),
    -(115 + 55/60 + 57/3600),
    -(116 +  4/60 + 45/3600),
    -(115 + 15/60 + 44/3600),
    -(115 + 36/60 + 15/3600),
    -(114 + 53/60 + 50/3600),
    -(114 + 51/60 + 36/3600),
    -(114 + 43/60 + 56/3600),
    -(114 + 44/60 + 56/3600)
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
ggplot(CYMA_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.10, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "CYMA-1 sampling locations"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(CYMA_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- CYMA_1_coords$site
colnames(geo_dist_km) <- CYMA_1_coords$site
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
# Table 3 lower triangle, row-wise
# species subset = localities 1-9 only
# negative FST values truncated to 0 afterward
# order: CLD, SFC, SPSS, CP, LS, CST, FDD, ED1, ED2
# -----------------------------
fst_vals <- c(

  # SFC
  0.033,

  # SPSS
  0.029, 0.007,

  # CP
  0.020, 0.016, 0.017,

  # LS
  0.050, 0.039, 0.037, 0.028,

  # CST
  0.040, 0.030, 0.030, 0.017, 0.029,

  # FDD
  0.030, 0.019, 0.015, 0.008, 0.022, 0.015,

  # ED1
  0.030, 0.007, 0.013, 0.008, 0.016, 0.014, 0.004,

  # ED2
  0.053, 0.021, 0.025, 0.019, 0.028, 0.027, 0.009, 0.008
)

CYMA_1_fst <- fill_sym_from_lower(
  pops = CYMA_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

CYMA_1_fst <- pmax(CYMA_1_fst, 0)
diag(CYMA_1_fst) <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(CYMA_1_fst)[row(CYMA_1_fst)[upper.tri(CYMA_1_fst)]],
  site2   = colnames(CYMA_1_fst)[col(CYMA_1_fst)[upper.tri(CYMA_1_fst)]],
  fst     = CYMA_1_fst[upper.tri(CYMA_1_fst)],
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
    title = "CYMA-1 isolation by distance"
  )

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(CYMA_1_fst), CYMA_1_coords$site))
stopifnot(identical(colnames(CYMA_1_fst), CYMA_1_coords$site))
stopifnot(isTRUE(all.equal(CYMA_1_fst, t(CYMA_1_fst))))
stopifnot(all(diag(CYMA_1_fst) == 0))
stopifnot(all(diag(geo_dist_km) == 0))

# -----------------------------
# 8) save RData
# -----------------------------
out_dir <- file.path(base_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  CYMA_1_fst,
  CYMA_1_coords,
  file = file.path(out_dir, "CYMA-1.RData")
)
