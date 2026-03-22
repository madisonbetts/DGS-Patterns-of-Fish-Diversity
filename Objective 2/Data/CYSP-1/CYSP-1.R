# -----------------------------
# CYSP-1 Ash Meadows pupfish complex
# Cyprinodon spp. (C. n. mionectes, C. n. pectoralis, C. diabolis)
# Martin & Wilcox 2004 Table 3 + spring coordinates
# -----------------------------

library(ggplot2)
library(geosphere)

# ------------------------------
# paths
# ------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/CYSP-1"

# -----------------------------
# 0) site coordinates
# order follows Table 3 exactly:
# F, R, KP, B, SS, M, NI, S, DH
# coordinates are approximate named spring coordinates from Ash Meadows sources
# -----------------------------
CYSP_1_coords <- data.frame(
  site = as.character(1:9),
  code = c("F", "R", "KP", "B", "SS", "M", "NI", "S", "DH"),
  taxon = c(
    rep("Cyprinodon nevadensis mionectes", 4),
    rep("Cyprinodon nevadensis pectoralis", 4),
    "Cyprinodon diabolis"
  ),
  location = c(
    "Fairbanks Spring, Ash Meadows, NV",
    "Rogers Spring, Ash Meadows, NV",
    "Kings Pool / Point of Rocks Springs, Ash Meadows, NV",
    "Big Spring, Ash Meadows, NV",
    "South Scruggs Spring, Ash Meadows, NV",
    "Marsh Spring, Ash Meadows, NV",
    "North Indian Spring, Ash Meadows, NV",
    "School Spring, Ash Meadows, NV",
    "Devils Hole, Ash Meadows, NV"
  ),
  lat = c(
    36.49051,
    36.47920,
    36.40130,
    36.37461,
    36.43167,
    36.42746,
    36.42551,
    36.42912,
    36.42580
  ),
  lon = c(
    -116.34119,
    -116.32530,
    -116.27280,
    -116.27440,
    -116.31139,
    -116.31504,
    -116.31421,
    -116.30560,
    -116.28700
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 1) map of sampling locations
# -----------------------------
ggplot(CYSP_1_coords, aes(x = lon, y = lat)) +
  geom_point(size = 3) +
  geom_text(aes(label = site), nudge_y = 0.004, size = 3.5) +
  coord_fixed() +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "CYSP-1 sampling locations"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(CYSP_1_coords[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- CYSP_1_coords$site
colnames(geo_dist_km) <- CYSP_1_coords$site
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
# order: F, R, KP, B, SS, M, NI, S, DH
# negative FST values truncated to 0 afterward
# -----------------------------
fst_vals <- c(

  # R
  0.09,

  # KP
  0.03, 0.13,

  # B
  0.00, 0.11, -0.01,

  # SS
  0.15, 0.14, 0.20, 0.18,

  # M
  0.05, 0.08, 0.06, 0.07, 0.11,

  # NI
  0.13, 0.13, 0.16, 0.16, 0.15, 0.09,

  # S
  0.28, 0.26, 0.28, 0.23, 0.19, 0.18, 0.30,

  # DH
  0.49, 0.50, 0.48, 0.37, 0.61, 0.49, 0.46, 0.67
)

CYSP_1_fst <- fill_sym_from_lower(
  pops = CYSP_1_coords$site,
  vals = fst_vals,
  diag_val = 0
)

CYSP_1_fst <- pmax(CYSP_1_fst, 0)
diag(CYSP_1_fst) <- 0

# -----------------------------
# 5) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(CYSP_1_fst)[row(CYSP_1_fst)[upper.tri(CYSP_1_fst)]],
  site2   = colnames(CYSP_1_fst)[col(CYSP_1_fst)[upper.tri(CYSP_1_fst)]],
  fst     = CYSP_1_fst[upper.tri(CYSP_1_fst)],
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
    title = "CYSP-1 isolation by distance"
  )

# -----------------------------
# 7) quick checks
# -----------------------------
stopifnot(identical(rownames(CYSP_1_fst), CYSP_1_coords$site))
stopifnot(identical(colnames(CYSP_1_fst), CYSP_1_coords$site))
stopifnot(isTRUE(all.equal(CYSP_1_fst, t(CYSP_1_fst))))
stopifnot(all(diag(CYSP_1_fst) == 0))
stopifnot(all(diag(geo_dist_km) == 0))

# -----------------------------
# 8) save RData
# -----------------------------
out_dir <- file.path(base_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

save(
  CYSP_1_fst,
  CYSP_1_coords,
  file = file.path(out_dir, "CYSP-1.RData")
)
