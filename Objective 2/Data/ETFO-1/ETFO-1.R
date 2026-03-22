# -----------------------------
# ETFO-1 fountain darter
# aggregated collection centroids + pairwise FST matrix + IBD plot
# Source study: Olsen et al. 2016 Conserv Genet 17:1393-1404
# -----------------------------

library(ggplot2)
library(geosphere)
library(sf)
library(tigris)

options(tigris_use_cache = TRUE)

# ------------------------------
# paths
# ------------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data"
save_dir <- file.path(base_dir, "ETFO-1")
out_dir  <- file.path(save_dir, "data")

if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# -----------------------------
# 0) site coordinates
# -----------------------------
# IMPORTANT:
# The paper reports 16 sampling sites, but pairwise FST only for 8 aggregated collections:
# CR1, CR2, CR3, CR4, SMR1, SMR2, SMR3, SMR4.
# To keep the FST matrix aligned with the published data, this script uses one centroid per
# aggregated collection.
#
# Coordinate logic used here:
# 1 = CR1: Landa Lake / Houston St / Liberty Ave / Landa Park centroid
# 2 = CR2: Elizabeth Ave / Old Channel
# 3 = CR3: Hinman Island / New-Old Channel centroid
# 4 = CR4: Garden St
# 5 = SMR1: Spring Lake / outflow dam centroid
# 6 = SMR2: Sewell Park / City Park / Rio Vista Park centroid
# 7 = SMR3: Cheatham St / I-35 centroid
# 8 = SMR4: Lower river / Cypress Tree / Todd Island centroid
#
# These are internet-estimated coordinates intended for Objective 2 IBD screening at the
# published collection scale, not exact capture GPS points.

ETFO_1_sites <- data.frame(
  site = as.character(1:8),
  collection = c("CR1", "CR2", "CR3", "CR4", "SMR1", "SMR2", "SMR3", "SMR4"),
  full_name = c(
    "Landa Lake / Houston St / Liberty Ave / Landa Park",
    "Elizabeth Ave / Old Channel",
    "Hinman Island / New-Old Channel",
    "Garden St",
    "Spring Lake / outflow dam",
    "Sewell Park / City Park / Rio Vista Park",
    "Cheatham St / I-35",
    "Lower river / Cypress Tree / Todd Island"
  ),
  lat = c(
    29.7102244,
    29.7082557,
    29.7087306,
    29.7020060,
    mean(c(29.8910424, 29.8901970)),
    mean(c(29.8880900, 29.8860200, 29.8784561)),
    mean(c(29.8765000, 29.8744444)),
    mean(c(29.8615000, 29.8574306))
  ),
  lon = c(
    -98.1302865,
    -98.1314163,
    mean(c(-98.1334417, -98.1253000)),
    -98.1190680,
    mean(c(-97.9311278, -97.9338450)),
    mean(c(-97.9340900, -97.9358000, -97.9339366)),
    mean(c(-97.9330000, -97.9311111)),
    mean(c(-97.9235000, -97.9206389))
  ),
  stringsAsFactors = FALSE
)

ETFO_1_sites

# -----------------------------
# 1) map of sampling collections
# -----------------------------
texas <- tigris::states(cb = TRUE, year = 2023)
texas <- texas[texas$STUSPS == "TX", ]

bbox_pad <- 0.03
xlim <- range(ETFO_1_sites$lon) + c(-bbox_pad, bbox_pad)
ylim <- range(ETFO_1_sites$lat) + c(-bbox_pad, bbox_pad)

site_points_sf <- st_as_sf(ETFO_1_sites, coords = c("lon", "lat"), crs = 4326)
site_labels_df <- ETFO_1_sites
site_labels_df$label <- paste(site_labels_df$site, site_labels_df$collection, sep = ": ")

ggplot() +
  geom_sf(data = texas, fill = "gray95", color = "gray40", linewidth = 0.4) +
  geom_sf(data = site_points_sf, size = 2.8) +
  geom_text(
    data = site_labels_df,
    aes(x = lon, y = lat, label = label),
    nudge_y = 0.003,
    size = 3.2
  ) +
  coord_sf(xlim = xlim, ylim = ylim, expand = FALSE) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ETFO-1 sampling collections"
  )

# -----------------------------
# 2) geographic distance matrix
# straight-line distance in km
# -----------------------------
geo_dist_km <- geosphere::distm(
  as.matrix(ETFO_1_sites[, c("lon", "lat")]),
  fun = geosphere::distHaversine
) / 1000

rownames(geo_dist_km) <- ETFO_1_sites$site
colnames(geo_dist_km) <- ETFO_1_sites$site

# -----------------------------
# 3) pairwise FST matrix
# Olsen et al. 2016 Table 2, collections in order:
# CR1, CR2, CR3, CR4, SMR1, SMR2, SMR3, SMR4
# Any FST < 0 forced to 0 per workflow
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

fst_vals <- c(
  0.0067,
  0.0029, 0.0028,
  0.0129, 0.0082, 0.0068,
  0.0114, 0.0113, 0.0099, 0.0142,
  0.0168, 0.0145, 0.0126, 0.0133, 0.0004,
  0.0209, 0.0155, 0.0119, 0.0196, 0.0031, 0.0017,
  0.0141, 0.0139, 0.0099, 0.0134, 0.0007, -0.0002, -0.0008
)

ETFO_1_fst <- fill_sym_from_lower(
  pops = ETFO_1_sites$site,
  vals = fst_vals,
  diag_val = 0
)

ETFO_1_fst[ETFO_1_fst < 0] <- 0

# -----------------------------
# 4) pairwise dataframe for IBD plot
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(ETFO_1_fst)[row(ETFO_1_fst)[upper.tri(ETFO_1_fst)]],
  site2   = colnames(ETFO_1_fst)[col(ETFO_1_fst)[upper.tri(ETFO_1_fst)]],
  fst     = ETFO_1_fst[upper.tri(ETFO_1_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# 5) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "ETFO-1 isolation by distance"
  )

# -----------------------------
# 6) quick checks
# -----------------------------
stopifnot(identical(rownames(ETFO_1_fst), ETFO_1_sites$site))
stopifnot(identical(colnames(ETFO_1_fst), ETFO_1_sites$site))
stopifnot(isTRUE(all.equal(ETFO_1_fst, t(ETFO_1_fst))))
stopifnot(all(ETFO_1_fst >= 0))

# -----------------------------
# 7) save RData
# -----------------------------
save(
  ETFO_1_fst,
  ETFO_1_sites,
  file = file.path(out_dir, "ETFO-1.RData")
)
