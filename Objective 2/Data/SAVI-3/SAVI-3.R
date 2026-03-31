#-------
# SAVI-3
#-------

# packages
library(dplyr)
library(ggplot2)
library(geosphere)
library(sf)
library(rnaturalearth)

# directories
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SAVI-3"
data_dir <- file.path(base_dir, "data")
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# 1) site coordinates
# -----------------------------
SAVI_3_coords <- data.frame(
  site = as.character(1:29),
  lat = c(
    46.7080, # 1 S_SL
    46.5030, # 2 S_LG
    46.6870, # 3 S_KA
    48.5630, # 4 S_NR
    48.7070, # 5 S_BS
    46.4860, # 6 S_SM
    44.3920, # 7 M_WR
    44.2610, # 8 M_FR
    45.9300, # 9 M_LB
    43.4230, # 10 M_MU
    43.6100, # 11 H_TR
    44.8550, # 12 H_MR
    42.5810, # 13 stC_CR
    42.2520, # 14 E_DE
    41.7070, # 15 E_MA
    41.3540, # 16 E_SA
    41.8260, # 17 E_RE
    41.9280, # 18 E_CI
    41.7560, # 19 E_OH
    41.8150, # 20 E_ZR
    42.3260, # 21 E_VB
    42.4460, # 22 E_SH
    42.5440, # 23 E_CC
    42.3680, # 24 E_BB
    42.8520, # 25 E_SC
    42.9000, # 26 E_ON
    44.1640, # 27 O_BQ
    44.0390, # 28 O_BR
    43.2300  # 29 O_OL
  ),
  lon = c(
    -92.1980, # 1 S_SL
    -89.5940, # 2 S_LG
    -90.8240, # 3 S_KA
    -88.2580, # 4 S_NR
    -88.4680, # 5 S_BS
    -84.3450, # 6 S_SM
    -88.7390, # 7 M_WR
    -88.0600, # 8 M_FR
    -86.9900, # 9 M_LB
    -85.6530, # 10 M_MU
    -84.2470, # 11 H_TR
    -79.8130, # 12 H_MR
    -82.8020, # 13 stC_CR
    -83.1690, # 14 E_DE
    -83.4850, # 15 E_MA
    -83.1060, # 16 E_SA
    -82.8820, # 17 E_RE
    -83.1370, # 18 E_CI
    -81.2760, # 19 E_OH
    -82.7770, # 20 E_ZR
    -79.5750, # 21 E_VB
    -79.5050, # 22 E_SH
    -79.0370, # 23 E_CC
    -78.9480, # 24 E_BB
    -78.8230, # 25 E_SC
    -79.6210, # 26 E_ON
    -77.3830, # 27 O_BQ
    -76.1140, # 28 O_BR
    -76.0050  # 29 O_OL
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) FST values
# lower triangle row-wise in figure order:
# O_BQ O_BR O_OL E_ON E_SC E_BB E_CC E_SH E_VB E_ZR E_OH E_CI
# E_RE E_SA E_MA E_DE stC_CR H_MR H_TR M_MU M_LB M_FR M_WR
# S_SM S_BS S_NR S_KA S_LG S_SL
# -----------------------------
fst_vals <- c(
  0.03,
  0.02, 0.02,
  0.02, 0.02, 0.02,
  0.04, 0.02, 0.02, 0.02,
  0.02, 0.02, 0.02, 0.05, 0.04,
  0.05, 0.06, 0.06, 0.06, 0.02, 0.03,
  0.07, 0.05, 0.06, 0.07, 0.06, 0.04, 0.02,
  0.04, 0.03, 0.03, 0.03, 0.03, 0.03, 0.05, 0.03,
  0.03, 0.03, 0.03, 0.03, 0.03, 0.05, 0.04, 0.06, 0.06,
  0.07, 0.07, 0.03, 0.08, 0.06, 0.07, 0.07, 0.07, 0.05, 0.03,
  0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.06, 0.04, 0.04, 0.04, 0.04,
  0.04, 0.04, 0.06, 0.05, 0.06, 0.07, 0.08, 0.07, 0.09, 0.06, 0.08, 0.08,
  0.07, 0.05, 0.03, 0.04, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03,
  0.03, 0.03, 0.05, 0.04, 0.06, 0.06, 0.06, 0.06, 0.08, 0.06, 0.07, 0.07, 0.07, 0.04,
  0.03, 0.00, 0.00, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.04, 0.03, 0.04,
  0.05, 0.06, 0.05, 0.07, 0.04, 0.06, 0.06, 0.05, 0.03, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01,
  0.01, 0.00, 0.00, 0.00, 0.04, 0.03, 0.05, 0.06, 0.06, 0.06, 0.08, 0.05, 0.07, 0.07, 0.06, 0.03, 0.01,
  0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.04, 0.03, 0.04, 0.05, 0.06, 0.05, 0.07, 0.04, 0.06,
  0.06, 0.05, 0.03, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.04, 0.03, 0.04, 0.05, 0.05, 0.05, 0.07,
  0.04, 0.06, 0.06, 0.05, 0.03, 0.01, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.04, 0.03, 0.04, 0.05, 0.05, 0.05, 0.06,
  0.04, 0.06, 0.06, 0.05, 0.03, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.04, 0.03, 0.04, 0.05, 0.05, 0.05, 0.06, 0.04, 0.06,
  0.06, 0.05, 0.03, 0.01, 0.00, 0.00, 0.00, 0.00, 0.00, 0.04, 0.03, 0.04, 0.05, 0.06, 0.05, 0.07, 0.04, 0.06, 0.06, 0.05, 0.03, 0.00,
  0.00, 0.00, 0.00, 0.00, 0.04, 0.03, 0.04, 0.05, 0.05, 0.05, 0.06, 0.04, 0.06, 0.06, 0.05, 0.03, 0.00, 0.00, 0.00, 0.00, 0.04, 0.03, 0.04,
  0.05, 0.05, 0.05, 0.06, 0.04, 0.06, 0.06, 0.05, 0.03, 0.00, 0.00, 0.00, 0.03, 0.02, 0.04, 0.04, 0.05, 0.05, 0.06, 0.04, 0.06, 0.06, 0.05, 0.02,
  0.00, 0.00, 0.04, 0.02, 0.04, 0.04, 0.05, 0.05, 0.06, 0.04, 0.06, 0.06, 0.05, 0.02, 0.00, 0.03, 0.02, 0.04, 0.04, 0.05, 0.05, 0.06, 0.04, 0.06, 0.06,
  0.05, 0.02, 0.00, 0.03, 0.01, 0.03, 0.04, 0.03, 0.04, 0.06, 0.04, 0.05, 0.05, 0.04, 0.02, 0.03, 0.04, 0.04, 0.04, 0.06, 0.07, 0.04, 0.07, 0.05, 0.05, 0.01,
  0.04, 0.04, 0.01, 0.05, 0.06, 0.04, 0.06, 0.05, 0.05, 0.02, 0.05, 0.05, 0.06, 0.07, 0.06, 0.08, 0.06, 0.06, 0.03, 0.03, 0.05, 0.06, 0.04, 0.06, 0.05, 0.05, 0.02,
  0.00, 0.05, 0.02, 0.03, 0.04, 0.03, 0.03, 0.06, 0.03, 0.03, 0.05, 0.04, 0.04, 0.05, 0.03, 0.05, 0.04, 0.03, 0.04, 0.07, 0.05, 0.04, 0.03, 0.06, 0.03, 0.03, 0.00, 0.03
)

stopifnot(length(fst_vals) == 406)

# -----------------------------
# 3) helper: build symmetric matrix
# -----------------------------
fill_sym_from_lower <- function(pops, vals, diag_val = 0) {
  n <- length(pops)
  expected <- n * (n - 1) / 2
  
  if (length(vals) != expected) {
    stop(
      paste0(
        "fst_vals has length ", length(vals),
        " but should be ", expected,
        " for ", n, " sites."
      )
    )
  }
  
  mat <- matrix(NA_real_, n, n)
  rownames(mat) <- pops
  colnames(mat) <- pops
  
  mat[lower.tri(mat)] <- vals
  mat <- t(mat)
  mat[lower.tri(mat)] <- vals
  diag(mat) <- diag_val
  mat[mat < 0] <- 0
  
  mat
}

# -----------------------------
# 4) build FST matrix
# -----------------------------
SAVI_3_fst <- fill_sym_from_lower(
  pops = SAVI_3_coords$site,
  vals = fst_vals,
  diag_val = 0
)

# -----------------------------
# 5) geographic distances
# -----------------------------
coords_mat <- as.matrix(SAVI_3_coords[, c("lon", "lat")])
geo_dist <- distm(coords_mat, fun = distHaversine) / 1000
rownames(geo_dist) <- SAVI_3_coords$site
colnames(geo_dist) <- SAVI_3_coords$site
diag(geo_dist) <- 0

# -----------------------------
# 6) flatten for IBD
# -----------------------------
get_lower <- function(mat) mat[lower.tri(mat)]

ibd_df <- data.frame(
  dist_km = get_lower(geo_dist),
  fst = get_lower(SAVI_3_fst)
)

# -----------------------------
# 7) basemap
# -----------------------------
world <- ne_countries(scale = "medium", returnclass = "sf")

bbox <- st_bbox(c(
  xmin = -93, xmax = -75,
  ymin = 41, ymax = 49
), crs = st_crs(4326))

world_crop <- st_crop(world, bbox)

# -----------------------------
# 8) QC map
# -----------------------------
ggplot() +
  geom_sf(data = world_crop, fill = "grey95", color = "grey70") +
  geom_point(
    data = SAVI_3_coords,
    aes(x = lon, y = lat),
    size = 2
  ) +
  geom_text(
    data = SAVI_3_coords,
    aes(x = lon, y = lat, label = site),
    size = 3,
    vjust = -0.8
  ) +
  coord_sf(xlim = c(-93, -75), ylim = c(41, 49)) +
  theme_classic(base_size = 12) +
  labs(
    x = "Longitude",
    y = "Latitude"
  )

# -----------------------------
# 9) IBD plot
# -----------------------------
ggplot(ibd_df, aes(x = log(dist_km), y = fst)) +
  geom_point(size = 1.5) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 12) +
  labs(
    x = "Geographic distance (km)",
    y = "FST"
  )


# -----------------------------
# 10) save outputs
# -----------------------------
save(
  SAVI_3_fst,
  SAVI_3_coords,
  file = file.path(data_dir, "SAVI-3.RData")
)


# -----------------------------
# 11) checks
# -----------------------------
stopifnot(nrow(SAVI_3_fst) == 29)
stopifnot(ncol(SAVI_3_fst) == 29)
stopifnot(identical(rownames(SAVI_3_fst), colnames(SAVI_3_fst)))
stopifnot(identical(rownames(SAVI_3_fst), SAVI_3_coords$site))
stopifnot(isTRUE(all.equal(SAVI_3_fst, t(SAVI_3_fst))))
stopifnot(all(diag(SAVI_3_fst) == 0))
stopifnot(all(diag(geo_dist) == 0))