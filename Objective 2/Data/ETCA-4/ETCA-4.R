############################
## ETCA-4
## Rainbow darter (Etheostoma caeruleum)
## Davis et al. 2015 Conserv Genet 16:167–179
## DOI: 10.1007/s10592-014-0649-1
##
## Inputs taken from paper:
## - Table 1: locality coordinates
## - Table 2: pairwise FST among 14 localities
##
## Workflow:
## 1) build coordinate dataframe in matrix order
## 2) build symmetrical FST matrix
## 3) set negative FST to 0
## 4) make site map
## 5) compute Euclidean geographic distances
## 6) make IBD plot
## 7) save only ETCA_4_fst and ETCA_4_coords to .RData
############################

# -----------------------------
# libraries
# -----------------------------
library(ggplot2)
library(maps)
library(geosphere)

# -----------------------------
# paths
# -----------------------------
base_dir <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/ETCA-4"
data_dir <- file.path(base_dir, "data")
dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# 1) coordinate dataframe
# matrix order follows Table 2:
# UI1 UI2 UI3 UI4 TK1 TK2 TK3 TK4 MQ1 MQ2 MQ3 CD1 CD2 CD3
# -----------------------------
ETCA_4_coords <- data.frame(
  site_id = 1:14,
  lat = c(
    43.421666, 43.489300, 43.402800, 43.490790,
    43.368110, 43.337150, 43.104022, 42.879481,
    42.658610, 42.547870, 42.619910,
    43.469588, 43.223680, 42.142530
  ),
  lon = c(
    -91.553500, -92.403333, -91.915550, -92.074530,
    -92.216790, -92.119710, -91.915420, -91.763618,
    -91.587500, -91.510150, -91.666950,
    -92.960773, -92.876650, -91.672300
  ),
  stringsAsFactors = FALSE
)

# -----------------------------
# 2) pairwise FST matrix from Table 2
# values below diagonal in paper
# negatives retained initially, then set to 0 below
# -----------------------------
n_sites <- nrow(ETCA_4_coords)

ETCA_4_fst <- matrix(NA_real_, n_sites, n_sites)
diag(ETCA_4_fst) <- 0

# fill lower triangle from Table 2
ETCA_4_fst[2,1]   <-  0.018
ETCA_4_fst[3,1]   <-  0.012
ETCA_4_fst[3,2]   <-  0.014
ETCA_4_fst[4,1]   <- -0.001
ETCA_4_fst[4,2]   <-  0.025
ETCA_4_fst[4,3]   <- -0.004
ETCA_4_fst[5,1]   <-  0.016
ETCA_4_fst[5,2]   <-  0.037
ETCA_4_fst[5,3]   <-  0.017
ETCA_4_fst[5,4]   <-  0.016
ETCA_4_fst[6,1]   <-  0.007
ETCA_4_fst[6,2]   <-  0.029
ETCA_4_fst[6,3]   <-  0.004
ETCA_4_fst[6,4]   <-  0.003
ETCA_4_fst[6,5]   <-  0.013
ETCA_4_fst[7,1]   <-  0.011
ETCA_4_fst[7,2]   <-  0.031
ETCA_4_fst[7,3]   <-  0.007
ETCA_4_fst[7,4]   <-  0.003
ETCA_4_fst[7,5]   <-  0.016
ETCA_4_fst[7,6]   <-  0.001
ETCA_4_fst[8,1]   <-  0.013
ETCA_4_fst[8,2]   <-  0.021
ETCA_4_fst[8,3]   <-  0.003
ETCA_4_fst[8,4]   <-  0.003
ETCA_4_fst[8,5]   <-  0.016
ETCA_4_fst[8,6]   <-  0.000
ETCA_4_fst[8,7]   <-  0.011
ETCA_4_fst[9,1]   <-  0.017
ETCA_4_fst[9,2]   <-  0.014
ETCA_4_fst[9,3]   <-  0.024
ETCA_4_fst[9,4]   <-  0.010
ETCA_4_fst[9,5]   <-  0.023
ETCA_4_fst[9,6]   <-  0.014
ETCA_4_fst[9,7]   <-  0.029
ETCA_4_fst[9,8]   <-  0.016
ETCA_4_fst[10,1]  <-  0.018
ETCA_4_fst[10,2]  <-  0.033
ETCA_4_fst[10,3]  <-  0.023
ETCA_4_fst[10,4]  <-  0.013
ETCA_4_fst[10,5]  <-  0.031
ETCA_4_fst[10,6]  <-  0.014
ETCA_4_fst[10,7]  <-  0.013
ETCA_4_fst[10,8]  <-  0.015
ETCA_4_fst[10,9]  <-  0.001
ETCA_4_fst[11,1]  <-  0.013
ETCA_4_fst[11,2]  <-  0.037
ETCA_4_fst[11,3]  <-  0.019
ETCA_4_fst[11,4]  <-  0.012
ETCA_4_fst[11,5]  <-  0.018
ETCA_4_fst[11,6]  <-  0.005
ETCA_4_fst[11,7]  <-  0.014
ETCA_4_fst[11,8]  <-  0.013
ETCA_4_fst[11,9]  <-  0.002
ETCA_4_fst[11,10] <- -0.003
ETCA_4_fst[12,1]  <-  0.005
ETCA_4_fst[12,2]  <-  0.029
ETCA_4_fst[12,3]  <-  0.010
ETCA_4_fst[12,4]  <-  0.007
ETCA_4_fst[12,5]  <-  0.026
ETCA_4_fst[12,6]  <-  0.011
ETCA_4_fst[12,7]  <-  0.002
ETCA_4_fst[12,8]  <-  0.008
ETCA_4_fst[12,9]  <-  0.020
ETCA_4_fst[12,10] <-  0.013
ETCA_4_fst[12,11] <-  0.024
ETCA_4_fst[13,1]  <-  0.007
ETCA_4_fst[13,2]  <-  0.020
ETCA_4_fst[13,3]  <-  0.011
ETCA_4_fst[13,4]  <-  0.012
ETCA_4_fst[13,5]  <-  0.026
ETCA_4_fst[13,6]  <-  0.014
ETCA_4_fst[13,7]  <-  0.006
ETCA_4_fst[13,8]  <-  0.013
ETCA_4_fst[13,9]  <-  0.018
ETCA_4_fst[13,10] <-  0.012
ETCA_4_fst[13,11] <-  0.026
ETCA_4_fst[13,12] <- -0.005
ETCA_4_fst[14,1]  <-  0.015
ETCA_4_fst[14,2]  <-  0.021
ETCA_4_fst[14,3]  <-  0.013
ETCA_4_fst[14,4]  <-  0.010
ETCA_4_fst[14,5]  <-  0.030
ETCA_4_fst[14,6]  <-  0.018
ETCA_4_fst[14,7]  <-  0.012
ETCA_4_fst[14,8]  <-  0.012
ETCA_4_fst[14,9]  <-  0.013
ETCA_4_fst[14,10] <-  0.015
ETCA_4_fst[14,11] <-  0.029
ETCA_4_fst[14,12] <-  0.004
ETCA_4_fst[14,13] <- -0.005

# make symmetrical
ETCA_4_fst[upper.tri(ETCA_4_fst)] <- t(ETCA_4_fst)[upper.tri(ETCA_4_fst)]

# set negatives to 0 per workflow
ETCA_4_fst[ETCA_4_fst < 0] <- 0

# use numeric site IDs as row/col names
rownames(ETCA_4_fst) <- ETCA_4_coords$site_id
colnames(ETCA_4_fst) <- ETCA_4_coords$site_id

# -----------------------------
# 3) site map
# -----------------------------
states_map <- map_data("state")

ggplot() +
  geom_polygon(
    data = states_map,
    aes(x = long, y = lat, group = group),
    fill = "grey95",
    color = "grey70",
    linewidth = 0.2
  ) +
  geom_point(
    data = ETCA_4_coords,
    aes(x = lon, y = lat),
    size = 3
  ) +
  geom_text(
    data = ETCA_4_coords,
    aes(x = lon, y = lat, label = site_id),
    nudge_y = 0.08,
    size = 3
  ) +
  coord_fixed(
    1.2,
    xlim = c(-94.5, -89.5),
    ylim = c(41.8, 44.0)
  ) +
  theme_classic(base_size = 12) +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "ETCA-4 sampling sites"
  )


# -----------------------------
# 4) geographic distance matrix
# -----------------------------
etca_4_geo_dist_m <- geosphere::distm(
  x = ETCA_4_coords[, c("lon", "lat")],
  fun = geosphere::distHaversine
)

etca_4_geo_dist_km <- etca_4_geo_dist_m / 1000
rownames(etca_4_geo_dist_km) <- ETCA_4_coords$site_id
colnames(etca_4_geo_dist_km) <- ETCA_4_coords$site_id

# -----------------------------
# 5) pairwise dataframe for IBD
# -----------------------------
etca_4_ibd_df <- data.frame(
  site1   = rownames(ETCA_4_fst)[row(ETCA_4_fst)[upper.tri(ETCA_4_fst)]],
  site2   = colnames(ETCA_4_fst)[col(ETCA_4_fst)[upper.tri(ETCA_4_fst)]],
  fst     = ETCA_4_fst[upper.tri(ETCA_4_fst)],
  dist_km = etca_4_geo_dist_km[upper.tri(etca_4_geo_dist_km)]
)

# -----------------------------
# 6) IBD plot
# -----------------------------
ggplot(etca_4_ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic(base_size = 12) +
  labs(
    x = "Euclidean distance (km)",
    y = expression(F[ST]),
    title = "ETCA-4 isolation by distance"
  )


# -----------------------------
# 7) save outputs
# save only fst and coords as requested
# -----------------------------
save(
  ETCA_4_fst,
  ETCA_4_coords,
  file = file.path(data_dir, "ETCA-4.RData")
)
