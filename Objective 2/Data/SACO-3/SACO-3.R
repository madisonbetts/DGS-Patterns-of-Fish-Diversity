
# -----------------------------
# SACO-3 Bull trout (Salvelinus confluentus)
# Whiteley et al. 2006
# Boise River basin, Idaho
#
# Cleaned workflow:
# - site coordinates are read directly from SACO-3-sites.csv
# - site 21 is removed because its pairwise FST values are unavailable
# - the FST matrix is built for the 20 fully observed sites
# - the matrix is forced to be symmetric
# - negative FST values are set to 0
#
# -----------------------------

library(geosphere)
library(ggplot2)
library(maps)

# -----------------------------
# paths
# -----------------------------
site_csv <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SACO-3/SACO-3-sites.csv"
out_rdata <- "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SACO-3/data/SACO-3.RData"

# -----------------------------
# 1) sites from csv
# -----------------------------
SACO_3_coords <- read.csv(site_csv, stringsAsFactors = FALSE)

names(SACO_3_coords)[names(SACO_3_coords) == "latitude"]  <- "lat"
names(SACO_3_coords)[names(SACO_3_coords) == "longitude"] <- "lon"

SACO_3_coords <- SACO_3_coords[order(SACO_3_coords$site), c("site", "stream", "lat", "lon")]
SACO_3_coords$site <- as.integer(SACO_3_coords$site)

# remove site 21 because pairwise FST values were not available in Table 4
SACO_3_coords <- SACO_3_coords[SACO_3_coords$site != 21, ]

stopifnot(nrow(SACO_3_coords) == 20)
stopifnot(identical(SACO_3_coords$site, 1:20))

# -----------------------------
# 2) pairwise FST matrix
# Table 4, above diagonal where available
# -----------------------------
SACO_3_fst <- matrix(NA_real_, nrow = 20, ncol = 20)
rownames(SACO_3_fst) <- colnames(SACO_3_fst) <- as.character(1:20)
diag(SACO_3_fst) <- 0

# row 1
SACO_3_fst[1,2]  <- 0.125
SACO_3_fst[1,3]  <- 0.179
SACO_3_fst[1,4]  <- 0.077
SACO_3_fst[1,5]  <- 0.065
SACO_3_fst[1,6]  <- 0.107
SACO_3_fst[1,7]  <- 0.071
SACO_3_fst[1,8]  <- 0.074
SACO_3_fst[1,9]  <- 0.084
SACO_3_fst[1,10] <- 0.052
SACO_3_fst[1,11] <- 0.098
SACO_3_fst[1,12] <- 0.134
SACO_3_fst[1,13] <- 0.101
SACO_3_fst[1,14] <- 0.175
SACO_3_fst[1,15] <- 0.096
SACO_3_fst[1,16] <- 0.077
SACO_3_fst[1,17] <- 0.201
SACO_3_fst[1,18] <- 0.052
SACO_3_fst[1,19] <- 0.091
SACO_3_fst[1,20] <- 0.120

# row 2
SACO_3_fst[2,3]  <- 0.013
SACO_3_fst[2,4]  <- 0.079
SACO_3_fst[2,5]  <- 0.060
SACO_3_fst[2,6]  <- 0.109
SACO_3_fst[2,7]  <- 0.004
SACO_3_fst[2,8]  <- 0.013
SACO_3_fst[2,9]  <- -0.012
SACO_3_fst[2,10] <- 0.006
SACO_3_fst[2,11] <- 0.015
SACO_3_fst[2,12] <- -0.006
SACO_3_fst[2,13] <- 0.052
SACO_3_fst[2,14] <- 0.047
SACO_3_fst[2,15] <- 0.152
SACO_3_fst[2,16] <- 0.032
SACO_3_fst[2,17] <- 0.205
SACO_3_fst[2,18] <- 0.053
SACO_3_fst[2,19] <- 0.073
SACO_3_fst[2,20] <- 0.010

# row 3
SACO_3_fst[3,4]  <- 0.164
SACO_3_fst[3,5]  <- 0.093
SACO_3_fst[3,6]  <- 0.210
SACO_3_fst[3,7]  <- 0.056
SACO_3_fst[3,8]  <- 0.064
SACO_3_fst[3,9]  <- 0.017
SACO_3_fst[3,10] <- 0.048
SACO_3_fst[3,11] <- 0.024
SACO_3_fst[3,12] <- 0.020
SACO_3_fst[3,13] <- 0.142
SACO_3_fst[3,14] <- 0.062
SACO_3_fst[3,15] <- 0.186
SACO_3_fst[3,16] <- 0.103
SACO_3_fst[3,17] <- 0.254
SACO_3_fst[3,18] <- 0.097
SACO_3_fst[3,19] <- 0.121
SACO_3_fst[3,20] <- 0.043

# row 4
SACO_3_fst[4,5]  <- 0.105
SACO_3_fst[4,6]  <- 0.028
SACO_3_fst[4,7]  <- 0.043
SACO_3_fst[4,8]  <- 0.045
SACO_3_fst[4,9]  <- 0.061
SACO_3_fst[4,10] <- 0.045
SACO_3_fst[4,11] <- 0.072
SACO_3_fst[4,12] <- 0.069
SACO_3_fst[4,13] <- 0.021
SACO_3_fst[4,14] <- 0.141
SACO_3_fst[4,15] <- 0.123
SACO_3_fst[4,16] <- 0.046
SACO_3_fst[4,17] <- 0.140
SACO_3_fst[4,18] <- 0.078
SACO_3_fst[4,19] <- 0.068
SACO_3_fst[4,20] <- 0.098

# row 5
SACO_3_fst[5,6]  <- 0.160
SACO_3_fst[5,7]  <- 0.073
SACO_3_fst[5,8]  <- 0.081
SACO_3_fst[5,9]  <- 0.063
SACO_3_fst[5,10] <- 0.059
SACO_3_fst[5,11] <- 0.069
SACO_3_fst[5,12] <- 0.096
SACO_3_fst[5,13] <- 0.105
SACO_3_fst[5,14] <- 0.179
SACO_3_fst[5,15] <- 0.193
SACO_3_fst[5,16] <- 0.083
SACO_3_fst[5,17] <- 0.310
SACO_3_fst[5,18] <- 0.059
SACO_3_fst[5,19] <- 0.152
SACO_3_fst[5,20] <- 0.093

# row 6
SACO_3_fst[6,7]  <- 0.064
SACO_3_fst[6,8]  <- 0.054
SACO_3_fst[6,9]  <- 0.084
SACO_3_fst[6,10] <- 0.057
SACO_3_fst[6,11] <- 0.097
SACO_3_fst[6,12] <- 0.090
SACO_3_fst[6,13] <- -0.014
SACO_3_fst[6,14] <- 0.175
SACO_3_fst[6,15] <- 0.145
SACO_3_fst[6,16] <- 0.045
SACO_3_fst[6,17] <- 0.153
SACO_3_fst[6,18] <- 0.091
SACO_3_fst[6,19] <- 0.072
SACO_3_fst[6,20] <- 0.115

# row 7
SACO_3_fst[7,8]  <- -0.012
SACO_3_fst[7,9]  <- -0.011
SACO_3_fst[7,10] <- -0.014
SACO_3_fst[7,11] <- 0.016
SACO_3_fst[7,12] <- 0.013
SACO_3_fst[7,13] <- 0.034
SACO_3_fst[7,14] <- 0.045
SACO_3_fst[7,15] <- 0.104
SACO_3_fst[7,16] <- 0.015
SACO_3_fst[7,17] <- 0.147
SACO_3_fst[7,18] <- 0.044
SACO_3_fst[7,19] <- 0.039
SACO_3_fst[7,20] <- 0.033

# row 8
SACO_3_fst[8,9]  <- -0.013
SACO_3_fst[8,10] <- -0.018
SACO_3_fst[8,11] <- 0.014
SACO_3_fst[8,12] <- 0.010
SACO_3_fst[8,13] <- 0.028
SACO_3_fst[8,14] <- 0.047
SACO_3_fst[8,15] <- 0.095
SACO_3_fst[8,16] <- 0.010
SACO_3_fst[8,17] <- 0.137
SACO_3_fst[8,18] <- 0.038
SACO_3_fst[8,19] <- 0.029
SACO_3_fst[8,20] <- 0.032

# row 9
SACO_3_fst[9,10] <- -0.020
SACO_3_fst[9,11] <- -0.003
SACO_3_fst[9,12] <- -0.014
SACO_3_fst[9,13] <- 0.042
SACO_3_fst[9,14] <- 0.014
SACO_3_fst[9,15] <- 0.083
SACO_3_fst[9,16] <- 0.015
SACO_3_fst[9,17] <- 0.130
SACO_3_fst[9,18] <- 0.021
SACO_3_fst[9,19] <- 0.023
SACO_3_fst[9,20] <- -0.005

# row 10
SACO_3_fst[10,11] <- 0.002
SACO_3_fst[10,12] <- 0.010
SACO_3_fst[10,13] <- 0.028
SACO_3_fst[10,14] <- 0.045
SACO_3_fst[10,15] <- 0.082
SACO_3_fst[10,16] <- 0.011
SACO_3_fst[10,17] <- 0.128
SACO_3_fst[10,18] <- 0.026
SACO_3_fst[10,19] <- 0.028
SACO_3_fst[10,20] <- 0.024

# row 11
SACO_3_fst[11,12] <- 0.002
SACO_3_fst[11,13] <- 0.058
SACO_3_fst[11,14] <- 0.079
SACO_3_fst[11,15] <- 0.119
SACO_3_fst[11,16] <- 0.048
SACO_3_fst[11,17] <- 0.165
SACO_3_fst[11,18] <- 0.056
SACO_3_fst[11,19] <- 0.057
SACO_3_fst[11,20] <- 0.045

# row 12
SACO_3_fst[12,13] <- 0.049
SACO_3_fst[12,14] <- 0.036
SACO_3_fst[12,15] <- 0.122
SACO_3_fst[12,16] <- 0.040
SACO_3_fst[12,17] <- 0.138
SACO_3_fst[12,18] <- 0.060
SACO_3_fst[12,19] <- 0.049
SACO_3_fst[12,20] <- 0.015

# row 13
SACO_3_fst[13,14] <- 0.132
SACO_3_fst[13,15] <- 0.141
SACO_3_fst[13,16] <- 0.022
SACO_3_fst[13,17] <- 0.165
SACO_3_fst[13,18] <- 0.064
SACO_3_fst[13,19] <- 0.062
SACO_3_fst[13,20] <- 0.073

# row 14
SACO_3_fst[14,15] <- 0.097
SACO_3_fst[14,16] <- 0.079
SACO_3_fst[14,17] <- 0.119
SACO_3_fst[14,18] <- 0.079
SACO_3_fst[14,19] <- 0.048
SACO_3_fst[14,20] <- 0.021

# row 15
SACO_3_fst[15,16] <- 0.078
SACO_3_fst[15,17] <- 0.051
SACO_3_fst[15,18] <- 0.041
SACO_3_fst[15,19] <- 0.027
SACO_3_fst[15,20] <- 0.092

# row 16
SACO_3_fst[16,17] <- 0.156
SACO_3_fst[16,18] <- 0.015
SACO_3_fst[16,19] <- 0.038
SACO_3_fst[16,20] <- 0.042

# row 17
SACO_3_fst[17,18] <- 0.135
SACO_3_fst[17,19] <- 0.043
SACO_3_fst[17,20] <- 0.129

# row 18
SACO_3_fst[18,19] <- 0.039
SACO_3_fst[18,20] <- 0.025

# row 19
SACO_3_fst[19,20] <- 0.043

# -----------------------------
# 3) make symmetric and clean
# -----------------------------
SACO_3_fst[lower.tri(SACO_3_fst)] <- t(SACO_3_fst)[lower.tri(SACO_3_fst)]
SACO_3_fst[SACO_3_fst < 0] <- 0
diag(SACO_3_fst) <- 0

stopifnot(identical(dim(SACO_3_fst), c(20L, 20L)))
stopifnot(isTRUE(all.equal(SACO_3_fst, t(SACO_3_fst), check.attributes = FALSE)))

# -----------------------------
# 4) geographic distances (km)
# -----------------------------
coords <- SACO_3_coords[, c("lon", "lat")]
geo_dist_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- SACO_3_coords$site

# -----------------------------
# 5) IBD dataframe
# keep only pairs with non-missing FST
# -----------------------------
ibd_df <- data.frame(
  site1   = rownames(SACO_3_fst)[row(SACO_3_fst)[upper.tri(SACO_3_fst)]],
  site2   = colnames(SACO_3_fst)[col(SACO_3_fst)[upper.tri(SACO_3_fst)]],
  fst     = SACO_3_fst[upper.tri(SACO_3_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

ibd_df <- ibd_df[!is.na(ibd_df$fst), ]

# -----------------------------
# 6) map of sampling sites
# -----------------------------
world_map <- map_data("world")

map_df <- SACO_3_coords
map_df$label <- as.character(map_df$site)

p_map <- ggplot() +
  geom_polygon(
    data = world_map,
    aes(x = long, y = lat, group = group),
    fill = "grey92",
    color = "grey55",
    linewidth = 0.2
  ) +
  geom_point(
    data = map_df,
    aes(x = lon, y = lat),
    color = "red3",
    size = 3
  ) +
  geom_text(
    data = map_df,
    aes(x = lon, y = lat, label = label),
    nudge_y = 0.03,
    size = 3.3
  ) +
  coord_fixed(
    ratio = 1.3,
    xlim = range(map_df$lon) + c(-2.5, 2.5),
    ylim = range(map_df$lat) + c(-1.5, 1.5)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "SACO-3 sampling sites"
  )

# -----------------------------
# 7) IBD plot
# -----------------------------
p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Geographic distance (km)",
    y = expression(F[ST]),
    title = "SACO-3 isolation by distance"
  )

# -----------------------------
# 8) save RData
# -----------------------------
SACO_3_coords$stream <- NULL # remove site names

save(
  SACO_3_fst,
  SACO_3_coords,
  file = out_rdata
)


SACO_3_coords$stream <- NULL

# -----------------------------
# 9) print outputs
# -----------------------------
print(SACO_3_coords)
print(round(SACO_3_fst, 4))
print(head(ibd_df))
print(p_map)
print(p_ibd)
