# -----------------------------
# SEAT-2 Creek chub
# Skalski et al. 2008
# Table 1 coordinates + Table A.2 pairwise FST
# -----------------------------

library(dplyr)
library(geosphere)
library(ggplot2)
library(maps)

# -----------------------------
# paths
# -----------------------------
setwd("/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SEAT-2")

# -----------------------------
# helper: degrees + decimal minutes to decimal degrees
# -----------------------------
dm_to_dd <- function(deg, min, west = FALSE) {
  x <- deg + min / 60
  if (west) x <- -x
  x
}

# -----------------------------
# sites (site, lat, lon ONLY)
# exact coordinates from Table 1
# sites 1-8 = Eno River tributaries
# sites 9-16 = Falls Lake tributaries
# -----------------------------
SEAT_2_sites <- data.frame(
  site = 1:16,
  lat = c(
    dm_to_dd(36, 5.09),
    dm_to_dd(36, 5.09),
    dm_to_dd(36, 5.17),
    dm_to_dd(36, 2.72),
    dm_to_dd(36, 3.09),
    dm_to_dd(36, 3.74),
    dm_to_dd(36, 3.16),
    dm_to_dd(36, 4.61),
    dm_to_dd(36, 2.52),
    dm_to_dd(36, 0.28),
    dm_to_dd(36, 0.92),
    dm_to_dd(36, 0.77),
    dm_to_dd(36, 1.40),
    dm_to_dd(36, 3.13),
    dm_to_dd(35, 58.25),
    dm_to_dd(35, 59.26)
  ),
  lon = c(
    dm_to_dd(79, 2.24, west = TRUE),
    dm_to_dd(79, 2.04, west = TRUE),
    dm_to_dd(79, 1.22, west = TRUE),
    dm_to_dd(78, 59.30, west = TRUE),
    dm_to_dd(78, 59.11, west = TRUE),
    dm_to_dd(78, 57.82, west = TRUE),
    dm_to_dd(78, 55.55, west = TRUE),
    dm_to_dd(78, 53.75, west = TRUE),
    dm_to_dd(78, 40.07, west = TRUE),
    dm_to_dd(78, 40.60, west = TRUE),
    dm_to_dd(78, 38.60, west = TRUE),
    dm_to_dd(78, 37.46, west = TRUE),
    dm_to_dd(78, 37.30, west = TRUE),
    dm_to_dd(78, 34.79, west = TRUE),
    dm_to_dd(78, 39.99, west = TRUE),
    dm_to_dd(78, 37.48, west = TRUE)
  )
)

# -----------------------------
# FST matrix (Table A.2 EXACT)
# entries are pairwise FST averaged over 17 loci
# -----------------------------
SEAT_2_fst <- matrix(0, 16, 16)
rownames(SEAT_2_fst) <- colnames(SEAT_2_fst) <- 1:16

# fill lower triangle row by row from Table A.2
SEAT_2_fst[2, 1] <- 0.062

SEAT_2_fst[3, 1] <- 0.095
SEAT_2_fst[3, 2] <- 0.091

SEAT_2_fst[4, 1] <- 0.063
SEAT_2_fst[4, 2] <- 0.071
SEAT_2_fst[4, 3] <- 0.079

SEAT_2_fst[5, 1] <- 0.050
SEAT_2_fst[5, 2] <- 0.059
SEAT_2_fst[5, 3] <- 0.067
SEAT_2_fst[5, 4] <- 0.043

SEAT_2_fst[6, 1] <- 0.040
SEAT_2_fst[6, 2] <- 0.065
SEAT_2_fst[6, 3] <- 0.077
SEAT_2_fst[6, 4] <- 0.045
SEAT_2_fst[6, 5] <- 0.027

SEAT_2_fst[7, 1] <- 0.097
SEAT_2_fst[7, 2] <- 0.071
SEAT_2_fst[7, 3] <- 0.097
SEAT_2_fst[7, 4] <- 0.080
SEAT_2_fst[7, 5] <- 0.052
SEAT_2_fst[7, 6] <- 0.073

SEAT_2_fst[8, 1] <- 0.095
SEAT_2_fst[8, 2] <- 0.085
SEAT_2_fst[8, 3] <- 0.131
SEAT_2_fst[8, 4] <- 0.103
SEAT_2_fst[8, 5] <- 0.072
SEAT_2_fst[8, 6] <- 0.074
SEAT_2_fst[8, 7] <- 0.055

SEAT_2_fst[9, 1] <- 0.122
SEAT_2_fst[9, 2] <- 0.138
SEAT_2_fst[9, 3] <- 0.131
SEAT_2_fst[9, 4] <- 0.097
SEAT_2_fst[9, 5] <- 0.110
SEAT_2_fst[9, 6] <- 0.099
SEAT_2_fst[9, 7] <- 0.121
SEAT_2_fst[9, 8] <- 0.099

SEAT_2_fst[10, 1] <- 0.092
SEAT_2_fst[10, 2] <- 0.101
SEAT_2_fst[10, 3] <- 0.114
SEAT_2_fst[10, 4] <- 0.063
SEAT_2_fst[10, 5] <- 0.083
SEAT_2_fst[10, 6] <- 0.059
SEAT_2_fst[10, 7] <- 0.102
SEAT_2_fst[10, 8] <- 0.106
SEAT_2_fst[10, 9] <- 0.073

SEAT_2_fst[11, 1] <- 0.126
SEAT_2_fst[11, 2] <- 0.133
SEAT_2_fst[11, 3] <- 0.135
SEAT_2_fst[11, 4] <- 0.112
SEAT_2_fst[11, 5] <- 0.121
SEAT_2_fst[11, 6] <- 0.099
SEAT_2_fst[11, 7] <- 0.137
SEAT_2_fst[11, 8] <- 0.139
SEAT_2_fst[11, 9] <- 0.116
SEAT_2_fst[11, 10] <- 0.095

SEAT_2_fst[12, 1] <- 0.138
SEAT_2_fst[12, 2] <- 0.139
SEAT_2_fst[12, 3] <- 0.111
SEAT_2_fst[12, 4] <- 0.112
SEAT_2_fst[12, 5] <- 0.136
SEAT_2_fst[12, 6] <- 0.106
SEAT_2_fst[12, 7] <- 0.151
SEAT_2_fst[12, 8] <- 0.150
SEAT_2_fst[12, 9] <- 0.091
SEAT_2_fst[12, 10] <- 0.091
SEAT_2_fst[12, 11] <- 0.084

SEAT_2_fst[13, 1] <- 0.113
SEAT_2_fst[13, 2] <- 0.109
SEAT_2_fst[13, 3] <- 0.098
SEAT_2_fst[13, 4] <- 0.067
SEAT_2_fst[13, 5] <- 0.089
SEAT_2_fst[13, 6] <- 0.075
SEAT_2_fst[13, 7] <- 0.107
SEAT_2_fst[13, 8] <- 0.100
SEAT_2_fst[13, 9] <- 0.059
SEAT_2_fst[13, 10] <- 0.062
SEAT_2_fst[13, 11] <- 0.046
SEAT_2_fst[13, 12] <- 0.049

SEAT_2_fst[14, 1] <- 0.091
SEAT_2_fst[14, 2] <- 0.107
SEAT_2_fst[14, 3] <- 0.101
SEAT_2_fst[14, 4] <- 0.083
SEAT_2_fst[14, 5] <- 0.106
SEAT_2_fst[14, 6] <- 0.080
SEAT_2_fst[14, 7] <- 0.119
SEAT_2_fst[14, 8] <- 0.094
SEAT_2_fst[14, 9] <- 0.055
SEAT_2_fst[14, 10] <- 0.067
SEAT_2_fst[14, 11] <- 0.051
SEAT_2_fst[14, 12] <- 0.058
SEAT_2_fst[14, 13] <- 0.024

SEAT_2_fst[15, 1] <- 0.115
SEAT_2_fst[15, 2] <- 0.116
SEAT_2_fst[15, 3] <- 0.142
SEAT_2_fst[15, 4] <- 0.079
SEAT_2_fst[15, 5] <- 0.101
SEAT_2_fst[15, 6] <- 0.085
SEAT_2_fst[15, 7] <- 0.123
SEAT_2_fst[15, 8] <- 0.112
SEAT_2_fst[15, 9] <- 0.073
SEAT_2_fst[15, 10] <- 0.060
SEAT_2_fst[15, 11] <- 0.042
SEAT_2_fst[15, 12] <- 0.095
SEAT_2_fst[15, 13] <- 0.041
SEAT_2_fst[15, 14] <- 0.051

SEAT_2_fst[16, 1] <- 0.133
SEAT_2_fst[16, 2] <- 0.139
SEAT_2_fst[16, 3] <- 0.175
SEAT_2_fst[16, 4] <- 0.119
SEAT_2_fst[16, 5] <- 0.134
SEAT_2_fst[16, 6] <- 0.104
SEAT_2_fst[16, 7] <- 0.144
SEAT_2_fst[16, 8] <- 0.128
SEAT_2_fst[16, 9] <- 0.075
SEAT_2_fst[16, 10] <- 0.086
SEAT_2_fst[16, 11] <- 0.061
SEAT_2_fst[16, 12] <- 0.083
SEAT_2_fst[16, 13] <- 0.071
SEAT_2_fst[16, 14] <- 0.061
SEAT_2_fst[16, 15] <- 0.048

# mirror + clean
SEAT_2_fst[upper.tri(SEAT_2_fst)] <- t(SEAT_2_fst)[upper.tri(SEAT_2_fst)]
SEAT_2_fst[SEAT_2_fst < 0] <- 0
diag(SEAT_2_fst) <- 0

# -----------------------------
# pairwise geographic distances for IBD
# note: paper used shortest in-water path distance;
# here Euclidean distance is used for workflow consistency
# -----------------------------
coords <- SEAT_2_sites %>% select(lon, lat)
geo_dist_km <- geosphere::distm(coords, fun = geosphere::distHaversine) / 1000
rownames(geo_dist_km) <- colnames(geo_dist_km) <- SEAT_2_sites$site

ibd_df <- data.frame(
  site1   = rownames(SEAT_2_fst)[row(SEAT_2_fst)[upper.tri(SEAT_2_fst)]],
  site2   = colnames(SEAT_2_fst)[col(SEAT_2_fst)[upper.tri(SEAT_2_fst)]],
  fst     = SEAT_2_fst[upper.tri(SEAT_2_fst)],
  dist_km = geo_dist_km[upper.tri(geo_dist_km)]
)

# -----------------------------
# map of sampling sites
# -----------------------------
world_map <- map_data("world")
usa_map   <- subset(world_map, region %in% c("USA"))

map_df <- SEAT_2_sites %>%
  mutate(label = as.character(site))

x_pad <- 0.12
y_pad <- 0.12

p_map <- ggplot() +
  geom_polygon(
    data = usa_map,
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
    nudge_y = 0.015,
    size = 3.5
  ) +
  coord_fixed(
    ratio = 1.3,
    xlim = range(map_df$lon) + c(-x_pad, x_pad),
    ylim = range(map_df$lat) + c(-y_pad, y_pad)
  ) +
  theme_classic() +
  labs(
    x = "Longitude",
    y = "Latitude",
    title = "SEAT-2 sampling sites"
  )

# -----------------------------
# IBD plot
# -----------------------------
p_ibd <- ggplot(ibd_df, aes(x = dist_km, y = fst)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  labs(
    x = "Geographic distance (km)",
    y = expression(F[ST]),
    title = "SEAT-2 isolation by distance"
  )

# -----------------------------
# save RData
# -----------------------------
save(
  SEAT_2_fst,
  SEAT_2_sites,
  file = "/Users/johnmccall/Library/CloudStorage/OneDrive-TheOhioStateUniversity/Spring_2026/Landgen_DGS/DGS-Patterns-of-Fish-Diversity/DGS-Patterns-of-Fish-Diversity/Objective 2/Data/SEAT-2/data/SEAT-2.RData"
)

# -----------------------------
# print outputs
# -----------------------------
print(SEAT_2_sites)
print(round(SEAT_2_fst, 4))
print(head(ibd_df))
print(p_map)
print(p_ibd)
